# 연구 방향 및 논문화 전략

## 1. 현재 상태 진단

### 보유 자산

| 생성기 | 설계 공간 | Ground State 보장 | SA-hard | 구별 불가능 |
|--------|----------|------------------|---------|-----------|
| `qubo_wishart.py` | Ising → QUBO | 수학적 보장 ($W^T t = 0$) | O | X |
| `qubo_zero_expectation.py` | QUBO 직접 | LP 최적화 페널티 | X | O |
| `qubo_hard_mode.py` | Ising → QUBO | 조건부 | X (SA 100%) | X |

### 핵심 딜레마

SA-hard와 구별 불가능, 두 성질이 근본적으로 충돌한다:

| 성질 | 요구 조건 |
|------|----------|
| SA-hard | J 행렬이 low-rank → Q 계수들이 **서로 상관** |
| 구별 불가능 | Q 계수가 **독립적이고 무작위** |

- **Wishart**: SA가 실패하지만(alpha < 1.0에서 성공률 0%), Q 행렬 분석으로 planted 문제임이 탐지 가능
- **Zero-Expectation**: Q 행렬이 랜덤과 구별 불가능하지만, SA가 100% 성공

벤치마크로서 "SA가 잘 푸는 문제"는 솔버 평가 가치가 약하고, "구별 가능한 문제"는 부정행위(Q 구조 분석 후 전용 알고리즘으로 풀기)에 취약하다.

---

## 2. 논문화를 위한 강점과 약점

### 강점

1. **명확한 문제 정의**: 정답이 보장된 QUBO 벤치마크 — 양자 우위 입증의 전제 조건
2. **체계적 실험**: Alpha sweep, N scaling, num_reads 분석, Hard Mode 비교 완비
3. **상전이 프로파일**: alpha_c ≈ 0.95 (N=100)에서 1차 상전이 실험적 확인
4. **두 성질의 독립성 분석**: 구별 불가능 vs SA-hard의 충돌 관계 정리

### 약점

1. **Wishart 구현 자체는 새롭지 않음**: Hamze et al. (2020) 재구현
2. **실제 QPU 실험 없음**: SA만으로는 "양자 어닐링 벤치마크" 논문이라 하기 부족
3. **Zero-Expectation의 이론적 검증 부족**: "구별 불가능" 주장에 대한 formal test (KS test, chi-square 등) 미비
4. **기존 벤치마크와의 비교 없음**: Chimera-native planted, random K-SAT 등과 비교 부재

---

## 3. 연구 방향

### 방향 A: 벤치마크 스위트 논문 (가장 현실적)

> **"QUBO Benchmark Suite with Verified Ground States for Quantum Annealing Evaluation"**

난이도 스펙트럼 전체를 제공하는 통합 벤치마크 스위트로 포지셔닝:

```
Easy ◄──────────────────────────────► Hard
Zero-Expectation    Hard Mode    Wishart(α=0.7)
(SA 100%)           (SA 100%)    (SA 0%)
구별 불가능          구별 가능      구별 가능
```

**핵심 주장**: 솔버 평가는 단일 난이도가 아니라 **난이도 프로파일** 전체에서 해야 한다. 쉬운 문제도 못 푸는 솔버와 어려운 문제까지 푸는 솔버를 구분하는 것이 벤치마크의 가치.

**추가 작업**:

| 우선순위 | 항목 | 내용 |
|---------|------|------|
| 필수 | D-Wave QPU 실험 | SA 실패 + QPU 성공 시연 |
| 필수 | 데이터셋 공개 | 재현성 + 커뮤니티 기여 |
| 높음 | 기존 벤치마크 비교 | Chimera-native, random K-SAT 등 |
| 중간 | QAOA 실험 | 다양한 솔버 유용성 시연 |

**타겟 저널**: *Scientific Data*, *Quantum Science and Technology*

**장점**: 이미 세 가지 생성기 + 실험 데이터 보유. QPU 실험 추가만으로 논문화 가능.

---

### 방향 B: SA-hard + 구별 불가능 결합의 정량적 한계 분석 (novelty 높음)

> **"Toward Indistinguishable Hard QUBO: Quantifying the Trade-off Between Computational Hardness and Statistical Detectability"**

"결합이 안 된다"가 아니라 **"어디까지 결합 가능한가"를 정량적으로 분석**.

**실험 설계**:

1. Wishart Q에 Zero-Expectation noise를 비율 λ로 혼합:
   - $Q_{mixed} = (1 - \lambda) \cdot Q_{Wishart} + \lambda \cdot Q_{ZeroExp}$
2. 3축 트레이드오프 측정:
   - λ vs SA 성공률
   - λ vs 통계적 탐지율 (KS test, chi-square)
   - λ vs ground state 보존율
3. 파레토 곡선 도출:
   - "구별 불가능성을 X% 포기하면 SA 성공률이 Y%로 떨어진다"

**추가 작업**:

| 우선순위 | 항목 | 내용 |
|---------|------|------|
| 필수 | 혼합 모델 구현 | Wishart + ZeroExp 비율 조절 |
| 필수 | 통계적 검정 구현 | KS test, chi-square, 고유값 분석 |
| 필수 | 파레토 곡선 실험 | λ sweep × N sweep |
| 높음 | 이론적 하한 분석 | 정보이론적 관점에서 결합 불가능성 경계 |

**타겟 저널**: *Physical Review E*, *Journal of Statistical Mechanics*

**장점**: 열린 연구 문제에 대한 실험적 경계 제시. Novelty 명확.
**단점**: 시간 소요 큼. 이론적 분석 난이도 높음.

---

### 방향 C: 솔버 간 비대칭 난이도 연구 (발견 의존적)

> **"Solver-Dependent Hardness of Planted QUBO: When SA Succeeds but QPU Fails (or Vice Versa)"**

Zero-Expectation QUBO가 SA에는 쉽지만 **다른 솔버에는 어려울 수 있다**는 가설 검증:

- SA-easy이면서 QPU-hard인 문제가 존재하는가?
- SA-easy이면서 QAOA-hard인 문제가 존재하는가?
- 반대 방향(SA-hard, QPU-easy)은 Wishart에서 이미 기대되는 시나리오

**실험 설계**:

| 생성기 | SA | QPU | QAOA | 해석 |
|--------|----|----|------|------|
| Zero-Expectation | Easy | ? | ? | 구별 불가능 문제의 양자 난이도 |
| Wishart (α=0.7) | Hard | ? | ? | 양자 우위 후보 |
| Hard Mode | Easy | ? | ? | 기준선 |

**추가 작업**:

| 우선순위 | 항목 | 내용 |
|---------|------|------|
| 필수 | D-Wave QPU 실험 | 세 생성기 모두 QPU에서 실험 |
| 필수 | QAOA 실험 | Qiskit으로 회로 깊이별 성능 측정 |
| 높음 | 난이도 메트릭 통일 | TTS(99%) 기준으로 솔버 간 비교 |

**타겟 저널**: *Nature Physics* (양자 우위 입증 시), *Physical Review Letters*

**장점**: 발견이 나오면 임팩트 매우 높음. 양자 우위 직접 증거.
**단점**: QPU 접근 필요. 결과가 "모든 솔버가 쉽게 풂"이면 negative result.

---

### 방향 D: Zero-Expectation 조건 하에서 난이도 올리기 (가장 어렵지만 임팩트 최대)

> **"Computationally Hard QUBO Indistinguishable from Random: A New Construction"**

E[q_ij] = 0을 유지하면서 SA가 실패하는 QUBO를 구성하는 새로운 방법 개발.

**가능한 접근**:

| 방법 | 아이디어 | 기댓값 0 유지 | SA-hard 가능성 |
|------|---------|-------------|--------------|
| LP 자유도 활용 | 페널티 비율의 남은 자유도로 local minima 심기 | O | 미지수 |
| 고차 상관 주입 | 기댓값은 0이되 분산/공분산에 구조 삽입 | O | 가능 |
| Sparse density | density를 낮춰 탐색 공간 복잡화 | O | 낮음 |
| 조건부 Wishart | Wishart 구성 후 기댓값 0으로 후보정 | 부분적 | 감소 |

**타겟 저널**: *Physical Review Letters*, *Journal of the ACM* (이론적 증명 포함 시)

**장점**: 성공 시 해당 분야의 열린 문제를 해결하는 것. 최고 임팩트.
**단점**: 불가능할 수도 있음. 연구 리스크 높음.

---

## 4. 방향별 비교

| 기준 | A: 스위트 | B: 결합 한계 | C: 비대칭 난이도 | D: 새 구성법 |
|------|----------|------------|----------------|------------|
| 현실성 | **높음** | 중간 | 중간 | 낮음 |
| Novelty | 중간 | **높음** | 발견 의존 | **최고** |
| 추가 작업량 | **적음** | 많음 | 중간 (QPU 필요) | 매우 많음 |
| 리스크 | **낮음** | 낮음 | 중간 | 높음 |
| 임팩트 | 중간 | 중간~높음 | **높음** (발견 시) | **최고** |
| QPU 필요 | O | X | **O (필수)** | X |

---

## 5. 권장 로드맵

### Phase 1: 방향 A (벤치마크 스위트) — 단기 목표

1. Zero-Expectation에 대한 formal 통계 검정 추가 (KS test, chi-square)
2. 데이터셋 공개 준비 (다양한 N, alpha, 생성기별)
3. QPU 실험 (D-Wave Leap 클라우드 이용 가능)
4. 논문 작성 → 투고

### Phase 2: 방향 B (결합 한계 분석) — 중기 목표

1. Phase 1 결과를 기반으로 혼합 모델 실험
2. 파레토 곡선 도출
3. 이론적 하한 분석 시도
4. 후속 논문 작성

### Phase 3: 방향 C 또는 D — 장기 목표

- QPU 결과에 따라 C (양자 우위 입증) 또는 D (새 구성법) 선택
