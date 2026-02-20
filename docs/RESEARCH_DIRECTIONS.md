# 연구 방향 및 논문화 전략

## 1. 현재 상태 진단

### 보유 자산 (7개 방법론)

| 생성기 | GS 보장 | SA-hard | 보조변수 | 구별 불가능 | 난이도 조절 | 상태 |
|--------|:-------:|:-------:|:--------:|:----------:|:----------:|:----:|
| Wishart | 수학적 | **O** (alpha<0.95) | 없음 | X (low-rank) | alpha | 완성 |
| McEliece | 조건부 (M 의존) | 미측정 | 있음 (대량) | 미분석 | m, t | 완성 |
| Quiet Planting | 조건부 (field 필요) | 중간 | 있음 | **O** (alpha<3.86) | alpha, field | 완성 |
| Posiform | **수학적 (유일)** | X (SA-trivial) | 없음 | 미분석 | X | 완성 |
| Posiform Hardened | **수학적 (유일)** | 조절 가능 | 없음 | X (block-diagonal) | alpha, coeff_type | 완성 |
| Zero Expectation | 수학적 | X (SA-trivial) | 없음 | **O** (E[q]=0) | density | 완성 |
| Hard Mode | 조건부 | X (SA-trivial) | 없음 | X (backbone) | noise_ratio | 완성 |

### 난이도 스펙트럼

```
SA-trivial ◄─────────────────────────────────────────────► SA-hard
Posiform   ZeroExp   HardMode   Hardened(α=0.1)   Quiet(f=0.5)   Wishart(α=0.7)
(100%)     (100%)    (100%)     (조절가능)          (N↑→0%)        (N≥100→0%)
                                                                    │
                                     McEliece ────── 미측정 (이론: O(2^k) PT)
```

### 핵심 딜레마: 3중 트레이드오프

SA-hard, 구별 불가능, 수학적 GS 보장 — 이 세 가지를 동시에 달성하는 방법론은 아직 없다:

| 성질 | 요구 조건 | 충돌 요인 |
|------|----------|----------|
| SA-hard | J 행렬이 low-rank → Q 계수들이 **서로 상관** | 구별 불가능과 충돌 |
| 구별 불가능 | Q 계수가 **독립적이고 무작위** | SA-hard와 충돌 |
| 수학적 GS 보장 | 구조적 제약 (양수 계수, 직교 투영 등) | 완전 무작위와 충돌 |

현재 가장 가까운 조합:

| 방법론 | SA-hard | 구별 불가능 | GS 보장 |
|--------|:-------:|:----------:|:-------:|
| Wishart | **O** | X | **O** |
| Quiet Planting | 중간 | **O** | 조건부 |
| Posiform Hardened | 조절 가능 | X | **O** |
| Zero Expectation | X | **O** | **O** |

---

## 2. 방법론별 미완 과제

### 2.1 Posiform Hardened — SA 난이도 체계적 검증

**현재 상태**: 생성기 완성, SA 실험 프레임워크 완성, 실제 실험 미수행.

**필요 작업**:

| 우선순위 | 항목 | 내용 | 예상 시간 |
|:--------:|------|------|:---------:|
| **필수** | Sweep 전이 실험 | N=500, sweep 1~10000 × (lin2/lin20) × (α=0.01/0.1) | 수 시간 |
| **필수** | N 스케일링 실험 | N=50~2000, sweeps=500 고정, 4개 config | 수 시간 |
| **필수** | Hardened vs Plain 비교 | 동일 N/sweep에서 hardened가 plain보다 어려운지 | 1시간 |
| 높음 | α 임계점 탐색 | α를 세밀하게 sweep하여 easy-hard 전이점 측정 | 2시간 |
| 중간 | Block-diagonal 구조 위장 | Kernighan-Lin 분할로 구조 은닉 가능성 검토 | 연구 필요 |

**기대 결과**: α=0.01에서 SA가 실패하고, α 증가에 따라 S-curve 형태의 전이가 관찰될 것. 논문 Fig 7-8과 정량 비교 가능.

### 2.2 McEliece — Eq.14 Exact Decomposition + QPU 실험

**현재 상태**: 페널티 M 강화 완료 (이론적 하한 기반 + 사후 검증). m=3에서 160/160, m=4에서 10/10 brute force 검증 성공.

**필요 작업**:

| 우선순위 | 항목 | 내용 | 예상 시간 |
|:--------:|------|------|:---------:|
| **필수** | Eq.14 exact decomposition | 페널티 M 의존성 완전 제거. 논문 Section 2의 identity 구현 | 연구 필요 |
| 높음 | p-local PT 실험 재현 | 논문의 핵심: p-local Ising에서 O(2^k) parallel tempering 스케일링 | QPU 필요 |
| 높음 | SA 스케일링 실험 | m 증가에 따른 SA 성공률 변화 측정 (m=3,4,5,6) | 수 시간 |
| 중간 | Stern 알고리즘 비교 | 고전 최적 공격 복잡도와의 비교 | 연구 필요 |
| 낮음 | 대규모 실험 (m≥5) | 현재 m=4까지만 실험. m=5 이상에서의 QUBO 크기/품질 평가 | 수 시간 |

**핵심 이슈**: Eq.14는 p-body 항을 보조변수 없이 2-body로 분해하는 exact identity. 현재 Rosenberg 방식은 보조변수 + 페널티 M이 필요하여 QUBO가 크고, M 부족 시 GS가 깨질 수 있음. Eq.14 구현이 되면 McEliece의 가치가 크게 올라감.

### 2.3 Quiet Planting — Field 최적화 + 대규모 실험

**현재 상태**: SA 스케일링 완료 (field=0.5: N=100→100%, N=500→20%, N=1000→0%).

**필요 작업**:

| 우선순위 | 항목 | 내용 |
|:--------:|------|------|
| 높음 | Field 강도 최적화 | field vs SA 성공률 트레이드오프 정밀 측정 |
| 높음 | Alpha 상전이 실험 | alpha=3.0~5.0에서 SA 성공률 변화 (SAT 상전이 근처) |
| 중간 | 보조변수 오버헤드 분석 | QUBO 크기 = n(1+alpha)의 실질적 영향 정량화 |

### 2.4 Wishart — QPU 실험 + 대규모 스케일링

**현재 상태**: SA 실험 완성 (alpha sweep, N scaling, 에너지 랜드스케이프 분석).

**필요 작업**:

| 우선순위 | 항목 | 내용 |
|:--------:|------|------|
| **필수** | D-Wave QPU 실험 | SA-hard가 QPU에서도 hard인지 확인 |
| 높음 | N=1000+ 스케일링 | alpha_c(N)의 N 의존성 측정 |
| 중간 | TTS(99%) 메트릭 | Time-To-Solution으로 솔버 간 공정 비교 |

---

## 3. 연구 방향

### 방향 A: 7-방법론 통합 벤치마크 스위트 논문 (가장 현실적)

> **"QUBO Benchmark Suite with Verified Ground States: A Comprehensive Toolkit for Quantum Annealing Evaluation"**

7개 방법론을 난이도/보장/은닉성 축으로 체계화한 통합 벤치마크 스위트.

**핵심 기여**:
1. 난이도 스펙트럼 전체를 커버하는 최초의 통합 벤치마크
2. 각 방법론의 SA 난이도를 동일 조건에서 정량 비교
3. 사용 시나리오별 최적 생성기 권장 가이드

**필요 추가 실험**:

| 우선순위 | 항목 | 내용 | 담당 방법론 |
|:--------:|------|------|:----------:|
| **필수** | Hardened SA 실험 | sweep 전이 + N 스케일링 + vs Plain | Posiform Hardened |
| **필수** | McEliece SA 실험 | m별 SA 성공률 스케일링 | McEliece |
| **필수** | 통합 비교 실험 | 7개 방법론 동일 N/sweeps에서 비교 | 전체 |
| 높음 | D-Wave QPU 실험 | 최소 Wishart + Hardened + Posiform | Wishart, Hardened |
| 높음 | 데이터셋 공개 | 재현성 + 커뮤니티 기여 | 전체 |
| 중간 | QAOA 실험 | 회로 깊이별 성능 (Qiskit) | 전체 |

**논문 구조 (예상)**:
```
1. Introduction — QUBO 벤치마크의 필요성
2. Background — QUBO, Ising, SA, QA
3. Methods — 7개 방법론 각각의 원리와 GS 보장
4. Experiments
   4.1 SA 난이도 비교 (7-way, 동일 조건)
   4.2 N 스케일링 (easy→hard 순서로 실패 시점)
   4.3 QPU 실험 (Wishart, Hardened 중심)
5. Discussion — 사용 시나리오, 트레이드오프
6. Conclusion + 데이터셋 공개
```

**타겟 저널**: *Scientific Data*, *Quantum Science and Technology*

**현실성**: **높음** — 대부분 실험 프레임워크가 이미 구현됨. SA 실험 실행 + QPU 접근이 주요 추가 작업.

---

### 방향 B: SA-hard + 구별 불가능 + GS 보장 결합의 정량적 한계 분석

> **"Toward Indistinguishable Hard QUBO: Quantifying the Three-Way Trade-off Between Hardness, Detectability, and Ground State Guarantee"**

3중 트레이드오프의 정량적 경계를 실험적으로 탐색.

**실험 설계**:

1. **Wishart + Zero-Expectation 혼합 모델**:
   ```
   Q_mixed = (1 - λ) × Q_Wishart + λ × Q_ZeroExp
   ```
   - λ sweep (0.0 ~ 1.0) × N sweep
   - 3축 동시 측정: SA 성공률 / 통계적 탐지율 (KS test) / GS 보존율

2. **Quiet Planting + Wishart 구조 비교**:
   - 둘 다 planted model이지만 은닉성 메커니즘이 다름
   - alpha < 3.86 (Quiet, 구별 불가)과 alpha < 0.95 (Wishart, SA-hard)의 교차점 탐색

3. **Posiform Hardened의 구별 불가능성 개선**:
   - Block-diagonal 구조를 위장하는 방법 연구
   - Random inter-subgraph coupling 추가 시 GS 보장이 유지되는 범위 측정

**파레토 곡선 도출**:
```
축 1: SA 성공률 (%) — 낮을수록 hard
축 2: 통계적 탐지율 (%) — 낮을수록 은닉
축 3: GS 보존율 (%) — 높을수록 신뢰

→ "구별 불가능성을 X% 포기하면 SA 성공률이 Y%로 떨어진다"
→ "SA-hard를 Z% 포기하면 GS 보장 신뢰도가 W%까지 올라간다"
```

**추가 작업**:

| 우선순위 | 항목 | 내용 |
|:--------:|------|------|
| 필수 | 혼합 모델 구현 | Wishart + ZeroExp 비율 조절 생성기 |
| 필수 | 통계적 검정 구현 | KS test, chi-square, 고유값 분포 분석 |
| 필수 | 파레토 곡선 실험 | λ sweep × N sweep |
| 높음 | 이론적 하한 분석 | 정보이론적 관점에서 결합 불가능성 경계 |

**타겟 저널**: *Physical Review E*, *Journal of Statistical Mechanics*

**현실성**: **중간** — 새로운 생성기 구현 + 대규모 실험 필요. 이론적 분석은 난이도 높음.

---

### 방향 C: 솔버 간 비대칭 난이도 연구

> **"Solver-Dependent Hardness of Planted QUBO: When SA Succeeds but QPU Fails (or Vice Versa)"**

같은 QUBO가 솔버마다 다른 난이도를 보이는 현상 연구.

**실험 설계**:

| 생성기 | SA | QPU | QAOA | 예상 시나리오 |
|--------|:--:|:---:|:----:|:----------:|
| Wishart (α=0.7) | **Hard** | ? | ? | QPU-easy라면 양자 우위 |
| Posiform Hardened (α=0.01) | 조절 가능 | ? | ? | QPU vs SA 교차점 |
| Quiet Planting (field=0.5) | 중간 | ? | ? | 은닉 문제의 양자 난이도 |
| Zero-Expectation | Easy | ? | ? | 구별 불가능 문제의 양자 난이도 |
| McEliece | 미측정 | ? | ? | 암호학적 난이도의 양자 한계 |

**핵심 질문**:
- SA-hard + QPU-easy 문제가 존재하는가? (= 양자 우위 증거)
- SA-easy + QPU-hard 문제가 존재하는가? (= 양자 열위 사례)
- Posiform Hardened에서 α를 조절하면 QPU-SA 교차점이 보이는가?

**타겟 저널**: *Nature Physics* (양자 우위 입증 시), *Physical Review Letters*

**현실성**: **중간** — QPU 접근이 필수. D-Wave Leap 클라우드로 가능하지만 비용 발생. 결과가 negative일 수도 있음.

---

### 방향 D: Posiform Hardened 확장 — 구별 불가능한 Hard QUBO

> **"Computationally Hard QUBO Indistinguishable from Random: Hiding Structure in Hardened Posiform"**

Posiform Hardened의 block-diagonal 구조를 위장하여 구별 불가능성을 확보하는 새로운 구성법.

**현재 한계**: Random QUBO R_i가 disjoint subgraph에만 정의되어 Q 행렬에 block-diagonal 구조가 노출됨.

**가능한 접근**:

| 방법 | 아이디어 | GS 보장 | 구별 불가능 | 난이도 유지 |
|------|---------|:-------:|:----------:|:----------:|
| Inter-subgraph noise | subgraph 간 작은 랜덤 coupling 추가 | 약화 가능 | 개선 | 유지 |
| Overlapping subgraphs | disjoint 대신 겹치는 subgraph 사용 | 연구 필요 | 개선 | 연구 필요 |
| Global random QUBO + Posiform | 단일 global random Q에 posiform overlay | **보장** (논문과 동일) | **가능** | 연구 필요 |
| Wishart base + Posiform overlay | Wishart Q에 posiform을 overlay | 연구 필요 | low-rank 문제 | **높음** |

**가장 유망**: Global random QUBO + Posiform — 현재 Hardened는 disjoint subgraph 위의 random QUBO를 사용하지만, 전체 변수에 걸친 단일 random QUBO를 base로 사용하고 그 위에 posiform을 overlay하면 block-diagonal 구조가 사라짐. 단, subproblem의 exact GS를 알아야 하므로 전체 random QUBO의 GS를 구해야 하는 문제가 있음 (현재는 subgraph 분할로 brute force 가능).

**타겟 저널**: *Physical Review Letters*, *npj Unconventional Computing*

**현실성**: **낮음** — 새로운 이론적 구성이 필요. 성공 시 임팩트 매우 높음.

---

### 방향 E: McEliece Eq.14 Exact Decomposition + 암호학적 벤치마크

> **"Cryptographically Hard QUBO Without Auxiliary Variables: Implementing Exact Degree Reduction for McEliece-Based Benchmarks"**

McEliece의 핵심 미구현 기능 (Eq.14 exact decomposition)을 완성하여 보조변수 없는 암호학적 QUBO를 실현.

**현재 상태**: Rosenberg reduction → 보조변수 O(kN)개 → QUBO가 매우 큼 → 실용성 제한

**Eq.14 구현 시 개선**:

| 항목 | 현재 (Rosenberg) | Eq.14 구현 후 |
|------|:----------------:|:------------:|
| 보조변수 | O(kN)개 | **0개** |
| QUBO 크기 | k + O(kN) | **k** |
| GS 보장 | 조건부 (M 의존) | **수학적** |
| 확장성 | m=4까지 실용적 | m=6+ 가능 |

**추가 작업**:

| 우선순위 | 항목 | 내용 |
|:--------:|------|------|
| 필수 | Eq.14 identity 이해 | 논문 Section 2의 exact decomposition 수학적 분석 |
| 필수 | 구현 + 검증 | brute force로 모든 가능한 입력에서 등가성 확인 |
| 높음 | p-local PT 실험 | O(2^k) 스케일링 재현 |
| 높음 | Wishart와 비교 | 동일 크기에서 SA/QPU 성공률 비교 |

**타겟 저널**: *Quantum Science and Technology*, *IEEE Transactions on Quantum Engineering*

**현실성**: **중간** — Eq.14의 수학적 구현이 핵심 난이도. 성공 시 McEliece가 현재 프로젝트에서 가장 가치 있는 방법론이 될 수 있음.

---

## 4. 방향별 비교

| 기준 | A: 통합 스위트 | B: 3중 트레이드오프 | C: 비대칭 난이도 | D: 구별 불가능 Hard | E: McEliece Exact |
|------|:-----------:|:----------------:|:--------------:|:----------------:|:----------------:|
| **현실성** | **높음** | 중간 | 중간 | 낮음 | 중간 |
| **Novelty** | 중간 | **높음** | 발견 의존 | **최고** | 높음 |
| **추가 작업량** | **적음** | 많음 | 중간 | 매우 많음 | 많음 |
| **리스크** | **낮음** | 낮음 | 중간 | 높음 | 중간 |
| **임팩트** | 중간 | 중간~높음 | **높음** (발견 시) | **최고** | 높음 |
| **QPU 필요** | O (권장) | X | **O (필수)** | X | O (권장) |
| **코드 변경** | 실험 실행만 | 새 생성기 | 실험 실행만 | 새 구성법 | 기존 코드 수정 |

---

## 5. 권장 로드맵

### Phase 1: 즉시 실행 가능한 실험 (1~2주)

**목표**: 방향 A 논문의 데이터 수집

1. **Posiform Hardened SA 실험 실행**
   ```bash
   python3 posiform_hardened/run_experiments.py sweep     # Sweep 전이
   python3 posiform_hardened/run_experiments.py scaling   # N 스케일링 + vs Plain
   ```

2. **McEliece SA 스케일링 실험**
   - m=3,4에서 SA 성공률 측정 (현재 GS 검증만 완료)
   - N 스케일링: k=4~20 범위

3. **7-way 통합 비교 실험**
   - 동일 조건 (N=100,500,1000 / sweeps=1000 / reads=100)
   - 7개 방법론 전부 SA로 비교

### Phase 2: 논문 작성 + QPU 실험 (1~2개월)

**목표**: 방향 A 논문 완성

1. D-Wave Leap 클라우드로 QPU 실험 (Wishart + Hardened 중심)
2. 데이터셋 공개 준비 (GitHub + Zenodo)
3. 논문 작성 + 투고

### Phase 3: 심화 연구 (3~6개월)

**목표**: 방향 B 또는 E

- Phase 1 결과에서 Posiform Hardened이 SA-hard임이 확인되면 → **방향 D** (구별 불가능화)
- McEliece Eq.14 이해에 진전이 있으면 → **방향 E** (exact decomposition)
- QPU 결과가 흥미로우면 → **방향 C** (양자 우위 증거)
- 이론적 관심이 있으면 → **방향 B** (3중 트레이드오프)

### Phase 4: 장기 목표

- 방향 D: 구별 불가능 + SA-hard + GS 보장 — 열린 문제 해결
- 방향 C: 양자 우위 직접 증거 — 최고 임팩트

---

## 6. 참고 문헌

1. Hamze, F., et al. "Wishart planted ensemble." *Physical Review E*, 101, 052102, 2020.
2. Mandra, S., et al. "McEliece-based cryptographic QUBO." 2025.
3. Krzakala, F. & Zdeborova, L. "Hiding quiet solutions in random CSP." *PRL*, 102, 238701, 2009.
4. Hahn, G., Pelofske, E. & Djidjev, H. "Using 2-SAT to generate QUBO instances." *QCE*, 2023.
5. Pelofske, E., Hahn, G. & Djidjev, H. "Increasing the hardness of posiform planting." *npj Unconventional Computing*, 2025.
6. Boros, E. & Hammer, P. L. "Pseudo-Boolean optimization." *DAM*, 123, 155-225, 2002.
