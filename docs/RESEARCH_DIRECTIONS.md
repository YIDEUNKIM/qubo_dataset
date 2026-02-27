# 연구 방향 및 논문화 전략

## 1. 현재 상태 진단

### 보유 자산 (6개 방법론)

| 생성기              | GS 보장 | SA-hard | 보조변수 | 구별 불가능 | 주요 탐지 경로 | 난이도 조절 | 상태 |
|------------------|:-------:|:-------:|:--------:|:----------:|:---:|:----------:|:----:|
| Wishart          | 수학적 (유한정밀도 제외) | **O** (alpha<0.95) | 없음 | X | low-rank, 행간 상관 | alpha | 완성 |
| McEliece         | 조건부 (M 의존) | **O** (m=4,t=2) | 있음 (대량) | 미분석 | — | m, t | 완성 |
| Quiet Planting   | 조건부 (field 필요) | 중간 | 있음 | **O** (alpha<3.86) | 없음 | alpha, field | 완성 |
| Posiform         | **수학적 (유일)** | X (SA-trivial) | 없음 | X | sparse(15%), 대각편향, 부호패턴 | X | 완성 |
| Posiform Hardened | **수학적 (유일)** | 조절 가능 | 없음 | X | block-diagonal, 대각편향, kurtosis | alpha, coeff_type | 완성 |
| Zero Expectation | 수학적 | X (SA-trivial) | 없음 | **대각만 가능** | 대각편향 (SNR~√n), 편향=갭 동치 | density | 완성 |

### SA 난이도 분류 기준

**표준 실험 조건**: num_sweeps = 10N, num_reads = 200, 랜덤 target으로 10회 이상 반복.

| 분류 | 조작적 정의 | TTS 스케일링 |
|------|------------|-------------|
| **SA-hard** | p(N=500) < 10% | TTS(99%) ~ exp(cN), c > 0 |
| **SA-moderate** | 10% ≤ p(N=500) < 90% | 지수적이나 느린 감소 |
| **SA-easy** | p(N=1000) ≥ 90% | TTS(99%) ~ poly(N) |

여기서 p(N)은 N-bit QUBO에서 SA가 ground state를 찾을 확률. TTS(99%)는 99% 확률로 GS를 찾는 데 필요한 총 SA 시간으로, TTS = τ × num_sweeps × ⌈log(1-0.99) / log(1-p_s)⌉.

### QUBO 벤치마크 적용성

생성된 QUBO를 벤치마크로 신뢰할 수 있는가 — 주장하는 ground state가 실제로 맞는가:

| 방법론 | GS 신뢰도 | 근거 | 벤치마크 적합성 |
|--------|:---------:|------|:--------------:|
| Posiform | **확정** | Tarjan SCC 유일성 증명 | **O** |
| Hardened Posiform | **확정** | subgraph brute force + posiform 보장 | **O** |
| Zero Expectation | **확정** | LP 최적화 구조적 보장 | **O** |
| Wishart | 주의 필요 | 이론: W^Tt=0. 실제: 부동소수점으로 W^Tt≈0, 큰 N에서 GS 이탈 가능 | **조건부** |
| Quiet Planting | 조건부 | field=0이면 모든 SAT 해가 동일 에너지 (축퇴). field>0 필요 | **조건부** |
| McEliece | 조건부 | 페널티 M 충분해야 GS 보장. M 부족 시 보조변수 해 오염 | **조건부** |

### 난이도 스펙트럼 (SA 실측 데이터 기반)

```
표준 조건 (sweeps=10N, reads=200):

SA-trivial ◄──────────────────────────────────────────────────────────► SA-hard
Posiform    ZeroExp    Hard(α=0.1)    Hard(α=0.01)    Quiet(f=0.5)    Wishart(α=0.7)
100%@N=1K   100%@1K    ~100%@N=500    ~90%@N=500      0%@N=500        0%@N=200
                                                          │
                                            McEliece(m=4,t=2): 0%@k=8 (total_vars=33)
```

> 전체 방법론의 정량적 비교 결과는 **[METHODOLOGY_COMPARISON.md](METHODOLOGY_COMPARISON.md)**를 참조.

### 핵심 딜레마: 3중 트레이드오프

SA-hard, 구별 불가능, 수학적 GS 보장 — 이 세 가지를 동시에 달성하는 방법론은 아직 없다:

| 성질 | 요구 조건 | 충돌 요인 |
|------|----------|----------|
| SA-hard | J 행렬이 low-rank → Q 계수들이 **서로 상관** | 구별 불가능과 충돌 |
| 구별 불가능 | Q 계수가 **독립적이고 무작위** | SA-hard와 충돌 |
| 수학적 GS 보장 | 구조적 제약 (양수 계수, 직교 투영 등) | 완전 무작위와 충돌 |

현재 가장 가까운 조합:

| 방법론 | SA-hard | 구별 불가능 | GS 보장 | 취약점 |
|--------|:-------:|:----------:|:-------:|:---:|
| Wishart | **O** | X (low-rank, 상관) | 조건부 | 고유값 분석으로 즉시 탐지 |
| McEliece | **O** | 미분석 | 조건부 | M 부족 시 GS 이탈 |
| Quiet Planting | 중간 | **O** (α<3.86) | 조건부 (field 필요) | field 존재 자체가 단서? |
| Posiform | X | X | **O** | sparse(15%), 부호패턴(2/3) |
| Posiform Hardened | 조절 가능 | X | **O** | block-diagonal(σ비 250:1) |
| Zero Expectation | X | **대각만 가능** (off-diag 완벽) | **O** | 대각편향(SNR~√n), **편향=갭 동치** → 제거 원리적 불가 |

> **결론**: 6개 방법론 중 3가지를 동시에 만족하는 것은 없음. Zero Expectation이 2가지(구별 불가능 + GS 보장)를 가장 강하게 만족하나 SA-trivial. 심층 분석 결과, ZeroExp의 대각 편향은 에너지 갭과 동일 물리량이라 제거가 원리적으로 불가능함이 확인됨 (노이즈 주입 최적 σ에서도 SNR 5%만 감소, Ising-derived는 비대각 100% 노출로 상보적). Quiet Planting이 SA 중간 + 구별 불가능의 유일한 조합이나 GS가 조건부.

---

## 2. 방법론별 미완 과제

### 2.1 Posiform Hardened — ~~SA 난이도 체계적 검증~~ 구별 불가능성 개선

**현재 상태**: SA 실험 완료, 구별 불가능성 분석 완료.

**완료된 작업**:
- ~~Sweep 전이 실험~~ **완료**: N=500에서 S-curve 관찰, lin2 α=0.01이 가장 어려움
- ~~N 스케일링 실험~~ **완료**: N=1000에서 lin2 α=0.01 → 40% 성공률
- ~~Hardened vs Plain 비교~~ **완료**: Hardened가 Plain보다 확실히 어려움
- ~~구별 불가능성 분석~~ **완료** (2026-02-25): block-diagonal, 대각편향, kurtosis 모두 탐지 가능

**발견된 탐지 경로** (N=100, 30회 반복):
1. **Block-diagonal 구조**: off-block std=0.004 vs on-block std=1.0 (KS p=0)
2. **대각 편향**: E[Q_ii|b=0]=+0.36, E[Q_ii|b=1]=-0.29 (lin2 α=0.01)
3. **Kurtosis**: 4.44 (bimodal) vs -2.0 (uniform random)

**남은 작업**:

| 우선순위 | 항목 | 내용 |
|:--------:|------|------|
| 높음 | Block-diagonal 구조 위장 | Global random QUBO base 사용 시 GS 보장 가능 범위 연구 |
| 중간 | Inter-subgraph noise 실험 | subgraph 간 랜덤 coupling 추가 시 GS/난이도 영향 측정 |
| 중간 | Kurtosis 정규화 | coefficient 분포를 Gaussian으로 변환하여 kurtosis 탐지 회피 가능성 |

### 2.2 McEliece — Eq.14 Exact Decomposition + QPU 실험

**현재 상태**: SA 실험 프레임워크 완성 (`test_mceliece.py`). m=3,4에서 SA 실험 완료.

**SA 실험 결과 요약**:
- m=3, t=2: k=2, total_vars=2, SA 100% (trivial — 변수 2개)
- m=3, t=1: k=4, total_vars=7, SA ~60-80% (aux 3개로 적당히 어려움)
- m=4, t=2: k=8, total_vars=33, SA **~0%** (aux 25개 → SA-hard)
- m=4, t=1: k=11, total_vars=71, SA 0% (aux 60개 → 매우 어려움)
- **m≥5는 Rosenberg 차수축소의 지수적 비용으로 QUBO 생성 자체가 60초+ 소요 → 실험 불가**

**필요 작업**:

| 우선순위 | 항목 | 내용 | 예상 시간 |
|:--------:|------|------|:---------:|
| **필수** | Eq.14 exact decomposition | 페널티 M 의존성 완전 제거. 논문 Section 2의 identity 구현 | 연구 필요 |
| 높음 | p-local PT 실험 재현 | 논문의 핵심: p-local Ising에서 O(2^k) parallel tempering 스케일링 | QPU 필요 |
| ~~높음~~ | ~~SA 스케일링 실험~~ | ~~m 증가에 따른 SA 성공률 변화 측정~~ | ~~완료~~ |
| 중간 | Stern 알고리즘 비교 | 고전 최적 공격 복잡도와의 비교 | 연구 필요 |
| 낮음 | 대규모 실험 (m≥5) | m≥5는 Rosenberg 차수축소 O(2^w) 비용으로 현재 불가. Eq.14 구현 후 재시도 | Eq.14 후 |

**핵심 이슈**: Eq.14는 p-body 항을 보조변수 없이 2-body로 분해하는 exact identity. 현재 Rosenberg 방식은 보조변수 + 페널티 M이 필요하여 QUBO가 크고, M 부족 시 GS가 깨질 수 있음. Eq.14 구현이 되면 McEliece의 가치가 크게 올라감. 또한 m≥5의 확장성 문제도 해결됨.

### 2.3 Quiet Planting — Field 최적화 + 대규모 실험

**현재 상태**: SA 스케일링 완료 (표준 조건, field=0.5, alpha=4.2):
- N=100: 100%, N=200: 100%, N=300: 40%, N=500: 0%, N=750: 0%
- 구별 불가능성: **O** (alpha<3.86에서 수학적으로 보장, 추가 분석 불필요)

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

### 방향 A: 6-방법론 통합 벤치마크 스위트 논문 (가장 현실적)

> **"QUBO Benchmark Suite with Verified Ground States: A Comprehensive Toolkit for Quantum Annealing Evaluation"**

6개 방법론을 난이도/보장/은닉성 축으로 체계화한 통합 벤치마크 스위트.

**핵심 기여**:
1. 난이도 스펙트럼 전체를 커버하는 최초의 통합 벤치마크
2. 각 방법론의 SA 난이도를 동일 조건에서 정량 비교
3. **각 방법론의 구별 불가능성을 체계적으로 정량 분석** (sparsity, 대각편향, 부호패턴, block-diagonal, kurtosis 등)
4. 사용 시나리오별 최적 생성기 권장 가이드

**필요 추가 실험**:

| 우선순위 | 항목 | 내용 | 담당 방법론 | 상태 |
|:--------:|------|------|:----------:|:----:|
| ~~**필수**~~ | ~~Hardened SA 실험~~ | ~~sweep 전이 + N 스케일링 + vs Plain~~ | ~~Posiform Hardened~~ | **완료** |
| ~~**필수**~~ | ~~McEliece SA 실험~~ | ~~m별 SA 성공률 스케일링~~ | ~~McEliece~~ | **완료** |
| ~~**필수**~~ | ~~통합 비교 실험~~ | ~~6개 방법론 동일 조건에서 비교~~ | ~~전체~~ | **완료** |
| ~~**필수**~~ | ~~구별 불가능성 분석~~ | ~~Posiform, Hardened, ZeroExp 정량 분석~~ | ~~3개 방법론~~ | **완료** |
| 높음 | D-Wave QPU 실험 | 최소 Wishart + Hardened + Posiform | Wishart, Hardened | 미수행 |
| 높음 | 데이터셋 공개 | 재현성 + 커뮤니티 기여 | 전체 | 미수행 |
| 중간 | QAOA 실험 | 회로 깊이별 성능 (Qiskit) | 전체 | 미수행 |

**논문 구조 (예상)**:
```
1. Introduction — QUBO 벤치마크의 필요성
2. Background — QUBO, Ising, SA, QA
3. Methods — 6개 방법론 각각의 원리와 GS 보장
4. Experiments
   4.1 SA 난이도 비교 (6-way, 동일 조건)
   4.2 N 스케일링 (easy→hard 순서로 실패 시점)
   4.3 구별 불가능성 분석 (각 방법론의 탐지 경로 정량화)  ← NEW
   4.4 QPU 실험 (Wishart, Hardened 중심)
5. Discussion — 3중 트레이드오프, 사용 시나리오 권장
6. Conclusion + 데이터셋 공개
```

**타겟 저널**: *Scientific Data*, *Quantum Science and Technology*

**현실성**: **높음** — 대부분 실험 프레임워크가 이미 구현됨. SA 실험 실행 + QPU 접근이 주요 추가 작업.

---

### 방향 B: SA-hard + 구별 불가능 + GS 보장 결합의 정량적 한계 분석

> **"Toward Indistinguishable Hard QUBO: Quantifying the Three-Way Trade-off Between Hardness, Detectability, and Ground State Guarantee"**

3중 트레이드오프의 정량적 경계를 실험적으로 탐색.

**현재까지 확인된 트레이드오프 현황** (2026-02-25 업데이트):

| 방법론 | SA-hard | 구별 불가능 | GS 보장 | 주요 탐지 경로 |
|--------|:-------:|:----------:|:-------:|:---:|
| Wishart | **O** | X | 조건부 | low-rank, 행간 상관, Marchenko-Pastur |
| Quiet Planting | 중간 | **O** (α<3.86) | 조건부 | 없음 (수학적 보장) |
| Posiform | X | X | **O** | sparse(15%), 부호패턴(2/3), 대각편향(±4.5) |
| Posiform Hardened | 조절 가능 | X | **O** | block-diagonal(σ비: 250:1), kurtosis(4.4 vs -2.0) |
| Zero Expectation | X | **대각만** | **O** | 대각편향(SNR~√n), **편향=갭 동치** (제거 불가) |
| McEliece | **O** | 미분석 | 조건부 | — |

**핵심 관찰**: 6개 방법론 모두 3가지를 동시에 만족하지 못함. 각 방법론의 약점이 정량적으로 명확해짐:
- Posiform 계열: GS 보장은 강하지만 Q 행렬 구조에서 target이 노출됨
- Wishart/QP: 은닉성은 있지만 GS 보장이 불완전하거나 SA-hard가 조건부
- Zero Expectation: off-diagonal 은닉 + GS 보장이지만 SA-trivial이고 대각이 노출됨

**Zero Expectation 심층 분석 결과** (`analyze_diagonal_deep.py`):
- 대각 편향 ≡ 에너지 갭: E[ΔE_k] = E[Q_kk], 대각=0 시 GS 100% 파괴
- 노이즈 주입: 최적 σ에서 SNR 4.8%만 감소 (편향 크기=갭 크기 → 은닉 원리적 불가)
- Ising-derived: 대각 SNR=0.105 (은닉)이지만 비대각 부호 탐지 100% (상보적)
- 하이브리드 λ=0.5 minimax: max_SNR=0.83 (17% 개선) — 여전히 높은 탐지율
- **결론**: 양의 페널티 하에서 off-diag=0이면 모든 GS 정보가 대각에 집중되어 제거 불가

**실험 설계**:

1. **Wishart + Zero-Expectation 혼합 모델**:
   ```
   Q_mixed = (1 - λ) × Q_Wishart + λ × Q_ZeroExp
   ```
   - λ sweep (0.0 ~ 1.0) × N sweep
   - 3축 동시 측정: SA 성공률 / 통계적 탐지율 (KS test) / GS 보존율

   > **참고**: ZeroOffDiag + Ising-derived 하이브리드(Q = (1-λ)Q_zero + λ·Q_ising)는 이미 분석 완료 → λ=0.5 minimax에서 max_SNR=0.83 (편향 채널 간 재분배만 가능, 총 정보량 보존). Wishart와의 혼합은 상관 구조 도입으로 다른 효과 기대.

2. **Quiet Planting + Wishart 구조 비교**:
   - 둘 다 planted model이지만 은닉성 메커니즘이 다름
   - alpha < 3.86 (Quiet, 구별 불가)과 alpha < 0.95 (Wishart, SA-hard)의 교차점 탐색

3. **Posiform Hardened의 구별 불가능성 개선** (방향 D와 연계):
   - Block-diagonal 구조 위장: global random QUBO base 또는 inter-subgraph coupling
   - 대각 편향 감소: posiform scale α 감소 시 대각 편향도 감소하는지 측정
   - Kurtosis 정규화: coefficient 분포 변환으로 bimodal → unimodal

4. **Posiform 부호패턴 취약점 분석** (신규):
   - Q_ij 부호가 target[i]==target[j]를 2/3 확률로 노출
   - 부호 패턴만으로 target을 복원하는 공격 알고리즘 구현 및 성공률 측정
   - 이 취약점이 Hardened Posiform에서도 존재하는지 확인

**파레토 곡선 도출**:
```
축 1: SA 성공률 (%) — 낮을수록 hard
축 2: 통계적 탐지율 (%) — 낮을수록 은닉
축 3: GS 보존율 (%) — 높을수록 신뢰

→ "구별 불가능성을 X% 포기하면 SA 성공률이 Y%로 떨어진다"
→ "SA-hard를 Z% 포기하면 GS 보장 신뢰도가 W%까지 올라간다"
```

**추가 작업**:

| 우선순위 | 항목 | 내용 | 상태 |
|:--------:|------|------|:----:|
| ~~필수~~ | ~~각 방법론 탐지 경로 정량화~~ | ~~sparsity, 부호패턴, 대각편향, block-diagonal 등~~ | **완료** |
| ~~필수~~ | ~~ZeroExp 대각 편향 심층 분석~~ | ~~편향=갭 동치, 노이즈 주입, Ising 비교, 하이브리드~~ | **완료** (max_SNR=0.83) |
| 필수 | 혼합 모델 구현 | Wishart + ZeroExp 비율 조절 생성기 | 미수행 |
| 필수 | 파레토 곡선 실험 | λ sweep × N sweep | 미수행 |
| 높음 | Posiform 부호패턴 공격 | target[i]==target[j] 복원 알고리즘 | 미수행 |
| 높음 | 이론적 하한 분석 | 정보이론적 관점에서 결합 불가능성 경계 | 미수행 |

**타겟 저널**: *Physical Review E*, *Journal of Statistical Mechanics*

**현실성**: **중간→높음** — 구별 불가능성 분석이 완료되어 3축 중 2축(SA 난이도, 탐지율)의 데이터가 이미 수집됨. 혼합 모델 구현과 파레토 곡선만 추가하면 논문화 가능.

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

Posiform Hardened의 구조적 취약점을 해결하여 구별 불가능성을 확보하는 새로운 구성법.

**정량화된 취약점** (2026-02-25 분석 결과):

| 탐지 경로 | 탐지 강도 | 탐지 방법 | 해결 난이도 |
|-----------|:---------:|----------|:-----------:|
| Block-diagonal 구조 | **치명적** | off-block σ=0.004 vs on-block σ=1.0 → KS p=0 | 중간 (global base로 해결) |
| 대각 편향 | 높음 | E[Q_ii\|b=0]=+0.36 vs E[Q_ii\|b=1]=-0.29 (t-test p≈0) | 어려움 (posiform의 본질적 속성) |
| Kurtosis 이상 | 중간 | 4.44 (bimodal) vs -2.0 (uniform) | 쉬움 (분포 변환) |
| 부호 패턴 | **미확인** | Posiform에서 2/3 확률로 target pair 노출 → Hardened에서도? | 분석 필요 |

**가능한 접근**:

| 방법 | 아이디어 | GS 보장 | 구별 불가능 | 난이도 유지 | 해결하는 취약점 |
|------|---------|:-------:|:----------:|:----------:|:---:|
| Global random QUBO + Posiform | 단일 global random Q에 posiform overlay | **보장** | **가능** | 연구 필요 | block-diagonal |
| Gaussian coefficient | lin2/lin20 대신 N(0,1)에서 sampling | 보장 | 개선 | 유지 | kurtosis |
| Inter-subgraph noise | subgraph 간 랜덤 coupling 추가 | 약화 가능 | 개선 | 유지 | block-diagonal (부분) |
| Wishart base + Posiform overlay | Wishart Q에 posiform overlay | 연구 필요 | low-rank 문제 | **높음** | block-diagonal |
| 대각 보정 | 대각 편향을 상쇄하는 추가 항 | **불가** (편향=갭 동치) | 미미 (SNR 5%↓) | 약화 가능 | 대각 편향 |

**가장 유망**: Global random QUBO + Posiform — 현재 Hardened는 disjoint subgraph 위의 random QUBO를 사용하지만, 전체 변수에 걸친 단일 random QUBO를 base로 사용하고 그 위에 posiform을 overlay하면 block-diagonal 구조가 사라짐. 단, subproblem의 exact GS를 알아야 하므로 전체 random QUBO의 GS를 구해야 하는 문제가 있음 (현재는 subgraph 분할로 brute force 가능).

**추가 유망**: Gaussian coefficient 변환 — kurtosis 탐지를 가장 쉽게 해결. lin2의 {-1,+1}을 N(0,1)로 교체 시 kurtosis가 정상화되고, Q 분포가 random QUBO에 더 가까워짐. GS 보장은 유지됨 (posiform overlay 크기만 충분하면).

**타겟 저널**: *Physical Review Letters*, *npj Unconventional Computing*

**현실성**: **낮음→중간** — block-diagonal은 global base로 해결 가능하나, 대각 편향은 posiform의 본질적 속성이라 근본적 해결이 어려움. 부분적 개선이라도 논문화 가능.

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
| **현실성** | **높음** | **중간→높음** | 중간 | 낮음→중간 | 중간 |
| **Novelty** | 중간→높음 | **높음** | 발견 의존 | **최고** | 높음 |
| **추가 작업량** | **적음** (QPU만) | 중간 (혼합모델) | 중간 | 많음 | 많음 |
| **리스크** | **낮음** | **낮음** | 중간 | 중간 | 중간 |
| **임팩트** | 중간→높음 | 높음 | **높음** (발견 시) | **최고** | 높음 |
| **QPU 필요** | O (권장) | X | **O (필수)** | X | O (권장) |
| **코드 변경** | 실험 실행만 | 혼합 생성기 | 실험 실행만 | 새 구성법 | 기존 코드 수정 |
| **기존 데이터 활용** | **높음** | **높음** | 낮음 | 중간 | 낮음 |

> **업데이트 (2026-02-25)**: 구별 불가능성 정량 분석 완료로 방향 A와 B의 현실성이 모두 상승. 특히 방향 A는 구별 불가능성 분석이 새로운 핵심 기여(Section 4.3)로 추가되어 novelty가 증가. 방향 B는 6개 방법론의 탐지 경로가 모두 정량화되어 파레토 곡선의 기초 데이터가 확보됨.

---

## 5. 권장 로드맵

### Phase 1: 즉시 실행 가능한 실험 — **완료** (2026-02-25)

**목표**: 방향 A 논문의 데이터 수집

1. ~~**Posiform Hardened SA 실험 실행**~~ **완료** (2026-02-20)
   - Sweep 전이 실험: 논문 4개 핵심 발견 모두 재현 확인
   - 결과: [METHODOLOGY_COMPARISON.md](METHODOLOGY_COMPARISON.md) Section 5

2. ~~**McEliece SA 스케일링 실험**~~ **완료** (2026-02-20)
   - `test_mceliece.py`: m-scaling, t-sweep, sweep 전이, 6-way 비교 실험 프레임워크 완성
   - m=3,4에서 SA 실험 완료: m=4,t=2(k=8, aux=25)에서 SA-hard 확인
   - **제한**: m≥5는 Rosenberg 차수축소 비용으로 QUBO 생성 불가 → Eq.14 필요

3. ~~**6-way 통합 비교 실험**~~ **완료** (2026-02-20)
   - 동일 조건 (N=50,100,200,500,1000 / sweeps=1000 / reads=50)
   - 6개 방법론 SA 비교 (McEliece 포함 — n_bits=8에서 비교)
   - 결과: [METHODOLOGY_COMPARISON.md](METHODOLOGY_COMPARISON.md)

4. ~~**Quiet Planting SA 스케일링 실험**~~ **완료** (2026-02-25)
   - 표준 조건 (sweeps=10*total_vars, reads=200): N=100(100%) → N=300(40%) → N=500(0%)
   - METHODOLOGY_COMPARISON 조건 (sweeps=1000, reads=50): N=500(0%), N=1000(0%)

5. ~~**구별 불가능성 정량 분석**~~ **완료** (2026-02-25)
   - Posiform: sparse(15%), 부호패턴(2/3), 대각편향(±4.5) → **구별 가능**
   - Hardened Posiform: block-diagonal(σ비 250:1), kurtosis(4.4), 대각편향 → **구별 가능**
   - Zero Expectation: off-diagonal E[q_ij]=0 (완벽 은닉), 대각편향(SNR~√n) → **대각만 구별**

6. ~~**Zero Expectation 대각 편향 심층 분석**~~ **완료** (2026-02-25)
   - 대각 편향 ≡ 에너지 갭 동치 증명: E[ΔE_k] = E[Q_kk], knockout 시 GS 100% 파괴
   - 노이즈 주입 한계: 최적 σ≈31.5(N=50)에서 SNR 4.8%만 감소 → 은닉 원리적 불가
   - Ising-derived 상보성: 대각 SNR=0.105 vs 비대각 부호탐지=100% (편향 채널 간 상보)
   - 하이브리드 minimax: λ=0.5에서 max_SNR=0.83 (순수 ZeroOffDiag 대비 17% 개선, 여전히 높음)
   - 결론: 양의 페널티 + off-diag=0 → GS 정보가 대각에 집중되어 제거 원리적 불가

### Phase 2: 논문 작성 + QPU 실험 (1~2개월)

**목표**: 방향 A 논문 완성 (구별 불가능성 분석을 핵심 기여로 추가)

1. **논문 작성** — SA 데이터 + 구별 불가능성 분석이 모두 준비됨
   - Section 4.1~4.2: SA 비교/스케일링 (데이터 완료)
   - Section 4.3: 구별 불가능성 분석 (**신규 핵심 기여**, 데이터 완료)
   - Section 5: 3중 트레이드오프 논의 (데이터 완료)
2. D-Wave Leap 클라우드로 QPU 실험 (Wishart + Hardened 중심) → Section 4.4
3. 데이터셋 공개 준비 (GitHub + Zenodo)
4. 논문 투고 → *Scientific Data* 또는 *Quantum Science and Technology*

### Phase 3: 심화 연구 (3~6개월)

**목표**: 방향 B (3중 트레이드오프 파레토 곡선) — 기초 데이터가 가장 많이 축적됨

**우선 순위 업데이트** (2026-02-25):
- **방향 B** (최우선): 6개 방법론의 탐지 경로가 모두 정량화되어 파레토 곡선의 기초 데이터 확보. 혼합 모델만 추가하면 논문화 가능
- **방향 D**: Hardened의 block-diagonal 취약점이 정량화되어 해결 방향이 구체화됨. global base 접근법이 유망
- **방향 E**: McEliece Eq.14 이해에 진전이 있으면
- **방향 C**: QPU 결과가 흥미로우면 (Phase 2에서 QPU 실험 후 판단)

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
