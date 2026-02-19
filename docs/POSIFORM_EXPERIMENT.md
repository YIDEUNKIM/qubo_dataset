# Posiform Planting (Planted 2-SAT → QUBO) — 실험 보고서

## 1. 배경 및 동기

### 기존 방법론과의 위치

| 생성기 | SA-Hard? | 보조변수 | Ground state 유일성 | 구별 불가능성 |
|--------|----------|---------|-------------------|-------------|
| Zero-Expectation | SA 100% 성공 | 없음 | 확률적 | E[q_ij]=0 보장 |
| Hard Mode | SA 100% 성공 | 없음 | 보장 | backbone 구조 노출 |
| Wishart | SA 실패 (alpha~0.7) | 없음 | Z2 대칭 | low-rank 구조 노출 |
| Quiet Planting | field 의존적 | **m개** (QUBO=5.2n) | 축퇴 (field 필요) | alpha < 3.86 보장 |
| **Posiform** | **SA 100% 성공** | **없음** (QUBO=n) | **수학적 보장** | 미분석 |

### 목표

**Posiform Planting** (Hahn, Pelofske, Djidjev, 2023)은 planted 2-SAT → posiform → QUBO 파이프라인으로, **보조변수 없이** QUBO 크기 = n을 유지하면서 target이 **수학적으로 유일한 ground state**임을 보장하는 방법론이다. 이 실험은 Posiform QUBO의 SA 난이도를 대규모(N=1000)까지 체계적으로 측정한다.

---

## 2. 알고리즘 요약

### 2.1 Planted 2-SAT 생성

1. Target x*로부터 2-SAT clause를 반복 생성 (랜덤 변수쌍 → wrong tuple 배제)
2. **Tarjan SCC** 기반 2-SAT solver로 satisfiability + uniqueness 검사
3. Target이 유일한 해가 될 때까지 2단계로 clause 추가:
   - **Phase 1**: 랜덤 clause 추가 (대부분의 변수 고정)
   - **Phase 2**: 아직 flippable한 변수에 대해 집중 clause 추가

### 2.2 Posiform → QUBO 변환

각 clause에 랜덤 양수 계수 b ~ Uniform(coeff_range) 부여 후 posiform 항을 QUBO로 변환:

| Wrong tuple | Posiform 항 | QUBO 변환 |
|:-----------:|:-----------:|:----------|
| (0, 0) | b(1-x_i)(1-x_j) | Q_ii -= b, Q_jj -= b, Q_ij += b, const += b |
| (0, 1) | b(1-x_i)x_j | Q_jj += b, Q_ij -= b |
| (1, 0) | bx_i(1-x_j) | Q_ii += b, Q_ij -= b |
| (1, 1) | bx_ix_j | Q_ij += b |

**핵심**: P(x*) = 0 (모든 posiform 항이 target에서 0), 다른 모든 x에서 P(x) > 0 → **유일한 ground state 보장**.

### 2.3 Quiet Planting과의 구조적 차이

| 특성 | Quiet Planting (3-SAT) | Posiform (2-SAT) |
|------|:---------------------:|:----------------:|
| SAT 유형 | 3-SAT (NP-complete) | **2-SAT (P)** |
| 보조변수 | m개 (QUBO 크기 = n+m) | **없음** (QUBO 크기 = n) |
| Ground state 유일성 | 축퇴 (field 필요) | **수학적 보장** |
| Planted field 필요 | 예 (축퇴 해소) | **불필요** |
| 상전이 구조 | alpha로 연속 제어 | 없음 (clause 수는 유일성에 의해 결정) |

---

## 3. 실험 설정

### 3.1 실험 파라미터

| 파라미터 | 값 |
|---------|-----|
| 문제 크기 (N) | 10, 20, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000 |
| 반복 횟수 (num_runs) | 10 (각 N에서 독립 QUBO 생성) |
| SA 시도 횟수 (num_reads) | **1** (최소 설정, 난이도 상한 측정) |
| SA sweeps | max(1000, 10 * N) |
| 계수 범위 (coeff_range) | (1.0, 3.0) |
| Solver | `neal.SimulatedAnnealingSampler` |
| 시드 | random.seed(42), np.random.seed(42) |

### 3.2 측정 지표

- **성공률**: SA가 target과 정확히 일치하는 해를 찾은 비율 (EXACT)
- **에너지비 (E_ratio)**: best_found_energy / target_energy (1.0이면 최적)
- **해밍거리**: target과 SA 해 사이의 비트 차이
- **생성 시간**: 2-SAT 생성 + posiform→QUBO 변환 시간

---

## 4. 실험 결과

### 4.1 N 스케일링 (num_reads=1)

| N | QUBO 크기 | Sweeps | 성공률 | E비 | 해밍거리 | 평균 생성시간 |
|---:|:---------:|:------:|:------:|:---:|:-------:|:-----------:|
| 10 | 10 | 1,000 | **10/10 (100%)** | 1.0000 | 0.0 | 0.00s |
| 20 | 20 | 1,000 | **10/10 (100%)** | 1.0000 | 0.0 | 0.01s |
| 50 | 50 | 1,000 | **10/10 (100%)** | 1.0000 | 0.0 | 0.11s |
| 100 | 100 | 1,000 | **10/10 (100%)** | 1.0000 | 0.0 | 0.57s |
| 150 | 150 | 1,500 | **10/10 (100%)** | 1.0000 | 0.0 | 1.54s |
| 200 | 200 | 2,000 | **10/10 (100%)** | 1.0000 | 0.0 | 2.35s |
| 300 | 300 | 3,000 | **10/10 (100%)** | 1.0000 | 0.0 | 6.60s |
| 400 | 400 | 4,000 | **10/10 (100%)** | 1.0000 | 0.0 | 13.46s |
| 500 | 500 | 5,000 | **10/10 (100%)** | 1.0000 | 0.0 | 22.27s |
| 600 | 600 | 6,000 | **10/10 (100%)** | 1.0000 | 0.0 | 29.89s |
| 700 | 700 | 7,000 | **10/10 (100%)** | 1.0000 | 0.0 | 45.13s |
| 800 | 800 | 8,000 | **10/10 (100%)** | 1.0000 | 0.0 | 63.60s |
| 900 | 900 | 9,000 | **10/10 (100%)** | 1.0000 | 0.0 | 91.49s |
| 1000 | 1000 | 10,000 | **10/10 (100%)** | 1.0000 | 0.0 | 108.24s |

**총 140/140 (100%) EXACT 성공. 단 1회 SA 시도로 N=1000까지 완벽 해결.**

### 4.2 생성 시간 스케일링

```
N=100:   ~0.6s
N=200:   ~2.4s    (4x)
N=300:   ~6.6s    (2.8x)
N=400:   ~13.5s   (2.0x)
N=500:   ~22.3s   (1.7x)
N=700:   ~45.1s
N=1000:  ~108.2s
```

생성 시간은 대략 **O(N^2) ~ O(N^2.5)** 수준. 유일성 검사 (`check_2sat_uniqueness`)가 O(N) * O(2-SAT solve) = O(N^2)의 비용을 가지므로 이론과 일치.

---

## 5. 분석 및 해석

### 5.1 SA에 대한 난이도가 극히 낮은 이유

Posiform QUBO가 SA에 대해 사실상 trivial한 이유를 구조적으로 분석하면:

#### (1) Smooth 에너지 지형 (Metastable 상태 없음)

2-SAT 기반 posiform은 **모든 계수가 양수**이므로:
- P(x*) = 0 (최소값)
- 각 clause 위반 시 정확히 하나의 posiform 항만큼 에너지 증가
- **에너지 장벽이 낮고, local minima 구조가 단순함**

Wishart의 경우 alpha_c 근처에서 ground state와 에너지가 ~92% 비슷하면서 해밍거리 ~N/2인 metastable 상태가 존재하여 SA가 빠져나올 수 없지만, Posiform에는 이러한 구조가 없다.

#### (2) 보조변수 부재로 인한 작은 탐색 공간

| 방법론 | QUBO 변수 수 | N=100 기준 |
|--------|:----------:|:---------:|
| Posiform | n | 100 |
| Quiet Planting (alpha=4.2) | n(1+alpha) | **520** |
| Wishart | n | 100 |

Quiet Planting은 동일 N에서 5.2배 큰 탐색 공간을 가짐. 반면 Posiform은 원래 문제 크기 그대로.

#### (3) 2-SAT의 다항 시간 풀이 가능성

2-SAT 자체가 P(다항 시간)에 풀리는 문제이므로, 이로부터 유도된 QUBO도 구조적으로 어려울 수 없다. SA의 single-bit-flip dynamics가 2-SAT의 함축 전파(unit propagation)와 유사하게 작동하여 효율적으로 ground state에 도달.

### 5.2 다른 방법론과의 정량적 비교

| 방법론 | N=100 성공률 | N=500 성공률 | N=1000 성공률 |
|--------|:----------:|:----------:|:-----------:|
| **Posiform** (num_reads=1) | **100%** | **100%** | **100%** |
| Quiet Planting (field=0.5) | 100% | 20% | 0% |
| Wishart (alpha=0.7) | ~10% | ~0% | ~0% |
| Zero-Expectation | ~100% | ~100% | ~100% |
| Hard Mode | ~100% | ~100% | ~100% |

Posiform은 Zero-Expectation, Hard Mode와 함께 "SA-easy" 카테고리에 속한다. 다만 **유일한 ground state를 수학적으로 보장**한다는 점에서 검증 용도로 차별화됨.

### 5.3 벤치마크로서의 가치

**장점**:
- 유일한 ground state가 수학적으로 보장 → **정답 검증이 명확**
- 보조변수 없음 → QUBO 크기가 원래 문제와 동일 → **공정한 비교 가능**
- SA에 대해 쉬움 → **솔버 정확성 검증(sanity check)**에 적합
- 생성이 효율적 (N=1000도 ~2분)

**한계**:
- SA에 대해 너무 쉬워서 **난이도 벤치마크로는 부적합**
- 연속적 난이도 파라미터(alpha 등)가 없어 **난이도 조절 불가**
- 2-SAT 기반이므로 **상전이 구조가 없음** (easy-hard-easy 프로파일 없음)

### 5.4 권장 사용 시나리오

| 시나리오 | Posiform | Wishart | Quiet Planting |
|---------|:--------:|:-------:|:--------------:|
| 솔버 정확성 검증 | **최적** | - | - |
| SA 한계 측정 | - | **최적** | 적합 |
| 양자 우위 벤치마크 | - | **최적** | 적합 |
| 통계적 은닉성 연구 | - | - | **최적** |
| 대규모 문제 (N>1000) | **적합** (QUBO=n) | 적합 (QUBO=n) | 비용 높음 (QUBO=5.2n) |

---

## 6. 결론

1. **Posiform QUBO는 SA에 대해 trivially easy**: num_reads=1 (단 1회 SA 시도)로 N=1000까지 140/140 (100%) 성공. 에너지비 = 1.0000, 해밍거리 = 0을 모든 실험에서 달성.

2. **이론적 근거**: 2-SAT의 다항 시간 풀이 가능성, posiform의 smooth 에너지 지형, 보조변수 부재로 인한 작은 탐색 공간이 복합적으로 작용.

3. **벤치마크 역할 분담**: Posiform은 솔버 정확성 검증용, Wishart/Quiet Planting은 난이도 벤치마크용으로 상호보완적 사용이 권장됨.

4. **병목은 생성 시간**: SA 풀이는 trivial하지만, 유일성 검사의 O(N^2) 비용으로 인해 N=1000에서 ~108초 소요. N > 2000 이상에서는 생성 자체가 실용적 한계에 도달할 수 있음.

---

## 7. 실험 환경

- **OS**: Linux 5.15.167.4 (WSL2)
- **Python**: 3.x
- **SA Solver**: `neal.SimulatedAnnealingSampler` (D-Wave Ocean SDK)
- **실험 일자**: 2026-02-19

---

## 8. 참고문헌

1. **Hahn, G., Pelofske, E., & Djidjev, H.** "Using 2-SAT to generate QUBO instances with known optimal solutions." *2023 IEEE International Conference on Quantum Computing and Engineering (QCE)*, 2023.
2. **Tarjan, R. E.** "Depth-first search and linear graph algorithms." *SIAM Journal on Computing*, 1(2), 146-160, 1972.
3. **Aspvall, B., Plass, M. F., & Tarjan, R. E.** "A linear-time algorithm for testing the truth of certain quantified Boolean formulas." *Information Processing Letters*, 8(3), 121-123, 1979.
4. **Boros, E. & Hammer, P. L.** "Pseudo-Boolean optimization." *Discrete Applied Mathematics*, 123(1-3), 155-225, 2002.
5. **Hamze, F., et al.** "Wishart planted ensemble: A tunably rugged pairwise Ising model with a first-order phase transition." *Physical Review E*, 101, 052102, 2020.
6. **Krzakala, F. & Zdeborova, L.** "Hiding quiet solutions in random constraint satisfaction problems." *Physical Review Letters*, 102(23), 238701, 2009.
