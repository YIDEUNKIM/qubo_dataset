# Wishart Planted Ensemble — 실험 보고서

## 1. 배경 및 동기

### 기존 Hard Mode의 한계
기존 `qubo_hard_mode.py`는 **Backbone(W=20.0) + Frustration(W=0.2)** 구조를 사용하여 QUBO 문제를 생성한다. 그러나 backbone과 noise 가중치의 비율이 100:1로, Simulated Annealing(SA)이 backbone 신호를 쉽게 추적하여 **N=100에서도 100% 성공률**을 보인다.

### 목표
SA가 **실제로 실패하는** 어려운 QUBO 벤치마크를 생성하기 위해, **Wishart Planted Ensemble** (Hamze et al., 2020)을 구현한다.

---

## 2. 알고리즘: Wishart Planted Ensemble

### 2.1 핵심 아이디어

임의의 목표 비트스트링 $t$가 주어졌을 때, $t$가 ground state임이 **수학적으로 보장**되면서도 SA가 풀기 어려운 Ising 모델을 구성한다.

난이도는 단일 파라미터 $\alpha = M/N$으로 제어된다.

### 2.2 구성 절차

**Step 1: Spin 변환**
$$t_i = \begin{cases} +1 & \text{if } b_i = 1 \\ -1 & \text{if } b_i = 0 \end{cases}$$

**Step 2: 직교 Gaussian 행렬 생성**

$M = \lfloor \alpha N \rfloor$개의 $N$차원 Gaussian 벡터 $\{w_\mu\}_{\mu=1}^{M}$를 생성한 뒤, 각 벡터를 $t$에 직교 투영한다:

$$w_\mu \leftarrow w_\mu - \frac{w_\mu \cdot t}{\|t\|^2} t$$

이로써 $W^T t = 0$이 보장된다.

**Step 3: Ising Coupling 행렬 계산**

$$J = -\frac{1}{N} W W^T$$

대각 성분은 제거하고, 상삼각 부분만 $J_{ij}$ ($i < j$)로 추출한다.

**Step 4: Ground State 보장 증명**

Ising 에너지: $H(s) = -\sum_{i<j} J_{ij} s_i s_j$

$t$에서의 에너지 기여:
$$-J_{ij} t_i t_j = \frac{1}{N} (W W^T)_{ij} t_i t_j$$

전체 에너지:
$$H(t) = \frac{1}{N} \sum_{i<j} (WW^T)_{ij} t_i t_j = \frac{1}{2N} t^T W W^T t - \text{diag terms}$$

$W^T t = 0$이므로 $t^T W W^T t = \|W^T t\|^2 = 0$.

따라서 $t$는 이차형식 기여가 0이며, 다른 임의의 spin 배열 $s$에서는:
$$H(s) \propto -\|W^T s\|^2 / N \leq 0 \text{ (대각 제거 후)}$$

$t$가 에너지 최소 상태가 됨.

**Step 5: QUBO 변환**

$$s_i = 2x_i - 1$$을 대입하여 Ising → QUBO:

$$q_{ij} = -4J_{ij}, \quad q_{ii} += 2\sum_{j \neq i} J_{ij}$$

### 2.3 난이도 메커니즘

$\alpha$가 제어하는 것은 **에너지 랜드스케이프의 복잡성**:

- **낮은 $\alpha$ (< 0.8)**: $J$ 행렬이 low-rank → 에너지 랜드스케이프가 평탄하고 수많은 near-degenerate local minima 존재 → SA가 metastable 상태에 갇힘
- **높은 $\alpha$ (≥ 1.0)**: $J$ 행렬이 full-rank에 가까움 → 명확한 global minimum → SA가 쉽게 탐색
- **상전이 구간 ($\alpha \approx 0.8-0.95$)**: 1차 상전이가 발생하며, 난이도가 급격히 변함

---

## 3. 구현

### 3.1 파일 구조

| 파일 | 역할 |
|------|------|
| `qubo_wishart.py` | Wishart Planted Ensemble QUBO 생성기 |
| `test_wishart.py` | SA 실험 프레임워크 (alpha sweep, scaling, 비교) |

### 3.2 주요 함수

```python
# qubo_wishart.py
create_wishart_ising(n, alpha, target, seed)  # Ising 모델 생성
ising_to_qubo(J_dict, h_dict, n)              # Ising → QUBO 변환
create_qubo_wishart(target, alpha, seed)       # 메인 진입점
verify_ground_state(Q, target)                 # 통계적 검증
verify_brute_force(Q, target, n)              # 완전 검증 (n ≤ 20)

# test_wishart.py
run_wishart_experiment(...)                    # Alpha sweep 실험
run_scaling_experiment(...)                    # N 스케일링 실험
run_comparison_experiment(...)                 # Wishart vs Hard Mode 비교
hardness_metrics(...)                          # TTS, 에너지비 등 메트릭
```

### 3.3 CLI 사용법

```bash
# QUBO 생성
python3 qubo_wishart.py 10110 0.7        # 타겟 10110, alpha=0.7
python3 qubo_wishart.py 10110 0.7 42     # seed 지정

# Alpha sweep 실험
python3 test_wishart.py 100 10           # N=100, 10 runs per alpha

# Scaling 실험
python3 test_wishart.py --scaling 0.7    # alpha=0.7 고정, N=[50,100,...,500]

# Wishart vs Hard Mode 비교
python3 test_wishart.py --compare

# Hardness metrics
python3 test_wishart.py --metrics
```

---

## 4. 실험 결과

### 4.1 Ground State 검증 (Brute Force)

| 조건 | 결과 |
|------|------|
| N=10, 30개 타겟 × 4개 alpha | **120/120 (100%) 성공** |

모든 경우에서 target이 정확히 ground state임을 확인.

### 4.2 Alpha Sweep — 거친 탐색 (N=100, num_reads=200, num_sweeps=1000)

| Alpha | M | SA 성공률 | 에너지 정확도 | 평균 해밍거리 |
|-------|---|----------|-------------|-------------|
| 0.3   | 30  | 0%   | 97.96% | 50.2 |
| 0.5   | 50  | 0%   | 95.35% | 52.6 |
| 0.7   | 70  | 0%   | 91.24% | 48.6 |
| 0.8   | 80  | 20%  | 93.60% | 56.2 |
| 1.0   | 100 | 100% | 100.0% | 40.0 |
| 1.5   | 150 | 100% | 100.0% | 60.0 |

**관찰**: $\alpha \approx 0.8 \sim 1.0$에서 상전이 발생. $\alpha = 0.7$에서 에너지 정확도가 최저(91.24%).

### 4.3 Alpha Sweep — 세밀한 탐색 (N=100, alpha=0.50~1.00, step=0.05)

| Alpha | M  | SA 성공률 | 에너지 정확도 | 평균 해밍거리 | 구간 |
|-------|----|----------|-------------|-------------|------|
| 0.50  | 50 | **0%**   | 94.58%      | 48.9        | Hard |
| 0.55  | 55 | **0%**   | 94.85%      | 48.9        | Hard |
| 0.60  | 60 | **0%**   | 93.72%      | 49.5        | Hard |
| 0.65  | 65 | 10%      | 93.90%      | 43.0        | Hard |
| 0.70  | 70 | **0%**   | 92.74%      | 48.3        | **최고 난이도** |
| 0.75  | 75 | 10%      | 92.20%      | 51.8        | Hard |
| 0.80  | 80 | 30%      | 93.17%      | 42.8        | 전이 시작 |
| 0.85  | 85 | 50%      | 95.30%      | 27.1        | 전이 구간 |
| 0.90  | 90 | 70%      | 96.30%      | 44.8        | 전이 구간 |
| 0.95  | 95 | **100%** | 100.0%      | 50.0        | Easy |
| 1.00  |100 | **100%** | 100.0%      | 50.0        | Easy |

**핵심 관찰**:

1. **Hard 구간 (alpha ≤ 0.75)**: SA 성공률 0~10%, 에너지 정확도 92~95% 수준에서 정체. alpha=0.70~0.75에서 에너지 정확도가 가장 낮음 (92.2~92.7%).

2. **상전이 구간 (alpha = 0.80~0.90)**: 성공률이 30% → 50% → 70%로 급격히 증가. 실패 시 에너지 gap도 점차 감소.

3. **Easy 구간 (alpha ≥ 0.95)**: SA가 100% ground state 도달. 임계점 $\alpha_c \approx 0.95$.

4. **에너지 정확도 최저점**: alpha=0.75에서 92.20%로 최저. SA가 ground state 대비 **~8% 에너지 gap**에 갇힘.

### 4.4 Scaling 실험 (alpha=0.7)

| N | SA 성공률 | 에너지 정확도 | 평균 해밍거리 | 해밍/N |
|---|----------|-------------|-------------|--------|
| 50  | **100%** | **100.0%** | 20.0 | 0.40 |
| 100 | 0%   | 92.02% | 51.2 | 0.51 |
| 150 | 0%   | 92.33% | 74.5 | 0.50 |
| 200 | 0%   | 91.49% | 100.0 | 0.50 |
| 300 | 0%   | 91.62% | 149.5 | 0.50 |
| 500 | 0%   | 91.29% | 238.2 | 0.48 |

**관찰**:
1. **N=50→100 사이에서 상전이**: N=50에서는 100% 성공, N=100부터 완전 실패
2. **에너지 정확도 ~91%로 포화**: N이 커져도 SA는 ground state 에너지의 약 91%에서 정체
3. **해밍거리 ≈ N/2**: SA가 찾은 해는 정답과 **랜덤 추측 수준**으로 다름

### 4.5 Wishart vs Hard Mode 비교 (N=100)

| 방식 | 파라미터 | SA 성공률 |
|------|---------|----------|
| **Wishart** | alpha=0.7 | **10%** |
| Hard Mode | noise=0.1 | **100%** |

Wishart가 Hard Mode 대비 **10배 이상 어려운** QUBO를 생성함.

---

## 5. 분석 및 결론

### 5.1 Wishart Ensemble의 난이도 특성

1. **1차 상전이**: 세밀한 sweep 결과, 임계점 $\alpha_c \approx 0.95$ (N=100 기준). alpha=0.90에서 70%, alpha=0.95에서 100%로 급격히 전환. 성공률이 0% → 100%로 바뀌는 구간이 alpha 0.15 폭(0.80~0.95) 안에 집중됨.

2. **에너지 정확도 최저점**: alpha=0.70~0.75 구간에서 SA의 에너지 정확도가 92.2~92.7%로 최저. Ground state 대비 **~8% 에너지 gap**에서 metastable 상태에 갇힘.

3. **Metastable Trapping의 N-독립성**: Scaling 실험에서 N=100~500까지 에너지 정확도가 91~92%로 일정. 이 gap은 문제 크기에 무관한 **구조적 한계**임.

4. **해밍거리 N/2**: 실패 시 SA가 찾은 해는 정답과 해밍거리 ≈ N/2 (랜덤 수준). 에너지만 근사할 뿐, 해의 구조는 전혀 다름.

### 5.2 기존 방식 대비 개선

| 비교 항목 | Hard Mode | Wishart |
|----------|-----------|---------|
| 난이도 제어 | noise_ratio (간접) | alpha (직접, 이론적) |
| SA 실패 유도 | 불가능 (100% 성공) | 가능 (alpha<1.0에서 실패) |
| Ground state 보장 | 조건부 (W_STRONG >> W_WEAK) | 수학적 보장 (W^T t = 0) |
| 상전이 프로파일 | 없음 | 명확한 1차 상전이 |

### 5.3 벤치마크 활용 권장

- **쉬운 문제**: alpha=1.5 (SA 검증용)
- **중간 난이도**: alpha=0.85~0.95 (상전이 구간, 부분 성공)
- **어려운 문제**: alpha=0.5~0.7 (SA 한계 테스트)
- **스케일링 테스트**: alpha=0.7 고정, N=50~500 (상전이 크기 의존성 확인)

---

## 6. 의존성

- `numpy`: 행렬 연산 (Wishart 구성)
- `neal`: D-Wave Simulated Annealing Sampler
- `dimod`: D-Wave 에코시스템 (neal 의존)
