# Wishart Planted Ensemble QUBO 생성기

## 개요

**Hamze et al. (2020)**에 기반한 SA-hard planted QUBO 생성 방법론. 직교 Gaussian 투영으로 목표 비트스트링이 ground state임을 수학적으로 보장하면서, alpha = M/N 단일 파라미터로 SA 난이도를 정밀하게 제어한다.

## 이론적 배경

### Wishart 행렬과 Planted Ensemble

**Wishart 행렬**은 W^T W 형태의 행렬로, 통계역학에서 무작위 공분산 구조를 모델링하는 데 사용된다. Planted Ensemble은 이 구조를 활용하여 **알려진 ground state를 가지면서도 SA에게 어려운** Ising 모델을 구성한다.

### 핵심 아이디어

1. 목표 비트스트링 t를 spin 벡터 (+1/-1)로 변환
2. t에 **직교**하는 M개의 N차원 Gaussian 벡터 W를 생성 (M = alpha * N)
3. Ising coupling: J = -(1/N) W W^T (대각 제거)
4. **W^T t = 0**이므로 t가 ground state임이 **수학적으로 보장**됨

### Ground State 보장 증명

Ising 에너지: H(s) = -(1/2N) sum_{i!=j} J_ij s_i s_j = (1/2N) ||Ws||^2 + const

W^T t = 0 이므로 H(t) = const (최소값). 다른 어떤 상태 s != t에 대해서는 Ws != 0이므로 H(s) > H(t).

### 유한 정밀도에서의 Ground State 보장 한계

위 증명은 **무한 정밀도(infinite precision)** 산술을 가정한다. 실제 구현에서는 부동소수점 연산의 유한 정밀도로 인해 ground state 보장이 약해질 수 있다:

1. **직교 투영 오차**: W^T t = 0을 만들기 위한 Gram-Schmidt 투영에서 수치 오차 발생. IEEE 754 double precision에서 W^T t ≈ O(epsilon_mach * sqrt(N)) 수준의 잔차가 남음.

2. **에너지 갭 축소**: 유한 정밀도 오차는 planted state와 다른 상태 사이의 에너지 갭을 줄이거나, 극단적인 경우 뒤집을 수 있다. 특히 **alpha가 작을 때** (M << N) 에너지 갭 자체가 작아서 수치 오차에 더 취약.

3. **실용적 영향**: Hamze et al. (2020)은 논문에서 이 문제를 인지하고, 실제 실험에서는 planted state가 ground state와 매우 높은 확률로 일치하지만 **수학적 보장은 유한 정밀도에서는 성립하지 않음**을 명시. 특히 양자 어닐러의 제한된 정밀도(D-Wave의 ~5 bit precision)에서는 이 문제가 더 심각할 수 있다.

4. **Posiform과의 차이**: Posiform 방법론은 P(x*) = 0이 정수 계수에서도 정확히 성립하므로, 유한 정밀도 문제가 없다. 반면 Wishart는 연속값 행렬 연산에 의존하므로 본질적으로 수치 오차에 노출.

### 상전이 (Phase Transition)

alpha = M/N이 증가하면 에너지 지형이 변한다:

| alpha 범위 | 에너지 지형 | SA 성공률 |
|:----------:|:----------:|:---------:|
| < 0.5 | 매끈함 (단일 바닥) | ~100% |
| 0.5 ~ 0.8 | **Metastable 상태 출현** | 급격히 감소 |
| ~ 0.95 | **1차 상전이** (alpha_c) | ~0% |
| > 1.5 | 정보가 충분해져 다시 쉬워짐 | 회복 |

이 **easy-hard-easy** 프로파일이 Wishart의 핵심 특성이다.

### SA가 실패하는 이유: Metastable 상태

alpha_c 근처에서 ground state와 에너지가 비슷한 (~92%) metastable 상태들이 존재:
- 에너지는 거의 같지만, 해밍거리는 ~N/2 (완전히 다른 위치)
- SA가 탈출하려면 수십 개의 비트를 **동시에** 뒤집어야 함
- SA는 한 번에 1개만 뒤집는 local search → 에너지 장벽을 넘을 수 없음

반면, **양자 어닐링(Quantum Annealing)**은 양자 터널링으로 장벽을 관통할 수 있어 SA 대비 우위가 기대됨.

## 구현 방식

### 알고리즘

```
Input: target (비트스트링), alpha (난이도), seed (재현성)

1. Spin 벡터 생성
   t[i] = +1 if target[i]=='1' else -1

2. Gaussian 벡터 생성
   W = N x M 행렬, 각 원소 ~ N(0,1)

3. 직교 투영 (각 열 mu에 대해)
   W[:, mu] -= (W[:, mu] . t / ||t||^2) * t
   → W^T t = 0 보장

4. Ising coupling 계산
   J = -(1/N) * W @ W^T

5. 대각 제거 + 상삼각 추출
   J_dict[(i,j)] = J_matrix[i,j]  (i < j)

6. Ising → QUBO 변환
   Q_ij = -4 * J_ij
   Q_ii += 2 * J_ij, Q_jj += 2 * J_ij
```

### 핵심 파라미터

| 파라미터 | 범위 | 설명 |
|---------|------|------|
| `alpha` | 0.3 ~ 1.5 | M/N 비율. **0.5~0.8이 가장 어려움** |
| `seed` | int | 재현성을 위한 난수 시드 |

### 주요 함수

| 함수 | 설명 |
|------|------|
| `create_wishart_ising(n, alpha, target, seed)` | Wishart Ising 모델 생성 (J, h) |
| `ising_to_qubo(J_dict, h_dict, n)` | Ising → QUBO 변환 |
| `create_qubo_wishart(target, alpha, seed)` | **메인 진입점** |
| `verify_ground_state(Q, target, num_random_samples)` | 통계적 검증 (single-flip + 랜덤 샘플) |
| `verify_brute_force(Q, target, n)` | N <= 20 완전 탐색 검증 |

## SA 난이도 특성

### 실험 결과 (N=100)

| alpha | SA 성공률 | 평균 에너지 비 | 평균 해밍거리 |
|:-----:|:---------:|:------------:|:------------:|
| 0.3 | 100% | 1.0000 | 0.0 |
| 0.5 | ~80% | 0.9997 | ~5 |
| 0.7 | ~10% | 0.9630 | ~43 |
| 0.8 | ~0% | 0.9319 | ~47 |
| 1.0 | ~10% | 0.9876 | ~20 |
| 1.5 | ~70% | 0.9993 | ~3 |

### N 스케일링 (alpha=0.7)

```
N=50:  ~70% 성공
N=100: ~10% 성공
N=150: ~0%  성공
N=200: ~0%  성공 (에너지 비 ~0.93)
```

성공률이 N에 따라 **지수적으로 감소**하며, 이는 1차 상전이의 특성.

## 구별 가능성 분석

Wishart QUBO는 Q 행렬 분석으로 **구별 가능**:

| 탐지 방법 | 설명 |
|----------|------|
| **Low-rank 구조** | J = -(1/N)WW^T는 rank-M. Q의 유효 rank ≈ M << N |
| **대각 편향** | Q_ii ∝ sum_j J_ij → target 비트 방향으로 편향 |
| **행간 상관** | W가 공유되므로 Q의 서로 다른 행이 상관됨 |

> 구별 가능하지만, 구별할 수 있다는 것과 풀 수 있다는 것은 별개 문제.

## 파일 구성

| 파일 | 역할 |
|------|------|
| `qubo_wishart.py` | 생성기 |
| `test_wishart.py` | SA 실험 (alpha sweep, N scaling, Wishart vs Hard Mode, hardness metrics) |

## 참고 문헌

### 핵심 논문

1. **Hamze, F., Raymond, J., Pattison, C. A., Biswas, K., & Katzgraber, H. G.** "Wishart planted ensemble: A tunably rugged pairwise Ising model with a first-order phase transition." *Physical Review E*, 101, 052102, 2020. [arXiv:2009.11217](https://arxiv.org/abs/2009.11217)
   - 직교 Gaussian 투영 기반 planted Ising 구성
   - alpha 파라미터에 의한 1차 상전이
   - alpha_c 근처에서 SA 실패 (metastable trapping)

### 관련 논문

2. **Barthel, W., & Hartmann, A. K.** "Clustering analysis of the ground-state structure of the vertex-cover problem." *Physical Review E*, 70(6), 066120, 2004. — Planted ensemble의 일반적 프레임워크.

3. **Krzakala, F., & Zdeborova, L.** "Hiding quiet solutions in random constraint satisfaction problems." *Physical Review Letters*, 102(23), 238701, 2009. — Quiet planting 개념의 원조.

### 양자 어닐링 벤치마크

4. **Hen, I., et al.** "Probing for quantum speedup in spin-glass problems with planted solutions." *Physical Review A*, 92(4), 042325, 2015.

5. **King, A. D., et al.** "Quantum critical dynamics in a 5,000-qubit programmable spin glass." *Nature*, 617, 61-66, 2023.

## 사용법

```bash
# 기본 생성 (target="10110", alpha=0.7)
python3 wishart/qubo_wishart.py 10110 0.7

# 시드 지정
python3 wishart/qubo_wishart.py 10110 0.7 42

# 길이로 랜덤 목표 생성
python3 wishart/qubo_wishart.py 100 0.8

# Alpha sweep 실험
python3 wishart/test_wishart.py 100 10

# N 스케일링 실험
python3 wishart/test_wishart.py --scaling 0.7

# Wishart vs Hard Mode 비교
python3 wishart/test_wishart.py --compare

# Hardness metrics (TTS, spectral gap)
python3 wishart/test_wishart.py --metrics
```
