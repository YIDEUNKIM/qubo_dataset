# Hard Mode QUBO 생성기

## 개요

**Backbone(강한 신호) + Frustration(약한 노이즈)** 구조로 SA 솔버에게 "적당히 어려운" QUBO를 생성하는 방법론. Ising 모델의 Star 그래프 백본이 ground state를 보장하고, 랜덤 약한 간선이 frustration을 추가한다.

## 이론적 배경

### Spin Glass와 Frustration

스핀 글래스(spin glass)에서 **frustration**이란, 인접한 스핀들이 동시에 모든 결합을 만족시킬 수 없는 상태를 말한다. Frustration이 많을수록 에너지 지형이 복잡해지고, 솔버가 ground state를 찾기 어려워진다.

### 이 방법론의 아이디어

1. **Backbone (Star Graph)**: 모든 노드가 중심 노드(node 0)에 강한 결합(W_STRONG=20.0)으로 연결. 이 결합만으로 ground state가 유일하게 결정됨.
2. **Noise Edges**: 나머지 간선에 약한 결합(W_WEAK=0.2)을 추가. `noise_ratio` 확률로 frustration(목표 반대 방향) 적용.

### Ground State 보장 증명

Star 그래프에서 모든 노드 i는 중심 0과 연결: J(0,i) = W_STRONG * s_0 * s_i.

노드 i를 뒤집으면 에너지 변화:
- 백본 기여: +2 * W_STRONG (항상 양)
- 노이즈 기여: 최대 N * density * W_WEAK (expectation)

안전 조건: **W_STRONG > N * density * W_WEAK**

N=300, density=1.0 기준: 300 * 0.1 * 0.2 = 6.0 << 20.0 (W_STRONG). 따라서 어떤 single-bit flip도 에너지를 증가시키며, cluster flip에 대해서도 100:1 비율로 안전.

## 구현 방식

### Ising → QUBO 변환

1. **Ising 결합**: J_ij = W * sign * s_i * s_j
   - W = W_STRONG (백본) 또는 W_WEAK (노이즈)
   - sign = +1 (정렬) 또는 -1 (frustration)

2. **QUBO 변환**: s_i = 2x_i - 1 치환
   ```
   H = -J_ij * s_i * s_j
     = -J_ij * (4x_ix_j - 2x_i - 2x_j + 1)
   Q_ij += -4 * J_ij
   Q_ii += 2 * J_ij
   Q_jj += 2 * J_ij
   ```

### 핵심 파라미터

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `W_STRONG` | 20.0 | 백본(Star) 결합 강도 |
| `W_WEAK` | 0.2 | 노이즈 결합 강도 |
| `noise_ratio` | 0.1 | Frustration 확률 (0.0~1.0) |
| `density` | 1.0 | 비-백본 간선 생성 확률 |

### 알고리즘

```
1. Star 그래프 구축 (중심 node 0 ↔ 모든 node i)
   - J(0,i) = W_STRONG * s_0 * s_i  (target 방향 정렬)
   - QUBO로 변환

2. 나머지 간선 (i,j) 순회 (i > 0, j > i):
   - density 확률로 간선 추가
   - noise_ratio 확률로 frustration (sign = -1)
   - J(i,j) = W_WEAK * sign * s_i * s_j
   - QUBO로 변환
```

### 주요 함수

| 함수 | 설명 |
|------|------|
| `create_qubo_hard(target, density, base_range, noise_ratio)` | 메인 생성 함수 |

## SA 난이도 특성

**SA에 대해 쉬운 문제** (SA 성공률 ~100%).

원인: Star 그래프 백본의 강한 결합(W_STRONG=20.0)이 약한 노이즈(W_WEAK=0.2)를 압도. SA가 에너지 지형의 전역적 방향을 쉽게 파악하여 ground state에 수렴.

### Wishart와의 비교

| 특성 | Hard Mode | Wishart |
|------|-----------|---------|
| SA 성공률 (N=100) | ~100% | ~0% (alpha=0.7) |
| 에너지 장벽 구조 | 얕고 넓음 | 깊고 좁음 (metastable) |
| 구별 가능성 | O (Star 구조 탐지 가능) | O (low-rank 탐지 가능) |
| 난이도 조절 | noise_ratio | alpha |

## 파일 구성

| 파일 | 역할 |
|------|------|
| `qubo_hard_mode.py` | 생성기 |

## 참고 문헌

### 이론적 배경

1. **Spin Glass 이론**: Edwards, S. F. & Anderson, P. W. "Theory of spin glasses." *Journal of Physics F: Metal Physics*, 5(5), 965, 1975.
2. **Frustration**: Toulouse, G. "Theory of the frustration effect in spin glasses." *Communications on Physics*, 2, 115-119, 1977.

### 관련 방법론

3. **Wishart Planted Ensemble**: Hamze, F., et al. "Wishart planted ensemble: A tunably rugged pairwise Ising model with a first-order phase transition." *Physical Review E*, 2020. [arXiv:2009.11217](https://arxiv.org/abs/2009.11217) — 보다 정교한 난이도 조절을 제공하는 대안.

### Ising-QUBO 변환

4. 표준 치환: s_i = 2x_i - 1.

## 사용법

```bash
# 기본 실행 (target="11001", noise=0.1)
python3 hard_mode/qubo_hard_mode.py

# 목표 지정 + noise ratio
python3 hard_mode/qubo_hard_mode.py 11001 0.1

# 길이로 랜덤 목표 생성
python3 hard_mode/qubo_hard_mode.py 100 0.2
```
