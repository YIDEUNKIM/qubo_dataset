# Hardened Posiform Planting QUBO Generator

## 개요

**Pelofske, Hahn, Djidjev (2024)** 기반: 기존 Posiform Planting의 SA-trivial 문제를 해결하기 위해, **Random QUBO subproblem + Posiform overlay**를 결합하여 난이도를 높인 방법론.

기존 Posiform은 모든 계수가 양수(smooth landscape)여서 SA가 100% 성공했다. Hardened Posiform은 random QUBO의 rugged landscape 위에 posiform을 겹쳐서 ground state는 보장하면서 난이도를 올린다.

**기반 논문**: Pelofske, Hahn, Djidjev, "Increasing the Hardness of Posiform Planting Using Random QUBOs for Programmable Quantum Annealer Benchmarking", 2024.

## Posiform과의 차이

| 특성 | Posiform (Plain) | Posiform Hardened |
|------|:---:|:---:|
| **구조** | Posiform만 사용 | Random QUBO + Posiform overlay |
| **에너지 지형** | Smooth (양수 계수만) | Rugged (random QUBO가 지배) |
| **SA 난이도** | Easy (N=1000도 100%) | 조절 가능 (alpha로 제어) |
| **GS 보장** | 수학적 (유일, P(x*)=0) | 수학적 (유일, 논문 Section 2.2) |
| **보조변수** | 없음 | 없음 |
| **QUBO 크기** | n | n |
| **난이도 파라미터** | 없음 (항상 easy) | alpha (posiform_scale), coeff_type |
| **기반 논문** | Hahn et al. 2023 | Pelofske et al. 2024 |

## 방법론

### 파이프라인

```
1. 변수 분할           2. Random QUBO 생성       3. Subproblem 풀기
n개 변수를 k개의     각 subgraph에 discrete    각 R_i의 ground state를
disjoint subgraph    coefficient random QUBO   brute force로 정확 계산
으로 분할             R_i 생성
    |                     |                        |
    v                     v                        v
[x_0..x_14]           R_1: {-1,+1} 계수       gs_1 = "01011..."
[x_15..x_29]          R_2: {-1,+1} 계수       gs_2 = "10100..."
...                   ...                      ...

4. Target 조합         5. Posiform 생성          6. 결합
subproblem 해를       target으로 planted       Q_final = Sigma R_i + alpha * P
concatenate           2-SAT posiform P 생성
    |                     |                        |
    v                     v                        v
target = gs_1+gs_2+..  P(target) = 0           alpha 작을수록 어려움
```

### Ground State 유일성 증명 (논문 Section 2.2)

임의의 X != X* 에 대해:

```
Q(X*) = Sigma R_i(X*) + alpha * P(X*)
      = Sigma R_i(X*) + 0          (P(X*) = 0, posiform 성질)

Q(X)  = Sigma R_i(X) + alpha * P(X)
      >= Sigma R_i(X*) + alpha * P(X)   (X*가 각 R_i의 최적해이므로)
      > Sigma R_i(X*)                   (P(X) > 0, X != X*에서 posiform)
      = Q(X*)
```

따라서 X*가 **유일한** ground state.

### 난이도 조절 원리

- `posiform_scale` (alpha): P의 기여도를 결정
  - alpha 크면 → P가 지배 → smooth landscape → SA-easy
  - alpha 작으면 → R_i가 지배 → rugged landscape → SA-hard
  - **alpha = 0.01이 가장 어려움** (논문 실험)
- `coeff_type`: random QUBO 계수 집합
  - `lin2`: {-1, +1} — 이산 계수 (가장 rugged)
  - `lin20`: {-1.0, -0.9, ..., 0.9, 1.0} — 21단계 (상대적으로 smooth)
- `subgraph_size`: 각 subgraph 크기
  - 크면 → subproblem이 더 어려워짐 (brute force 비용 증가)
  - 기본값 15 (2^15 = 32,768 상태로 brute force 가능)

## 파라미터

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `n` | 30 | 총 변수 수 |
| `max_subgraph_size` | 15 | 각 subgraph 최대 크기 (brute force 한계 ~23) |
| `coeff_type` | 'lin2' | 계수 집합: 'lin2' ({-1,+1}) 또는 'lin20' (21단계) |
| `posiform_scale` | 0.1 | alpha: posiform 스케일링 (작을수록 어려움) |
| `posiform_coeff_range` | (1.0, 1.0) | posiform clause 계수 범위 |
| `seed` | None | 난수 시드 |

## 한계점

### 1. SA 난이도 검증 미완

논문에서는 D-Wave QPU에서 hardened가 plain보다 어렵다고 보고하지만, 본 구현에서 SA 스케일링 실험은 아직 체계적으로 수행하지 않았다.

### 2. Subgraph 크기 제한

각 subgraph의 ground state를 brute force로 계산하므로, subgraph_size <= 23 제약이 있다. 더 큰 subgraph는 SA나 다른 솔버로 approximate할 수 있지만 GS 보장이 약해진다.

### 3. 구조적 정보 노출

Random QUBO R_i는 disjoint subgraph 위에만 정의되므로, Q 행렬의 block-diagonal 구조가 솔버에 힌트를 줄 수 있다. 논문에서는 hardware graph의 Kernighan-Lin 분할을 사용하여 이를 완화한다.

### 4. Posiform 기여도 vs 난이도 트레이드오프

alpha가 너무 작으면 posiform의 GS 강제력이 약해져서, SA가 random QUBO의 local minimum에 갇혀 posiform이 의도한 target을 찾지 못할 수 있다. 반면 alpha가 크면 posiform이 지배하여 easy해진다.

## 사용법

```bash
# 기본 생성 (N=30, lin2, alpha=0.1)
python3 posiform_hardened/qubo_posiform_hardened.py

# 파라미터 지정 (N, coeff_type, alpha, seed)
python3 posiform_hardened/qubo_posiform_hardened.py 50 lin2 0.01 42

# lin20 계수로 생성
python3 posiform_hardened/qubo_posiform_hardened.py 100 lin20 0.1

# SA 실험: sweep 전이
python3 posiform_hardened/test_posiform_hardened.py --sweep

# SA 실험: N 스케일링
python3 posiform_hardened/test_posiform_hardened.py --scaling

# SA 실험: Hardened vs Plain 비교
python3 posiform_hardened/test_posiform_hardened.py --compare
```

## 파일 구조

```
posiform_hardened/
  qubo_posiform_hardened.py      # 생성기 (메인 코드)
  test_posiform_hardened.py      # SA 실험 프레임워크
  run_experiments.py             # 실험 스크립트
  README.md                      # 이 문서
  papers/
    Pelofske2022_QA_Boolean_Tensor.pdf
  results/
    hardened_*.txt               # 생성된 QUBO 파일
```

## 참고 문헌

1. **Pelofske, E., Hahn, G., & Djidjev, H.** "Increasing the Hardness of Posiform Planting Using Random QUBOs for Programmable Quantum Annealer Benchmarking." [arXiv:2411.03626](https://arxiv.org/abs/2411.03626), 2024. Published in [npj Unconventional Computing](https://www.nature.com/articles/s44335-025-00032-6), 2025.
2. **Hahn, G., Pelofske, E., & Djidjev, H.** "Using 2-SAT to generate QUBO instances with known optimal solutions." QCE 2023. (기반이 되는 plain posiform)
