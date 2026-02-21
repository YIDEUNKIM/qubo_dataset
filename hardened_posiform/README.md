# Hardened Posiform Planting QUBO 생성기

## 개요

**Pelofske, Hahn, Djidjev (2024)**에 기반한 hardened QUBO 생성 방법론. 기본 posiform planting (Hahn et al. 2023)이 SA에 trivially easy한 문제를 만드는 한계를 극복하기 위해, **discrete coefficient random QUBO**를 추가하여 에너지 지형을 복잡하게 만든다. Ground state의 유일성은 수학적으로 보장된다.

## 이론적 배경

### 기본 Posiform의 한계

기본 posiform planting은 모든 계수가 양수인 smooth한 에너지 지형을 만들어, SA가 N=1000에서도 단 1회 시도로 100% 성공한다. 이는 벤치마크로서 부적합하다.

### Hardening 전략

논문의 핵심 아이디어: **random QUBO의 rugged 에너지 지형**과 **posiform의 ground state 보장**을 결합.

```
Q_final = Σ R_i + α × P

R_i: 각 subgraph의 discrete coefficient random QUBO
P:   posiform planted QUBO (target이 유일한 ground state)
α:   posiform scaling coefficient (작을수록 어려움)
```

### Ground State 유일성 증명 (논문 Section 2.2)

1. `sub_i(X*)` minimizes `R_i` (by construction — brute force로 최적해 계산)
2. `X*` uniquely minimizes `P` (posiform guarantee — 2-SAT 유일성)
3. 따라서: `Q(X*) = Σ R_i(X*) + α·P(X*) < Σ R_i(X̂) + α·P(X̂)` for any `X̂ ≠ X*`

직관적으로: random QUBO 부분은 target에서 이미 최솟값이고, posiform 부분은 target에서만 0이므로, 어떤 다른 상태도 더 낮은 에너지를 가질 수 없다.

### 난이도 조절 파라미터

| 파라미터 | 값 | 효과 |
|---------|-----|------|
| `posiform_scale` (α) | 0.01 (hardest), 0.1 (easier) | α가 작을수록 random QUBO가 지배 → rugged landscape |
| `coeff_type` | `lin2` ({-1,+1}), `lin20` ({-1,...,+1}) | lin2가 더 어려움 (이산적, 더 많은 frustration) |
| `max_subgraph_size` | ≤23 (default 15) | 각 random subproblem 크기 (brute force 가능해야 함) |

## 알고리즘

```
Input: n (변수 수), max_subgraph_size, coeff_type, posiform_scale (α), seed

Step 1: 변수 분할
  - n개 변수를 max_subgraph_size 이하의 disjoint 그룹으로 분할
  - 논문: Kernighan-Lin recursive bisection (hardware graph용)
  - 구현: complete graph이므로 순차적 동일 크기 분할

Step 2: 각 subgraph에 random QUBO 생성
  - discrete coefficient set에서 무작위 선택
  - lin2: {-1, +1}
  - lin20: {-1.0, -0.9, ..., 0.9, 1.0} (21개 값)
  - complete graph 내 모든 (linear + quadratic) 항에 할당

Step 3: 각 subproblem의 ground state 계산
  - brute force 완전 탐색 (2^k, k ≤ 23)
  - 사전 매핑 (var_to_shift) + bit 연산으로 최적화

Step 4: Ground state 결합 → target bitstring
  - 모든 subproblem 최적해를 concatenate

Step 5: Posiform planted QUBO 생성
  - target을 유일한 해로 갖는 2-SAT → posiform → QUBO

Step 6: 결합
  - Q_final = Σ R_i + α × P_posiform
```

## 구현 방식

### 주요 함수

| 함수 | 설명 |
|------|------|
| `partition_variables(n, max_subgraph_size)` | n개 변수를 disjoint 그룹으로 분할 |
| `generate_random_qubo(variables, coeff_type, seed)` | subgraph에 discrete coefficient QUBO 생성 |
| `solve_qubo_brute_force(Q, variables)` | Brute force 최적해 계산 (≤23 변수) |
| `create_qubo_hardened_posiform(n, ...)` | **메인 진입점** |
| `verify_hardened_ground_state(Q, target, n)` | Ground state 검증 (brute force / statistical) |

### 반환값

```python
Q, info = create_qubo_hardened_posiform(
    n, max_subgraph_size=15, coeff_type='lin2',
    posiform_scale=0.1, seed=42
)
# Q: QUBO 딕셔너리 {(i,j): weight}
# info: {
#   'n', 'target', 'num_partitions', 'partition_sizes',
#   'coeff_type', 'posiform_scale', 'posiform_is_unique',
#   'posiform_num_clauses', 'random_total_energy',
#   'random_total_degenerate', 'target_energy', 'num_qubo_terms'
# }
```

## SA 실험 결과

### 실험 1: Sweep 전이 (N=500)

SA sweep 수에 따른 per-sample ground-state sampling rate. 논문 (Pelofske 2024) Fig 7-8 재현.

**실험 A** (8 instances, 100 reads/sweep):

| Config | 1 | 5 | 10 | 50 | 100 | 500 | 1000 | 5000 |
|--------|--:|--:|---:|---:|----:|----:|-----:|-----:|
| lin2, α=0.01 | 0.0% | 8.8% | 4.5% | 4.9% | 4.6% | 13.2% | 18.9% | 29.8% |
| lin2, α=0.1 | 0.0% | 16.9% | 25.2% | 64.1% | 86.5% | 99.5% | 99.9% | 100.0% |
| lin20, α=0.01 | 0.0% | 0.0% | 0.0% | 2.6% | 8.4% | 34.8% | 50.0% | 73.1% |
| lin20, α=0.1 | 0.0% | 76.4% | 83.2% | 96.4% | 99.2% | 100.0% | 100.0% | 100.0% |

**실험 B** (5 instances, 20 reads/sweep — METHODOLOGY_COMPARISON 기준):

| Config | 10 sw | 50 sw | 100 sw | 500 sw | 1000 sw | 5000 sw |
|--------|:-----:|:-----:|:------:|:------:|:-------:|:-------:|
| lin2, α=0.01 | 2.0% | 1.0% | 4.0% | 6.0% | 12.0% | 24.0% |
| lin2, α=0.1 | 24.0% | 59.0% | 83.0% | 100% | 100% | 100% |
| lin20, α=0.01 | 0.0% | 2.0% | 12.0% | 31.0% | 52.0% | 68.0% |
| lin20, α=0.1 | 87.0% | 98.0% | 99.0% | 100% | 100% | 100% |

### 논문 핵심 발견 재현 확인

| 논문 발견 | 코드 결과 | 일치 |
|-----------|-----------|:---:|
| α=0.01이 α=0.1보다 어려움 | lin2: 24% vs 100% (5000sw) | O |
| lin2가 lin20보다 어려움 | α=0.01에서 lin2(24%) < lin20(68%) | O |
| Sweep 수 증가 → S-curve 전이 | 모든 config에서 단조 증가 | O |
| α=0.01은 매우 많은 sweep 필요 | lin2,α=0.01: 5000sw에서도 24% | O |

### 핵심 관찰

1. **α=0.01이 α=0.1보다 훨씬 어려움**: lin2,α=0.01은 5000 sweeps에서도 24~30%에 그침
2. **lin2가 lin20보다 어려움**: 이산적 계수 {-1,+1}이 더 많은 frustration 생성
3. **비단조적 행동**: lin2,α=0.01에서 sweeps=5(8.8%) → sweeps=10(4.5%) 구간에서 성공률 하락 관찰 — 논문 Section 3.3과 일치. SA가 중간 sweeps에서 local minimum에 갇히는 현상.
4. **난이도 순서**: lin2,α=0.01 (hardest) > lin20,α=0.01 > lin2,α=0.1 > lin20,α=0.1 (easiest)

### N 스케일링 (SA 성공률)

num_reads=50, num_sweeps=1000:

| Config | N=50 | N=100 | N=200 | N=500 | N=1000 |
|--------|:----:|:-----:|:-----:|:-----:|:------:|
| **lin2, α=0.01** | 100% | 100% | 100% | 90% | **40%** |
| **lin2, α=0.1** | 100% | 100% | 100% | 100% | 100% |
| Plain Posiform | 100% | 100% | 100% | 100% | 100% |

Hardened posiform (α=0.01)은 N=500부터 SA 성공률이 감소하기 시작하며, N=1000에서 40%로 SA-moderate 영역에 진입한다.

### 다른 방법론과의 위치

| 방법론 | SA 난이도 | GS 보장 | 난이도 조절 | N=500 성공률 |
|--------|:---------:|:-------:|:----------:|:----------:|
| Plain Posiform | trivially easy | **수학적 (유일)** | X | 100% |
| **Hardened (α=0.1)** | easy~moderate | **수학적 (유일)** | **O** | 100% |
| **Hardened (α=0.01)** | **moderate** | **수학적 (유일)** | **O** | **90%** |
| Zero Expectation | trivially easy | 수학적 | X | 100% |
| Quiet Planting (f=0.5) | medium | 조건부 | O | — |
| Wishart (α=0.7) | **hard** | 수학적 (유한정밀도) | O | 0% |
| McEliece (m=4,t=2) | **hard** | 조건부 | O | — |

> 전체 방법론의 정량적 비교: [`docs/METHODOLOGY_COMPARISON.md`](../docs/METHODOLOGY_COMPARISON.md)

## 파일 구성

| 파일 | 역할 |
|------|------|
| `qubo_posiform_hardened.py` | 생성기 (partition → random QUBO → brute force solve → posiform → glue) |
| `test_posiform_hardened.py` | SA 실험 (sweep 전이, N-scaling, hardened vs plain 비교) |
| `papers/` | 참고 논문 PDF |
| `results/` | 생성된 QUBO 파일 |

## 참고 문헌

### 핵심 논문

1. **Pelofske, E., Hahn, G., & Djidjev, H.** "Increasing the Hardness of Posiform Planting Using Random QUBOs for Programmable Quantum Annealer Benchmarking." *arXiv:2411.03626*, 2024.
   - Random QUBO + posiform planting 결합으로 SA 난이도 증가
   - discrete coefficient sets (lin2, lin20)
   - posiform scaling coefficient α로 난이도 제어

### 기반 논문

2. **Hahn, G., Pelofske, E., & Djidjev, H.** "Using 2-SAT to generate QUBO instances with known optimal solutions." *IEEE QCE*, 2023.
   - 기본 posiform planting 방법론
   - 2-SAT → posiform → QUBO 파이프라인

## 사용법

```bash
# 기본 생성 (n=30, lin2, α=0.1)
python3 hardened_posiform/qubo_posiform_hardened.py

# N, coeff_type, α, seed 지정
python3 hardened_posiform/qubo_posiform_hardened.py 500 lin2 0.01 42
python3 hardened_posiform/qubo_posiform_hardened.py 100 lin20 0.1

# Sweep 전이 실험
python3 hardened_posiform/test_posiform_hardened.py --sweep 10

# N-Scaling 실험
python3 hardened_posiform/test_posiform_hardened.py --scaling 10

# Hardened vs Plain 비교
python3 hardened_posiform/test_posiform_hardened.py --compare 10
```
