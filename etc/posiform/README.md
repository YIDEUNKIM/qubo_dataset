# Posiform Planting QUBO 생성기

## 개요

**Hahn, Pelofske, Djidjev (2023)**에 기반한 planted QUBO 생성 방법론. Planted 2-SAT 인스턴스를 posiform(양의 계수만 갖는 다항식)으로 변환한 후 QUBO로 변환한다. **보조변수가 필요 없으므로 QUBO 크기가 원래 문제와 동일(= n)**하며, target이 **수학적으로 유일한 ground state**임이 보장된다.

## 이론적 배경

### Posiform이란?

**Posiform**은 이진 변수 x_i와 그 부정 x̄_i = 1 - x_i의 곱으로 구성된 다항식으로, **모든 계수가 양수**인 형태:

$$P(x) = \sum_k c_k \cdot \prod_{i \in S_k} l_i, \quad c_k > 0$$

여기서 l_i는 x_i 또는 x̄_i. 핵심 성질: **P(x) >= 0** (모든 x에 대해).

### 2-SAT과 Posiform의 관계

2-SAT clause (l_1 OR l_2)가 위반되는 유일한 할당은 **두 리터럴 모두 False**인 경우. 이를 posiform으로 표현:

| Clause | 위반 조건 | Posiform 항 |
|--------|----------|------------|
| (x_i OR x_j) | x_i=0, x_j=0 | b * x̄_i * x̄_j |
| (x_i OR x̄_j) | x_i=0, x_j=1 | b * x̄_i * x_j |
| (x̄_i OR x_j) | x_i=1, x_j=0 | b * x_i * x̄_j |
| (x̄_i OR x̄_j) | x_i=1, x_j=1 | b * x_i * x_j |

target x\*가 모든 clause를 만족하면 → 모든 posiform 항이 x\*에서 0 → P(x\*) = 0 (최솟값).
target이 2-SAT의 **유일한** 해이면 → 다른 모든 할당에서 P(x) > 0 → x\*가 **유일한** ground state.

### Posiform → QUBO 변환

x̄_i = 1 - x_i를 대입하면 표준 QUBO 형태가 된다:

| Wrong tuple | Posiform | QUBO 변환 |
|:-----------:|:--------:|:----------|
| (0, 0) | b(1-x_i)(1-x_j) | Q_ii -= b, Q_jj -= b, Q_ij += b, const += b |
| (0, 1) | b(1-x_i)x_j | Q_jj += b, Q_ij -= b |
| (1, 0) | bx_i(1-x_j) | Q_ii += b, Q_ij -= b |
| (1, 1) | bx_ix_j | Q_ij += b |

### 왜 Quiet Planting과 다른가?

| 특성 | Quiet Planting (3-SAT) | Posiform (2-SAT) |
|------|:---------------------:|:----------------:|
| SAT 유형 | 3-SAT | 2-SAT |
| 보조변수 | M개 (QUBO 크기 = n+m) | **없음** (QUBO 크기 = n) |
| Ground state 유일성 | 보장 안 됨 (축퇴 존재) | **수학적으로 보장** |
| Field 필요 | 예 (축퇴 해소) | **불필요** |
| SAT 풀이 복잡도 | NP-complete | **P** (다항 시간) |

## 구현 방식

### 2-SAT 유일 해 생성 (Tarjan SCC 기반)

2-SAT의 해는 **함축 그래프(implication graph)**의 강연결 성분(SCC) 분석으로 다항 시간에 구할 수 있다.

**리터럴 인코딩**: 2*i = x_i, 2*i+1 = NOT x_i

**함축 규칙**: clause (l_1 OR l_2)은 두 개의 함축을 생성:
- NOT l_1 → l_2
- NOT l_2 → l_1

**UNSAT 판정**: x_i와 NOT x_i가 같은 SCC에 있으면 모순.

**유일성 검사**: 각 변수를 반대로 강제한 후 SAT인지 확인. 모든 변수에서 UNSAT이면 유일한 해.

### 2단계 Clause 생성

```
Phase 1: 랜덤 clause 추가
  - 무작위 변수 쌍 (i, j) 선택
  - target이 만족하는 clause 생성 (wrong tuple 무작위 선택)
  - 주기적으로 유일성 검사
  - 유일해가 되면 종료

Phase 2: Targeted clause 추가 (Phase 1에서 유일성 미달 시)
  - 아직 "flippable"한 변수에 대해 집중적으로 clause 추가
  - flippable 변수 = 반대값으로 강제해도 SAT인 변수
  - flippable 변수와 다른 변수의 쌍으로 clause 추가
```

### 전체 파이프라인

```
1. Planted 2-SAT 생성 (create_planted_2sat)
   - target이 유일한 해가 될 때까지 clause 추가
   - Tarjan SCC로 해 구성, 유일성 검사

2. Posiform 구성 (posiform_to_qubo)
   - 각 clause에 랜덤 양의 계수 b ~ Uniform(coeff_range) 부여
   - wrong tuple → posiform 항 → x̄=1-x 대입 → QUBO 변환

3. 검증
   - P(target) = 0 확인
   - Brute force (N<=20) 또는 통계적 검증
```

### 핵심 파라미터

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `coeff_range` | (1.0, 3.0) | Posiform 계수 범위 (lo, hi). 넓을수록 에너지 지형 변화 |
| `max_clauses_factor` | 10 | 최대 clause 수 = factor * n |
| `seed` | None | 재현성을 위한 난수 시드 |

### 주요 함수

| 함수 | 설명 |
|------|------|
| `tarjan_scc(adj, n_nodes)` | 반복적(iterative) Tarjan SCC 알고리즘 |
| `build_implication_graph(n, clauses)` | 2-SAT 함축 그래프 구축 |
| `solve_2sat(n, clauses)` | Tarjan SCC 기반 2-SAT 풀이 |
| `check_2sat_uniqueness(n, clauses, solution)` | 해의 유일성 검사 |
| `create_planted_2sat(target, max_clauses_factor, seed)` | Planted 2-SAT 생성 |
| `posiform_to_qubo(n, clauses, coeff_range, seed)` | Posiform → QUBO 변환 |
| `create_qubo_posiform(target, coeff_range, ...)` | **메인 진입점** |

### 반환값

```python
Q, info = create_qubo_posiform(target, coeff_range=(1.0, 3.0))
# Q: QUBO 딕셔너리 {(i,j): weight}
# info: {'n', 'num_clauses', 'is_unique', 'coeff_range', 'target_energy', ...}
```

## SA 난이도 특성

### N 스케일링 실험 결과 (num_reads=1, num_runs=10)

**SA 단 1회 시도로 N=1000까지 100% 성공.**

| N | QUBO 크기 | Sweeps | 성공률 | E비 | 해밍거리 | 생성시간 |
|---:|:---------:|:------:|:------:|:---:|:-------:|:-------:|
| 10 | 10 | 1,000 | 100% (10/10) | 1.0000 | 0.0 | 0.00s |
| 20 | 20 | 1,000 | 100% (10/10) | 1.0000 | 0.0 | 0.01s |
| 50 | 50 | 1,000 | 100% (10/10) | 1.0000 | 0.0 | 0.11s |
| 100 | 100 | 1,000 | 100% (10/10) | 1.0000 | 0.0 | 0.57s |
| 200 | 200 | 2,000 | 100% (10/10) | 1.0000 | 0.0 | 2.35s |
| 300 | 300 | 3,000 | 100% (10/10) | 1.0000 | 0.0 | 6.60s |
| 500 | 500 | 5,000 | 100% (10/10) | 1.0000 | 0.0 | 22.27s |
| 700 | 700 | 7,000 | 100% (10/10) | 1.0000 | 0.0 | 45.13s |
| 1000 | 1000 | 10,000 | 100% (10/10) | 1.0000 | 0.0 | 108.24s |

> 자세한 분석: [`docs/POSIFORM_EXPERIMENT.md`](../docs/POSIFORM_EXPERIMENT.md)

### SA-Easy인 이유

1. **Smooth 에너지 지형**: 모든 posiform 계수가 양수 → metastable 상태 없음
2. **작은 탐색 공간**: 보조변수 없이 QUBO 크기 = n (Quiet Planting은 5.2n)
3. **2-SAT의 다항 시간 풀이 가능성**: 원래 문제가 P에 속함

### 다른 방법론과의 비교

| 방법론 | N=50 | N=100 | N=200 | N=500 | N=1000 |
|--------|:----:|:-----:|:-----:|:-----:|:------:|
| **Posiform** (num_reads=1) | **100%** | **100%** | **100%** | **100%** | **100%** |
| Zero Expectation | 100% | 100% | 100% | 100% | 100% |
| Hard Posiform (lin2, α=0.1) | 100% | 100% | 100% | 100% | 100% |
| Hard Posiform (lin2, α=0.01) | 100% | 100% | 100% | 90% | **40%** |
| Quiet Planting (field=0.5) | 100% | 60% | 0% | — | — |
| Wishart (alpha=0.7) | 70% | 0% | 0% | 0% | 0% |
| McEliece (m=4, t=2) | — | — | — | — | — |

> McEliece는 N별 비교 불가 (k=8 고정, total_vars=33~46). 별도 실험에서 SA ~0%. 전체 비교: [`docs/METHODOLOGY_COMPARISON.md`](../docs/METHODOLOGY_COMPARISON.md)

### 구별 불가능성 분석

N=100, 30회 반복 실험 결과. Posiform QUBO는 순수 random QUBO와 **다수의 경로로 쉽게 구별 가능**하다.

#### 1. 밀도 (Sparsity) — 즉시 탐지

| | Posiform | Random QUBO |
|---|:---:|:---:|
| **Density** (비영 항 / 전체) | **15.1%** | **100%** |

Posiform은 2-SAT clause에 등장하는 변수 쌍만 Q_ij가 비영이다. 나머지 85%의 항이 0이므로, Q 행렬의 sparsity 패턴만으로 즉시 구별된다. Random QUBO는 모든 항이 비영.

#### 2. 비대각 부호 편향 — target pair 노출 (핵심!)

각 clause는 하나의 wrong tuple을 배제하며, Q_ij 기여의 부호는 wrong tuple에 의해 결정된다:

| Target pair (b_i, b_j) | E[Q_ij] | 양수 비율 | 해석 |
|:---:|:---:|:---:|---|
| (0, 0) | **-0.70** | 33% | 음수 편향 |
| (0, 1) | **+0.72** | 67% | 양수 편향 |
| (1, 0) | **+0.70** | 67% | 양수 편향 |
| (1, 1) | **-0.70** | 33% | 음수 편향 |

**규칙: target[i] == target[j]이면 Q_ij < 0, 다르면 Q_ij > 0 (확률 2/3).**

t-test (same vs diff): stat=-50.14, **p ≈ 0**. Q_ij의 부호 패턴에서 target[i]와 target[j]의 동일 여부가 즉시 드러난다.

**이론적 근거**: target pair (t_i, t_j)에 대해 3개의 wrong tuple 중:
- same bits (00 또는 11): wrong tuple 중 2개가 다른 비트 → Q_ij에 음수 기여 2/3
- diff bits (01 또는 10): wrong tuple 중 2개가 같은 비트 → Q_ij에 양수 기여 2/3

#### 3. 대각 편향 — target 비트 직접 노출

| | E[Q_ii\|b=0] | E[Q_ii\|b=1] | 편향 |
|---|:---:|:---:|:---:|
| **Posiform** | **+4.47** | **-4.55** | **-9.02** |

t-test: p ≈ 0. 대각 항의 부호만으로 각 변수의 target 값을 직접 추측 가능.

#### 4. 행 합 편향

| | E[row\_sum\|b=0] | E[row\_sum\|b=1] | t-test p |
|---|:---:|:---:|:---:|
| **Posiform** | **+4.87** | **-4.71** | ≈ 0 |

#### 5. 구별 불가능한 속성

| 속성 | 구별 가능? | 비고 |
|------|:---:|------|
| **밀도 (sparsity)** | **O (trivial)** | 15% vs 100% |
| **비대각 부호 패턴** | **O** | target pair에 의존 |
| **대각 편향** | **O** | target bit 직접 노출 |
| **행 합 편향** | **O** | target bit 직접 노출 |
| Rank 구조 | X | 둘 다 full-rank (100/100) |
| 비대각 양수/음수 비율 | X | 50.2% / 49.8% ≈ 50/50 |

**종합**: Posiform QUBO는 밀도, 대각 편향, 비대각 부호 패턴 등 **다수의 경로로 쉽게 탐지 가능**하다. 이는 설계 의도와 일치한다 — Posiform의 목적은 통계적 은닉이 아니라 **유일한 GS의 수학적 보장 + 보조변수 없는 컴팩트 QUBO 생성**이다.

### 장점과 한계

**장점**:
- 보조변수 없음 → QUBO 크기 = n (매우 컴팩트)
- 유일한 ground state 수학적 보장
- 2-SAT이므로 planted 인스턴스 생성이 효율적
- SA 정확성 검증(sanity check)에 최적

**한계**:
- SA에 대해 trivially easy → 난이도 벤치마크로 부적합
- 연속적 난이도 파라미터(alpha 등)가 없어 난이도 조절 불가
- 2-SAT 기반이므로 상전이 구조가 없음 (easy-hard-easy 프로파일 없음)
- 대규모(N>1000)에서 유일성 검사의 O(n^2) 비용으로 생성 시간 증가

## 파일 구성

| 파일 | 역할 |
|------|------|
| `qubo_posiform.py` | 생성기 (Tarjan SCC + 2-SAT solver + posiform 변환) |
| `test_posiform.py` | SA 실험 (N scaling, coeff sweep, 비교) |

## 참고 문헌

### 핵심 논문

1. **Hahn, G., Pelofske, E., & Djidjev, H.** "Using 2-SAT to generate QUBO instances with known optimal solutions." *2023 IEEE International Conference on Quantum Computing and Engineering (QCE)*, 2023.
   - Planted 2-SAT → posiform → QUBO 파이프라인
   - 보조변수 없이 N-bit QUBO 생성
   - 유일한 ground state 수학적 보장

### 2-SAT 알고리즘

2. **Tarjan, R. E.** "Depth-first search and linear graph algorithms." *SIAM Journal on Computing*, 1(2), 146-160, 1972.
   - 강연결 성분(SCC) 분해 알고리즘
   - 2-SAT 풀이의 기반

3. **Aspvall, B., Plass, M. F., & Tarjan, R. E.** "A linear-time algorithm for testing the truth of certain quantified Boolean formulas." *Information Processing Letters*, 8(3), 121-123, 1979.
   - 함축 그래프 기반 2-SAT 선형 시간 알고리즘

### Posiform 이론

4. **Boros, E. & Hammer, P. L.** "Pseudo-Boolean optimization." *Discrete Applied Mathematics*, 123(1-3), 155-225, 2002.
   - Posiform의 정의와 QUBO 변환 이론

5. **Rosenberg, I. G.** "Reduction of bivalent maximization to the quadratic case." *Cahiers du Centre d'Etudes de Recherche Operationnelle*, 17, 71-74, 1975.
   - 고차 의사불 함수의 2차화 (posiform은 이미 2차이므로 Rosenberg 불필요)

### QUBO 벤치마크

6. **Pelofske, E., Hahn, G., & Djidjev, H.** "Quantum annealing algorithms for Boolean tensor networks." *Scientific Reports*, 12, 8539, 2022. — 관련 planted QUBO 벤치마크 연구.

## 사용법

```bash
# 기본 생성 (target="10110", coeff_range=(1.0, 3.0))
python3 posiform/qubo_posiform.py 10110

# 시드 지정
python3 posiform/qubo_posiform.py 10110 42

# 길이로 랜덤 목표 생성
python3 posiform/qubo_posiform.py 100

# N 스케일링 실험
python3 posiform/test_posiform.py --scaling 10

# 계수 범위 sweep 실험
python3 posiform/test_posiform.py --coeff 10

# 4-way 비교 (Posiform vs Quiet vs Wishart vs ZeroExp)
python3 posiform/test_posiform.py --compare
```
