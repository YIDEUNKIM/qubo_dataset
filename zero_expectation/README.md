# Zero-Expectation QUBO 생성기

## 개요

**랜덤 QUBO와 통계적으로 구별 불가능한** planted QUBO를 생성하는 방법론. Q 행렬의 비대각 계수가 **E[q_ij] = 0** 을 만족하도록, 각 큐빗 쌍의 3가지 오답 상태에 부여하는 페널티 비율을 **선형 프로그래밍(LP)**으로 최적화한다.

## 이론적 배경

### QUBO와 Planted 문제

QUBO: E(x) = x^T Q x, x_i in {0,1}을 최소화하는 문제.

**Planted QUBO**: 목표 비트스트링 x\*가 주어졌을 때, x\*가 ground state인 Q를 구성하되, Q의 계수 분포가 랜덤 QUBO와 구별 불가능하도록 만드는 것이 핵심 질문이다.

### 페널티 기반 QUBO 구성법

각 큐빗 쌍 (i, j)에서 목표 상태 (b_i, b_j) 외의 3가지 오답 상태에 독립적인 양의 페널티 r > 0를 부여한다. 오답 상태 (s_i, s_j)에 대한 페널티 r은 Q 행렬에 다음과 같이 분배된다:

| 오답 상태 (s_i, s_j) | 페널티 함수 | Q_ii 기여 | Q_jj 기여 | Q_ij 기여 |
|:-------------------:|:----------:|:---------:|:---------:|:---------:|
| (0, 0) | r(1-x_i)(1-x_j) | -r | -r | +r |
| (0, 1) | r(1-x_i)x_j | 0 | +r | -r |
| (1, 0) | rx_i(1-x_j) | +r | 0 | -r |
| (1, 1) | rx_ix_j | 0 | 0 | +r |

목표 상태 (b_i, b_j)에서는 모든 페널티 함수가 0이므로, target에서의 에너지 기여 = 0. 오답 상태에서는 해당 페널티 > 0이므로 에너지 증가 → **target이 ground state**.

### LP 최적화: E[q_ij] = 0 조건

Q 행렬의 비대각 항 q_ij는 3개 페널티의 부호 합으로 결정된다. 목표 쌍이 4가지 ((0,0), (0,1), (1,0), (1,1))이고 각각 동일한 확률로 등장한다고 가정하면, 전체 기댓값 E[q_ij] = 0 조건은 각 목표 쌍의 페널티 비율에 대한 선형 제약이 된다.

이를 LP로 풀면 다음 **최적 비율**이 도출된다:

```
Target (0,0): r(0,1) = 1.00, r(1,0) = 1.00, r(1,1) = 1.65
Target (0,1): r(0,0) = 2.00, r(1,0) = 1.00, r(1,1) = 1.68
Target (1,0): r(0,0) = 2.00, r(0,1) = 1.00, r(1,1) = 1.68
Target (1,1): r(0,0) = 1.00, r(0,1) = 3.00, r(1,0) = 3.00
```

각 오답에 부여하는 실제 페널티: `r = base * ratio`, base ~ Uniform(1, 3).

### 대각 편향과 비대각 편향의 동시 달성 불가능성

`test_diagonal_zero.py`에서 수학적으로 증명:

- **E[Q_ij] = 0** (비대각 무편향)과 **E[Q_ii|b=0] = E[Q_ii|b=1]** (대각 무편향)을 **동시에** 달성하는 것은, 양의 페널티 제약 하에서 **불가능**하다.
- 증명: 동시 달성 조건을 전개하면 r_2 + r_6 + r_7 + r_11 = 0이 필요하나, 모두 양수이므로 모순.

따라서 여러 모델이 서로 다른 트레이드오프를 제공한다.

## 구현 방식

### Strategy Pattern

`PenaltyModel` 추상 기본 클래스(ABC)를 통해 다양한 페널티 전략을 플러그인 방식으로 교체 가능:

```python
class PenaltyModel(ABC):
    @abstractmethod
    def get_ratios(self, target_pair: tuple) -> dict:
        """penalty_state -> ratio 매핑 반환"""
        pass
```

### 제공 모델 4종

#### 1. DefaultZeroExpectationModel (기본)

LP로 최적화된 비율. **E[q_ij] = 0** 보장 (비대각 기댓값 0).

```python
self._ratios_table = {
    (0, 0): {(0, 1): 1.00, (1, 0): 1.00, (1, 1): 1.65},
    (0, 1): {(0, 0): 2.00, (1, 0): 1.00, (1, 1): 1.68},
    (1, 0): {(0, 0): 2.00, (0, 1): 1.00, (1, 1): 1.68},
    (1, 1): {(0, 0): 1.00, (0, 1): 3.00, (1, 0): 3.00},
}
```

- 장점: 비대각 항의 부호/크기로 target pair를 추론할 수 없음
- 단점: 대각 항 E[Q_ii|b=0] != E[Q_ii|b=1] (대각 편향 존재)

#### 2. ZeroOffDiagonalModel

개별 q_ij 기댓값 0. **double-flip ratio = single-flip 두 ratio의 합**.

```python
(0, 0): {(0, 1): 1.0, (1, 0): 1.0, (1, 1): 2.0},  # 2.0 = 1.0 + 1.0
(0, 1): {(0, 0): 1.0, (1, 0): 2.0, (1, 1): 1.0},
(1, 0): {(0, 0): 1.0, (0, 1): 2.0, (1, 1): 1.0},
(1, 1): {(0, 0): 2.0, (0, 1): 1.0, (1, 0): 1.0},
```

- 조건 도출: q_ij 부호 패턴이 (0,0)→+r, (0,1)→-r, (1,0)→-r, (1,1)→+r이므로, E[q_ij] = -R(01) - R(10) + R(11) = 0 → R(11) = R(01) + R(10)
- 장점: |q_ij|의 multiset이 항상 {1, 1, 2}로 target pair와 **완전히 독립** (모든 모멘트 동일)
- 단점: 대각 편향 존재. 다만 탐지 SNR이 sqrt(n)으로, off-diagonal 누출(SNR ~ n)보다 유리

#### 3. BalancedModel

Minimax 균형: max(|off_bias|, |diag_bias|) = 5/3 ≈ 1.67로 **이론적 하한** 달성.

```python
(0, 0): {(0, 1): 1.0,  (1, 0): 1.0,  (1, 1): 2.0},
(0, 1): {(0, 0): 7/6,  (1, 0): 1.0,  (1, 1): 1.0},
(1, 0): {(0, 0): 7/6,  (0, 1): 1.0,  (1, 1): 1.0},
(1, 1): {(0, 0): 1.0,  (0, 1): 4/3,  (1, 0): 4/3},
```

- 비대각 편향과 대각 편향을 동시에 최소화하는 LP minimax 풀이
- 양쪽 편향 모두 5/3으로 동일하게 억제

#### 4. SimpleUniformModel (기준선)

모든 오답에 동일 페널티 (ratio = 1.0). 편향 최적화 없음. 비교 실험용.

### Ising-Derived 모드

`create_qubo_ising_derived(target)`: Ising 모델에서 유도된 QUBO.

```
J_ij = alpha * s_i * s_j,  alpha ~ Uniform(1, 3)
s_i = 2*x_i - 1 치환으로 QUBO 변환:
  Q_ij = -4 * J_ij
  Q_ii += 2 * J_ij, Q_jj += 2 * J_ij
```

- **E[Row Sum] = 0** 보장 (Ising 구조에 의해 행 합 기댓값이 대칭)
- Ground state가 target과 일치함이 수학적으로 보장 (J_ij * s_i * s_j > 0)

### QUBO 생성 파이프라인

```
1. 모든 큐빗 쌍 (i, j) 순회  (i < j)
   ├── density 확률로 상호작용 추가 여부 결정
   ├── target pair (b_i, b_j) 확인
   ├── PenaltyModel.get_ratios(target_pair) → 3개 오답 비율
   └── 각 오답 상태에 대해:
       ├── r = Uniform(1, 3) * ratio
       └── 페널티 함수 전개 → Q_ii, Q_jj, Q_ij에 기여 누적

2. Q 행렬 반환 (상삼각 딕셔너리 형태)
```

### 핵심 파라미터

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `density` | 1.0 | 큐빗 쌍 상호작용 추가 확률. 1.0이면 모든 쌍 사용 |
| `base_range` | (1, 3) | 기본 페널티 r의 균등분포 범위 |
| `model` | DefaultZeroExpectationModel | 페널티 비율 전략 |
| `balance_rows` | False | True이면 Ising-derived 모드 사용 |

### 주요 함수

| 함수 | 설명 |
|------|------|
| `create_qubo_precise(target, density, model)` | **메인 진입점**. 밀도와 모델 지정 |
| `create_qubo_ising_derived(target, density)` | Ising 모델 기반 행-균형 QUBO |
| `solve_brute_force(Q, n)` | N <= 20일 때 완전 탐색 검증 |
| `batch_test(num_tests, n_bits)` | 배치 정확도 테스트 |
| `large_scale_analysis(n_tests, n_bits)` | 대규모 계수 분포 분석 |
| `compare_with_random_qubo(n_bits, n_tests)` | 생성 QUBO vs 순수 랜덤 QUBO 비교 |

### 반환값

```python
Q = create_qubo_precise(target, density=1.0)
# Q: QUBO 딕셔너리 {(i,j): weight} (i <= j, 상삼각)
```

## SA 난이도 특성

### SA에 대해 Trivially 쉬운 문제

num_reads=1 (단 1회 SA 시도), num_sweeps=1000에서:

| N | 성공률 | 에너지비 | 해밍거리 | 시간 |
|---:|:------:|:-------:|:-------:|:----:|
| 10 | **100%** (10/10) | 1.0000 | 0 | 0.00s |
| 20 | **100%** (10/10) | 1.0000 | 0 | 0.00s |
| 50 | **100%** (10/10) | 1.0000 | 0 | 0.00s |
| 100 | **100%** (10/10) | 1.0000 | 0 | 0.01s |
| 200 | **100%** (10/10) | 1.0000 | 0 | 0.03s |
| 300 | **100%** (10/10) | 1.0000 | 0 | 0.09s |
| 500 | **100%** (10/10) | 1.0000 | 0 | 0.23s |
| 1000 | **100%** (10/10) | 1.0000 | 0 | 0.95s |

N=1000에서도 100%로, 7개 방법론 중 가장 안정적으로 SA-trivial.

### 원인 분석: 왜 SA가 즉시 푸는가

#### 1. 로컬 미니마가 1개뿐 (= target 자체)

N=16에서 전수 탐색(2^16 = 65,536개 상태) 결과:

| 방식 | 로컬 미니마 수 | 비율 |
|------|:------------:|:----:|
| **Zero-Expectation** | **1개** | 0.0015% |
| Wishart (alpha=0.7) | **224개** | 0.3418% |

Zero-Expectation은 에너지 랜드스케이프 전체에 **함정(trap)이 하나도 없다**. SA가 아무 상태에서 출발해도 1-bit flip으로 에너지가 낮아지는 방향이 항상 존재하고, 그 방향을 따라가면 반드시 target에 도달한다.

#### 2. Frustration이 0

N=20에서 각 pair (i,j)의 off-diagonal 커플링이 target 상태와 일관되는지 분석:

| 방식 | Frustrated pairs |
|------|:---------------:|
| **Zero-Expectation** | **0 / 190** (0%) |
| Wishart (alpha=0.7) | 25 / 190 (~13%) |

**Frustration = 서로 다른 pair가 모순되는 상태를 원하는 것.** Zero-Expectation은 모든 pair가 독립적으로 동일한 target 방향을 가리킨다. 변수 간 충돌이 전혀 없으므로 로컬 미니마가 생길 수 없다.

Wishart는 J_ij = -(1/N) sum_mu W_imu W_jmu로 구성되어, 공유 벡터 W를 통해 커플링이 상관되고, 일부 pair가 서로 모순되는 방향을 선호한다.

#### 3. 에너지가 해밍거리에 선형 비례 (Funnel 구조)

N=16 전수 탐색, 해밍거리별 에너지 프로파일:

**Zero-Expectation:**

| 해밍거리 d | 평균 에너지 | 표준편차 | 신호/노이즈 |
|:---------:|:---------:|:------:|:---------:|
| 0 (target) | -109.41 | 0.00 | - |
| 1 | -72.38 | 8.90 | - |
| 2 | -37.34 | 10.11 | - |
| 3 | -4.29 | 9.60 | - |

에너지 증가 기울기: **~37/bit flip**, 셸 내 표준편차: **~9** → **신호(37) >> 노이즈(9)**. SA가 어디서든 "target 방향이 어딘지" 명확히 감지.

**Wishart (alpha=0.7):**

| 해밍거리 d | 평균 에너지 | 표준편차 | 신호/노이즈 |
|:---------:|:---------:|:------:|:---------:|
| 0 (target) | -2.43 | 0.00 | - |
| 1 | -1.16 | 0.49 | - |
| 2 | -0.06 | 1.04 | - |
| 3 | +0.87 | 1.50 | - |

에너지 증가 기울기: **~1.3/bit flip**, 셸 내 표준편차: **~1.5** → **신호(1.3) < 노이즈(1.5)**. SA가 방향을 잡을 수 없음. d=3에서 이미 d=0보다 낮은 에너지 상태가 존재.

#### 4. 에너지 갭이 거대함

N=20 전수 탐색:

| 방식 | Ground (E0) | 1st Excited (E1) | 절대 갭 | 상대 갭 |
|------|:----------:|:---------------:|:------:|:------:|
| **Zero-Expectation** | -329.22 | -290.30 | **38.92** | **11.8%** |
| Wishart (alpha=0.7) | -4.71 | -4.71 | **~0** | **~0%** |

#### 5. Single-bit flip 에너지 장벽

N=20, target으로부터 1비트 뒤집었을 때:

| 방식 | 최소 Delta E | 평균 Delta E | 최대 Delta E |
|------|:----------:|:----------:|:----------:|
| **Zero-Expectation** | **+39.67** | +54.21 | +93.25 |
| Wishart (alpha=0.7) | +0.61 | +1.21 | +2.37 |

ZeroExp의 최소 delta(39.67)가 Wishart의 최대 delta(2.37)보다 **17배** 큼.

### 구조적 원인 요약

```
Zero-Expectation (깔때기):           Wishart (울퉁불퉁):

에너지                                에너지
  |                                    |
  |\                                   |\   /\   /\
  | \                                  | \_/  \_/  \
  |  \                                 |  ^    ^    \___  <- ground state
  |   \                                | trap  trap
  |    \___  <- ground state           |
  +-------- 상태 공간                   +-------- 상태 공간

모든 경로가 바닥으로 향함             여러 함정이 존재, SA 탈출 불가
```

**핵심**: 페널티 구조에서 모든 pair가 독립적으로 target 방향을 가리키므로 frustration = 0, local minima = 1. Q 행렬의 1차 통계량(평균, 분산)을 위장하는 것은 에너지 **함수** 구조(상관 구조, frustration, local minima)를 위장하는 것과 전혀 다르다.

### 다른 방법론과의 비교

| 방법론 | SA 난이도 | GS 보장 | 구별 불가능 | 난이도 조절 | N=500 성공률 |
|--------|:---------:|:-------:|:----------:|:----------:|:----------:|
| **Zero Expectation** | **trivial** | 수학적 | **E[q_ij]=0** | X | **100%** |
| Posiform | trivial | 수학적 (유일) | 미분석 | X | 100% |
| Hard Posiform (α=0.01) | moderate | 수학적 (유일) | X | **O** | 90% |
| Quiet Planting (f=0.5) | medium | 조건부 | **O** (α<3.86) | **O** | — |
| Wishart (α=0.7) | **hard** | 수학적 (유한정밀도) | X | **O** | 0% |
| McEliece (m=4,t=2) | **hard** | 조건부 | 미분석 | **O** | — |

> 전체 비교: [`docs/METHODOLOGY_COMPARISON.md`](../docs/METHODOLOGY_COMPARISON.md)

### 벤치마크로서의 위치

**장점**:
- E[q_ij] = 0으로 Q 행렬 비대각 항이 랜덤 QUBO와 통계적으로 유사
- LP 최적화된 비율이 수학적으로 보장됨
- Strategy Pattern으로 다양한 페널티 전략 실험 가능
- 솔버 정확성 검증(sanity check)에 적합 — 이 문제를 못 풀면 솔버 구현에 문제가 있는 것

**한계**:
- SA가 num_reads=1로도 N=500을 즉시 풀어버림 → 난이도 벤치마크로 부적합
- 비대각 무편향(E[q_ij]=0)과 대각 무편향을 동시에 달성 불가 → 완전한 통계적 위장 아님
- 솔버 간 성능 변별력이 전혀 없음

## 파일 구성

| 파일 | 역할 |
|------|------|
| `qubo_zero_expectation.py` | 생성기 (Strategy Pattern + PenaltyModel 4종 + Ising-derived) |
| `test_zero_expectation.py` | SA N-scaling 실험 |
| `test_diagonal_zero.py` | 대각 편향 분석 + 불가능성 증명 + 대안 모델 탐색 |
| `analyze_q_structure.py` | Q 행렬 구조 비교 (ZeroExp vs Wishart: 고유값, rank, 행간 상관) |
| `results/` | 생성된 QUBO 파일 |

## 참고 문헌

### 직접 참고

1. **프로젝트 내부 PDF**: `기댓값0만들기_2트.pdf` — Zero-expectation QUBO 구성의 LP 최적화 방법론. 페널티 비율 도출 과정, 비대각 기댓값 0 조건의 수학적 유도.
2. **프로젝트 내부 PDF**: `읽어봐줘!.pdf` — QUBO 벤치마크 데이터셋 생성 동기 및 요구사항.

### QUBO 일반

3. **Kochenberger, G. A., et al.** "The unconstrained binary quadratic programming problem: a survey." *Journal of Combinatorial Optimization*, 28(1), 58-81, 2014.
4. **Glover, F., Kochenberger, G., & Du, Y.** "Quantum Bridge Analytics I: a tutorial on formulating and using QUBO models." *4OR*, 17(4), 335-371, 2019.

### 의사 불 최적화 (Pseudo-Boolean Optimization)

5. **Boros, E. & Hammer, P. L.** "Pseudo-Boolean optimization." *Discrete Applied Mathematics*, 123(1-3), 155-225, 2002. — 페널티 함수의 이론적 기반. 이진 변수 다항식의 최적화 프레임워크.

### Ising-QUBO 변환

6. 표준 치환: s_i = 2x_i - 1. Ising H = -sum J_ij s_i s_j → QUBO Q_ij = -4J_ij, Q_ii += 2J_ij.

## 사용법

```bash
# 기본 실행 (기본 목표 해)
python3 zero_expectation/qubo_zero_expectation.py

# 이진 목표 지정
python3 zero_expectation/qubo_zero_expectation.py 10110

# 랜덤 목표 (길이 50)
python3 zero_expectation/qubo_zero_expectation.py 50

# Ising-derived 행-균형 모드
python3 zero_expectation/qubo_zero_expectation.py 10110 balance

# SA 스케일링 실험
python3 zero_expectation/test_zero_expectation.py 10,20,50,100 10

# 대각 편향 분석 + 불가능성 증명
python3 zero_expectation/test_diagonal_zero.py

# Q 행렬 구조 분석 (ZeroExp vs Wishart)
python3 zero_expectation/analyze_q_structure.py
```
