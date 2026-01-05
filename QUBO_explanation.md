# Quantum Annealing & QUBO

## 1. 양자 어닐링 (Quantum Annealing)

양자 어닐링(Quantum Annealing)은 **조합 최적화 문제(Combinatorial Optimization Problem)**의 전역 최적해(Global Minimum)를 찾기 위한 메타휴리스틱 알고리즘입니다.

-   **기본 원리**: 금속을 고온으로 가열했다가 천천히 식히면(Annealing) 내부 입자들이 에너지가 가장 낮은 안정된 상태(결정 구조)를 찾아가는 물리적 현상에서 착안했습니다.
-   **양자 터널링 (Quantum Tunneling)**: 고전적인 어닐링(Simulated Annealing)은 에너지 장벽을 넘기 위해 열적 요동(Thermal Fluctuation)을 이용하지만, 양자 어닐링은 **양자 터널링 효과**를 이용해 에너지 장벽을 뚫고 지나가 더 효율적으로 최적해를 찾을 수 있습니다.

## 2. QUBO (Quadratic Unconstrained Binary Optimization)

QUBO는 **제약 조건이 없는 2차 이진 최적화** 문제를 의미합니다. 양자 어닐링 머신(예: D-Wave)이 풀 수 있는 문제의 표준 형태입니다.

### 수식 정의
목적 함수(Objective Function)는 다음과 같이 정의됩니다.

$$ E(x) = x^T Q x = \sum_{i} Q_{ii} x_i + \sum_{i < j} Q_{ij} x_i x_j $$

여기서:
-   $x$: 이진 변수 벡터 ($x_i \in \{0, 1\}$)
-   $Q$: 실수 대칭 행렬 (상삼각 행렬로 표현하기도 함)
-   $Q_{ii}$: 선형 항(Linear term)의 계수 (Bias)
-   $Q_{ij}$: 2차 항(Quadratic term)의 계수 (Coupler strength)

우리의 목표는 에너지 $E(x)$를 **최소화**하는 $x$를 찾는 것입니다.

---

## 3. 특정 해를 갖는 QUBO 문제 만들기

우리의 목표는 $x = b_0 b_1 \dots b_{n-1}$ (여기서 $b_i \in \{0, 1\}$)가 정답(최소 에너지 상태)이 되도록 $Q$ 행렬을 설계하는 것입니다.

### 아이디어
각 비트 $x_i$가 독립적으로 목표값 $b_i$와 같을 때 에너지가 최소가 되도록 만들면 됩니다. 즉, 각 비트 $i$에 대해 다음 함수를 최소화합니다.

$$ f_i(x_i) = (x_i - b_i)^2 $$

이 함수는 $x_i = b_i$일 때 0(최소)이고, 다를 때 1이 됩니다.

### 수식 전개
$$ (x_i - b_i)^2 = x_i^2 - 2x_i b_i + b_i^2 $$

이진 변수의 성질($x_i \in \{0, 1\}$)에 의해 $x_i^2 = x_i$입니다. 따라서,

$$ = x_i - 2x_i b_i + b_i^2 = x_i(1 - 2b_i) + b_i^2 $$

여기서 $b_i^2$은 상수이므로 최적화 위치에 영향을 주지 않습니다. 따라서 우리가 최소화해야 할 실제 목적 함수 부분은 다음과 같습니다.

$$ \text{minimize } \sum_i x_i (1 - 2b_i) $$

### Q 행렬 구성
위 식은 $x_i$와 $x_j$의 곱(상호작용) 항이 없고, 오직 $x_i$ 항만 존재합니다. 즉, $Q$ 행렬은 **대각 행렬(Diagonal Matrix)**이 됩니다.

-   **대각 성분 ($Q_{ii}$)**: $1 - 2b_i$
    -   만약 $b_i = 1$이면, $Q_{ii} = 1 - 2 = -1$. ($-x_i$를 최소화하려면 $x_i=1$이어야 함)
    -   만약 $b_i = 0$이면, $Q_{ii} = 1 - 0 = 1$. ($x_i$를 최소화하려면 $x_i=0$이어야 함)
-   **비대각 성분 ($Q_{ij}$)**: $0$ (변수 간 상호작용 없음)

### 결론
원하는 해가 $b$일 때, QUBO 행렬 $Q$는 다음과 같이 설정합니다.
$$ Q_{ii} = \begin{cases} -1 & \text{if } b_i = 1 \\ 1 & \text{if } b_i = 0 \end{cases} $$
$$ Q_{ij} = 0 \quad (\text{for } i \neq j) $$

## 4. [심화] 복잡한 QUBO 문제 만들기 (상호작용 추가)

위의 기본 구성은 비트 간의 상호작용($Q_{ij}$)이 없는 단순한 형태입니다. 하지만 실제 문제나 더 복잡한 에너지 랜드스케이프를 만들기 위해, **정답은 유지하면서 변수 간의 상호작용을 추가**할 수 있습니다. `qubo_generator.py` 코드는 이 방식을 사용합니다.

### 기본 아이디어
기본 목적 함수 $\sum (x_i - b_i)^2$에, 정답 상태에서 에너지가 0이 되는 "벌점 항(Penalty Term)"을 추가합니다.

#### 1. 두 비트의 목표 값이 같은 경우 ($b_i = b_j$)
두 변수 $x_i, x_j$가 같은 값을 가질 때 에너지가 낮아지도록 합니다.
$$ P_{equal} = C(x_i - x_j)^2 $$
- $x_i = x_j$이면 0 (벌점 없음)
- $x_i \ne x_j$이면 $C$ (벌점 부과)

전개하면:
$$ C(x_i^2 - 2x_i x_j + x_j^2) = C(x_i - 2x_i x_j + x_j) $$
- $Q_{ii}$에 $+C$ 추가
- $Q_{jj}$에 $+C$ 추가
- $Q_{ij}$에 $-2C$ 추가

#### 2. 두 비트의 목표 값이 다른 경우 ($b_i \neq b_j$)
두 변수 $x_i, x_j$가 다른 값을 가질 때 에너지가 낮아지도록 합니다.
$$ P_{diff} = C(x_i + x_j - 1)^2 $$
- $x_i \ne x_j$ (하나는 0, 하나는 1)이면 합이 1 $\rightarrow$ 0 (벌점 없음)
- $x_i = x_j$ (둘 다 0이거나 둘 다 1)이면 합이 0 또는 2 $\rightarrow$ $C$ (벌점 부과)

전개하면:
$$ C(x_i^2 + x_j^2 + 1 + 2x_i x_j - 2x_i - 2x_j) $$
$$ = C(2x_i x_j - x_i - x_j) + \text{const} $$
- $Q_{ii}$에 $-C$ 추가
- $Q_{jj}$에 $-C$ 추가
- $Q_{ij}$에 $+2C$ 추가

이러한 상호작용 항들을 무작위로 추가하면, **정답(Global Optimum)은 변하지 않지만** $Q$ 행렬이 0이 아닌 비대각 성분을 많이 가지는(Dense) 복잡한 문제가 됩니다.

---

## 5. [고급] 기댓값 0 QUBO 생성 및 전략 패턴 (`qubo_zero_expectation.py`)

`qubo_zero_expectation.py`는 단순한 상호작용 추가를 넘어, 생성된 문제가 **통계적으로 편향되지 않도록(Zero Expectation)** 정밀하게 설계된 페널티 비율을 사용합니다. 이를 위해 소프트웨어 공학적 설계인 **전략 패턴(Strategy Pattern)**을 도입했습니다.

### 5.1 전략 패턴 (Strategy Pattern) 적용
코드의 유연성과 유지보수성을 높이기 위해, 페널티 비율 계산 로직을 별도의 클래스로 분리했습니다.

-   **`PenaltyModel` (Abstract Base Class)**: 모든 페널티 전략이 따라야 할 인터페이스를 정의합니다. `get_ratios(target_pair)` 메서드를 통해 상황별 페널티 비율을 동적으로 반환합니다.
-   **`DefaultZeroExpectationModel`**: 선형 프로그래밍(LP)을 통해 최적화된 "기댓값 0" 비율을 제공하는 기본 구현체입니다.
    -   예: 정답이 (0,0)일 때, (0,1) 오답에는 1.0배, (1,1) 오답에는 1.65배의 페널티를 부과합니다.
-   **`SimpleUniformModel`**: 모든 오답에 대해 균등한(1.0) 페널티를 부과하는 단순 모델입니다. (비교 실험용)

### 5.2 의존성 주입 (Dependency Injection)
`create_qubo_precise` 함수는 특정 모델에 강결합되지 않고, `PenaltyModel` 타입의 객체를 인자로 받습니다. 이를 통해 사용자는 언제든지 자신만의 페널티 전략을 정의하여 주입할 수 있습니다.

```python
# 기본 기댓값 0 모델 사용
Q_standard = create_qubo_precise(target)

# 단순 균등 모델 사용
Q_simple = create_qubo_precise(target, model=SimpleUniformModel())
```
