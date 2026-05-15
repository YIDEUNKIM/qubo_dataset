# Truth Table QUBO 생성기

## 개요

임의의 **진리표**(모든 비트스트링 → 에너지 매핑)로부터 QUBO를 구성한다. 모든 $2^n$개 비트스트링에 에너지를 직접 지정할 수 있어, 에너지 랜드스케이프를 **완전히 제어**하는 벤치마크 생성이 가능하다. 두 가지 모드 제공:

- **정확 모드 (Exact)**: Möbius 변환 + Rosenberg 차수축소. 에너지 100% 일치, 보조변수 $O(n \cdot 2^n)$ 폭발.
- **근사 모드 (Approx)**: 제약 최소제곱(QP). 보조변수 0개, ground state 보장, 에너지 근사.

```
입력:  {000→3, 001→4, 010→4, 011→5, 100→3, 101→5, 110→2, 111→1}
출력:  Q 행렬 (QUBO), ground state = 111 (energy=1) 보장
```

## 이론적 배경

### 정확 모드: Möbius 변환 + Rosenberg 차수축소

#### Step 1: Möbius 변환 (진리표 → 다중선형 다항식)

임의의 함수 $f: \{0,1\}^n \to \mathbb{R}$은 **유일한** 다중선형 다항식으로 표현 가능:

$$f(x) = \sum_{S \subseteq [n]} c_S \prod_{i \in S} x_i$$

계수 $c_S$는 포함-배제 원리(Möbius 역변환)로 계산:

$$c_S = \sum_{T \subseteq S} (-1)^{|S|-|T|} f(\mathbf{1}_T)$$

Fast Möbius 알고리즘으로 $O(n \cdot 2^n)$에 계산. 구현: `mobius_transform()`

**n=3 예시:**

```
진리표: {000→3, 001→4, 010→4, 011→5, 100→3, 101→5, 110→2, 111→1}

f(x) = 3 + x₂ + x₃ + x₁x₃ - 2x₁x₂ - 3x₁x₂x₃
         ↑상수   ↑선형      ↑이차         ↑고차(3차)
```

#### Step 2: 항 분류

| 차수 | 처리 | QUBO 배치 |
|------|------|-----------|
| 0 (상수) | 에너지 오프셋 | Q에 미포함 |
| 1 (선형) | $Q_{ii} += c_i$ | 대각 |
| 2 (이차) | $Q_{ij} += c_{ij}$ | 비대각 |
| ≥3 (고차) | Rosenberg 차수축소 | 보조변수 도입 |

구현: `classify_terms()`

#### Step 3: Rosenberg 차수축소

$k$차 단항식 $c \cdot x_{i_1} x_{i_2} \cdots x_{i_k}$에 대해:

1. 보조변수 $y = x_{i_1} \cdot x_{i_2}$ 도입
2. 패널티 $M(x_{i_1}x_{i_2} - 2x_{i_1}y - 2x_{i_2}y + 3y)$ 추가 → $y = x_{i_1}x_{i_2}$ 강제
3. 원래 항을 $c \cdot y \cdot x_{i_3} \cdots x_{i_k}$로 치환
4. 여전히 차수 > 2이면 반복 (체이닝)

$k$차 단항식 1개당 보조변수 $k-2$개 필요.

**3가지 차수축소 전략:**

| 전략 | 설명 | 구현 |
|------|------|------|
| `original` | 매번 새 보조변수 생성 | `rosenberg_reduce()` |
| `cache` | 동일 쌍 `(x_a, x_b)` 재활용 | `rosenberg_reduce_reuse()` |
| `greedy` (기본) | 빈도 기반 쌍 선택 + 재활용 | `rosenberg_reduce_greedy()` |

**보조변수 절감 효과:**

| n | original | cache | greedy | greedy 절감률 |
|:-:|:--------:|:-----:|:------:|:------------:|
| 5 | 23 | 9 | 5 | 78% |
| 7 | 201 | 35 | 15 | 92.5% |
| 8 | 522 | 66 | 22 | 95.8% |

#### Step 4: 패널티 강도

```
M = max(tt_range, max_y(S_y)) + 1

tt_range = max(E) - min(E)                      ← 진리표 에너지 범위
S_y = Σ |coeff| for reduced terms involving y    ← 보조변수 y 관련 축소 항 계수 합
```

두 조건 모두 만족해야 함:
1. M > tt_range: 진리표 에너지 범위 초과
2. M > max(S_y): Möbius 변환 고차항 계수 합 초과 (진리표 범위보다 클 수 있음)

조건 2가 없으면 SA가 보조변수 제약을 위반하여 더 낮은 에너지를 찾을 수 있음. 구현: `compute_penalty_strength()`

### 근사 모드: 제약 최소제곱 (QP)

보조변수 없이 $n$-변수 QUBO를 직접 구한다. 고차항을 2차화하는 대신, 진리표 에너지에 **가장 가까운 2차 다항식**을 최소제곱으로 직접 fitting한다.

$$\min_Q \sum_x \left( E_Q(x) - E_{\text{truth}}(x) \right)^2 \quad \text{s.t.} \quad E_Q(\text{target}) + \varepsilon \leq E_Q(x) \;\; \forall x \neq \text{target}$$

- **자유 파라미터**: $p = n(n+1)/2$개 (Q 행렬 upper triangular)
- **제약조건**: $2^n - 1$개 (target이 모든 다른 상태보다 $\varepsilon$ 이상 낮은 에너지)

이하 Step 1~3에서 이 문제를 어떻게 푸는지 설명한다.

#### Step 1: Feature Matrix 구성

QUBO 에너지의 일반식:

$$E_Q(x) = \sum_{i} q_{ii} x_i + \sum_{i \lt j} q_{ij} x_i x_j$$

여기서 $x_i \in \{0,1\}$은 비트스트링의 각 비트(주어진 값), $q_{ij}$는 Q 행렬 성분(미지수)이다. 미지수는 총 $p = n(n+1)/2$개. **2차까지만 사용하므로** $x_i x_j x_k$ 같은 3차 이상의 항이 없고, 보조변수가 필요 없다.

각 비트스트링에 비트값을 대입하면, $q_{ij}$에 대한 연립방정식이 된다. QUBO는 상수항이 없으므로($E_Q(\mathbf{0}) = 0$ 항상), 진리표에서 `offset = E_truth(target)`을 빼서 target 에너지가 0이 되도록 조정한 $\mathbf{e}_{adj}$를 우변에 놓는다.

**n=3 예시**: 정확 모드에서 Möbius 변환하면 3차항 $-3 x_0 x_1 x_2$가 발생하여 보조변수가 필요한 진리표를 사용한다.

```
진리표: {000→0, 100→1, 010→2, 110→1, 001→2, 101→3, 011→3, 111→-1}
target = 111 (에너지 -1, 최소)
offset = -1,  e_adj = [1, 2, 3, 2, 3, 4, 4, 0]

정확 모드: f(x) = x₀ + 2x₁ + 2x₂ - 2x₀x₁ - x₁x₂ - 3x₀x₁x₂
                                                       ↑ 3차항 → 보조변수 필요!
근사 모드: E_Q(x) = q₀₀x₀ + q₀₁x₀x₁ + q₀₂x₀x₂ + q₁₁x₁ + q₁₂x₁x₂ + q₂₂x₂
                     ↑ 2차까지만 → 3차항 무시 → 보조변수 불필요
```

$E_Q(x)$에 미지수 6개($q_{00}, q_{01}, q_{02}, q_{11}, q_{12}, q_{22}$). 모든 $2^3 = 8$개 비트스트링에 비트값을 대입하면:

```
           q00  q01  q02  q11  q12  q22   e_adj
x=000: [   0    0    0    0    0    0 ]     1
x=100: [   1    0    0    0    0    0 ]     2
x=010: [   0    0    0    1    0    0 ]     3
x=110: [   1    1    0    1    0    0 ]     2
x=001: [   0    0    0    0    0    1 ]     3
x=101: [   1    0    1    0    0    1 ]     4
x=011: [   0    0    0    1    1    1 ]     4
x=111: [   1    1    1    1    1    1 ]     0
```

이것이 $A \mathbf{q} = \mathbf{e}_{adj}$ (방정식 8개, 미지수 6개)이다.

**일반화**: $n$-변수 QUBO에서 $A \in \mathbb{R}^{2^n \times p}$. 각 행은 비트스트링 $x$를 대입했을 때 각 $q_{ij}$에 곱해지는 계수($x_i$ 또는 $x_i x_j$)를 나열한 것이다.

#### Step 2: 정규방정식 (무제약 최소제곱)

$A$는 $2^n$행 $\times$ $p$열이고, $n \geq 2$이면 $2^n \gt p$이므로 방정식이 미지수보다 많은 **과결정 시스템**이다. 정확한 해는 일반적으로 없고, 잔차 제곱합을 최소화하는 최소제곱 해를 구한다:

$$\mathbf{q}^{*} = (A^T A)^{-1} A^T \mathbf{e}_{adj}$$

rank$(A) = p$이면 $(A^T A)$가 정칙이므로 **해가 유일**하다.

**n=3 예시 계속**: 방정식 8개, 미지수 6개 → 과결정. 정규방정식을 풀면:

```
q00=2.29, q01=-3.86, q02=-1.86, q11=3.29, q12=-2.86, q22=3.29

검증 (E_Q = A·q vs e_adj):
  x=000: E_Q= 0.00, e_adj=1 → 잔차=-1.00  (상수항 없음 → 불가피)
  x=100: E_Q= 2.29, e_adj=2 → 잔차=+0.29
  x=010: E_Q= 3.29, e_adj=3 → 잔차=+0.29
  x=110: E_Q= 1.71, e_adj=2 → 잔차=-0.29
  x=001: E_Q= 3.29, e_adj=3 → 잔차=+0.29
  x=101: E_Q= 3.71, e_adj=4 → 잔차=-0.29
  x=011: E_Q= 3.71, e_adj=4 → 잔차=-0.29
  x=111: E_Q= 0.29, e_adj=0 → 잔차=+0.29  ← target도 오차 있음
```

잔차 발생 원인: 원래 진리표에 3차 성분($-3 x_0 x_1 x_2$)이 있는데, 2차까지만으로 표현하니 근사 오차가 생긴다. 이 오차가 보조변수 대신 지불하는 **대가**다.

구현: `np.linalg.lstsq(A, e_adj)`

#### Step 3: Ground State 제약 확인 및 SLSQP

무제약 해가 ground state를 보장하지 못할 수 있다. 이 경우 **제약 최적화**로 전환:

$$\min_{\mathbf{q}} \frac{1}{2} \mathbf{q}^T (A^T A) \mathbf{q} - (A^T \mathbf{e}_{adj})^T \mathbf{q}$$

$$\text{s.t.} \quad (A_x - A_{\text{target}}) \mathbf{q} \geq \varepsilon \quad \forall x \neq \text{target}$$

- 목적함수는 잔차 제곱합 $\frac{1}{2} \lVert A\mathbf{q} - \mathbf{e}_{adj} \rVert^2$과 동치 (상수 차이)
- 제약은 선형: $E_Q(x) - E_Q(\text{target}) \geq \varepsilon$
- **Iterative Cutting Plane**: 매 반복 위반+마진 부족 제약만 active set에 추가하여 SLSQP 호출 (최대 100회, 정체 감지 시 조기 종료)

**n=3 예시 계속**: Step 2의 해에서 $E_Q(111) = 0.29$, $E_Q(000) = 0$ → 갭 = $-0.29 \lt \varepsilon$ → target이 ground state가 아님 → **위반**. SLSQP가 파라미터를 미세 조정:

```
q00=2.34, q01=-4.01, q02=-2.01, q11=3.34, q12=-3.01, q22=3.34

  x=000: E_Q= 0.00   갭=0.01=ε  ✓  (등호: 경계에서 멈춤)
  x=100: E_Q= 2.34   갭=2.35    ✓
  x=010: E_Q= 3.34   갭=3.35    ✓
  x=110: E_Q= 1.67   갭=1.68    ✓
  x=001: E_Q= 3.34   갭=3.35    ✓
  x=101: E_Q= 3.67   갭=3.68    ✓
  x=011: E_Q= 3.67   갭=3.68    ✓
  x=111: E_Q=-0.01   ← target (유일한 최소)
```

구현: `scipy.optimize.minimize(method='SLSQP')`

**핵심 장점**: 보조변수 0개 → QUBO 크기 = n. 최적화 버전으로 n≤23까지 실용적.

#### 최적화 버전 (`create_qubo_approx_optimized`)

근사 모드의 최적화 버전. 동일한 Q 행렬을 생성하면서 메모리·연산을 개선하고, feasible fallback으로 수렴을 보장한다.

**최적화 내용:**

| 항목 | 원본 | 최적화 |
|------|------|--------|
| Phase 1 (lstsq) | A 생성 + lstsq: O(n⁴·2^n) | **A-free**: ATA 해석적 + ATe 스트리밍: O(n²·2^n) |
| Phase 2 메모리 | A ($2^n \times p$) 유지 | ATA ($p \times p$) + ATe ($p$) 만 유지 |
| Phase 2 SLSQP | 전체 active 제약 한번에 | **마진 포함 절단 평면** + 정체 감지 |
| Phase 2 수렴 | lstsq start만 (일부 seed 미수렴) | lstsq start → **feasible fallback** (수렴 보장) |
| 배치 에너지 | Python 이중 루프 | 벡터화 행렬 곱 (`_batched_energies`) |
| 순서 보존율 | O($4^n$) 쌍별 비교 | O($n \cdot 2^n$) Kendall tau |

**Feasible Fallback 메커니즘:**

lstsq start + SLSQP cutting plane이 정체(동일 위반 수 3회 반복)하면:
- **Feasible start**: 대각 Q 구성 ($q_{ii} = -\varepsilon$ if $x_i^* = 1$, $+\varepsilon$ otherwise) → target이 확정적 최소, 위반 0개

RMSE: lstsq 수렴 시 ~1.4, feasible fallback 시 ~2.9. 모든 seed에서 ground state 보장.

**Seed별 수렴 차이의 원인과 해결:**

lstsq start + SLSQP cutting plane만 사용하면 일부 seed에서 ground state를 보장하지 못한다. 원인은 SLSQP의 수치 정밀도 경계 문제:

1. SLSQP는 제약 $E_Q(x) - E_Q(\text{target}) \geq \varepsilon$을 내부적으로 ~$10^{-8}$ 정밀도로 만족시킴
2. 그러나 `_batched_energies`로 재평가하면 다른 부동소수점 연산 경로를 거치므로 gap이 미세하게 달라짐
3. Seed에 따라 에너지 랜드스케이프가 달라지고, 일부 제약이 정확히 $\varepsilon$ 경계에 위치
4. 재평가 시 gap이 $\varepsilon$ 위 → 수렴 / $\varepsilon$ 아래 → 위반 판정 → 이미 active set에 있으므로 새 제약 없음 → 무한 루프

Feasible fallback은 이 문제를 근본적으로 회피: 대각 Q만 사용하면 위반이 0개이므로 SLSQP 자체가 불필요하다.

**Ground State 유일성 증명:**

알고리즘의 세 실행 경로 모두 출력 Q가 target $x^*$를 유일한 ground state로 갖는다.

*경로 1–2 (lstsq feasible / SLSQP 수렴).* 수렴 판정과 Phase 3 검증이 동일한 결정적 함수 `_batched_energies`에 동일 입력 $\mathbf{q}$로 호출된다. 결정적 프로그램은 동일 입력에 동일 출력을 생성하므로 (IEEE 754-2019 §11 [15]), 수렴 시 $\min_{x \neq x^*} \text{gap}(x) \geq \varepsilon - \varepsilon \cdot 10^{-4} > 0$이면 Phase 3에서도 동일 값. $\square$

*경로 3 (feasible fallback).* 별도 증명이 필요하다.

**정리.** $q_{ii} = -\varepsilon \text{ if } x_i^* = 1$, $q_{ii} = +\varepsilon \text{ if } x_i^* = 0$, $q_{ij} = 0 \; (i \neq j)$, $\varepsilon > 0$이면:

$$\forall x \neq x^*: \quad E_Q(x) - E_Q(x^*) = \varepsilon \cdot d_H(x, x^*) \geq \varepsilon$$

*증명.* 비대각 성분이 0이므로 $E_Q(x) = \sum_i q_{ii} x_i$. 따라서:

$$E_Q(x) - E_Q(x^*) = \sum_i q_{ii}(x_i - x_i^*)$$

$x_i = x_i^*$인 비트의 기여는 0. $x_i \neq x_i^*$인 각 비트의 기여:

$$x_i^* = 1: \; (-\varepsilon)(0-1) = +\varepsilon \qquad x_i^* = 0: \; (+\varepsilon)(1-0) = +\varepsilon$$

따라서 $E_Q(x) - E_Q(x^*) = \varepsilon \cdot |\{i : x_i \neq x_i^*\}| = \varepsilon \cdot d_H(x, x^*)$. $x \neq x^* \Rightarrow d_H \geq 1$. $\square$

**IEEE 754 부동소수점 안전성.** 위 정리는 실수 산술 기준이다. 구현에서 `gs_verified = (min_gap > 0)` 판정이 부동소수점 오차로 뒤집히지 않음을 보인다.

`_batched_energies`는 $n$차원 내적 $\mathbf{b} \cdot \mathbf{q}_{\text{diag}}$를 계산한다 ($b_i \in \{0,1\}$, $q_i = \pm\varepsilon$). IEEE 754 배정밀도에서 $n$차원 내적의 전방향 오차 상한은 (Higham, 2002, Thm 3.1 [13]):

$$|\text{fl}(\mathbf{b} \cdot \mathbf{q}) - \mathbf{b} \cdot \mathbf{q}| \leq \gamma_n \sum_i |b_i q_i|, \quad \gamma_n = \frac{nu}{1-nu}, \; u = 2^{-53}$$

Gap 계산 오차:

$$|\text{gap}_{\text{fl}} - \text{gap}_{\text{exact}}| \leq 2\gamma_n \cdot n\varepsilon + \tfrac{1}{2}\text{ulp}(\text{gap})$$

$n = 23$, $\varepsilon = 0.01$일 때: $\gamma_{23} \approx 2.6 \times 10^{-15}$이므로

$$\text{오차} \leq 2 \times 2.6 \times 10^{-15} \times 0.23 \approx 1.2 \times 10^{-15}$$

실제 gap $\geq \varepsilon = 10^{-2}$이므로 **13자릿수 여유**. 부동소수점 오차가 부호를 뒤집을 수 없다. $\square$

**검증 결과**: n=3~23, 80개 테스트 전부 ground state 보장 (gs_verified=True).

**성능 비교:**

| n | t_orig(s) | t_opt(s) | speedup |
|:-:|:---------:|:--------:|:-------:|
| 6 | 0.08 | 0.01 | **5.5x** |
| 8 | 1.6 | 0.9 | **1.8x** |
| 18 | (불가) | ~30 | - |
| 20 | (불가) | ~50 | - |
| 22 | (불가) | ~120 | - |

n≥13에서 원본은 A 행렬 구성($2^n \times p$) 자체가 병목. 최적화 버전만 실용적.

### 왜 두 모드가 필요한가?

|                    | 정확 original | 정확 greedy | 근사 (QP)         |
|--------------------|:------------:|:-----------:|:-----------------:|
| 보조변수            | $O(n \cdot 2^n)$ | **78~96% 절감** | **0개**       |
| QUBO 크기 (n=8)    | 530변수      | **30변수**   | **8변수**         |
| 에너지 정확도       | 100% 일치    | 100% 일치    | 근사 (RMSE ~1.4/~2.9) |
| Ground state       | 보장          | 보장         | **보장** (QP 제약) |
| SA 성공률 (n=5)    | 10.70%       | 18.30%      | **68.00%**        |
| 실용 한계           | n ≤ 4~5     | **n ≤ 7~8** | **n ≤ 23** (최적화) |

- **n = 3**: greedy 정확 모드 권장 (SA 성공률 최고, 에너지 100% 일치)
- **n ≥ 4**: approx 근사 모드 권장 (안정적 50~70% 성공률, QUBO 크기 = n)

## 에너지 함수 프리셋

### Random Landscape

```python
E(target) = 0
E(x != target) = uniform(0.1, 5.0)
```

target이 유일한 ground state이고, 나머지는 0.1~5.0 범위에서 균등 분포. 구현: `preset_random_landscape()`

### Multi-Valley

```python
targets[0]: global minimum (energy = 0)
targets[1:]: local minima (energy = gap)
나머지: Hamming 거리 비례 에너지 + 장벽
```

local minima 수 ↑ = SA가 계곡에 갇힘. 양자 tunneling 벤치마크. 구현: `preset_multi_valley()`

## 구현 방식

### 전체 파이프라인

```
정확 모드:
  진리표 → Möbius 변환 → 항 분류 → Rosenberg 차수축소 → QUBO 조립 → 전수검사 검증

근사 모드 (원본):
  진리표 → Feature matrix A 구성 → lstsq → 제약 위반 시 SLSQP → Q dict 생성 → 검증

근사 모드 (최적화):
  진리표 → A-free (ATA 해석적 + ATe 스트리밍) → lstsq → SLSQP cutting plane (정체 시 feasible fallback) → 검증
```

### 핵심 파라미터

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `truth_table` | (필수) | dict / list / callable. 비트스트링 → 에너지 매핑 |
| `n` | 자동 추론 | 변수 개수 |
| `epsilon` | 0.01 | 근사 모드의 ground state 에너지 갭 하한 |
| `barrier_height` | 5.0 | Multi-Valley 프리셋의 장벽 높이 |
| `seed` | None | 재현성을 위한 난수 시드 |
| `reduce_strategy` | `'greedy'` | Rosenberg 전략 (`'original'`/`'cache'`/`'greedy'`) |
| `verbose` | True | 진행 상황 출력 여부 |

### 주요 함수

| 함수 | 설명 |
|------|------|
| `mobius_transform(truth_table, n)` | 진리표 → 다중선형 다항식 계수 (fast Möbius, $O(n \cdot 2^n)$) |
| `classify_terms(coefficients)` | 계수를 상수/선형/이차/고차로 분류 |
| `rosenberg_reduce(higher_order, n)` | 고차항 → 2차화 (매번 새 보조변수) |
| `rosenberg_reduce_reuse(higher_order, n)` | 동일 쌍 보조변수 재활용 |
| `rosenberg_reduce_greedy(higher_order, n)` | 빈도 기반 쌍 선택 + 재활용 (기본값) |
| `compute_penalty_strength(truth_table, n)` | 최소 패널티 강도 M 계산 |
| `assemble_qubo(...)` | 선형 + 이차 + 축소된 항 + 패널티 → Q dict 조립 |
| `compute_aux_values(x_orig, aux_info)` | 원래 변수 → 올바른 보조변수 값 계산 |
| `verify_qubo(Q, truth_table, n, aux_info, offset)` | 모든 $2^n$ 비트스트링 전수검사 |
| `preset_random_landscape(n, target, seed)` | Random Landscape 에너지 함수 생성 |
| `preset_multi_valley(n, targets, gap, ...)` | Multi-Valley 에너지 함수 생성 |
| `create_qubo_approx(truth_table, n, epsilon)` | **근사 모드 진입점** |
| `create_qubo_approx_optimized(truth_table, n, epsilon, use_cplex)` | **근사 모드 (최적화)**: A-free + lstsq start + feasible fallback + Kendall tau. n≤23 |
| `create_qubo_truthtable(truth_table, n, seed, reduce_strategy)` | **정확 모드 진입점** (strategy: original/cache/greedy) |

### 반환값

```python
# 정확 모드
Q, info = create_qubo_truthtable(truth_table)
# Q: QUBO 딕셔너리 {(i,j): weight}
# info: {'n_original', 'n_aux', 'n_total', 'offset', 'penalty_M',
#        'ground_state', 'n_higher_order', 'aux_info', 'target',
#        'reduce_strategy'}

# 근사 모드
Q, info = create_qubo_approx(truth_table)
# Q: QUBO 딕셔너리 {(i,j): weight}
# info: {'n_original', 'n_aux'(=0), 'n_total', 'offset', 'rmse',
#        'max_error', 'energy_gap', 'order_preservation', 'gs_verified',
#        'ground_state', 'target', 'aux_info'}
```

## SA 난이도 특성

### 공통 SA 설정

```
솔버: neal.SimulatedAnnealingSampler (D-Wave)
num_reads = 100 (인스턴스당 SA 샘플 수)
num_sweeps = 5000
instances = 10 (각 설정마다 랜덤 target으로 10회 반복)
성공률 = GS 찾은 read 수 / 전체 read 수 (10 × 100 = 1000 samples)
프리셋: preset_random_landscape (E(target)=0, 나머지=uniform(0.1, 5.0))
```

### 정확 모드: N-Scaling (greedy)

| n | QUBO 크기 | 보조변수 | SA 성공률 |
|:-:|:---------:|:-------:|:---------:|
| 3 | 4 | 1 | 76.60% |
| 4 | 6 | 2 | 43.20% |
| 5 | 10 | 5 | 16.30% |
| 6 | 18 | 12 | 8.00% |
| 7 | 22 | 15 | 5.10% |

M 보정 후 n=6~7에서도 소폭 성공 (8.0%, 5.1%). 정확 모드 실용 한계 = n ≤ 6~7 (greedy + 보정된 M).

### 정확 모드: 3-전략 비교 (10 runs)

#### 보조변수 수

| n | original | cache | greedy | greedy 절감률 |
|:-:|:--------:|:-----:|:------:|:------------:|
| 3 | 1 | 1 | 1 | 0.0% |
| 4 | 6 | 4 | 2 | 66.7% |
| 5 | 23 | 9 | 5 | 78.3% |
| 6 | 72 | 18 | 12 | 83.3% |
| 7 | 201 | 35 | 15 | 92.5% |
| 8 | 522 | 66 | 22 | 95.8% |

#### SA 성공률

| n | original | cache | greedy |
|:-:|:--------:|:-----:|:------:|
| 3 | 76.60% | 80.10% | 77.30% |
| 4 | 39.60% | 40.90% | **44.00%** |
| 5 | 10.70% | 8.00% | **18.30%** |
| 6 | 3.20% | 4.40% | **10.00%** |
| 7 | 3.00% | 1.40% | **3.50%** |
| 8 | 0.00% | 0.60% | **2.10%** |

### Greedy 확장 스케일링 (exact greedy vs approx, 10 runs)

| n | Greedy QUBO | Greedy Rate | Approx QUBO | Approx Rate |
|:-:|:-----------:|:-----------:|:-----------:|:-----------:|
| 3 | 4 | **76.60%** | 3 | 65.10% |
| 4 | 6 | 38.20% | 4 | **70.00%** |
| 5 | 10 | 23.20% | 5 | **68.00%** |
| 6 | 18 | 7.50% | 6 | **50.80%** |
| 7 | 22 | 3.30% | 7 | **59.30%** |
| 8 | 30 | 1.20% | 8 | **54.50%** |
| 9 | 46 | 0.70% | 9 | **50.20%** |
| 10 | 78 | 0.20% | 10 | **59.10%** |

교차점은 n=3~4. n=3에서만 greedy 우세, n≥4에서 approx가 안정적으로 50~70% 유지. M 보정 후 greedy도 n=8~10에서 소폭 성공 (0.2~1.2%).

### Sweep 전이 (n=8, greedy, 20 instances)

| Sweeps | GS Rate |
|:------:|:-------:|
| 50 | 0.65% |
| 100 | 1.00% |
| 200 | 1.05% |
| 500 | 1.75% |
| 1000 | 1.75% |
| 2000 | 1.95% |
| 5000 | 2.15% |
| 10000 | 2.10% |

M 보정 후 SA가 소폭 성공하지만, sweep 200배 증가(50→10000)에도 0.65% → 2.10%로 거의 포화. Rosenberg 패널티가 SA에 근본적으로 어려운 에너지 랜드스케이프를 생성.

### 다른 방법론과의 비교 (n=5)

```
성공률 기준 주의:
  - Truth Table (Exact/Approx): GS 찾은 read 수 / 전체 read 수 (1000 samples)
  - Wishart, ZeroExp, Posiform: best sample 기준 (10 runs 중 GS 찾은 run 수)
```

| 방법론 | SA 성공률 |
|--------|:---------:|
| Exact-Rosenberg | 22.90% |
| Wishart (α=0.7) | 40.00% |
| Approx-Valley | 54.00% |
| Approx-Random | 67.60% |
| ZeroExpectation | 100.00% |
| Posiform | 100.00% |

> 자세한 실험 결과: [`TRUTHTABLE_EXPERIMENT.md`](TRUTHTABLE_EXPERIMENT.md)

### 벤치마크 활용 방안

| 활용 | 설명 |
|------|------|
| **에너지 랜드스케이프 제어** | 모든 $2^n$개 비트스트링의 에너지를 직접 설계 가능 |
| **난이도 조절** | local minima 수, 장벽 높이 등을 파라미터로 제어 |
| **솔버 정밀 평가** | ground state뿐 아니라 전체 에너지 순위를 알고 있으므로 "근사 품질" 평가 가능 |
| **소규모 양자 프로세서** | n ≤ 23에서 에너지 지형을 정확히 알고 있는 벤치마크 |

## 한계

### 정확 모드

1. **보조변수 지수 폭발**: n=8 → greedy로도 30변수 (original: 530변수). SA 실용 한계 greedy로 n ≤ 7~8.
2. **변수 효율**: 30변수 QUBO로 8비트 문제만 표현. Posiform은 같은 크기로 30비트 표현.
3. **Hardness 출처가 인위적**: 난이도가 Rosenberg 패널티 구조에서 발생.
4. **Sweep 포화**: n=8 greedy에서 sweep 200배 증가에도 0.65% → 2.10%로 거의 포화.

### 근사 모드

1. **에너지 정확도 상실**: RMSE ~1.4 (lstsq 수렴) / ~2.9 (feasible fallback), 에너지 순서 보존율 ~50%. 전체 에너지 스펙트럼은 근사.
2. **SLSQP 생성 시간**: 최적화 버전 기준 n=20에서 ~50s, n=22에서 ~120s. n+1마다 약 2배 증가.
3. **진리표 크기 한계**: 입력이 $2^n$개이므로 $n \leq \sim 23$에서만 실용적.

### 공통

진리표 기반 접근 자체가 $n \leq \sim 23$으로 제한됨 ($2^n$ 열거). 구조화된 에너지 함수는 대부분 저차(≤2차)가 되어 이 방법론의 차별점이 약화됨.

## 파일 구성

| 파일 | 역할 |
|------|------|
| `qubo_truthtable.py` | 생성기 (정확: Möbius + Rosenberg / 근사: QP / 근사 최적화) |
| `test_truthtable.py` | SA 실험 (valley sweep, N-scaling, 비교, 전략, sweep) |
| `test_approx_comparison.py` | 근사 원본 vs 최적화 비교 검증 (Unit + E2E + 성능) |
| `results/` | 생성된 QUBO 파일 (edge-list CSV) |

## 참고 문헌

### Möbius 변환

1. **Rota, G.-C.** "On the foundations of combinatorial theory I. Theory of Möbius functions." *Zeitschrift für Wahrscheinlichkeitstheorie und verwandte Gebiete*, 2(4), 340-368, 1964.
   - 포함-배제 원리의 일반화, Möbius 역변환 이론.

### Rosenberg 차수축소

2. **Rosenberg, I. G.** "Reduction of bivalent maximization to the quadratic case." *Cahiers du Centre d'Etudes de Recherche Operationnelle*, 17, 71-74, 1975.
   - 고차 의사불 함수의 2차화. 보조변수 $y = x_i x_j$ 도입 + 패널티 항.

### Pseudo-Boolean 최적화

3. **Boros, E. & Hammer, P. L.** "Pseudo-Boolean optimization." *Discrete Applied Mathematics*, 123(1-3), 155-225, 2002.
   - 다중선형 다항식과 QUBO의 관계, quadratization 이론.

### Quadratization 서베이

4. **Dattani, N.** "Quadratization in discrete optimization and quantum mechanics." *arXiv:1901.04405*, 2019.
   - 70+ quadratization 기법 서베이. NTR-KZFD, PTR-BCR 등. [`papers/`](papers/Dattani2019_Quadratization_Survey.pdf)

5. **Boros, E. & Gruber, A.** "On quadratization of pseudo-Boolean functions." *arXiv:1404.6538*, 2014.
   - Term-wise quadratization 기법 기초. [`papers/`](papers/Boros_Gruber2014_Quadratization_PseudoBoolean.pdf)

6. **Anthony, M., Boros, E., Crama, Y. & Gruber, A.** "Quadratization of symmetric pseudo-Boolean functions." *Discrete Applied Mathematics*, 203, 1-12, 2016.
   - 대칭 의사불 함수의 보조변수 하한 증명. [`papers/`](papers/Anthony_Boros_Crama_Gruber2014_Symmetric_Quadratization.pdf)

### 근사 이차화

보조변수 없이 고차 다항식을 2차(QUBO)로 근사하는 방법론. 본 모듈의 근사 모드(최소제곱 + 제약 QP)와 관련된 핵심 수식이 포함된 논문들:

8. **Dragoi, S., Baiardi, A. & Egger, D. J.** "Approximate quadratization of high-order Hamiltonians for combinatorial quantum optimization." *arXiv:2505.04700*, 2025.
   - 고차 해밀토니안의 보조변수 없는 근사 2차화. [`papers/`](papers/Dragoi2025_Approximate_Quadratization_QAOA.pdf)
   - **Section III-A (p.3-4)**: 두 가지 근사 이차화 방법을 수식으로 정리:
     - Hypergraph clique expansion: 고차 항의 가중치를 l2-norm 최소화로 2차 그래프에 배분 (식 7-8)
     - Variational quadratic Hamiltonian: 2차 해밀토니안 파라미터를 원래 고차 에너지 기준으로 동시 최적화 (식 10-12)
   - 본 모듈의 근사 모드와 가장 유사: 둘 다 "고차→2차 투영 시 오차 제곱합 최소화" 원리

10. **Zheng, G. & Krikidis, I.** "Constrained higher-order binary optimization for wireless communications systems using Ising machines." *arXiv:2509.20092*, 2025.
    - Taylor 전개 + augmented Lagrangian으로 고차→2차 근사. [`papers/`](papers/Zheng_Krikidis2025_Constrained_HOBO.pdf)
    - **Section III-A (p.5)**: 고차 다항식 f(x)를 현재 점 x0에서 2차 Taylor 전개로 근사하는 반복 알고리즘 (식 16-17, Algorithm 2). 보조변수 없이 고차→2차 근사하는 수식을 가장 직관적으로 설명

7. **Gabor, T., Rosenfeld, M. L., Linnhoff-Popien, C. & Feld, S.** "How to approximate any objective function via quadratic unconstrained binary optimization." *IEEE Q-SANER*, 2022. arXiv:2204.11035.
   - 임의 목적함수 → 다항식 근사(Taylor/Lagrange/Fourier) → QUBO 변환 툴킷. [`papers/`](papers/Gabor2022_Approximate_Any_Objective_QUBO.pdf)
   - **Section II-A**: 6가지 함수→다항식 변환 기법 (Lagrange, Spline, 단순화, 단조함수, Fourier, Taylor)
   - **Section II-B**: 다항식→QUBO 변환 절차 (이진 인코딩 + Rosenberg 차수축소)
   - 근사 이차화 수식 자체보다는 전체 변환 파이프라인의 범용 툴킷에 초점

9. **Nakada, H. & Tanaka, S.** "Systematic and efficient construction of QUBO forms for high-order and dense interactions." *arXiv:2506.08448*, 2025.
   - ReLU 기저 기반 범용 근사→QUBO 파이프라인. [`papers/`](papers/Nakada_Tanaka2025_ReLU_QUBO.pdf)
   - **Section 3**: ReLU 함수의 Legendre 변환을 이용한 2차화 (식 2-3) + 이산화 방법 (식 4-5). ML 회귀 모델을 QUBO로 변환하는 체계적 방법론

### 수치 해석 및 최적화

13. **Higham, N. J.** *Accuracy and Stability of Numerical Algorithms*. 2nd ed., SIAM, 2002.
    - 부동소수점 오차 분석의 표준 교과서. **Theorem 3.1** (p.63): $n$차원 내적의 전방향 오차 상한 $|\text{fl}(\mathbf{x}^T \mathbf{y}) - \mathbf{x}^T \mathbf{y}| \leq \gamma_n \sum |x_i y_i|$. Feasible fallback의 gap 부호 안전성 증명에 사용.

14. **Boyd, S. & Vandenberghe, L.** *Convex Optimization*. Cambridge University Press, 2004.
    - 볼록 QP 이론 (Chapter 4.4) 및 절단 평면법 (Chapter 6.4, Localization methods). $A^T A$가 양정치이면 목적함수가 순볼록 → 전역 최적해 유일. 근사 모드의 SLSQP cutting plane 수렴 근거.

15. **IEEE 754-2019.** "IEEE Standard for Floating-Point Arithmetic." IEEE, 2019.
    - 배정밀도 부동소수점 연산의 결정성 규정. 동일 입력·동일 연산 → 동일 출력 보장 (§11, Reproducibility). 경로 1–2의 `_batched_energies` 일관성 근거.

16. **Kraft, D.** "A Software Package for Sequential Quadratic Programming." DFVLR-FB 88-28, 1988.
    - `scipy.optimize.minimize(method='SLSQP')`의 기반 알고리즘. 등식·부등식 제약 하의 비선형 최적화. 제약 만족 정밀도 ~$10^{-8}$ 특성의 출처.

### 관련 QUBO 구성

11. **Mandal, A., Roy, A., Upadhyay, S. & Ushijima-Mwesigwa, H.** "Compressed quadratization of higher order binary optimization problems." *arXiv:2001.00658*, 2020.
    - Ising 공간 차수축소 시 보조변수 2개 필요 증명. [`papers/`](papers/Mandal2020_Compressed_Quadratization.pdf)

12. **Verma, A. & Lewis, M.** "Goal seeking quadratic unconstrained binary optimization." *arXiv:2103.12951*, 2021.
    - 목표 에너지값에 근접하는 이진 벡터 탐색. [`papers/`](papers/Verma_Lewis2021_Goal_Seeking_QUBO.pdf)

## 사용법

```bash
# 정확 모드 (진리표 직접 입력)
python3 truthtable/qubo_truthtable.py '{"000":3,"001":4,"010":4,"011":5,"100":3,"101":5,"110":2,"111":1}'

# 근사 모드 (--approx)
python3 truthtable/qubo_truthtable.py '{"000":3,"001":4,"010":4,"011":5,"100":3,"101":5,"110":2,"111":1}' --approx

# 프리셋: Random Landscape (n=10, target=1011001100)
python3 truthtable/qubo_truthtable.py --preset random 10 1011001100 --approx --seed 42

# 프리셋: Multi-Valley (n=6, 2개 계곡)
python3 truthtable/qubo_truthtable.py --preset valley 6 101010,010101 --approx --seed 42

# SA 실험: Multi-Valley Sweep
python3 truthtable/test_truthtable.py --valley 10

# SA 실험: N-Scaling
python3 truthtable/test_truthtable.py --scaling 10

# SA 실험: 7-way 비교
python3 truthtable/test_truthtable.py --compare 10

# SA 실험: 차수축소 전략 비교 (original/cache/greedy)
python3 truthtable/test_truthtable.py --strategy 10
python3 truthtable/test_truthtable.py --strategy 10 3,4,5,6

# SA 실험: Greedy 확장 스케일링 (exact greedy vs approx)
python3 truthtable/test_truthtable.py --greedy-scaling 10
python3 truthtable/test_truthtable.py --greedy-scaling 10 3,4,5,6,7,8

# SA 실험: Sweep 전이 (S-curve, n=8)
python3 truthtable/test_truthtable.py --sweep 20

# 전략 지정 생성
python3 truthtable/qubo_truthtable.py --preset random 8 10110011 --strategy greedy
python3 truthtable/qubo_truthtable.py --preset random 8 10110011 --strategy original
```
