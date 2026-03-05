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

보조변수 없이 $n$-변수 QUBO를 직접 구한다.

$$\min_Q \sum_x \left(E_Q(x) - E_{\text{truth}}(x)\right)^2 \quad \text{s.t.} \quad E_Q(\text{target}) + \varepsilon \leq E_Q(x) \;\; \forall x \neq \text{target}$$

- **자유 파라미터**: $n(n+1)/2$개 (Q 행렬 upper triangular)
- **제약조건**: $2^n - 1$개 (target이 모든 다른 상태보다 낮은 에너지)
- **풀이**: 무제약 최소제곱 → 위반 시 SLSQP (iterative cutting plane)

구현: `create_qubo_approx()`

**핵심 장점**: 보조변수 0개 → QUBO 크기 = n. n=16까지 실용적.

### 왜 두 모드가 필요한가?

|                    | 정확 original | 정확 greedy | 근사 (QP)         |
|--------------------|:------------:|:-----------:|:-----------------:|
| 보조변수            | $O(n \cdot 2^n)$ | **78~96% 절감** | **0개**       |
| QUBO 크기 (n=8)    | 530변수      | **30변수**   | **8변수**         |
| 에너지 정확도       | 100% 일치    | 100% 일치    | 근사 (RMSE ~0.7)  |
| Ground state       | 보장          | 보장         | **보장** (QP 제약) |
| SA 성공률 (n=5)    | 10.70%       | 18.30%      | **68.00%**        |
| 실용 한계           | n ≤ 4~5     | **n ≤ 7~8** | **n ≤ 16**        |

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

근사 모드:
  진리표 → Feature matrix 구성 → 최소제곱 → 제약 위반 시 SLSQP → Q dict 생성 → 검증
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
| **소규모 양자 프로세서** | n ≤ 16에서 에너지 지형을 정확히 알고 있는 벤치마크 |

## 한계

### 정확 모드

1. **보조변수 지수 폭발**: n=8 → greedy로도 30변수 (original: 530변수). SA 실용 한계 greedy로 n ≤ 7~8.
2. **변수 효율**: 30변수 QUBO로 8비트 문제만 표현. Posiform은 같은 크기로 30비트 표현.
3. **Hardness 출처가 인위적**: 난이도가 Rosenberg 패널티 구조에서 발생.
4. **Sweep 포화**: n=8 greedy에서 sweep 200배 증가에도 0.65% → 2.10%로 거의 포화.

### 근사 모드

1. **에너지 정확도 상실**: RMSE ~0.7, 에너지 순서 보존율 ~55%. 전체 에너지 스펙트럼은 근사.
2. **SLSQP 생성 시간**: n=16에서 ~20분. n+1마다 약 3~4배 증가. n=17+는 진리표 열거($2^n$) 자체가 병목.
3. **진리표 크기 한계**: 입력이 $2^n$개이므로 $n \leq \sim 20$에서만 진리표 열거 가능.

### 공통

진리표 기반 접근 자체가 $n \leq \sim 20$으로 제한됨 ($2^n$ 열거). 구조화된 에너지 함수는 대부분 저차(≤2차)가 되어 이 방법론의 차별점이 약화됨.

## 파일 구성

| 파일 | 역할 |
|------|------|
| `qubo_truthtable.py` | 생성기 (정확: Möbius + Rosenberg / 근사: QP) |
| `test_truthtable.py` | SA 실험 (valley sweep, N-scaling, 비교, 전략, sweep) |
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

7. **Dragoi, S., Baiardi, A. & Egger, D. J.** "Approximate quadratization of high-order Hamiltonians for combinatorial quantum optimization." *arXiv:2505.04700*, 2025.
   - 고차 해밀토니안의 보조변수 없는 근사 2차화. [`papers/`](papers/Dragoi2025_Approximate_Quadratization_QAOA.pdf)

8. **Nakada, H. & Tanaka, S.** "Systematic and efficient construction of QUBO forms for high-order and dense interactions." *arXiv:2506.08448*, 2025.
   - ReLU 기저 기반 범용 근사→QUBO 파이프라인. [`papers/`](papers/Nakada_Tanaka2025_ReLU_QUBO.pdf)

9. **Zheng, G. & Krikidis, I.** "Constrained higher-order binary optimization for wireless communications systems using Ising machines." *arXiv:2509.20092*, 2025.
   - Taylor 전개 + augmented Lagrangian으로 고차→2차 근사. [`papers/`](papers/Zheng_Krikidis2025_Constrained_HOBO.pdf)

### 관련 QUBO 구성

10. **Mandal, A., Roy, A., Upadhyay, S. & Ushijima-Mwesigwa, H.** "Compressed quadratization of higher order binary optimization problems." *arXiv:2001.00658*, 2020.
    - Ising 공간 차수축소 시 보조변수 2개 필요 증명. [`papers/`](papers/Mandal2020_Compressed_Quadratization.pdf)

11. **Verma, A. & Lewis, M.** "Goal seeking quadratic unconstrained binary optimization." *arXiv:2103.12951*, 2021.
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
