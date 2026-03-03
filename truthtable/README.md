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

$k$차 단항식 1개당 보조변수 $k-2$개 필요. 구현: `rosenberg_reduce()`

#### Step 4: 패널티 강도

```
M = max(truth_table) - min(truth_table) + 1
```

보조변수 위반 시 에너지가 반드시 진리표 최대값을 초과 → ground state 보장. 구현: `compute_penalty_strength()`

### 근사 모드: 제약 최소제곱 (QP)

보조변수 없이 $n$-변수 QUBO를 직접 구한다.

$$\min_Q \sum_x \left(E_Q(x) - E_{\text{truth}}(x)\right)^2 \quad \text{s.t.} \quad E_Q(\text{target}) + \varepsilon \leq E_Q(x) \;\; \forall x \neq \text{target}$$

- **자유 파라미터**: $n(n+1)/2$개 (Q 행렬 upper triangular)
- **제약조건**: $2^n - 1$개 (target이 모든 다른 상태보다 낮은 에너지)
- **풀이**: 무제약 최소제곱 → 위반 시 SLSQP (iterative cutting plane)

구현: `create_qubo_approx()`

**Iterative Cutting Plane**: $n \geq 12$에서 $2^n - 1$개 제약 전부를 넣으면 SLSQP가 느려짐. 매 반복마다 위반/근접 제약만 active set에 추가하여 효율적으로 풀이.

**핵심 장점**: 보조변수 0개 → QUBO 크기 = n. n=14까지 실용적.

### 왜 두 모드가 필요한가?

|                    | 정확 모드 (Rosenberg) | 근사 모드 (QP)         |
|--------------------|---------------------|----------------------|
| 보조변수            | $O(n \cdot 2^n)$ 폭발 | **0개**               |
| QUBO 크기 (n=8)    | 530변수             | **8변수**             |
| 에너지 정확도       | 100% 일치           | 근사 (RMSE ~0.7)      |
| Ground state       | 보장                | **보장** (QP 제약)     |
| SA 성공률 (n=8)    | 0%                  | **100%**              |
| 실용 한계           | n ≤ 5~6            | **n ≤ 14+**           |
| 생성 시간 (n=14)   | 불가                | ~100s (SLSQP 병목)    |
| 에너지 순서 보존율   | 100%               | ~55%                  |

정확 모드는 **수학적 정확성**이 요구되는 소규모 검증에, 근사 모드는 **실용적 벤치마크 생성**에 적합.

## 에너지 함수 프리셋

### Energy Gap

```python
E(target) = 0
E(x != target) = gap + |N(0, noise_scale)|
```

gap ↓ = 난이도 ↑. 양자 어닐링의 minimum spectral gap과 직결. 구현: `preset_energy_gap()`

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
| `gap` | 2.0 | Energy Gap 프리셋의 갭 크기 |
| `noise_scale` | 1.0 | Energy Gap 프리셋의 노이즈 스케일 |
| `barrier_height` | 5.0 | Multi-Valley 프리셋의 장벽 높이 |
| `seed` | None | 재현성을 위한 난수 시드 |
| `verbose` | True | 진행 상황 출력 여부 |

### 주요 함수

| 함수 | 설명 |
|------|------|
| `mobius_transform(truth_table, n)` | 진리표 → 다중선형 다항식 계수 (fast Möbius, $O(n \cdot 2^n)$) |
| `classify_terms(coefficients)` | 계수를 상수/선형/이차/고차로 분류 |
| `rosenberg_reduce(higher_order, n)` | 고차항 → 보조변수 + 패널티 → 2차화 (체이닝) |
| `compute_penalty_strength(truth_table, n)` | 최소 패널티 강도 M 계산 |
| `assemble_qubo(...)` | 선형 + 이차 + 축소된 항 + 패널티 → Q dict 조립 |
| `compute_aux_values(x_orig, aux_info)` | 원래 변수 → 올바른 보조변수 값 계산 |
| `verify_qubo(Q, truth_table, n, aux_info, offset)` | 모든 $2^n$ 비트스트링 전수검사 |
| `preset_energy_gap(n, target, gap, ...)` | Energy Gap 에너지 함수 생성 |
| `preset_multi_valley(n, targets, gap, ...)` | Multi-Valley 에너지 함수 생성 |
| `create_qubo_approx(truth_table, n, epsilon)` | **근사 모드 진입점** |
| `create_qubo_truthtable(truth_table, n, seed)` | **정확 모드 진입점** |

### 반환값

```python
# 정확 모드
Q, info = create_qubo_truthtable(truth_table)
# Q: QUBO 딕셔너리 {(i,j): weight}
# info: {'n_original', 'n_aux', 'n_total', 'offset', 'penalty_M',
#        'ground_state', 'n_higher_order', 'aux_info', 'target'}

# 근사 모드
Q, info = create_qubo_approx(truth_table)
# Q: QUBO 딕셔너리 {(i,j): weight}
# info: {'n_original', 'n_aux'(=0), 'n_total', 'offset', 'rmse',
#        'max_error', 'energy_gap', 'order_preservation', 'gs_verified',
#        'ground_state', 'target', 'aux_info'}
```

### 확장성 분석

| n | 진리표 크기 | 정확 모드 QUBO 크기 | 정확 모드 보조변수 | 근사 모드 QUBO 크기 |
|---|------------|:------------------:|:-----------------:|:------------------:|
| 3 | 8 | 4 | 1 | 3 |
| 5 | 32 | 28 | 23 | 5 |
| 8 | 256 | ~530 | ~522 | 8 |
| 10 | 1,024 | ~3,094 | ~3,084 | 10 |
| 14 | 16,384 | 불가 | 불가 | 14 |

## SA 난이도 특성

### 공통 SA 설정

```
솔버: neal.SimulatedAnnealingSampler (D-Wave)
num_reads = 100 (인스턴스당 SA 샘플 수)
num_sweeps = 5000
instances = 10 (각 설정마다 랜덤 target으로 10회 반복)
성공률 = GS 찾은 read 수 / 전체 read 수 (10 × 100 = 1000 samples)
```

### 정확 모드: Energy Gap Sweep

n=5, 프리셋: `preset_energy_gap(gap=각 값, noise_scale=1.0)`

| Gap | GS Rate | Avg Hamming | QUBO 크기 |
|:---:|:-------:|:-----------:|:---------:|
| 0.1 | 6.60% | 2.65 | 28변수 |
| 0.5 | 5.80% | 2.53 | 28변수 |
| 1.0 | 9.90% | 2.40 | 28변수 |
| 2.0 | 12.70% | 2.46 | 28변수 |
| 5.0 | 17.20% | 2.42 | 28변수 |
| 10.0 | 21.40% | 2.03 | 28변수 |

정확 모드의 낮은 성공률은 에너지 갭이 아니라 **Rosenberg 보조변수 폭발**로 인한 탐색 공간 증가가 원인.

### 정확 모드: N-Scaling

sizes=[3,4,5,6,7], 프리셋: `preset_energy_gap(gap=2.0, noise_scale=1.0)`

| n | QUBO 크기 | 보조변수 | SA 성공률 |
|:-:|:---------:|:-------:|:---------:|
| 3 | 4 | 1 | 98.9% |
| 4 | 10 | 6 | 51.6% |
| 5 | 28 | 23 | 13.6% |
| 6 | 78 | 72 | 0.6% |
| 7 | 208 | 201 | 0.0% |

n=7 이상에서 SA 성공률 0%. QUBO 크기가 기하급수적으로 폭발.

### 근사 모드: Energy Gap Sweep

n=5, 프리셋: `preset_energy_gap(gap=각 값, noise_scale=1.0)`

| Gap | GS Rate | Avg Hamming | QUBO 크기 |
|:---:|:-------:|:-----------:|:---------:|
| 0.1 | 76.30% | 0.55 | 5변수 |
| 0.5 | 68.70% | 0.95 | 5변수 |
| 1.0 | 65.50% | 0.93 | 5변수 |
| 2.0 | 65.80% | 0.84 | 5변수 |
| 5.0 | 65.70% | 1.00 | 5변수 |
| 10.0 | 53.90% | 1.15 | 5변수 |

정확 모드 대비 **5~12배 개선**. QUBO 크기 = n으로 탐색 공간이 작음.

### 근사 모드: N-Scaling

프리셋: `preset_energy_gap(gap=2.0, noise_scale=1.0)`, 단일 인스턴스

| n | QUBO 크기 | SA 성공률 | 생성 시간 | GS 보장 |
|:-:|:---------:|:---------:|:---------:|:-------:|
| 3 | 3 | 71.9% | <0.1s | YES |
| 5 | 5 | 53.4% | <0.1s | YES |
| 8 | 8 | 100% | 0.2s | YES |
| 10 | 10 | 67% | 0.3s | YES |
| 12 | 12 | 55% | 11.7s | YES |
| 14 | 14 | 52% | 103s | YES |

보조변수 없이 n=14까지 실용적. SA 성공률은 n이 커져도 50% 이상 유지.

### 다른 방법론과의 비교 (n=5)

```
성공률 기준 주의:
  - Truth Table (Exact/Approx): GS 찾은 read 수 / 전체 read 수 (1000 samples)
  - Wishart, ZeroExp, Posiform: best sample 기준 (10 runs 중 GS 찾은 run 수)
```

| 방법론 | SA 성공률 |
|--------|:---------:|
| Exact-Rosenberg | 13.80% |
| Approx-Gap | 50.00% |
| Wishart (α=0.7) | 50.00% |
| Approx-Valley | 56.10% |
| ZeroExpectation | 100.00% |
| Posiform | 100.00% |

> 자세한 실험 결과: [`docs/TRUTHTABLE_EXPERIMENT.md`](../docs/TRUTHTABLE_EXPERIMENT.md)

### 벤치마크 활용 방안

| 활용 | 설명 |
|------|------|
| **에너지 랜드스케이프 제어** | 모든 $2^n$개 비트스트링의 에너지를 직접 설계 가능 |
| **난이도 조절** | gap, local minima 수, 장벽 높이 등을 파라미터로 제어 |
| **솔버 정밀 평가** | ground state뿐 아니라 전체 에너지 순위를 알고 있으므로 "근사 품질" 평가 가능 |
| **소규모 양자 프로세서** | n ≤ 14에서 에너지 지형을 정확히 알고 있는 벤치마크 |

## 한계

### 정확 모드

1. **보조변수 지수 폭발**: n=8 → QUBO 530변수. SA 실용 한계 n ≤ 5~6.
2. **변수 효율**: 28변수 QUBO로 5비트 문제만 표현. Posiform은 같은 크기로 28비트 표현.
3. **Hardness 출처가 인위적**: 난이도가 Rosenberg 패널티 구조에서 발생. 에너지 갭 제어 효과가 묻힘.

### 근사 모드

1. **에너지 정확도 상실**: RMSE ~0.7, 에너지 순서 보존율 ~55%. 전체 에너지 스펙트럼은 근사.
2. **SLSQP 생성 시간**: n=14에서 ~100초. n=15+ 는 진리표 열거($2^n$) 자체가 병목.
3. **진리표 크기 한계**: 입력이 $2^n$개이므로 $n \leq \sim 20$에서만 진리표 열거 가능.

### 공통

진리표 기반 접근 자체가 $n \leq \sim 20$으로 제한됨 ($2^n$ 열거). 구조화된 에너지 함수는 대부분 저차(≤2차)가 되어 이 방법론의 차별점이 약화됨.

## 파일 구성

| 파일 | 역할 |
|------|------|
| `qubo_truthtable.py` | 생성기 (정확: Möbius + Rosenberg / 근사: QP) |
| `test_truthtable.py` | SA 실험 (gap sweep, valley sweep, N-scaling, 7-way 비교) |
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
   - Term-wise quadratization 기법 기초. 다중 분할 및 공통 부분 기반 집약적 접근. [`papers/`](papers/Boros_Gruber2014_Quadratization_PseudoBoolean.pdf)

6. **Anthony, M., Boros, E., Crama, Y. & Gruber, A.** "Quadratization of symmetric pseudo-Boolean functions." *Discrete Applied Mathematics*, 203, 1-12, 2016. (*arXiv:1404.6535*)
   - 대칭 의사불 함수의 보조변수 하한 증명. 정확 2차화의 한계 동기. [`papers/`](papers/Anthony_Boros_Crama_Gruber2014_Symmetric_Quadratization.pdf)

### 근사 이차화 (보조변수 없는 접근)

7. **Dragoi, S., Baiardi, A. & Egger, D. J.** "Approximate quadratization of high-order Hamiltonians for combinatorial quantum optimization." *arXiv:2505.04700*, 2025.
   - 고차 해밀토니안의 보조변수 없는 근사 2차화 (QAOA 대상). 본 프로젝트와 가장 유사한 접근. [`papers/`](papers/Dragoi2025_Approximate_Quadratization_QAOA.pdf)

8. **Nakada, H. & Tanaka, S.** "Systematic and efficient construction of QUBO forms for high-order and dense interactions." *arXiv:2506.08448*, 2025.
   - ReLU 기저 기반 범용 근사→QUBO 파이프라인. 진리표 최소제곱과 대안적 접근. [`papers/`](papers/Nakada_Tanaka2025_ReLU_QUBO.pdf)

9. **Zheng, G. & Krikidis, I.** "Constrained higher-order binary optimization for wireless communications systems using Ising machines." *arXiv:2509.20092*, 2025.
   - Taylor 전개 + augmented Lagrangian으로 고차→2차 근사. 반복적 Lagrangian 구조가 본 프로젝트의 cutting plane + SLSQP와 유사. [`papers/`](papers/Zheng_Krikidis2025_Constrained_HOBO.pdf)

### 관련 QUBO 구성

10. **Mandal, A., Roy, A., Upadhyay, S. & Ushijima-Mwesigwa, H.** "Compressed quadratization of higher order binary optimization problems." *arXiv:2001.00658*, 2020.
    - Ising 공간 차수축소 시 보조변수 2개 필요 증명. 압축 표현으로 보조변수 절감. [`papers/`](papers/Mandal2020_Compressed_Quadratization.pdf)

11. **Verma, A. & Lewis, M.** "Goal seeking quadratic unconstrained binary optimization." *arXiv:2103.12951*, 2021.
    - 목표 에너지값에 근접하는 이진 벡터 탐색. 최소제곱 편차 프레임워크가 본 프로젝트와 개념적으로 유사. [`papers/`](papers/Verma_Lewis2021_Goal_Seeking_QUBO.pdf)

## 사용법

```bash
# 정확 모드 (진리표 직접 입력)
python3 truthtable/qubo_truthtable.py '{"000":3,"001":4,"010":4,"011":5,"100":3,"101":5,"110":2,"111":1}'

# 근사 모드 (--approx)
python3 truthtable/qubo_truthtable.py '{"000":3,"001":4,"010":4,"011":5,"100":3,"101":5,"110":2,"111":1}' --approx

# 프리셋: Energy Gap (n=10, target=1011001100, gap=2.0)
python3 truthtable/qubo_truthtable.py --preset gap 10 1011001100 --approx --seed 42

# 프리셋: Multi-Valley (n=6, 2개 계곡)
python3 truthtable/qubo_truthtable.py --preset valley 6 101010,010101 --approx --seed 42

# SA 실험: Energy Gap Sweep
python3 truthtable/test_truthtable.py --gap 10

# SA 실험: Multi-Valley Sweep
python3 truthtable/test_truthtable.py --valley 10

# SA 실험: N-Scaling
python3 truthtable/test_truthtable.py --scaling 10

# SA 실험: 7-way 비교
python3 truthtable/test_truthtable.py --compare 10
```
