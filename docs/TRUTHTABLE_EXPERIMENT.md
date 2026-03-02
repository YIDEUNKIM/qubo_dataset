# README: Truth Table QUBO

## 개요

임의의 진리표(모든 비트스트링 → 에너지 매핑)로부터 **수학적으로 정확한 QUBO**를 구성한다.

```
입력:  {000→3, 001→4, 010→4, 011→5, 100→3, 101→5, 110→2, 111→1}
출력:  Q 행렬 (QUBO), ground state = 111 (energy=1) 보장
```

## 수학적 기반

### Step 1: Möbius 변환 (진리표 → 다중선형 다항식)

임의의 함수 f: {0,1}^n → R은 유일한 다중선형 다항식으로 표현 가능:

```
f(x) = Σ_{S ⊆ [n]} c_S · Π_{i∈S} x_i
```

계수 c_S는 포함-배제 원리(Möbius 역변환)로 계산:

```
c_S = Σ_{T ⊆ S} (-1)^{|S|-|T|} · f(1_T)
```

Fast Möbius 알고리즘으로 O(n · 2^n)에 계산. 구현: `mobius_transform()`

**n=3 예시 검증:**

```
f(x) = 3 + x₂ + x₃ + x₁x₃ - 2x₁x₂ - 3x₁x₂x₃
         ↑상수   ↑선형      ↑이차         ↑고차(3차)
```

### Step 2: 항 분류

| 차수 | 처리 | QUBO 배치 |
|------|------|-----------|
| 0 (상수) | 에너지 오프셋으로 기록 | Q에 포함 안 함 |
| 1 (선형) | c_i · x_i → Q[i,i] += c_i | 대각 |
| 2 (이차) | c_ij · x_i x_j → Q[i,j] += c_ij | 비대각 |
| ≥3 (고차) | **Rosenberg 차수축소 필요** | 보조변수 도입 |

구현: `classify_terms()`

### Step 3: Rosenberg 차수축소

k차 단항식 c · x_{i1} · x_{i2} · ... · x_{ik}에 대해:

1. 보조변수 y = x_{i1} · x_{i2} 도입
2. 패널티 M · (x_{i1}·x_{i2} - 2·x_{i1}·y - 2·x_{i2}·y + 3·y) 추가
   - y = x_{i1}·x_{i2}일 때 패널티 = 0
   - y ≠ x_{i1}·x_{i2}일 때 패널티 ≥ M
3. 원래 항을 c · y · x_{i3} · ... · x_{ik}로 치환
4. 여전히 차수 > 2이면 반복 (체이닝)

**필요 보조변수**: k차 단항식 1개당 k-2개

구현: `rosenberg_reduce()`

### Step 4: 패널티 강도 결정

```
M = max(truth_table) - min(truth_table) + 1
```

보조변수 위반 시 에너지가 진리표 최대값을 초과 → ground state 보장.

구현: `compute_penalty_strength()`

### Step 5: QUBO 조립

원래 2차 이하 항 + 축소된 이차항 + Rosenberg 패널티 → Q dict 생성.

구현: `assemble_qubo()`

### Step 6: 전수검사 검증

모든 2^n 비트스트링에 대해 (원래 변수 + 최적 보조변수)의 QUBO 에너지가 진리표 에너지와 일치하는지 확인.

구현: `verify_qubo()`

## 에너지 함수 프리셋

### Energy Gap (프리셋 4)

```python
E(target) = 0
E(x ≠ target) = gap + |N(0, noise_scale)|
```

gap이 작을수록 난이도 ↑. 양자 어닐링의 minimum spectral gap과 직결.

구현: `preset_energy_gap(n, target, gap, noise_scale, seed)`

### Multi-Valley (프리셋 2)

```python
targets[0]: global minimum (energy = 0)
targets[1:]: local minima (energy = gap)
나머지: Hamming 거리 비례 에너지 + 장벽
```

local minima 수 증가 → SA가 계곡에 갇힘. 양자 tunneling 벤치마크.

구현: `preset_multi_valley(n, targets, gap, barrier_height, seed)`

## 사용법

```bash
# JSON 진리표 직접 입력
python truthtable/qubo_truthtable.py '{"000":3,"001":4,"010":4,"011":5,"100":3,"101":5,"110":2,"111":1}'

# Energy Gap 프리셋
python truthtable/qubo_truthtable.py --preset gap 5 10110 --gap 2.0 --seed 42

# Multi-Valley 프리셋
python truthtable/qubo_truthtable.py --preset valley 5 10110,01001 --barrier 5.0 --seed 42

# SA 실험
python truthtable/test_truthtable.py --gap 10        # gap sweep
python truthtable/test_truthtable.py --valley 10     # valley sweep
python truthtable/test_truthtable.py --scaling 10    # N-scaling
python truthtable/test_truthtable.py --compare 10    # 5-way 비교
```

## 파일 구조

```
truthtable/
├── qubo_truthtable.py   ← 생성기 (Möbius + Rosenberg + 프리셋)
├── test_truthtable.py   ← SA 실험 (gap/valley/scaling/비교)
└── results/             ← 생성된 QUBO (edge-list CSV)
```

## 파이프라인 복잡도

| Step | 내용 | 복잡도 |
|------|------|--------|
| 1 | Möbius 변환 | O(n · 2^n) |
| 2 | 항 분류 | O(2^n) |
| 3 | Rosenberg 차수축소 | O(Σ C(n,k)(k-2)) |
| 4 | 패널티 강도 계산 | O(2^n) |
| 5 | QUBO 조립 | O(n + aux) |
| 6 | 전수검사 검증 | O(2^n) |

---

## 실험 결과

### 실험 1: Energy Gap Sweep (n=5)

gap 크기별 SA 성공률. gap ↓ = 난이도 ↑.

```
설정: n=5, QUBO=28변수, instances=10, reads=100, sweeps=5000

Gap    | GS Rate  | Avg Hamming
-------|----------|------------
0.1    |   6.60%  |  2.65
0.5    |   5.80%  |  2.53
1.0    |   9.90%  |  2.40
2.0    |  12.70%  |  2.46
5.0    |  17.20%  |  2.42
10.0   |  21.40%  |  2.03
```

**관찰**: gap ↑ → SA 성공률 ↑ 확인. 난이도 제어 파라미터로 유효.

### 실험 2: Multi-Valley Sweep (n=5)

local minima 수 증가에 따른 SA 성공률 변화.

```
설정: n=5, QUBO=28변수, instances=10, reads=100, sweeps=5000

Valleys | GS Rate  | Avg Hamming
--------|----------|------------
1       |  17.30%  |  1.69
2       |  18.20%  |  1.81
3       |  16.00%  |  2.05
4       |  14.10%  |  2.16
```

**관찰**: valley ↑ → SA 성공률 소폭 ↓, Hamming ↑. n=5에서는 차이 미미하나 경향 존재.

### 실험 3: N-Scaling

n 증가에 따른 보조변수 폭발 + SA 성능 급락.

```
설정: gap=2.0, runs=10, reads=100, sweeps=5000

n  | QUBO 크기 | 보조변수 | SA 성공률
---|-----------|---------|----------
3  |         4 |       1 |   98.9%
4  |        10 |       6 |   51.6%
5  |        28 |      23 |   13.6%
6  |        78 |      72 |    0.6%
7  |       208 |     201 |    0.0%
```

**QUBO 크기 이론값** (worst case: 모든 고차 계수 ≠ 0):

| n | 고차항 수 | 보조변수 (최대) | 총 QUBO 크기 |
|---|----------|---------------|-------------|
| 3 | 1 | 1 | 4 |
| 5 | 16 | 26 | ~31 |
| 8 | 219 | 502 | ~510 |
| 10 | 968 | 3,084 | ~3,094 |
| 12 | 4,017 | ~14,000 | ~14,012 |
| 15 | ~32,600 | ~130,000 | ~130,015 |

### 실험 4: 5-way 비교 (n=5)

동일 n에서 기존 방법론과 SA 성공률 비교.

```
설정: n=5, runs=10, reads=100, sweeps=5000

방법론              | SA 성공률
--------------------|----------
TruthTable-Gap      |  13.20%
TruthTable-Valley   |  20.80%
Wishart (α=0.7)     |  60.00%
ZeroExpectation     | 100.00%
Posiform            | 100.00%
```

**관찰**: 같은 n=5에서 Truth Table QUBO가 가장 어려움. 단, hardness의 원인은 에너지 랜드스케이프가 아니라 보조변수 구조.

---

## 한계

### 1. 보조변수 지수 폭발 (핵심 한계)

n개 원래 변수의 진리표 표현에 최대 Σ C(n,k)(k-2) ≈ O(n · 2^n)개 보조변수 필요.

```
n=5  → 원래 5변수,  QUBO 28변수    (5.6배)
n=8  → 원래 8변수,  QUBO 530변수   (66배)
n=10 → 원래 10변수, QUBO 3094변수  (309배)
```

SA 실용 한계: **n ≤ 5~6** (QUBO ~80변수)

### 2. 변수 효율 문제

28변수 QUBO로 얻는 정보는 5비트 ground state뿐. 동일 QUBO 크기에서 다른 방법론은 훨씬 큰 문제를 표현 가능.

| 방법론 | 28변수 QUBO → 원래 문제 크기 |
|--------|---------------------------|
| Truth Table | n=5 (보조변수 23개) |
| Posiform | n=28 (보조변수 없음) |
| Wishart | n=28 (보조변수 없음) |
| Hardened Posiform | n=28 (보조변수 없음) |

### 3. Hardness 출처가 인위적

난이도는 에너지 랜드스케이프 설계가 아니라 Rosenberg 패널티 구조에서 발생. 패널티 M이 실제 에너지 스케일보다 크기 때문에 솔버는 "패널티 위반 않는 영역"을 먼저 탐색해야 함. 의도한 "에너지 갭에 의한 난이도 제어"와는 다른 종류의 hardness.

### 4. 진리표 확장성 한계

2^n개 에너지를 모두 지정해야 하므로 입력 자체가 n ≤ ~20에서만 가능. 구조화된 에너지 함수(예: Hamming 거리 기반)는 대부분 저차(≤2차)가 되어 보조변수 불필요 → 이 방법론의 존재 의의 소멸.

### 5. 에너지 갭 제어 효과가 보조변수 noise에 묻힘

gap=0.1 vs gap=10.0의 SA 성공률 차이가 6.6% vs 21.4%에 불과. 난이도 대부분이 보조변수 구조에서 오기 때문에 gap 파라미터의 실질적 제어력 제한적.

## 결론

Truth Table QUBO는 **수학적으로 완전히 정확**하고 (전수검사 100% 통과) **유효한 QUBO**로서 모든 솔버에 적용 가능.

그러나 벤치마크로서의 실용적 가치는 제한적:
- 변수 효율이 나쁨 (5비트 문제에 28변수 QUBO)
- 난이도 제어가 보조변수 구조에 묻힘
- 기존 방법론 대비 같은 QUBO 크기에서 표현 가능한 문제 규모가 작음

**잠재적 활용**: 소규모(n ≤ 5~6) 양자 프로세서의 정밀 벤치마크에서, "모든 에너지 준위를 알고 있는" 특성을 살려 솔버의 근사 품질(approximation ratio)을 정밀 평가하는 용도.
