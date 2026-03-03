# README: Truth Table QUBO

## 개요

임의의 진리표(모든 비트스트링 → 에너지 매핑)로부터 QUBO를 구성한다. 두 가지 모드 제공:

- **정확 모드 (Exact)**: Möbius 변환 + Rosenberg 차수축소. 에너지 100% 일치, 보조변수 O(n·2^n) 폭발.
- **근사 모드 (Approx)**: 제약 최소제곱(QP). 보조변수 0개, ground state 보장, 에너지 근사.

```
입력:  {000→3, 001→4, 010→4, 011→5, 100→3, 101→5, 110→2, 111→1}
출력:  Q 행렬 (QUBO), ground state = 111 (energy=1) 보장
```

## 수학적 기반

### 정확 모드: Möbius 변환 + Rosenberg 차수축소

#### Step 1: Möbius 변환 (진리표 → 다중선형 다항식)

임의의 함수 f: {0,1}^n → R은 유일한 다중선형 다항식으로 표현 가능:

```
f(x) = Σ_{S ⊆ [n]} c_S · Π_{i∈S} x_i
c_S = Σ_{T ⊆ S} (-1)^{|S|-|T|} · f(1_T)
```

Fast Möbius 알고리즘 O(n · 2^n). 구현: `mobius_transform()`

**n=3 예시:**

```
f(x) = 3 + x₂ + x₃ + x₁x₃ - 2x₁x₂ - 3x₁x₂x₃
         ↑상수   ↑선형      ↑이차         ↑고차(3차)
```

#### Step 2: 항 분류

| 차수 | 처리 | QUBO 배치 |
|------|------|-----------|
| 0 (상수) | 에너지 오프셋 | Q에 미포함 |
| 1 (선형) | Q[i,i] += c_i | 대각 |
| 2 (이차) | Q[i,j] += c_ij | 비대각 |
| ≥3 (고차) | Rosenberg 차수축소 | 보조변수 도입 |

#### Step 3: Rosenberg 차수축소

k차 단항식에 대해 보조변수 y = x_{i1}·x_{i2} 도입, 패널티 M·(x_{i1}x_{i2} - 2x_{i1}y - 2x_{i2}y + 3y) 추가. 차수 > 2이면 반복 (체이닝). 단항식당 보조변수 k-2개.

#### Step 4: 패널티 강도

```
M = max(truth_table) - min(truth_table) + 1
```

### 근사 모드: 제약 최소제곱 (QP)

보조변수 없이 n-변수 QUBO를 직접 구한다.

```
min_Q  Σ_x (E_Q(x) - E_truth(x))^2      ← 진리표에 최대한 가깝게
s.t.   E_Q(target) + ε ≤ E_Q(x)          ← ground state 보장 (∀x ≠ target)
```

- **자유 파라미터**: n(n+1)/2개 (Q 행렬 upper triangular)
- **제약조건**: 2^n - 1개 (target이 모든 다른 상태보다 낮은 에너지)
- **풀이**: 무제약 최소제곱 → 위반 시 SLSQP (iterative cutting plane)

구현: `create_qubo_approx()`

**핵심 장점**: 보조변수 0개 → QUBO 크기 = n. n=14까지 실용적.

## 에너지 함수 프리셋

### Energy Gap

```python
E(target) = 0
E(x ≠ target) = gap + |N(0, noise_scale)|
```

gap ↓ = 난이도 ↑. 양자 어닐링 minimum spectral gap과 직결.

### Multi-Valley

```python
targets[0]: global minimum (energy = 0)
targets[1:]: local minima (energy = gap)
나머지: Hamming 거리 비례 에너지 + 장벽
```

local minima 수 ↑ = SA가 계곡에 갇힘. 양자 tunneling 벤치마크.

## 사용법

```bash
# 정확 모드
python truthtable/qubo_truthtable.py '{"000":3,"001":4,"010":4,"011":5,"100":3,"101":5,"110":2,"111":1}'

# 근사 모드 (--approx)
python truthtable/qubo_truthtable.py '{"000":3,"001":4,"010":4,"011":5,"100":3,"101":5,"110":2,"111":1}' --approx

# 프리셋 + 근사 모드
python truthtable/qubo_truthtable.py --preset gap 10 1011001100 --approx --seed 42
python truthtable/qubo_truthtable.py --preset valley 6 101010,010101 --approx --seed 42

# SA 실험
python truthtable/test_truthtable.py --gap 10
python truthtable/test_truthtable.py --valley 10
python truthtable/test_truthtable.py --scaling 10
python truthtable/test_truthtable.py --compare 10
```

## 파일 구조

```
truthtable/
├── qubo_truthtable.py   ← 생성기 (정확: Möbius+Rosenberg / 근사: QP)
├── test_truthtable.py   ← SA 실험 (gap/valley/scaling/비교)
└── results/             ← 생성된 QUBO (edge-list CSV)
```

---

## 실험 결과: 정확 모드 (Rosenberg)

### 실험 1: Energy Gap Sweep

에너지 갭 크기별 SA 성공률. gap이 작을수록 ground state 찾기 어려움.

```
설정:
  n = 5
  instances = 10 (인스턴스당 랜덤 target 생성)
  num_reads = 100 (SA 샘플 수/인스턴스)
  num_sweeps = 5000
  프리셋: preset_energy_gap(gap=각 값, noise_scale=1.0)
  성공률 = GS 찾은 read 수 / 전체 read 수 (10 instances × 100 reads = 1000 samples)
```

```
Gap    | GS Rate  | Avg Hamming | QUBO 크기
-------|----------|-------------|----------
0.1    |   6.60%  |  2.65       | 28변수
0.5    |   5.80%  |  2.53       | 28변수
1.0    |   9.90%  |  2.40       | 28변수
2.0    |  12.70%  |  2.46       | 28변수
5.0    |  17.20%  |  2.42       | 28변수
10.0   |  21.40%  |  2.03       | 28변수
```

### 실험 2: Multi-Valley Sweep

Local minima 수 증가에 따른 SA 성공률 변화.

```
설정:
  n = 5
  instances = 10
  num_reads = 100, num_sweeps = 5000
  프리셋: preset_multi_valley(gap=0.5, barrier_height=5.0)
  targets: 서로 Hamming 거리 ≥ n//3인 랜덤 bitstring
  성공률 = GS 찾은 read 수 / 전체 read 수
```

```
Valleys | GS Rate  | Avg Hamming
--------|----------|------------
1       |  17.30%  |  1.69
2       |  18.20%  |  1.81
3       |  16.00%  |  2.05
4       |  14.10%  |  2.16
```

### 실험 3: N-Scaling

n 증가에 따른 QUBO 크기 폭발 및 SA 성공률 변화.

```
설정:
  sizes = [3, 4, 5, 6, 7]
  runs = 10 (각 n마다 랜덤 target 10회)
  num_reads = 100, num_sweeps = 5000
  프리셋: preset_energy_gap(gap=2.0, noise_scale=1.0)
  성공률 = GS 찾은 read 수 / 전체 read 수
```

```
n  | QUBO 크기 | 보조변수 | SA 성공률
---|-----------|---------|----------
3  |         4 |       1 |   98.9%
4  |        10 |       6 |   51.6%
5  |        28 |      23 |   13.6%
6  |        78 |      72 |    0.6%
7  |       208 |     201 |    0.0%
```

### 실험 4: 5-way 비교

```
설정:
  n = 5, runs = 10
  num_reads = 100, num_sweeps = 5000
  Truth Table: preset_energy_gap(gap=1.0)
  성공률 기준 차이 주의:
    - Exact-Rosenberg: GS 찾은 read 수 / 전체 read 수 (1000 samples)
    - Wishart, ZeroExp, Posiform: best sample 기준 (10 runs 중 GS 찾은 run 수)
```

```
방법론              | SA 성공률
--------------------|----------
Exact-Rosenberg     |  13.20%
Wishart (α=0.7)     |  60.00%
ZeroExpectation     | 100.00%
Posiform            | 100.00%
```

---

## 실험 결과: 근사 모드 (QP)

### 실험 1: Energy Gap Sweep

정확 모드와 동일 조건에서 근사 모드로 실행.

```
설정:
  n = 5
  instances = 10 (인스턴스당 랜덤 target 생성)
  num_reads = 100, num_sweeps = 5000
  프리셋: preset_energy_gap(gap=각 값, noise_scale=1.0)
  성공률 = GS 찾은 read 수 / 전체 read 수 (1000 samples)
```

```
Gap    | GS Rate  | Avg Hamming | QUBO 크기
-------|----------|-------------|----------
0.1    |  76.30%  |  0.55       | 5변수
0.5    |  68.70%  |  0.95       | 5변수
1.0    |  65.50%  |  0.93       | 5변수
2.0    |  65.80%  |  0.84       | 5변수
5.0    |  65.70%  |  1.00       | 5변수
10.0   |  53.90%  |  1.15       | 5변수
```

### 실험 2: Multi-Valley Sweep

```
설정:
  n = 5, instances = 10
  num_reads = 100, num_sweeps = 5000
  프리셋: preset_multi_valley(gap=0.5, barrier_height=5.0)
  targets: 서로 Hamming 거리 ≥ n//3인 랜덤 bitstring
  성공률 = GS 찾은 read 수 / 전체 read 수
```

```
Valleys | GS Rate  | Avg Hamming
--------|----------|------------
1       |  59.60%  |  0.98
2       |  60.30%  |  1.04
3       |  65.70%  |  0.88
4       |  66.90%  |  0.83
```

### 실험 3: N-Scaling

```
설정:
  sizes = [3, 4, 5, 6, 7]
  runs = 10, num_reads = 100, num_sweeps = 5000
  프리셋: preset_energy_gap(gap=2.0, noise_scale=1.0)
  성공률 = GS 찾은 read 수 / 전체 read 수
```

```
n  | QUBO 크기 | 보조변수 | SA 성공률
---|-----------|---------|----------
3  |         3 |       0 |   71.9%
4  |         4 |       0 |   60.9%
5  |         5 |       0 |   53.4%
6  |         6 |       0 |   59.9%
7  |         7 |       0 |   62.2%
```

### 실험 4: 비교

```
설정:
  n = 5, runs = 10
  num_reads = 100, num_sweeps = 5000
  Approx-Gap: preset_energy_gap(gap=1.0)
  Approx-Valley: preset_multi_valley(gap=0.5, 2개 계곡)
  Exact-Rosenberg: preset_energy_gap(gap=1.0)
  성공률 기준 차이 주의:
    - Approx-Gap/Valley, Exact-Rosenberg: GS 찾은 read 수 / 전체 read 수 (1000 samples)
    - Wishart, ZeroExp, Posiform: best sample 기준 (10 runs 중 GS 찾은 run 수)
```

```
방법론              | SA 성공률
--------------------|----------
Exact-Rosenberg     |  13.80%
Approx-Gap          |  50.00%
Wishart (α=0.7)     |  50.00%
Approx-Valley       |  56.10%
ZeroExpectation     | 100.00%
Posiform            | 100.00%
```

---

## 정확 vs 근사 비교

### Energy Gap Sweep (n=5)

```
Gap  | Exact (Rosenberg) | Approx (QP) | 배율
-----|-------------------|-------------|-----
0.1  |   6.60%           |  76.30%     | 11.6x
0.5  |   5.80%           |  68.70%     | 11.8x
1.0  |   9.90%           |  65.50%     |  6.6x
2.0  |  12.70%           |  65.80%     |  5.2x
5.0  |  17.20%           |  65.70%     |  3.8x
10.0 |  21.40%           |  53.90%     |  2.5x
```

### N-Scaling

```
n | Exact (QUBO크기, SA)    | Approx (QUBO크기, SA)
--|------------------------|---------------------
3 | 4변수,  98.9%          | 3변수,  71.9%
4 | 10변수, 51.6%          | 4변수,  60.9%
5 | 28변수, 13.6%          | 5변수,  53.4%
6 | 78변수,  0.6%          | 6변수,  59.9%
7 | 208변수, 0.0%          | 7변수,  62.2%
```

### 대규모 근사 모드 테스트

```
설정:
  별도 스크립트로 실행 (test_truthtable.py 외부)
  각 n에 대해 단일 인스턴스 (랜덤 target)
  프리셋: preset_energy_gap(gap=2.0, noise_scale=1.0)
  num_reads = 100, num_sweeps = 5000
  성공률 = GS 찾은 read 수 / 전체 read 수
```

```
n   | QUBO 크기 | SA 성공률 | 생성 시간 | GS 보장
----|-----------|---------|----------|--------
8   |  8변수    |  100%   |  0.2s    | YES
10  | 10변수    |   67%   |  0.3s    | YES
12  | 12변수    |   55%   | 11.7s    | YES
14  | 14변수    |   52%   |  103s    | YES
```

### 모드 비교 요약

|                    | 정확 모드 (Rosenberg) | 근사 모드 (QP)         |
|--------------------|---------------------|----------------------|
| 보조변수            | O(n·2^n) 폭발       | **0개**               |
| QUBO 크기 (n=8)    | 530변수             | **8변수**             |
| 에너지 정확도       | 100% 일치           | 근사 (RMSE ~0.7)      |
| Ground state       | 보장                | **보장** (QP 제약)     |
| SA 성공률 (n=8)    | 0%                  | **100%**              |
| 실용 한계           | n ≤ 5~6            | **n ≤ 14+**           |
| 생성 시간 (n=14)   | —                   | ~100s (SLSQP 병목)    |
| 에너지 순서 보존율   | 100%               | ~55%                  |

---

## 한계

### 정확 모드

1. **보조변수 지수 폭발**: n=8 → QUBO 530변수. SA 실용 한계 n ≤ 5~6.
2. **변수 효율**: 28변수 QUBO로 5비트 문제만 표현. Posiform은 같은 크기로 28비트 표현.
3. **Hardness 출처가 인위적**: 난이도가 Rosenberg 패널티 구조에서 발생. 에너지 갭 제어 효과가 묻힘.

### 근사 모드

1. **에너지 정확도 상실**: RMSE ~0.7, 에너지 순서 보존율 ~55%. 전체 에너지 스펙트럼은 근사.
2. **SLSQP 생성 시간**: n=14에서 ~100초. n=15+ 는 진리표 열거(2^n) 자체가 병목.
3. **진리표 크기 한계**: 입력이 2^n개이므로 n ≤ ~20에서만 진리표 열거 가능.

### 공통

진리표 기반 접근 자체가 n ≤ ~20으로 제한됨 (2^n 열거). 구조화된 에너지 함수는 대부분 저차(≤2차)가 되어 이 방법론의 차별점이 약화됨.

## 결론

- **정확 모드**: 수학적으로 완벽하나 보조변수 폭발로 실용성 제한 (n ≤ 5~6).
- **근사 모드**: 보조변수 0개로 n=14까지 실용적. Ground state 보장. SA 성공률 정확 모드 대비 5~12배 개선. Wishart(50%)과 동급 난이도.
- **벤치마크 활용**: 소규모 양자 프로세서에서 "에너지 랜드스케이프를 알고 있는" 벤치마크로 활용 가능. 근사 모드 권장.
