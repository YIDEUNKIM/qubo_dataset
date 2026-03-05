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
|------|------|-----------
| 0 (상수) | 에너지 오프셋 | Q에 미포함 |
| 1 (선형) | Q[i,i] += c_i | 대각 |
| 2 (이차) | Q[i,j] += c_ij | 비대각 |
| ≥3 (고차) | Rosenberg 차수축소 | 보조변수 도입 |

#### Step 3: Rosenberg 차수축소

k차 단항식에 대해 보조변수 y = x_{i1}·x_{i2} 도입, 패널티 M·(x_{i1}x_{i2} - 2x_{i1}y - 2x_{i2}y + 3y) 추가. 차수 > 2이면 반복 (체이닝). 단항식당 보조변수 k-2개.

#### Step 4: 패널티 강도

```
M = max(tt_range, max_y(S_y)) + 1

tt_range = max(E) - min(E)                  ← 진리표 에너지 범위
S_y = Σ |coeff| for reduced terms involving y  ← 보조변수 y 관련 축소 항 계수 합
```

두 조건 모두 만족해야 함:
1. M > tt_range: 진리표 에너지 범위 초과 (패널티 위반이 에너지 이득보다 큼)
2. M > max(S_y): 축소된 항의 계수 합 초과 (Möbius 변환 고차항이 진리표 범위보다 클 수 있음)

조건 2가 없으면 SA가 보조변수 제약을 위반하여 더 낮은 에너지를 찾을 수 있음.

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

**핵심 장점**: 보조변수 0개 → QUBO 크기 = n. n=16까지 실용적.

## 에너지 함수 프리셋

### Random Landscape

```python
E(target) = 0
E(x ≠ target) = uniform(0.1, 5.0)
```

target이 유일한 ground state이고, 나머지는 0.1~5.0 범위에서 균등 분포. 구현: `preset_random_landscape()`

### Multi-Valley

```python
targets[0]: global minimum (energy = 0)
targets[1:]: local minima (energy = gap)
나머지: Hamming 거리 비례 에너지 + 장벽
```

local minima 수 ↑ = SA가 계곡에 갇힘. 양자 tunneling 벤치마크. 구현: `preset_multi_valley()`

## 사용법

```bash
# 정확 모드
python truthtable/qubo_truthtable.py '{"000":3,"001":4,"010":4,"011":5,"100":3,"101":5,"110":2,"111":1}'

# 근사 모드 (--approx)
python truthtable/qubo_truthtable.py '{"000":3,"001":4,"010":4,"011":5,"100":3,"101":5,"110":2,"111":1}' --approx

# 프리셋 + 근사 모드
python truthtable/qubo_truthtable.py --preset random 10 1011001100 --approx --seed 42
python truthtable/qubo_truthtable.py --preset valley 6 101010,010101 --approx --seed 42

# SA 실험
python truthtable/test_truthtable.py --valley 10
python truthtable/test_truthtable.py --scaling 10
python truthtable/test_truthtable.py --compare 10
python truthtable/test_truthtable.py --strategy 10
python truthtable/test_truthtable.py --greedy-scaling 10
python truthtable/test_truthtable.py --sweep 20
```

## 파일 구조

```
truthtable/
├── qubo_truthtable.py   ← 생성기 (정확: Möbius+Rosenberg / 근사: QP)
├── test_truthtable.py   ← SA 실험 (valley/scaling/비교/전략/sweep)
└── results/             ← 생성된 QUBO (edge-list CSV)
```

---

## 실험 설정 총정리

### 공통 설정

| 항목 | 값 |
|------|-----|
| SA solver | neal.SimulatedAnnealingSampler (D-Wave) |
| num_reads | 100 |
| num_sweeps | 5000 (sweep 실험 제외) |
| 프리셋 | `preset_random_landscape` (E(target)=0, 나머지=uniform(0.1, 5.0)) |

### 실험별 설정

| # | 실험명 | 모드 | n | runs | 프리셋 | 변수 파라미터 | CLI |
|:-:|--------|:----:|:-:|:----:|--------|-------------|-----|
| 1 | Valley Sweep | Exact | 5 | 10 | `multi_valley(gap=0.5, barrier=5.0)` | valleys: 1,2,3,4 | `--valley 10` |
| 2 | N-Scaling | Exact greedy | 3~7 | 10 | `random_landscape` | n: 3,4,5,6,7 | `--scaling 10` |
| 3 | 7-way 비교 | Both | 5 | 10 | `random_landscape` / `multi_valley` | — | `--compare 10` |
| 4 | 3-전략 비교 | Exact | 3~8 | 10 | `random_landscape` | 전략: original/cache/greedy | `--strategy 10 3,4,5,6,7,8` |
| 5 | Greedy vs Approx | Both | 3~10 | 10 | `random_landscape` | n: 3~10 | `--greedy-scaling 10` |
| 6 | Sweep 전이 | Exact greedy | 8 | 20 | `random_landscape` | sweeps: 50~10000 | `--sweep 20` |

### 성공률 측정 기준

| 방법론 | 기준 |
|--------|------|
| Exact-Rosenberg / Approx | **전체 read 기준** (instances × 100 reads 중 GS 찾은 비율) |
| Wishart / ZeroExp / Posiform | **best sample 기준** (instance당 최적 1개만 확인) |

---

## 실험 결과

### 실험 1: Multi-Valley Sweep

Local minima 수 증가에 따른 SA 성공률 변화.

> 설정: n=5, instances=10, num_reads=100, num_sweeps=5000, 프리셋: preset_multi_valley(gap=0.5, barrier_height=5.0)
> targets: 서로 Hamming 거리 ≥ n//3인 랜덤 bitstring

| Valleys | GS Rate | Avg Hamming |
|:-------:|:-------:|:-----------:|
| 1 | 49.40% | 0.71 |
| 2 | 40.40% | 1.25 |
| 3 | 36.10% | 1.59 |
| 4 | 32.80% | 1.69 |

Valley 수 증가에 따라 성공률이 단조 감소 (49.4% → 32.8%). Hamming 거리도 증가 — SA가 local minima에 갇힘.

### 실험 2: N-Scaling (정확 모드, greedy)

n 증가에 따른 QUBO 크기 폭발 및 SA 성공률 변화.

> 설정: sizes=[3,4,5,6,7], runs=10, num_reads=100, num_sweeps=5000, 프리셋: preset_random_landscape

| n | QUBO 크기 | 보조변수 | SA 성공률 |
|:-:|:---------:|:-------:|:---------:|
| 3 | 4 | 1 | 76.60% |
| 4 | 6 | 2 | 43.20% |
| 5 | 10 | 5 | 16.30% |
| 6 | 18 | 12 | 8.00% |
| 7 | 22 | 15 | 5.10% |

n=6~7에서도 M 보정 후 SA가 소폭 성공 (8.0%, 5.1%). 정확 모드 실용 한계 = n ≤ 6~7 (greedy + 보정된 M).

### 실험 3: 7-way 비교

> 설정: n=5, runs=10, num_reads=100, num_sweeps=5000
> 성공률 기준 차이: Exact/Approx는 전체 read 기준, Wishart/ZeroExp/Posiform은 best sample 기준

| 방법론 | SA 성공률 |
|--------|:---------:|
| Exact-Rosenberg | 22.90% |
| Wishart (α=0.7) | 40.00% |
| Approx-Valley | 54.00% |
| Approx-Random | 67.60% |
| ZeroExpectation | 100.00% |
| Posiform | 100.00% |

### 실험 4: 차수축소 전략 비교

3가지 Rosenberg 차수축소 전략의 성능 비교.

> 설정: n=[3,4,5,6,7,8], 전략=original/cache/greedy, 10 runs, SA: reads=100, sweeps=5000
> 프리셋: preset_random_landscape

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

- **n=4~6**: greedy가 original 대비 1.1~3.1배 개선
- **n=7~8**: greedy가 유일하게 의미 있는 성공률 유지 (3.5%, 2.1%)
- **보조변수 절감**: greedy가 n=8에서 95.8% 절감 (522→22)

### 실험 5: Greedy 확장 스케일링 (exact greedy vs approx)

> 설정: n=[3..10], 모드=greedy/approx, 10 runs, SA: reads=100, sweeps=5000
> 프리셋: preset_random_landscape

| n | Greedy QUBO | Greedy Rate | Approx QUBO | Approx Rate | Greedy Gen | Approx Gen |
|:-:|:-----------:|:-----------:|:-----------:|:-----------:|:----------:|:----------:|
| 3 | 4 | **76.60%** | 3 | 65.10% | 0.000s | 0.010s |
| 4 | 6 | 38.20% | 4 | **70.00%** | 0.000s | 0.001s |
| 5 | 10 | 23.20% | 5 | **68.00%** | 0.000s | 0.001s |
| 6 | 18 | 7.50% | 6 | **50.80%** | 0.001s | 0.002s |
| 7 | 22 | 3.30% | 7 | **59.30%** | 0.002s | 0.471s |
| 8 | 30 | 1.20% | 8 | **54.50%** | 0.007s | 0.581s |
| 9 | 46 | 0.70% | 9 | **50.20%** | 0.025s | 1.817s |
| 10 | 78 | 0.20% | 10 | **59.10%** | 0.091s | 3.604s |

- **교차점**: n=3에서만 greedy 우세, n≥4부터 approx가 역전
- **Approx 안정성**: n=3~10 전 범위에서 50~70% 유지
- **Greedy 한계**: n=8부터 ~1% 이하
- **결론**: random landscape에서는 **n≥4부터 approx 모드 권장**

### 실험 6: Sweep 전이 (S-curve)

> 설정: n=8 (greedy: ~30변수 QUBO), 20 instances (동일 인스턴스 재사용), SA: reads=100
> 프리셋: preset_random_landscape

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

M 보정 후 SA가 소폭 성공하지만, sweep 200배 증가(50→10000)에도 0.65% → 2.10%로 거의 포화. Rosenberg 패널티 구조가 만드는 에너지 랜드스케이프는 SA에 근본적으로 어려움. 30변수 QUBO이지만 유효 자유도는 8개뿐이며 나머지 22개는 패널티에 의해 구속됨.

---

## 모드 비교 요약

|  | 정확 original | 정확 greedy | 근사 (QP) |
|--|:------------:|:-----------:|:---------:|
| 보조변수 | O(n·2^n) 폭발 | **78~96% 절감** | **0개** |
| QUBO 크기 (n=8) | 530변수 | **30변수** | **8변수** |
| 에너지 정확도 | 100% 일치 | 100% 일치 | 근사 (RMSE ~0.7) |
| Ground state | 보장 | 보장 | **보장** (QP 제약) |
| SA 성공률 (n=5) | 10.70% | 18.30% | **68.00%** |
| 실용 한계 | n ≤ 4~5 | **n ≤ 7~8** | **n ≤ 16** |

---

## 한계

### 정확 모드

1. **보조변수 지수 폭발**: n=8 → greedy로도 30변수 (original: 530변수). SA 실용 한계 greedy로 n ≤ 7~8.
2. **변수 효율**: 30변수 QUBO로 8비트 문제만 표현. Posiform은 같은 크기로 30비트 표현.
3. **Hardness 출처가 인위적**: 난이도가 Rosenberg 패널티 구조에서 발생.
4. **SA-hard 구조**: n=8에서 sweep 200배 증가에도 성공률 0.65% → 2.10%로 거의 포화.

### 근사 모드

1. **에너지 정확도 상실**: RMSE ~0.7, 에너지 순서 보존율 ~55%. 전체 에너지 스펙트럼은 근사.
2. **SLSQP 생성 시간**: n=14에서 ~100초. n=15+ 는 진리표 열거(2^n) 자체가 병목.
3. **진리표 크기 한계**: 입력이 2^n개이므로 n ≤ ~20에서만 진리표 열거 가능.

### 공통

진리표 기반 접근 자체가 n ≤ ~20으로 제한됨 (2^n 열거). 구조화된 에너지 함수는 대부분 저차(≤2차)가 되어 이 방법론의 차별점이 약화됨.

## 결론

- **정확 모드 + greedy**: original 대비 보조변수 78~96% 절감. M 보정 후 n=7~8에서도 소폭 성공 (3.5~5.1%).
- **근사 모드**: 보조변수 0개로 n=16까지 실용적. Ground state 보장. n≥4에서 안정적 50~70%.
- **전략 선택 기준**: n=3에서만 greedy 정확 모드 우세, n≥4에서는 approx 근사 모드 권장.
- **벤치마크 활용**: 소규모 양자 프로세서에서 "에너지 랜드스케이프를 알고 있는" 벤치마크로 활용 가능.
