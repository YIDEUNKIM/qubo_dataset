# Degree-Bounded Möbius QUBO — 실험 보고서

## 1. 배경 및 동기

### Truth Table QUBO의 두 가지 한계

1. **크기 제한**: 진리표 2^n → 단일 블록 n ≤ 12
2. **Rosenberg 비용**: exact 모드에서 고차항 차수축소 → 보조변수 폭발 → SA-hard (penalty 지배)

Concat 방식은 한계 1을 우회하지만, exact 모드의 Rosenberg 비용(한계 2)은 해결하지 못함.

### 아이디어: Degree-Bounded Möbius 계수

2^n 진리표를 거치지 않고, **degree ≤ d의 Möbius 계수를 직접 생성**:
- O(n^d) 시간 (2^n 불필요)
- 양수 1차항 + 제한된 음수 고차항 → ground state 수학적 보장
- d = 2면 Rosenberg 불필요 → 보조변수 0

---

## 2. 수학적 기반

### Theorem 1 (Ground State Guarantee)

다변량 다중선형 다항식 f(z) = Σ_{S⊆[n], |S|≥1} c_S · ∏_{i∈S} z_i  (z ∈ {0,1}^n) 에서:

**조건:**
- (C1) c\_{i} > 0  ∀i ∈ [n]  (양수 1차항)
- (C2) Σ\_{|S|≥2} |c_S| < min_i c\_{i}  (고차항 바운드)

**결론:** f(0) = 0 이고 f(z) > 0  ∀z ≠ 0  (z=0이 유일한 ground state)

**증명:**

z ≠ 0이면 T = supp(z) = {i : z_i = 1} ≠ ∅.

```
f(z) = Σ_{S⊆T, |S|≥1} c_S
     = Σ_{i∈T} c_{i}  +  Σ_{S⊆T, |S|≥2} c_S
     ≥ Σ_{i∈T} c_{i}  -  Σ_{S⊆T, |S|≥2} |c_S|
     ≥ min_i c_{i}     -  Σ_{|S|≥2} |c_S|       (∵ |T| ≥ 1)
     > 0                                          (∵ C2)    □
```

### Corollary: SA-Trivial Landscape (d = 2)

d = 2이고 모든 c_S ≥ 0 (|S| ≥ 1) 이면, 임의의 z에서 z_i를 1→0으로 flip할 때:

```
ΔE = -Σ_{S: i∈S} c_S · ∏_{j∈S\{i}} z_j ≤ 0
```

모든 항이 비음수이므로 ΔE ≤ 0. 따라서 **로컬 미니마가 z=0 외에 존재하지 않음**.

음수 고차항이 있어도, (C2) 제약으로 크기가 제한되므로 d=2에서는 SA가 greedy descent만으로 ground state를 찾을 수 있음.

### 변수 치환

Target t ∈ {0,1}^n에 대해 z_i = x_i ⊕ t_i로 치환:
- z_i = x_i (t_i = 0)
- z_i = 1 - x_i (t_i = 1)

z-공간에서 z = 0이 ground state ↔ x-공간에서 x = target이 ground state.

---

## 3. 알고리즘

### Pipeline

```
Step 1: z-공간 Möbius 계수 생성 (degree ≤ d, Theorem 1 조건 충족)
Step 2: z→x 변환 (다항식 치환, O(|terms| · 2^d))
Step 3: 항 분류 (constant / linear / quadratic / higher_order)
Step 4: Rosenberg 차수축소 (degree ≥ 3, greedy — d ≥ 3일 때만)
Step 5: 패널티 강도 계산 (분석적, 진리표 불필요)
Step 6: QUBO 조립
Step 7: (선택) Posiform hardening
```

### 시간 복잡도

| 단계 | 복잡도 | 비고 |
|------|--------|------|
| 계수 생성 | O(Σ_{k=1}^{d} C(n,k)) ≈ O(n^d/d!) | 고차항 후보 열거 |
| z→x 변환 | O(|terms| · 2^d) | 다항식 치환, d 고정 시 선형 |
| Rosenberg | O(|higher_order| · d) | d ≥ 3일 때만 |
| **전체** | **O(n^d)** | **d 고정 시 다항식** |

### 실용적 크기 한계

| d | 계수 수 (n=100) | 계수 수 (n=500) | Rosenberg | 실용 한계 |
|---|----------------|----------------|-----------|-----------|
| 2 | 5,050 | 125,250 | 불필요 | **n 무제한** |
| 3 | 166,750 | 20,833,750 | 필요 | n ≤ 100 |
| 4 | 4,087,975 | — | 필요 | n ≤ 50 |

d = 2가 Rosenberg 없이 O(n²)으로 생성 가능하여 가장 실용적.

---

## 4. 실험 결과

### 4.1 실험 1: Ground State 검증 (brute force)

n ≤ 15에서 2^n 완전 열거로 Theorem 1의 ground state 보장 검증.

| n | d | GS 정확 | 유일성 | 보조변수 |
|---|---|---------|--------|----------|
| 8 | 2 | **20/20** | **20/20** | 0 |
| 8 | 3 | **20/20** | **20/20** | 12 |
| 10 | 2 | **20/20** | **20/20** | 0 |
| 10 | 3 | **20/20** | **20/20** | 21 |
| 10 | 4 | **20/20** | **20/20** | 41 |
| 12 | 2 | **20/20** | **20/20** | 0 |
| 12 | 3 | **20/20** | **20/20** | 34 |
| 15 | 2 | **20/20** | **20/20** | 0 |
| 15 | 3 | **20/20** | **20/20** | 49 |

**모든 config에서 100% ground state + 100% 유일성.** Theorem 1이 완벽히 동작.

### 4.2 실험 2: d-sweep (n=20)

degree별 SA 성공률. density=0.5.

**설정**: n=20, runs=10, reads=100, sweeps=1000

| d | QUBO Size | 보조변수 | GS Rate |
|---|-----------|----------|---------|
| 2 | 20 | 0 | **100.00%** |
| 3 | 118 | 98 | **0.00%** |
| 4 | 210 | 190 | **0.00%** |
| 5 | 702 | 682 | **0.00%** |

**관찰:**

1. **d=2 → d=3에서 급격한 전이**: d=2(100%) → d=3(0%). QUBO 크기가 20→118(×5.9)으로 팽창하면서 SA가 완전히 실패.
2. **d≥3의 보조변수 폭발**: n=20에서도 d=3은 98개, d=5는 682개 보조변수. Rosenberg penalty가 원래 에너지를 완전히 지배.
3. **d=2의 SA-trivial 특성**: Corollary에 의해 로컬 미니마가 없으므로 SA가 100% 성공.

### 4.3 실험 3: n-scaling (d=2)

**설정**: d=2, runs=10, reads=100, sweeps=1000

| n | QUBO Size | GS Rate | 생성+SA 시간 |
|---|-----------|---------|-------------|
| 20 | 20 | **100.00%** | 0.1s |
| 50 | 50 | **100.00%** | 0.3s |
| 100 | 100 | **100.00%** | 0.6s |
| 200 | 200 | **100.00%** | 1.9s |
| 500 | 500 | **100.00%** | 9.8s |

**n=500에서도 100%.** d=2의 SA-trivial 특성으로 인해 n에 무관하게 항상 풀림.

### 4.4 실험 4: higher_ratio (γ) sweep (d=2, n=100)

γ 증가 → 음수 2차항 증가 → 난이도 변화 여부 확인.

| γ | GS Margin | GS Rate |
|---|-----------|---------|
| 0.3 | 0.8510 | **100.00%** |
| 0.5 | 0.7513 | **100.00%** |
| 0.8 | 0.6017 | **100.00%** |
| 0.95 | 0.5269 | **100.00%** |
| 0.99 | 0.5070 | **100.00%** |

**γ=0.99에서도 100%.** 음수 2차항의 크기가 제한되어 로컬 미니마를 충분히 만들지 못함.

### 4.5 실험 5: 방법론 비교 (N=50)

**설정**: N=50, runs=10, reads=100, sweeps=1000

| Method | GS Rate | Time |
|--------|---------|------|
| DegreeBounded(d=3) | **0.00%** | 12.9s |
| Hardened(α=0.01) | **85.60%** | 2.0s |
| DegreeBounded(d=2) | **100.00%** | 0.3s |
| DegreeBounded(d=2, α=0.1) | **100.00%** | 1.6s |
| Hardened(α=0.1) | **100.00%** | 2.0s |

---

## 5. 핵심 발견: 근본적 트레이드오프

### SA-Hardness vs Ground State Guarantee

```
d = 2:  GS 보장 ✓ + Rosenberg 불필요 ✓ → SA-trivial (벤치마크 가치 없음)
d ≥ 3:  GS 보장 ✓ + Rosenberg 필요    → SA-unsolvable (보조변수 지배)
```

이 트레이드오프는 **Theorem 1의 (C2) 조건에서 필연적으로 발생**:

- (C2)는 고차항 크기를 1차항 최소값 이하로 제한
- d=2에서 이 제약은 에너지 landscape를 "거의 양수"로 만들어 로컬 미니마 소멸
- d≥3에서 Rosenberg penalty M ≫ 원래 에너지 → SA가 penalty만 봄

**Hardened Posiform과의 차이:**
- Hardened는 random QUBO(로컬 미니마 풍부) + posiform(GS 보장) 결합
- DegreeBounded는 단일 다항식으로 GS 보장 + 난이도를 동시에 달성하려 함
- 단일 다항식에서 (C2) 조건과 SA-hardness는 구조적으로 양립 불가

### 결론

Degree-Bounded Möbius 접근법은:

1. **수학적으로 정확**: Theorem 1에 의한 ground state 보장이 100% 동작 (20/20 brute force 검증)
2. **시간 복잡도 우수**: O(n^d)로 2^n 불필요, n=500까지 즉시 생성 (d=2)
3. **SA 벤치마크로는 부적합**:
   - d=2: 로컬 미니마 없음 → 항상 100% → 솔버 평가 불가
   - d≥3: Rosenberg 비용 → 항상 0% → 솔버 평가 불가
4. **유용한 용도**: 정답이 보장된 QUBO 인스턴스 생성 (검증/테스트용), d=2 + posiform hardening 조합으로 난이도 조절 가능

---

## 6. 사용법

```bash
# 기본 생성 (d=3)
python3 degree_bounded/qubo_degree_bounded.py 10110 3 --seed 42

# d=2 (보조변수 없음, 대규모 가능)
python3 degree_bounded/qubo_degree_bounded.py 10110011 2

# 대규모 (n=100, d=2)
python3 degree_bounded/qubo_degree_bounded.py 10011...  2 --seed 42

# Posiform hardening
python3 degree_bounded/qubo_degree_bounded.py 10110 2 --harden 0.1

# 실험
python3 degree_bounded/test_degree_bounded.py --verify 20    # GS 검증
python3 degree_bounded/test_degree_bounded.py --d-sweep 10   # d별 SA
python3 degree_bounded/test_degree_bounded.py --scaling 10   # n-scaling
python3 degree_bounded/test_degree_bounded.py --compare 10   # 방법론 비교
```
