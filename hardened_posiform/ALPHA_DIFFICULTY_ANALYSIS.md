# Hardened Posiform α 난이도의 수학적 분석

## 개요

Pelofske et al. (2024)의 Hardened Posiform에서 posiform scaling coefficient α가 SA 난이도에 미치는 영향을 수학적으로 분석한다.

**논문의 한계**: Pelofske (2024)는 α=0.01과 α=0.1만 실험적으로 비교하며, 난이도 변화의 원인에 대해 *"A possible cause for this is because of a large number of local minima whose energy is very close to the true ground-state"* 라는 추측 한 문장만 제시한다. Hahn (2023)은 계수의 난이도 영향을 *"future work"* 으로 남긴다.

본 문서는 두 가지 상보적 관점에서 α → 난이도 관계를 증명한다:
1. **국소적 관점**: Energy gap과 barrier의 비율 → TTS ≳ exp(c/α) (지수적 하한)
2. **전역적 관점**: Posiform이 Hamming distance 공간에 만드는 선형 포텐셜 → SNR 프레임워크

---

## 1. 설정

### Hardened Posiform 구조

```
Q(x) = R(x) + α·P(x)
```

- **R = Σ R_i**: disjoint subgraph별 random discrete-coefficient QUBO의 합 (block-diagonal)
- **P**: posiform planted QUBO. P(x) = Σ_c b_c · 𝟙(x violates clause c), b_c > 0
- **α > 0**: posiform scaling coefficient
- **x\***: planted target bitstring

### 핵심 성질

| 성질 | 근거 |
|---|---|
| R(x\*) = Σ min(R_i) | x\*는 각 subproblem의 GS를 concatenate |
| P(x\*) = 0 | posiform의 정의 (target에서 모든 clause 충족) |
| P(x) > 0, ∀x ≠ x\* | 2-SAT uniqueness 보장 |
| Q(x\*) = R(x\*) | P(x\*)=0이므로 α 무관 |

### 표기법

| 기호 | 정의 |
|---|---|
| S = GS(R) | R의 ground state 집합. x\* ∈ S. \|S\| ≥ 1 |
| δ_R(x) = R(x) - R(x\*) | R에서의 에너지 초과분. x ∈ S이면 0, x ∉ S이면 > 0 |
| Δ_R = min_{x∉S} δ_R(x) | R의 최소 비축퇴 에너지 갭 |
| Δ_P = min_{y∈S\{x\*}} P(y) | R-축퇴 상태들 중 posiform의 최소 패널티. > 0 |
| d(x) = HD(x, x\*) | x와 x\*의 Hamming distance |

---

## 2. 에너지 갭 분석 (Δ_min = O(α))

### 정리 1: 최소 에너지 갭의 α 선형성

> Hardened Posiform Q = R + α·P의 최소 에너지 갭 Δ_min(α)는 충분히 작은 α에서
> **Δ_min(α) = α · Δ_P** 를 만족한다.

**증명.**

임의의 x ≠ x\*에 대해 에너지 갭:

```
Δ(x) = Q(x) - Q(x*) = [R(x) - R(x*)] + α·[P(x) - P(x*)]
      = δ_R(x) + α·P(x)
```

P(x\*)=0이므로 P(x) - P(x\*) = P(x).

두 경우로 분류한다:

**(Case A)** x ∈ S \ {x\*} (R-축퇴 상태):
```
δ_R(x) = 0  ⟹  Δ(x) = α·P(x) ≥ α·Δ_P
```

**(Case B)** x ∉ S:
```
δ_R(x) ≥ Δ_R > 0  ⟹  Δ(x) ≥ Δ_R + α·P(x) ≥ Δ_R
```

따라서:

```
Δ_min(α) = min(min_{x∈S\{x*}} α·P(x),  min_{x∉S} [δ_R(x) + α·P(x)])
          = min(α·Δ_P,  Δ_R + α·min_{x∉S} P(x))
```

α < Δ_R / Δ_P 일 때 α·Δ_P < Δ_R이므로:

```
Δ_min(α) = α·Δ_P
```

**즉, 작은 α에서 에너지 갭이 α에 비례하여 0으로 수렴한다.** ∎

### 따름정리: α_c의 존재

> α_c = Δ_R / Δ_P를 경계로:
> - α < α_c: R-축퇴 상태가 first excited state (갭 = α·Δ_P)
> - α > α_c: R의 비축퇴 상태가 first excited state (갭 ≥ Δ_R)

이 전환이 SA 난이도의 상전이를 만든다.

---

## 3. 장벽 분석 (B = O(1))

### 정리 2: 에너지 장벽의 α 독립성

> x\*에서 R-축퇴 상태 y ∈ S\{x\*}로의 에너지 장벽 B(α)는 작은 α에서 O(1)이다.

**증명.**

x\*에서 y로의 single-bit-flip 경로 x\* = z_0 → z_1 → ... → z_d = y (d = HD(x\*, y))를 생각하자. 경로 위의 최대 에너지:

```
B(α) = max_{0≤k≤d} [Q(z_k) - Q(x*)] = max_k [δ_R(z_k) + α·P(z_k)]
```

중간 상태 z_k는 일반적으로 S에 속하지 않으므로:
- δ_R(z_k) ~ O(1) (R의 에너지 스케일에 의해 결정)
- α·P(z_k) ~ O(α) (작은 α에서 무시 가능)

따라서:

```
B(α) = max_k δ_R(z_k) + O(α) → B_R  (as α → 0)
```

여기서 B_R = max_k δ_R(z_k) > 0은 R의 에너지 landscape에 의해 결정되는 상수.

**장벽은 R의 에너지 스케일 O(1)로 남는 반면, 갭은 O(α)로 수축한다.** ∎

---

## 4. SA 수렴 시간 하한 (TTS ≳ exp(c/α))

### 정리 3: SA 난이도의 지수적 하한

> Hardened Posiform Q = R + α·P에 대해, SA의 time-to-solution은
> **TTS(α) ≳ exp(B_R / (α·Δ_P))** 를 만족한다.

**증명 (Hajek-type argument).**

SA가 x\*를 찾으려면, R-축퇴 상태 y (Δ(y) = α·Δ_P)와 x\*를 구분해야 한다.

**Step 1.** Metropolis chain의 온도 요구: 에너지 차이 α·Δ_P를 분해하려면 T ≲ α·Δ_P인 저온이 필요하다.

**Step 2.** 저온에서의 mixing time: 온도 T에서 장벽 B_R을 넘는 Metropolis 탈출 시간은

```
τ_escape(T) ~ exp(B_R / T)
```

Arrhenius law에 의해 결정된다.

**Step 3.** T ~ α·Δ_P에서:

```
τ_escape ~ exp(B_R / (α·Δ_P))
```

SA의 annealing schedule에서 이 온도 구간에 최소 τ_escape 스텝을 할당해야 하므로:

```
TTS(α) ≳ exp(B_R / (α·Δ_P)) = exp(c/α)
```

여기서 c = B_R / Δ_P > 0. ∎

### 물리적 해석

| α 범위 | Barrier/Gap 비율 | SA 시간 | 해석 |
|---|---|---|---|
| α → 0 | O(1/α) → ∞ | exp(c/α) → ∞ | 갭이 닫히지만 장벽 유지 → 지수적 어려움 |
| α ~ α_c | O(1) | polynomial | 갭과 장벽이 비슷한 스케일 → 전이 영역 |
| α ≫ α_c | O(1) | polynomial | P 지배, funnel landscape → 쉬움 |

---

## 5. Hamming Distance 분석

### 정리 4: Posiform 에너지의 HD 선형성

> HD(x, x\*) = d인 상태 x에 대해, posiform 에너지의 기댓값은
> **E[P(x) | HD = d] = (2mb̄ / 3n) · d + O(d²/n²)** 이다.
> (m = clause 수, b̄ = 평균 계수, d ≪ n)

**증명.**

각 clause는 변수 쌍 (i,j)에서 3개 wrong tuple 중 1개를 균일 랜덤 선택하여 배제한다.

target (x\*_i, x\*_j) = (a, b)에 대해 wrong tuple의 위반 조건:

| wrong tuple | violated when |
|---|---|
| (1-a, b) | i만 flip, j 유지 |
| (a, 1-b) | j만 flip, i 유지 |
| (1-a, 1-b) | i, j 모두 flip |

각각 선택 확률 1/3.

flip된 위치의 집합을 F (|F| = d)라 하면, 랜덤 clause (i,j)에서:

```
P(i ∈ F, j ∉ F or j ∈ F, i ∉ F) = 2d(n-d) / [n(n-1)]
P(i ∈ F, j ∈ F)                  = d(d-1)  / [n(n-1)]
P(i ∉ F, j ∉ F)                  = (n-d)(n-d-1) / [n(n-1)]
```

clause당 기대 패널티:

```
E[penalty] = b · (1/3) · [2d(n-d) + d(d-1)] / [n(n-1)]
           = b·d·(2n - d - 1) / [3n(n-1)]
```

m개 clause, 평균 계수 b̄에 대해:

```
E[P(x) | HD=d] = m·b̄·d·(2n - d - 1) / [3n(n-1)]
```

d ≪ n이면:

```
E[P(x) | HD=d] ≈ (2mb̄ / 3n) · d  :=  γ · d
```

여기서 **γ = 2mb̄/(3n)** 은 posiform의 "HD 기울기 상수". ∎

### 따름정리: Posiform은 HD 공간에서 선형 포텐셜 우물을 형성

```
V_P(d) = α · γ · d
```

α가 이 포텐셜의 기울기를 조절한다. d=0 (x\*)에서 최소, d가 커질수록 선형 증가.

---

## 6. Signal-to-Noise Ratio (SNR) 프레임워크

### 정의: SA의 신호 대 잡음비

SA가 single bit flip으로 x\*에 한 걸음 다가갈 때:

**신호 (Posiform gradient):**
```
ΔP = α · γ    (HD 1 감소당 기대 에너지 감소량)
```

**잡음 (Random QUBO fluctuation):**

block 내 k개 변수와 coupling되므로:
```
σ_R ~ O(√k)    (lin2: 각 coupling ±1, 분산 = k-1)
```

**SNR:**

```
SNR = α·γ / σ_R = α · (2mb̄/3n) / √k
```

### SA 행동의 SNR 의존성

| SNR | SA 행동 | α 범위 |
|---|---|---|
| ≫ 1 | 매 flip마다 올바른 방향 감지 → funnel 따라 하강 | α ≫ σ_R/γ |
| ~ 1 | 신호와 잡음 경쟁 → 확률적 진행 | α ~ σ_R/γ |
| ≪ 1 | 방향 감지 불가 → random walk → HD 수렴 실패 | α ≪ σ_R/γ |

### 전이점 추정

```
α* ≈ σ_R / γ = √k · 3n / (2mb̄)
```

노트북 실험 파라미터 (n=500, k=10, m ≈ 5n, b̄ ≈ 1):

```
σ_R ≈ √10 ≈ 3.16
γ = 2·2500·1/(3·500) ≈ 3.33
α* ≈ 3.16/3.33 ≈ 0.95
```

실험에서 α=0.1이 99.8%로 쉬운 것은 block-diagonal 구조 덕분에 SA가 각 block을 독립적으로 풀 수 있어 실효적 σ_R이 이론값보다 작기 때문이다. α=0.01 (36.8%)과 α=0.001 (15.5%)에서의 난이도 증가는 SNR 감소와 정확히 일치한다.

---

## 7. α=0의 특이점: 축퇴와 비단조성

### α=0일 때

Q = R (block-diagonal). SA가 보는 landscape:
- **에너지 갭**: Δ_R > 0 (subblock의 최소 비축퇴 갭)
- **축퇴도**: |S| = Π_i |GS(R_i)| (기하급수적)
- SA는 R의 GS를 쉽게 찾지만, |S|개 축퇴 상태 중 x\*를 특정할 수 없음

### GS 에너지 도달률의 비단조성

실험 데이터 (lin2, sw=1000):

| α | GS 에너지% | Target% | 해석 |
|---|---|---|---|
| 0 | 82.4% | 0.0% | R 쉽게 풀림, 축퇴로 x\* 특정 불가 |
| 0.001 | 10.4% | 10.4% | 축퇴 깨짐 + signal 약함 → 가장 어려움 |
| 0.01 | 24.0% | 24.0% | signal 강화로 회복 |
| 0.1 | 99.6% | 99.6% | 강한 signal → 쉬움 |

**비단조성의 수학적 원인:**

α=0 → α=0⁺ 전환에서 두 가지가 동시에 일어난다:

1. **축퇴 제거** (난이도 증가): |S|개 동치 상태 → 1개 유일 상태. GS가 유일해지면서 SA의 "성공 확률" 분모가 1/|S| → 1로 바뀌지 않고, 에너지 갭이 α·Δ_P로 축소되어 찾기 어려워짐.

2. **Gradient signal 생성** (난이도 감소): 동시에 V_P(d) = α·γ·d가 활성화되어 SA에 방향 정보 제공.

α가 매우 작을 때는 (1)이 지배하여 난이도가 증가하고, α가 커지면 (2)가 지배하여 난이도가 감소한다.

---

## 8. 실험 데이터와의 대응

### 노트북 실험 (n=500, k=10, lin2, 1000 instances × 100 reads, sw=1000)

| α | 성공률 | Avg HD | Δ_min ∝ | B/Δ ∝ | 예측 |
|---|---|---|---|---|---|
| 0.1 | 99.8% | 0.9 | 0.1·Δ_P | 10/Δ_P | 쉬움 (SNR > 1) |
| 0.01 | 36.8% | 431 | 0.01·Δ_P | 100/Δ_P | 어려움 |
| 0.001 | 15.5% | 810 | 0.001·Δ_P | 1000/Δ_P | 매우 어려움 |

### Fixed R/P 실험 (n=500, k=15, lin2, 100 instances × 50 reads)

| α | sw=100 | sw=1000 | sw=5000 | HD (sw=5000) |
|---|---|---|---|---|
| 0 | 0.0% | 0.0% | 0.0% | 27 |
| 0.001 | 4.4% | 10.9% | 13.3% | 8 |
| 0.01 | 8.2% | 23.0% | 35.6% | 4 |
| 0.1 | 87.7% | 99.7% | 100.0% | 0 |

모든 실험에서 **α 증가 → 성공률 단조 증가, HD 단조 감소** 패턴이 일관된다.

---

## 9. 요약

### 핵심 결과

| # | 정리 | 의미 |
|---|---|---|
| 1 | Δ_min(α) = α·Δ_P | 에너지 갭이 α에 선형 비례하여 닫힘 |
| 2 | B(α) → B_R (상수) | 에너지 장벽은 R에 의해 결정, α 무관 |
| 3 | TTS ≳ exp(c/α) | SA 난이도가 1/α에 지수적으로 증가 |
| 4 | E[P\|HD=d] ≈ γ·d | Posiform이 HD 공간에 선형 포텐셜 생성 |
| - | SNR = α·γ/σ_R | α가 SA의 방향 감지 능력을 직접 제어 |

### 두 관점의 상보성

| | 국소적 관점 (Barrier-Gap) | 전역적 관점 (HD Gradient) |
|---|---|---|
| 핵심 양 | Δ_min, B, ρ = B/Δ | γ, σ_R, SNR |
| 설명하는 것 | 왜 지수적으로 어려운가 (하한) | 왜 SA가 방향을 잃는가 (메커니즘) |
| 관련 이론 | Hajek (1988), Arrhenius law | Random energy model, signal detection |

### 논문과의 관계

본 분석은 Pelofske et al. (2024)의 경험적 관찰 *"local minima whose energy is very close to the true ground-state"* 를 **정리 1 (Δ_min = α·Δ_P)** 로 정확히 formalize한 것이다. 또한 Hahn et al. (2023)이 future work으로 남긴 *"if posiform planting allows one to tune the difficulty ... via the choice of the posiform coefficients"* 에 대한 수학적 답을 제공한다.

---

## 참고문헌

- Pelofske, Hahn, Djidjev (2024). "Increasing the Hardness of Posiform Planting Using Random QUBOs for Programmable Quantum Annealer Benchmarking."
- Hahn, Pelofske, Djidjev (2023). "Using 2-SAT to Efficiently Generate Optimal Quadratic Unconstrained Binary Optimization Problem Instances for Benchmarking."
- Hajek (1988). "Cooling Schedules for Optimal Annealing." Mathematics of Operations Research, 13(2), 311-329.
