# Quiet Planting Field Signal 실험 결과

## 실험 목적

Hardened Posiform에서 확인한 **Signal-Hardness Dilemma**가 Quiet Planting에서도 동일하게 나타나는지 검증한다. 동일 패턴이 재현되면 이 현상이 특정 방법론의 특이성이 아닌 **planted problem 일반의 구조적 한계**임을 증명할 수 있다.

### 두 방법론의 대응 관계

| | Hardened Posiform | Quiet Planting |
|---|---|---|
| 기본 문제 | Σ R_i (random discrete QUBO) | SAT penalty QUBO (Rosenberg) |
| 유일성 장치 | α × P (posiform) | field_strength (planted field) |
| 난이도 파라미터 | α (posiform scale) | field_strength |
| 축퇴 조건 | α = 0 → subproblem 축퇴도의 곱 | field = 0 → 모든 SAT 해가 동일 에너지 |
| anti-signal | α < 0 | field < 0 |

### 축퇴의 원인

Quiet Planting에서 field=0일 때: 모든 SAT 해가 모든 clause를 만족하므로 penalty = 0. clause 가중치를 다르게 해도 (clause_weight_range) 모든 SAT 해에 대해 `Σ w_k × 0 = 0`이므로 축퇴를 깰 수 없다. planted field만이 선형 편향을 추가하여 target을 유일한 GS로 만든다.

## 실험 조건

- N = 100, alpha = 4.2 (QUBO size = 520)
- **field_strength: -1.0, -0.1, 0, 0.1, 0.5, 1.0** (6개 조건)
- SA: sweep 수 [1, 5, 10, 50, 100, 500, 1000, 5000], 인스턴스 10개, read 50개/sweep
- 에너지 분석: sweep=1000 고정
- 모든 field 조건이 **동일한 base QUBO**를 공유 (field만 변경)

## 실험 1: Target이 Local Minimum인지 검사

| field | min_flip_delta 범위 | local_min? |
|---|---|---|
| **-1.0** | 음수 | **X (전부)** |
| **-0.1** | 음수 | **X (전부)** |
| 0 | 0.0000 | O (전부) |
| 0.1 | 양수 | O (전부) |
| 0.5 | 양수 | O (전부) |
| 1.0 | 양수 | O (전부) |

Hardened Posiform과 동일: field < 0에서 target은 local minimum조차 아니다.

## 실험 2: Sweep 전이 — Target Sampling Rate (%)

| field | sw=1 | sw=5 | sw=10 | sw=50 | sw=100 | sw=500 | sw=1000 | sw=5000 |
|---|---|---|---|---|---|---|---|---|
| **-1.0** | 0% | 0% | 0% | 0% | 0% | 0% | 0% | **0%** |
| **-0.1** | 0% | 0% | 0% | 0% | 0% | 0% | 0% | **0%** |
| 0 | 0% | 0% | 0% | 0% | 0% | 0% | 0% | **0%** |
| 0.1 | 0% | 0% | 0% | 0% | 0% | 0% | 0% | **0%** |
| 0.5 | 0% | 0% | 0% | 0% | 0% | 1.2% | 1.8% | **11.2%** |
| 1.0 | 0% | 0% | 0% | 0% | 1.6% | 23.6% | 38.6% | **64.6%** |

field ≤ 0.1에서는 5000 sweeps에서도 target을 찾지 못한다. field=0.5부터 찾기 시작하고, field=1.0에서 64.6%.

## 실험 3: Hamming Distance — 핵심 결과 (sw=5000)

| field | Hamming med | Hamming avg | 비고 |
|---|---|---|---|
| **-1.0** | **81** | 80.0 | anti-signal → target에서 적극적으로 멀어짐 |
| **-0.1** | **42** | 42.2 | 약한 anti-signal |
| 0 | 30 | 30.7 | signal 없음 (축퇴) |
| 0.1 | 18 | 18.5 | 약한 signal (불충분) |
| 0.5 | 3 | 3.1 | signal → target 근처 |
| **1.0** | **0** | 0.5 | 강한 signal → target 도달 |

**field가 -1.0에서 +1.0으로 증가할수록 Hamming distance가 단조 감소한다.** Hardened Posiform의 α 실험과 정확히 동일한 패턴.

### Hardened Posiform과의 Hamming 비교

| 방향 | Hardened Posiform (N=500) | Quiet Planting (N=100) |
|---|---|---|
| anti-signal (강) | α=-0.1 → Ham=143 (N의 29%) | field=-1.0 → Ham=81 (N의 81%) |
| signal 없음 | α=0 → Ham=28 (N의 6%) | field=0 → Ham=30 (N의 30%) |
| signal (강) | α=0.1 → Ham=0 | field=1.0 → Ham=0 |

N에 대한 비율은 다르지만 (문제 구조가 다르므로), **단조 감소 패턴은 동일**하다.

## 실험 4: 에너지 분석 — Degeneracy vs Genuine Hardness (sweep=1000)

| field | GS에너지% | Target% | Hamming med | GS찾았지만 target 아닌 수 |
|---|---|---|---|---|
| **-1.0** | 0.0% | 0.0% | 81 | 0 |
| **-0.1** | 0.0% | 0.0% | 44 | 0 |
| **0** | **2.0%** | **0.0%** | 33 | **10** |
| 0.1 | 0.0% | 0.0% | 24 | 0 |
| 0.5 | 2.8% | 2.8% | 6 | 0 |
| 1.0 | 42.0% | 42.0% | 1 | 0 |

### 핵심 관찰

**field=0에서 "GS 에너지를 찾았지만 target이 아닌" 패턴이 재현된다.**

- GS 에너지 일치: 2.0% (10개) → SA가 최적 에너지에 도달
- Target 일치: 0.0% (0개) → 그 해가 target은 아님
- 차이: 10개 → 모든 SAT 해가 동일 에너지이므로 다른 SAT 해를 찾음

이는 Hardened Posiform α=0 (GS 에너지 80% vs target 0%)과 **질적으로 동일한 현상**이다. 비율의 차이 (80% vs 2%)는 문제 구조의 차이 (discrete coefficient plateau vs continuous SAT penalty)에 기인한다.

**field > 0에서 GS 에너지 일치 = Target 일치 (완벽한 대응):**

field=0.5: 2.8% = 2.8%, field=1.0: 42.0% = 42.0%. Hardened Posiform α > 0과 동일.

## 해석

### 1. Signal-Hardness Dilemma는 planted problem 일반의 구조적 한계이다

| 현상 | Hardened Posiform | Quiet Planting |
|---|---|---|
| signal=0 → 축퇴 | α=0: GS 80% but target 0% | field=0: GS 2% but target 0% |
| signal<0 → anti-signal | α<0: Hamming 증가 | field<0: Hamming 증가 |
| signal>0 → GS=target | α>0: GS%=target% | field>0: GS%=target% |
| Hamming 단조 감소 | α↑ → Hamming↓ | field↑ → Hamming↓ |

**네 가지 현상이 모두 동일하게 재현된다.** 이는 방법론에 무관한 일반적 원리이다.

### 2. 일반화된 딜레마의 구조

모든 planted QUBO 방법은 다음 구조를 공유한다:

```
Q = Q_base + s × Q_signal
```

| 항 | Hardened Posiform | Quiet Planting |
|---|---|---|
| Q_base | Σ R_i (random QUBO) | SAT penalty (Rosenberg) |
| Q_signal | P (posiform) | planted field (linear bias) |
| s | α (posiform scale) | field_strength |

이때:
- Q_base만으로는 target이 GS이지만 유일하지 않다 (축퇴)
- Q_signal이 축퇴를 깨서 target을 유일한 GS로 만든다
- **그러나 Q_signal은 동시에 솔버에게 target 방향의 gradient를 제공한다**
- s를 줄이면 gradient가 약해져 "어려워 보이지만", 실제로는 축퇴가 증가할 뿐
- s를 0으로 만들면 축퇴가 완전히 발생하여 벤치마크가 무효화된다

### 3. Wishart Planted Ensemble과의 차이점

Wishart 방법은 이 구조에서 다소 다르다:
- Q = W^T diag(λ) W 형태로, signal과 base가 명확히 분리되지 않음
- 하지만 α = M/N (projection 비율)이 유사한 역할을 수행
- α가 작으면 문제가 어려워지고, α가 크면 쉬워지는 동일 패턴

## 논문에서의 위치

이 실험은 논문의 **Section 6: The Signal-Hardness Dilemma** (일반화 섹션)에서 핵심 증거로 사용된다:

> "Hardened Posiform에서 발견한 signal-hardness dilemma는 Quiet Planting에서도 동일하게 나타난다. 두 방법론은 기본 문제 (random QUBO vs SAT penalty), 유일성 장치 (posiform vs planted field), 차수축소 여부 (없음 vs Rosenberg) 등 모든 세부 사항이 다르지만, signal parameter를 0으로 만들면 동일한 축퇴 현상이 발생하고, 음수로 만들면 동일한 anti-signal 효과가 나타난다. 이는 이 딜레마가 planted problem의 구조적 한계임을 강하게 시사한다."

## 실험 파일

- `quiet_planting/experiment_field_signal.py` — Field signal 실험 (sweep 전이 + 에너지 분석)
- `hardened_posiform/ALPHA_NEGATIVE_EXPERIMENT.md` — Hardened Posiform α 실험 (대응 실험)
- `hardened_posiform/ALPHA_ZERO_EXPERIMENT.md` — 최초 발견 (α=0 실험)
