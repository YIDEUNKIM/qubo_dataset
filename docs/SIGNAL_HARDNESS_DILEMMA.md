# Signal-Hardness Dilemma: Planted QUBO 벤치마크의 구조적 한계

## 개요

Planted QUBO 벤치마크에서 **답을 유일하게 만드는 장치는 필연적으로 솔버에게 gradient signal을 제공**한다. 이 "Signal-Hardness Dilemma"가 특정 방법론의 특이성이 아닌 planted problem 일반의 구조적 한계임을 네 차례의 실험으로 증명한다.

### 실험 일람

| # | 실험 | 방법론 | 핵심 발견 |
|---|---|---|---|
| 1 | α=0 실험 | Hardened Posiform | SA가 GS 에너지를 79% 찾지만 target은 0% → 축퇴 |
| 2 | α<0 실험 | Hardened Posiform | Hamming이 단조 감소 → posiform은 방향성 signal |
| 3 | α=0.001 실험 | Hardened Posiform | 비단조 GS 도달률 (80%→9%→19%→99%) → posiform의 이중 효과 |
| 4 | Field signal 실험 | Quiet Planting | 4가지 현상 모두 재현 → 방법론에 무관한 일반 원리 |

---

## 1. 통합 프레임워크

모든 planted QUBO 방법은 다음 구조를 공유한다:

```
Q = Q_base + s × Q_signal
```

| 항 | Hardened Posiform | Quiet Planting |
|---|---|---|
| Q_base | Σ R_i (random discrete QUBO) | SAT penalty QUBO (Rosenberg) |
| Q_signal | P (posiform planted QUBO) | planted field (linear bias) |
| s (signal parameter) | α (posiform scale) | field_strength |

이 구조에서:
- Q_base만으로도 target은 GS이지만 **유일하지 않다** (축퇴)
- Q_signal이 축퇴를 깨서 target을 유일한 GS로 만든다
- **그러나 Q_signal은 동시에 솔버에게 target 방향의 gradient를 제공한다**

---

## 2. 실험 1: α=0 — 축퇴의 발견 (Hardened Posiform)

**출처**: `hardened_posiform/ALPHA_ZERO_EXPERIMENT.md`

### 배경

Pelofske et al. (2024)의 Hardened Posiform: `Q_final = Σ R_i + α × P`. α를 줄이면 SA 성공률이 낮아지고, 논문은 이를 "난이도 증가"로 해석한다. α=0이면 posiform 기여가 완전히 사라진다.

### 실험 조건

- N=500, subgraph_size=15, coeff_type: lin2/lin20
- α: 0, 0.01, 0.1
- SA: sweep [1, 5, 10, 50, 100, 500, 1000, 5000], 인스턴스 10개, read 50개

### 핵심 결과: 에너지 분석 (sweep=1000)

| Config | GS 에너지 일치 | Target 일치 | 차이 |
|---|---|---|---|
| **lin2, α=0** | **79.0%** | **0.0%** | **395개** |
| lin2, α=0.01 | 21.0% | 21.0% | 0개 |
| lin2, α=0.1 | 99.8% | 99.8% | 0개 |
| **lin20, α=0** | **8.2%** | **0.2%** | **40개** |
| lin20, α=0.01 | 49.0% | 49.0% | 0개 |
| lin20, α=0.1 | 100% | 100% | 0개 |

### 발견

**α=0에서 SA는 GS 에너지를 79% 찾고 있다. 하지만 그 해가 target이 아니다.**

- 축퇴도: lin2 평균 ~4×10^14, lin20 평균 ~209
- α=0의 "성공률 0%"는 landscape 복잡성이 아닌 **축퇴(degeneracy)에 의한 실패**
- α>0에서는 GS 에너지 일치 = target 일치가 **정확히 동일** → posiform이 축퇴를 완전히 깸

### 함의

논문의 벤치마크 메트릭(target sampling rate)은 **genuine hardness와 degeneracy를 구분하지 못한다.**

---

## 3. 실험 2: α<0 — 방향성 signal의 증명 (Hardened Posiform)

**출처**: `hardened_posiform/ALPHA_NEGATIVE_EXPERIMENT.md`

### 동기

α=0 실험만으로는 posiform의 역할에 대해 두 가지 해석이 가능:

| 해석 | α=0 호환 | α<0 예측 |
|---|---|---|
| (A) 축퇴 제거 장치 | O | Hamming ≈ α=0 (방향이 없으므로) |
| (B) 방향성 signal | O | Hamming > α=0 (**반대 방향으로 작용**) |

### 실험 조건

- N=500, α: **-0.1, -0.01, 0, 0.01, 0.1**
- 나머지 α=0 실험과 동일

### 핵심 결과 1: Local Minimum 검사

| α | lin2 min_flip_delta | lin20 min_flip_delta | local_min? |
|---|---|---|---|
| -0.1 | -1.0 ~ -1.3 | -0.6 ~ -1.5 | **X (전부)** |
| -0.01 | -0.10 ~ -0.13 | -0.01 ~ -0.12 | **X (전부)** |
| 0 | 0.0000 | 0.0000 | O |
| 0.01 | +0.01 ~ +0.04 | +0.01 ~ +0.13 | O |
| 0.1 | +0.1 ~ +0.4 | +0.1 ~ +0.5 | O |

**α<0에서 target은 local minimum조차 아니다.** Pelofske et al. 논문의 uniqueness 증명이 α>0을 암묵적으로 전제함을 실험적으로 확인.

### 핵심 결과 2: Hamming Distance (sw=5000)

| α | lin2 med | lin20 med | 해석 |
|---|---|---|---|
| **-0.1** | **143** | **180** | anti-signal → target에서 적극적으로 멀어짐 |
| **-0.01** | **54** | **34** | 약한 anti-signal |
| 0 | 28 | 10 | signal 없음 (축퇴) |
| 0.01 | 4 | 0 | 약한 signal → target 근처 |
| **0.1** | **0** | **0** | 강한 signal → target 정확히 도달 |

**α가 -0.1에서 +0.1로 증가할수록 Hamming distance가 단조 감소.** 해석 (B)가 맞다: posiform은 방향성 gradient signal이다.

### 핵심 결과 3: 에너지 분해 (sweep=1000)

SA가 찾은 해의 에너지를 `Σ R_i` 성분과 `P` 성분으로 분리:

| Config | R gap | P gap | Q(SA)-Q(target) |
|---|---|---|---|
| lin2, α=-0.1 | +28.2 | +809 | **-18.73** |
| lin2, α=-0.01 | +1.3 | +349 | **-2.20** |

α=-0.1 해석: R 기준 +28.2 (target보다 나쁨), P 기준 +809 (target에서 멀어짐). 그러나 αP gap = -0.1 × 809 = **-80.9**이므로 전체 에너지는 target보다 낮다. SA는 posiform signal을 **읽고 반대로 따라가고 있다.**

### 결론

```
Q = Σ R_i + αP
      ↑         ↑
  target 방향    α>0: target 방향 (협력)
  (항상)         α<0: 반-target 방향 (경쟁)
```

"축퇴 제거 장치"라면 α의 부호에 따라 방향이 바뀌지 않는다. **방향 반전은 방향성 signal만이 가능한 행동이다.**

---

## 4. 실험 3: α=0.001 — Landscape 변형의 비단조성 (Hardened Posiform)

**출처**: `hardened_posiform/ALPHA_NEGATIVE_EXPERIMENT.md` (실험 5)

### 동기

α=0 (축퇴 지배) → α=0.01 (genuine hardness) 사이의 전환점 탐색. α=0.001이라는 극히 작은 signal이 축퇴를 깨기에 충분한가?

### 핵심 결과: 에너지 분석 (sweep=1000)

| Config | GS에너지% | Target% | 차이 |
|---|---|---|---|
| lin2, α=0 | **80.2%** | 0.0% | 401개 |
| **lin2, α=0.001** | **8.6%** | **8.6%** | **0개** |
| lin2, α=0.01 | 18.8% | 18.8% | 0개 |
| lin2, α=0.1 | 99.4% | 99.4% | 0개 |
| lin20, α=0 | **8.4%** | 0.4% | 40개 |
| **lin20, α=0.001** | **7.2%** | **7.2%** | **0개** |
| lin20, α=0.01 | 47.0% | 47.0% | 0개 |
| lin20, α=0.1 | 100% | 100% | 0개 |

### 핵심 발견: 축퇴 제거의 임계점

**α=0.001에서 이미 "GS 에너지 일치 = Target 일치"가 완벽하게 성립한다.** 0.001이라는 극히 작은 posiform만으로도 축퇴가 완전히 깨진다.

### 핵심 발견: GS 에너지 도달률의 비단조적 패턴

| α | GS 에너지 도달률 (lin2) | 해석 |
|---|---|---|
| 0 | **80.2%** | landscape 평탄 → 어디서든 GS 에너지에 도달 (but 축퇴) |
| 0.001 | **8.6%** | posiform이 plateau를 깨서 유일한 골짜기 생성 → 찾기 어려움 |
| 0.01 | 18.8% | signal 강화로 gradient 생겨 오히려 회복 |
| 0.1 | 99.4% | 강한 signal로 쉽게 도달 |

**80% → 9% → 19% → 99%** 의 비단조 패턴은 posiform의 이중 효과를 극적으로 보여준다:

1. **축퇴 제거**: α=0 → 0.001에서 plateau가 사라지면서 유일한 GS 생성
2. **Landscape 변형**: 동시에 평탄했던 에너지 landscape에 골짜기와 장벽이 생성

α가 아주 작을 때 (0.001) 축퇴는 깨졌지만 signal이 너무 약해서 SA가 유일한 골짜기를 찾기 어려움. α를 키우면 gradient signal이 강해져서 다시 쉬워진다.

---

## 5. 실험 4: Quiet Planting Field Signal — 일반화 (Quiet Planting)

**출처**: `quiet_planting/FIELD_SIGNAL_EXPERIMENT.md`

### 동기

Hardened Posiform에서의 발견이 방법론 특이적이 아닌 planted problem 일반의 현상인지 검증. Quiet Planting은 기본 문제(SAT penalty), 유일성 장치(planted field), 차수축소(Rosenberg) 등 모든 세부 사항이 Hardened Posiform과 다르다.

### 실험 조건

- N=100, alpha=4.2 (QUBO size=520)
- field_strength: **-1.0, -0.1, 0, 0.1, 0.5, 1.0**
- 모든 field 조건이 **동일한 base QUBO** 공유

### 축퇴의 원인

Quiet Planting에서 field=0일 때: 모든 SAT 해가 모든 clause를 만족하므로 penalty=0. clause 가중치를 다르게 해도 (clause_weight_range) 모든 SAT 해에 대해 `Σ w_k × 0 = 0`이므로 축퇴를 깰 수 없다. planted field만이 선형 편향을 추가하여 target을 유일한 GS로 만든다.

### 핵심 결과 1: Hamming Distance (sw=5000)

| field | Hamming med | 비고 |
|---|---|---|
| **-1.0** | **81** | anti-signal → target에서 적극적으로 멀어짐 |
| **-0.1** | **42** | 약한 anti-signal |
| 0 | 30 | signal 없음 (축퇴) |
| 0.1 | 18 | 약한 signal (불충분) |
| 0.5 | 3 | signal → target 근처 |
| **1.0** | **0** | 강한 signal → target 도달 |

**Hardened Posiform의 α 실험과 정확히 동일한 단조 감소 패턴.**

### 핵심 결과 2: 에너지 분석 (sweep=1000)

| field | GS에너지% | Target% | GS찾았지만 target 아닌 수 |
|---|---|---|---|
| -1.0 | 0.0% | 0.0% | 0 |
| -0.1 | 0.0% | 0.0% | 0 |
| **0** | **2.0%** | **0.0%** | **10** |
| 0.1 | 0.0% | 0.0% | 0 |
| 0.5 | 2.8% | 2.8% | 0 |
| 1.0 | 42.0% | 42.0% | 0 |

**field=0에서 "GS 에너지를 찾았지만 target이 아닌" 패턴이 재현된다.** Hardened Posiform α=0 (79% vs 0%)과 질적으로 동일한 현상.

### 네 가지 현상의 완벽한 대응

| 현상 | Hardened Posiform | Quiet Planting |
|---|---|---|
| signal=0 → 축퇴 | α=0: GS 80% but target 0% | field=0: GS 2% but target 0% |
| signal<0 → anti-signal | α<0: Hamming 증가 | field<0: Hamming 증가 |
| signal>0 → GS=target | α>0: GS%=target% | field>0: GS%=target% |
| Hamming 단조 감소 | α↑ → Hamming↓ | field↑ → Hamming↓ |

**네 가지 현상이 모두 동일하게 재현된다.** 이는 방법론에 무관한 일반적 원리이다.

---

## 6. Signal-Hardness Dilemma의 정식화

### 딜레마의 구조

모든 planted QUBO 벤치마크는 다음 딜레마에 직면한다:

```
Signal 증가  ──→  유일성 보장 ✓ + 솔버에게 gradient 제공 → 문제가 쉬워짐
Signal 감소  ──→  gradient 약화 → "어려워 보임" but 실제로는 축퇴 증가
Signal = 0   ──→  유일성 상실 → 벤치마크 무효화
Signal < 0   ──→  anti-signal → target에서 적극적으로 멀어짐
```

이를 수식으로 표현하면:

```
Q = Q_base + s × Q_signal     (s: signal parameter)

s > 0:  target = unique GS, SA에게 gradient signal 제공 (쉬움)
s → 0:  target ∈ degenerate GS, signal 소멸 (축퇴)
s = 0:  target ∈ exponentially many degenerate GS (벤치마크 무효)
s < 0:  target ≠ GS, anti-signal로 SA가 target에서 멀어짐
```

### Posiform의 이중 역할 (Dual Role)

Pelofske et al. (2024) 논문은 posiform의 역할을 "GS 유일성 보장"으로만 서술한다. 그러나 실험 결과는 두 가지 역할이 동시에 존재함을 보여준다:

| 역할 | 설명 | 증거 |
|---|---|---|
| (a) 축퇴 제거 | target을 유일한 GS로 만듦 | α=0.001에서 이미 완벽한 축퇴 제거 |
| (b) Gradient signal | SA에게 target 방향의 에너지 경사 제공 | α<0에서 방향 반전, Hamming 단조 감소 |

**역할 (b)는 논문에서 간과되었으며, 이것이 벤치마크의 난이도를 근본적으로 훼손한다.**

### 벤치마크 메트릭의 한계

기존 메트릭 "target sampling rate"의 문제:

| α | target sampling rate | 실체 |
|---|---|---|
| 0 | 0% (가장 어려움?) | SA가 GS 에너지 80% 찾음 — **축퇴**, not hardness |
| 0.001 | 9% | genuine hardness (축퇴 깨졌지만 signal 약함) |
| 0.01 | 19% | genuine hardness |
| 0.1 | 100% (가장 쉬움) | 강한 signal |

α=0의 0%와 α=0.001의 9%는 질적으로 완전히 다른 현상이지만, target sampling rate만으로는 구분할 수 없다. **에너지 기반 메트릭** (GS 에너지 도달률, GS 에너지 일치 vs target 일치 비교)이 필요하다.

---

## 7. Pelofske et al. (2024)에 대한 함의

### 원 논문의 주장과 본 연구의 반론

| 논문의 주장 | 본 연구의 반론 |
|---|---|
| α를 줄이면 난이도가 증가한다 | α를 줄이면 signal이 약해질 뿐, α=0에서는 축퇴 |
| Posiform은 GS 유일성 보장 장치이다 | Posiform은 동시에 gradient signal을 제공한다 |
| 작은 α에서 SA 실패 = 높은 난이도 | α=0에서 SA 실패 = 축퇴에 의한 실패 (에너지 최적화는 성공) |
| Uniqueness 증명은 일반적이다 | 증명은 α>0을 암묵적으로 전제한다 (α<0에서 target이 local min 아님) |

### 논문 증명의 암묵적 전제

Section 2.2의 uniqueness 증명:
```
Q(X*) = Σ R_i(X*) + αP(X*) < Σ R_i(X̂) + αP(X̂) = Q(X̂)
```

P(X\*) < P(X̂)에 α를 곱할 때 부등호 방향이 유지되려면 **α > 0**이어야 한다. 논문은 이를 명시하지 않지만, α<0 실험에서 target이 local minimum조차 아님이 확인되어 이 전제의 필요성이 실험적으로 입증되었다.

---

## 8. 요약: 실험별 증명 구조

```
실험 1 (α=0)
  → 발견: SA가 GS 에너지를 찾지만 target은 못 찾음
  → 증명: "성공률 0%"의 실체는 축퇴이지 landscape 복잡성이 아니다
  → 시사: Posiform이 축퇴 제거 이상의 역할을 한다

실험 2 (α<0)
  → 발견: Hamming이 α에 대해 단조 감소, α<0에서 target은 local min 아님
  → 증명: Posiform은 "축퇴 제거 장치"가 아닌 "방향성 gradient signal"이다
  → 시사: signal의 부호에 따라 방향이 반전됨 = 벡터적 성질

실험 3 (α=0.001)
  → 발견: α=0.001에서 축퇴 완전 제거, 그러나 GS 도달률 80%→9% 급락
  → 증명: Posiform의 이중 효과 — 축퇴 제거 + landscape 변형이 동시 발생
  → 시사: 비단조적 패턴이 딜레마의 정량적 구조를 드러냄

실험 4 (Quiet Planting field signal)
  → 발견: Hardened Posiform의 4가지 현상이 모두 재현
  → 증명: Signal-Hardness Dilemma는 특정 방법론이 아닌 planted problem 일반의 구조적 한계
  → 시사: Q = Q_base + s × Q_signal 구조를 공유하는 모든 방법에 적용
```

---

## 9. 논문 방향

### 제안: "The Signal-Hardness Dilemma in Planted QUBO Benchmarks"

| 섹션 | 내용 | 핵심 실험 |
|---|---|---|
| Introduction | Planted QUBO 벤치마크의 필요성과 기존 방법의 한계 | — |
| Background | Hardened Posiform, Quiet Planting 방법론 소개 | — |
| The α=0 Anomaly | "가장 어려운" 설정이 실제로는 축퇴일 뿐 | 실험 1 |
| Posiform as Directional Signal | α<0 실험으로 방향성 증명 | 실험 2 |
| Landscape Transformation | α=0.001의 비단조 패턴으로 이중 효과 규명 | 실험 3 |
| Generalization | Quiet Planting에서의 재현으로 일반성 증명 | 실험 4 |
| The Dilemma | Q = Q_base + s × Q_signal 구조의 근본적 한계 정식화 | 통합 |
| Implications | 벤치마크 메트릭 개선 제안, 향후 연구 방향 | — |

### 핵심 메시지

> Planted QUBO 벤치마크에서 답을 유일하게 만드는 장치(posiform/planted field)는 **단순한 축퇴 제거가 아니라 솔버에게 방향성 gradient signal을 제공**한다. 이 signal은 파라미터의 부호에 따라 방향이 반전되며, 크기에 따라 강도가 조절된다. Signal 없이 유일성을 보장할 수 없고, signal이 있으면 문제가 쉬워진다 — 이것이 planted QUBO 벤치마크의 근본적 딜레마이다.

---

## 실험 파일

### Hardened Posiform
- `hardened_posiform/experiment_alpha_zero.py` — α=0 sweep 전이
- `hardened_posiform/experiment_alpha_zero_energy.py` — α=0 에너지 분석
- `hardened_posiform/experiment_alpha_negative.py` — α<0 실험 (sweep 전이 + 에너지 분해)
- `hardened_posiform/experiment_alpha_fine.py` — α=0.001 fine-grained 분석
- `hardened_posiform/ALPHA_ZERO_EXPERIMENT.md` — α=0 결과 문서
- `hardened_posiform/ALPHA_NEGATIVE_EXPERIMENT.md` — α<0 및 α=0.001 결과 문서

### Quiet Planting
- `quiet_planting/experiment_field_signal.py` — Field signal 실험
- `quiet_planting/FIELD_SIGNAL_EXPERIMENT.md` — Field signal 결과 문서
