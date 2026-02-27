# Hardened Posiform α=0 실험 결과

## 실험 목적

Pelofske et al. (2024)의 Hardened Posiform 방법에서 posiform scaling coefficient α를 0으로 설정했을 때의 동작을 분석한다.

**핵심 질문**: α=0.01이 α=0.1보다 어렵다면, α=0은 어떤가?

## 배경

Hardened Posiform의 최종 QUBO:

```
Q_final = Σ R_i + α × P
```

- `Σ R_i`: disjoint subgraph별 random discrete-coefficient QUBO의 합
- `P`: posiform planted QUBO (target이 유일한 ground state)
- `α`: posiform scaling coefficient

α=0이면 `Q_final = Σ R_i`로, posiform 기여가 완전히 사라진다.

### Ground state 보장 여부

| α 값 | target이 ground state? | target이 **유일한** ground state? |
|---|---|---|
| α > 0 | O | O (posiform이 축퇴를 깸) |
| α = 0 | **O** | **X** (subproblem 축퇴도의 곱만큼 degenerate) |

α=0에서도 target은 ground state이다. 각 subproblem R_i의 ground state를 concatenate한 것이 target이므로, `Σ R_i(target) = Σ min(R_i)` = 전체 최솟값. 그러나 각 subproblem에 축퇴(동일 에너지의 다른 해)가 있으면, 전체 축퇴도 = Π(각 subproblem 축퇴도)로 기하급수적으로 증가한다.

## 실험 조건

- N = 500, subgraph_size = 15
- coeff_type: lin2 ({-1, +1}), lin20 ({-1.0, -0.9, ..., 0.9, 1.0})
- α: 0, 0.01, 0.1
- SA: sweep 수 [1, 5, 10, 50, 100, 500, 1000, 5000], 인스턴스 10개, read 50개/sweep
- 에너지 분석: sweep=1000 고정, SA가 찾은 해의 에너지 vs target 에너지 비교

## 실험 1: Sweep 전이 (Ground-State Sampling Rate %)

### lin2

| α | sw=1 | sw=5 | sw=10 | sw=50 | sw=100 | sw=500 | sw=1000 | sw=5000 |
|---|---|---|---|---|---|---|---|---|
| **0** | 0% | 0% | 0% | 0% | 0% | 0% | 0% | **0%** |
| 0.01 | 0% | 10.4% | 4.2% | 5.2% | 5.2% | 16.6% | 21.2% | 28.4% |
| 0.1 | 0% | 15.6% | 24.8% | 66.0% | 88.8% | 99.8% | 100% | 100% |

### lin20

| α | sw=1 | sw=5 | sw=10 | sw=50 | sw=100 | sw=500 | sw=1000 | sw=5000 |
|---|---|---|---|---|---|---|---|---|
| **0** | 0% | 0% | 0% | 0% | 0% | 0% | 0.4% | **1.0%** |
| 0.01 | 0% | 0% | 0% | 2.8% | 11.0% | 37.8% | 47.2% | 76.0% |
| 0.1 | 0% | 80.6% | 84.4% | 97.2% | 99.4% | 100% | 100% | 100% |

표면적으로는 **α=0이 압도적으로 어렵다**. 5000 sweeps에서도 lin2는 0%, lin20은 1%.

## 실험 2: 에너지 분석 — "어려움"의 실체

SA가 찾은 해의 에너지가 target 에너지와 같은지 비교 (sweep=1000 고정):

| Config | GS 에너지 일치 | Target bitstring 일치 | GS 찾았지만 target 아닌 수 |
|---|---|---|---|
| **lin2, α=0** | **395/500 (79.0%)** | **0/500 (0.0%)** | **395개** |
| lin2, α=0.01 | 105/500 (21.0%) | 105/500 (21.0%) | 0개 |
| lin2, α=0.1 | 499/500 (99.8%) | 499/500 (99.8%) | 0개 |
| **lin20, α=0** | **41/500 (8.2%)** | **1/500 (0.2%)** | **40개** |
| lin20, α=0.01 | 245/500 (49.0%) | 245/500 (49.0%) | 0개 |
| lin20, α=0.1 | 500/500 (100%) | 500/500 (100%) | 0개 |

### 핵심 관찰

**α=0에서 SA는 ground state 에너지를 79%나 찾고 있다.** 하지만 그 해가 target이 아니다.

- lin2, α=0: GS 에너지 일치 시 target과의 Hamming distance 중앙값 = **28** (500비트 중)
- 축퇴도: lin2 평균 ~4×10^14, lin20 평균 ~209

**α>0에서는 GS 에너지 일치 = target 일치가 정확히 같다.** posiform이 축퇴를 완전히 깨서, ground state를 찾으면 그것이 곧 target이다.

## 해석

### 1. α=0의 "어려움"은 에너지 landscape의 복잡성이 아닌 축퇴(degeneracy)

α=0에서 SA 성공률 0%는 에너지 최적화가 어려워서가 아니다. SA는 ground state 에너지를 79% 확률로 찾는다. "실패"의 원인은 ~10^14개의 degenerate ground state 중 특정한 하나(target)를 찍어야 하는 확률 문제일 뿐이다.

### 2. Posiform의 이중 역할 (dual role)

Pelofske et al. 논문은 posiform의 역할을 "ground state 유일성 보장"으로만 서술한다. 그러나 실험 결과는 posiform이 두 가지 역할을 동시에 수행함을 보여준다:

| 역할 | 설명 |
|---|---|
| (a) 축퇴 제거 | target을 유일한 GS로 만듦 (논문이 명시적으로 서술) |
| (b) Gradient signal 제공 | SA에게 target 방향으로의 에너지 경사를 제공 (논문이 간과) |

α를 줄이면 (a)는 유지되지만 (b)의 signal이 약해진다. α=0이면 (a)도 (b)도 사라진다.

### 3. "난이도"와 "축퇴"의 혼동

논문의 벤치마크 메트릭은 "SA가 planted target을 찾는 비율"이다. 이 메트릭 하에서:

- α=0.01: SA가 에너지 최적화를 21%만 성공 → 실제로 어렵다 (landscape 복잡성)
- α=0: SA가 에너지 최적화를 79% 성공하지만 target 매칭은 0% → 어려운 게 아니라 degenerate할 뿐

α=0.01의 21%는 genuine hardness이고, α=0의 0%는 degeneracy-induced failure이다. 벤치마크로서의 의미가 근본적으로 다르다.

### 4. lin2 vs lin20 차이의 원인

| | lin2 ({-1, +1}) | lin20 ({-1, -0.9, ..., 1}) |
|---|---|---|
| α=0 축퇴도 | ~4×10^14 | ~209 |
| α=0 GS에너지 일치율 | 79.0% | 8.2% |

lin2는 계수가 {-1, +1}뿐이어서 에너지 landscape에 거대한 평탄 영역(plateau)이 형성된다. 이 평탄 영역은:
- 축퇴도를 폭발적으로 증가시킴 (target 매칭을 불가능하게 만듦)
- 동시에 에너지 최적화를 쉽게 만듦 (어디에서든 최솟값에 도달 가능)

lin20은 계수가 21가지로 다양하여 plateau가 적고, 축퇴도가 낮은 대신 에너지 landscape이 더 울퉁불퉁하여 실질적으로 더 어렵다.

## 실험 파일

- `hardened_posiform/experiment_alpha_zero.py` — Sweep 전이 비교 (α=0 vs 0.01 vs 0.1)
- `hardened_posiform/experiment_alpha_zero_energy.py` — 에너지 분석 (GS 에너지 일치 vs target 일치)
