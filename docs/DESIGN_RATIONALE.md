# 설계 근거 및 개념 정리

## 1. 왜 Ising 기반으로 설계하는가?

### QUBO 직접 설계의 어려움

QUBO 에너지: $E(x) = \sum Q_{ii}x_i + \sum_{i<j} Q_{ij}x_ix_j$, $x_i \in \{0, 1\}$

이진 변수의 비대칭성 때문에 ground state 보장이 복잡하다:
- $x_i = 0$이면 관련 항이 모두 사라짐
- $x_i = 1$이면 모든 항이 활성화됨
- 대각/비대각 항의 상호작용이 비선형적으로 얽힘

### Ising 설계의 장점

Ising 에너지: $H(s) = -\sum J_{ij} s_i s_j$, $s_i \in \{-1, +1\}$

target spin $t$에서 에너지를 최소로 만드는 조건이 단순하다:

$$J_{ij} \cdot t_i \cdot t_j > 0$$

이것만 만족하면 ground state가 보장된다. 이후 $s_i = 2x_i - 1$ 치환으로 QUBO로 변환:

$$Q_{ij} = -4J_{ij}, \quad Q_{ii} += 2J_{ij}$$

### 파일별 설계 공간

| 파일 | 설계 공간 | 목적 |
|------|----------|------|
| `qubo_zero_expectation.py` | QUBO 직접 | 랜덤과 구별 불가능한 QUBO |
| `qubo_hard_mode.py` | Ising → QUBO | SA가 풀기 어려운 QUBO |
| `qubo_wishart.py` | Ising → QUBO | SA가 풀기 어려운 QUBO |
| `qubo_posiform.py` | 2-SAT → Posiform → QUBO | GS 유일성 수학적 보장 |
| `qubo_posiform_hardened.py` | Random QUBO + Posiform | GS 보장 + SA 난이도 조절 |
| `qubo_quiet_planted.py` | 3-SAT → Rosenberg → QUBO | 통계적 구별 불가능 |
| `qubo_mceliece.py` | McEliece → Ising → QUBO | 암호학적 hardness |

---

## 2. 구별 가능성 vs 계산 난이도

### 핵심: "구별 가능 ≠ 풀 수 있음"

Wishart QUBO는 Q 행렬을 분석하면 랜덤이 아님을 탐지할 수 있다 (구별 가능). 그러나 SA는 여전히 ground state를 찾지 못한다 (계산적으로 어려움). 이 두 개념은 독립적이다.

### SA가 하는 일

SA는 Q 행렬의 구조를 분석하지 않는다. 에너지 지형 위를 걸어다닐 뿐이다:

```
1. 랜덤 비트스트링에서 시작
2. 비트 하나 뒤집어봄
3. 에너지가 낮아지면 수락, 아니면 확률적으로 수락/거부
4. 2-3 반복
```

SA는 2^N개 상태의 에너지 지형을 local move로 탐색하는 것이지, Q 행렬의 고유값이나 rank를 보는 게 아니다.

### 구별 가능하다는 것

Q 행렬을 분석하면 planted 문제임을 **탐지**할 수 있다는 뜻이다:
- 고유값 분해 → rank-M 구조 노출
- 대각 항 통계 → target 비트 방향 편향 감지

그러나 "planted 문제라는 것을 아는 것"과 "정답이 무엇인지 아는 것"은 완전히 다른 문제이다.

비유:
- **구별 가능** = 이 금고에 보물이 있다는 걸 아는 것
- **풀 수 있음** = 금고의 비밀번호를 아는 것

### SA가 실패하는 이유

Wishart 에너지 지형에는 ground state 근처에 **metastable 상태**(가짜 바닥)가 많다:
- 에너지가 ground state의 ~92% 수준으로 비슷함
- 그러나 비트스트링은 ground state와 해밍거리 ≈ N/2 (완전히 다른 위치)
- SA가 metastable 상태에서 탈출하려면 수십 개의 비트를 **동시에** 뒤집어야 함
- SA는 한 번에 1개만 뒤집으므로 에너지 장벽을 넘을 수 없음

### Wishart에서 구별이 가능한 이유

Ising → QUBO 변환 시 구조가 노출된다:

```
Q_ii = 2 * Σ_j J_ij    ← target 비트에 따라 부호 편향
Q_ij = -4 * J_ij       ← J가 low-rank(rank M)이므로 상관 구조
```

공격자가 Q의 고유값 분해를 하면 rank-M 구조를 감지하고, 대각 항의 편향으로 target 비트를 추측할 수 있다.

---

## 3. 구별 불가능은 언제 필요한가?

### 필요한 경우: 공정한 블라인드 벤치마크

솔버 개발자가 **부정행위**하는 것을 방지해야 할 때:
1. Q 행렬 구조를 분석하여 "이건 Wishart planted 문제"라고 식별
2. Wishart 전용 디코딩 알고리즘으로 빠르게 풀기
3. 벤치마크에서 부당한 고득점

이를 방지하려면 문제가 랜덤 QUBO와 통계적으로 구별 불가능해야 한다 → `qubo_zero_expectation.py`

### 불필요한 경우: 순수 솔버 성능 테스트

SA, 양자 어닐링 등의 **범용 솔버 성능**만 측정할 때는 구별 가능 여부가 무관하다. 솔버가 Q 행렬의 구조 분석 없이 순수하게 에너지 지형을 탐색하기 때문이다 → `qubo_wishart.py`

---

## 4. 두 성질의 결합이 어려운 이유

| 성질 | 요구 조건 |
|------|----------|
| SA-hard | J 행렬이 low-rank 구조를 가져야 metastable 상태 생성 → Q 계수들이 **서로 상관** |
| 구별 불가능 | Q 행렬 계수가 독립적이고 무작위처럼 보여야 함 → 구조가 **없어야** 함 |

SA-hard하려면 구조가 필요하고, 구별 불가능하려면 구조가 없어야 한다. 근본적으로 충돌하는 요구사항이다.

### 가능한 타협안

| 방법 | SA-hard | 구별 불가능 | 한계 |
|------|---------|-----------|------|
| Wishart + zero-expectation noise | noise 약하면 유지 | 부분적 | noise 강하면 난이도 파괴 |
| 대각 항만 보정 | 유지 가능 | 대각만 은닉 | 비대각 상관 구조 노출 |
| Wishart + 밀도 조절 | 유지 | 부분적 | 근본적 해결 아님 |

완전한 결합은 열린 연구 문제이다.

---

## 5. SA 실패 문제를 QPU는 풀 수 있는가?

### SA가 실패하는 근본 원인

SA는 **한 번에 비트 1개만 뒤집는** local search이다. Wishart 에너지 지형의 metastable 상태에서 ground state로 이동하려면 수십 개의 비트를 동시에 뒤집어야 하는데, SA는 그 사이의 에너지 장벽을 넘을 수 없다.

### QPU(양자 어닐링)의 차이

양자 어닐링은 **양자 터널링(Quantum Tunneling)** 효과를 이용한다. 에너지 장벽을 "넘는" 것이 아니라 "뚫고 지나간다." 따라서 SA가 갇히는 metastable 상태를 QPU는 탈출할 수 있는 가능성이 있다.

| 솔버 | 장벽 극복 방법 | Wishart 문제 예상 |
|------|-------------|-----------------|
| SA (Simulated Annealing) | 열적 요동 (비트 1개 뒤집기) | 실패 (실험으로 확인) |
| QPU (양자 어닐링, D-Wave) | 양자 터널링 (장벽 관통) | **풀 수도 있음** |
| QAOA (변분 양자) | 양자 중첩 + 변분 최적화 | 회로 깊이에 따라 다름 |

### 이 벤치마크의 가치

SA는 실패하고 QPU는 성공하는 문제가 존재한다면, 이는 **양자 우위(Quantum Advantage)**의 직접적 증거가 된다.

Wishart planted ensemble이 좋은 벤치마크인 이유:
1. **Ground state가 수학적으로 보장됨** → 정답을 알고 있으므로 솔버 성능을 정확히 측정 가능
2. **SA가 실패함** → 고전 솔버의 한계가 명확
3. **난이도 조절 가능** → alpha 파라미터로 easy-hard 전이를 정밀 제어
4. **확장 가능** → N을 키우면 난이도가 증가하므로 스케일링 연구에 적합

이 벤치마크로 "QPU가 SA 대비 얼마나 더 잘 푸는가"를 정량적으로 측정할 수 있다.

---

## 6. 용도별 선택 가이드

> 6개 방법론 전체의 용도별 권장 가이드 및 SA 벤치마크 결과는 **[METHODOLOGY_COMPARISON.md](METHODOLOGY_COMPARISON.md)**를 참조.
