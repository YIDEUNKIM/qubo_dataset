# Quiet Planting (Planted 3-SAT → QUBO) — 실험 보고서

## 1. 배경 및 동기

### 기존 방법론의 한계

| 생성기 | SA-Hard? | 구별 불가능? | 비고 |
|--------|----------|-------------|------|
| Zero-Expectation | SA 100% 성공 | E[q_ij]=0 보장 | 통계적으로 은닉되지만 너무 쉬움 |
| Hard Mode | SA 100% 성공 | backbone 구조 노출 | W_STRONG/W_WEAK 비율로 탐지 가능 |
| Wishart | SA 실패 (alpha~0.7) | low-rank 구조 노출 | J = -(1/N)WW^T의 rank = M < N |

### 목표

**Quiet Planting** (Krzakala & Zdeborova, 2009)은 planted random 3-SAT에서 유래한 방법으로, 3-SAT의 상전이 현상을 활용하여 SA-hard 문제를 생성하면서 랜덤 3-SAT와 통계적으로 구별 불가능한 인스턴스를 만든다.

---

## 2. 알고리즘

### 2.1 Planted 3-SAT 생성

목표 비트스트링 `target`이 주어지면:
1. 3개 변수를 무작위 선택
2. 8가지 violating assignment 중 target과 다른 7가지에서 균일 샘플
3. m = alpha * n개 clause 반복 생성

이렇게 생성된 3-SAT 인스턴스는:
- target이 반드시 만족 (by construction)
- alpha < 3.86이면 랜덤 3-SAT와 통계적으로 구별 불가능

### 2.2 Rosenberg Reduction (3-SAT → QUBO)

각 clause의 위반 페널티는 cubic: `z1 * z2 * z3`

Rosenberg reduction으로 보조변수 y를 도입하여 quadratic화:

```
g(x, y) = y*z3 + z1*z2 - 2*z1*y - 2*z2*y + 3*y
```

여기서 `zi = xi` (vi=1) 또는 `(1-xi)` (vi=0)

**증명**: y = z1*z2일 때 g = z1*z2*z3 (원래 페널티). y != z1*z2이면 추가 penalty >= 1.

QUBO 크기 = n + m = n(1 + alpha). alpha=4.2이면 QUBO는 5.2n개 변수.

### 2.3 3-SAT 상전이

- **alpha < 3.86**: condensation threshold 이하. quiet planting 보장
- **alpha ≈ 4.27**: satisfiability threshold (SAT/UNSAT 경계)
- **alpha > 4.27**: 랜덤 3-SAT는 UNSAT이지만 planted는 여전히 SAT

---

## 3. 실험 결과

### 3.1 기본 구현 (field 없음) — 축퇴 문제 발견

**설정**: N=50, alpha sweep, num_runs=10, num_reads=200

| Alpha | EXACT | ENERGY_MATCH | FAIL | Avg Hamming |
|-------|-------|--------------|------|-------------|
| 2.0   | 0%    | 100%         | 0%   | 19.7        |
| 3.0   | 0%    | 100%         | 0%   | 18.2        |
| 3.5   | 0%    | 100%         | 0%   | 13.9        |
| 3.86  | 0%    | 100%         | 0%   | 15.4        |
| 4.0   | 0%    | 100%         | 0%   | 12.2        |
| 4.2   | 0%    | 100%         | 0%   | 13.6        |
| 4.5   | 0%    | 100%         | 0%   | 10.3        |
| 5.0   | 0%    | 100%         | 0%   | 4.4         |

**관찰**: SA가 모든 alpha에서 ground state 에너지를 완벽히 찾지만 (ENERGY_MATCH 100%), planted target 자체는 전혀 못 찾음 (EXACT 0%).

**원인 분석 — 대규모 축퇴(degeneracy)**:
- 3-SAT에는 exponentially many SAT 해가 존재 (특히 alpha < 4.27)
- 모든 SAT 해는 모든 clause를 만족 → penalty = 0 → **동일한 QUBO 에너지**
- SA는 아무 SAT 해나 찾으면 되므로, planted target이 아닌 다른 해를 반환
- alpha 증가 시 SAT 해 수 감소 → Hamming distance 감소 (19.7 → 4.4)

### 3.2 축퇴 해소 시도

#### Clause 가중치만 (실패)

각 clause에 랜덤 가중치 w_k ~ Uniform(1, 3) 적용.

| 설정 | EXACT | ENERGY_MATCH | Avg Hamming |
|------|-------|--------------|-------------|
| baseline (w=1) | 0/10 | 10/10 | 13.7 |
| **clause_weight (1-3)** | **0/10** | **10/10** | **11.3** |

**실패 원인**: 모든 SAT 해가 모든 clause를 만족하므로, 가중치가 뭐든 0 × w = 0. clause 가중치는 SAT 해 간 에너지 차이를 만들 수 없음.

#### Planted Field (성공)

각 변수 i에 target 방향으로 미세한 선형 편향 추가:
- target_i = 1 → Q[(i,i)] -= epsilon_i (x_i=1 선호)
- target_i = 0 → Q[(i,i)] += epsilon_i (x_i=0 선호)
- epsilon_i ~ Uniform(field/2, field*1.5)

| 설정 | EXACT | ENERGY_MATCH | Avg Hamming |
|------|-------|--------------|-------------|
| field=0.1 | 3/10 | 0/10 | 0.9 |
| **field=0.5** | **10/10** | **0/10** | **0.0** |
| field=1.0 | 10/10 | 0/10 | 0.0 |

**원리**: 모든 SAT 해가 clause penalty = 0이지만, target만이 field penalty도 0. 다른 SAT 해는 target과 다른 비트마다 +epsilon_i penalty를 받아 에너지가 높아짐.

### 3.3 Field 강도 vs 문제 크기 (alpha=4.2)

| N \ Field | 0.01 | 0.05 | 0.1  | 0.2  | 0.5   |
|-----------|------|------|------|------|-------|
| **50**    | 0%   | 0%   | 50%  | 90%  | 100%  |
| **100**   | 0%   | 0%   | 0%   | 0%   | 100%  |
| **200**   | 0%   | 0%   | 0%   | 0%   | 90%   |

**핵심 관찰**:
1. N이 커질수록 같은 field에서 성공률 급락 (N=50 field=0.1 → 50%, N=100 → 0%)
2. field=0.5에서도 N=200부터 실패 시작 (90%)
3. QUBO 크기가 n(1+alpha) = 5.2n이므로 변수 수가 빠르게 증가

### 3.4 대규모 N 스케일링 (field=0.5, alpha=4.2)

field=0.5로 고정하고 N을 100에서 1000까지 키운 실험:

| N | QUBO 크기 | EXACT | ENERGY | FAIL | Avg Hamming | Avg Time |
|---|-----------|-------|--------|------|-------------|----------|
| 100 | 520 | **10/10** | 0/10 | 0/10 | 0.0 | 5s |
| 300 | 1560 | **7/10** | 0/10 | 3/10 | 0.3 | 48s |
| 500 | 2600 | **2/10** | 0/10 | 8/10 | 1.6 | ~2.5h |
| 750 | 3900 | **0/10** | 0/10 | 10/10 | 2.6 | ~2.9h |
| 1000 | 5200 | **0/10** | 0/10 | 10/10 | 4.5 | 15min |

**핵심 발견**:

1. **SA 상전이 확인**: N=100(100%) → N=300(70%) → N=500(20%) → N=750+(0%). field=0.5에서 SA가 완전히 실패하는 임계 크기가 존재함.
2. **에너지도 못 찾음**: ENERGY_MATCH = 0. field가 축퇴를 깬 후에는 ground state가 유일해져서, SA가 에너지 근처에도 도달 못 함. Wishart와 유사한 "진짜 SA-hard" 패턴.
3. **Glassy landscape 증거**: Hamming distance가 2.6~4.5로 작음. SA가 target 근처까지는 접근하지만 마지막 몇 비트를 뒤집지 못하는 전형적 glassy 거동.
4. **시간 비용**: N=500~750에서 SA run당 ~3시간 소요 (num_sweeps = 10 * QUBO_size). N=1000은 생성+탐색 자체가 빨라짐 (탈출 포기?).

---

## 4. Planted Field의 한계와 트레이드오프

### 은닉성 vs 유일성

- **field가 크면**: target이 유일한 ground state, SA가 찾기 쉬움 → 벤치마크로 무의미
- **field가 작으면**: 축퇴가 깨지지 않아 target 식별 불가 → 정답 검증 불가
- **sweet spot**: field가 축퇴는 깨되, N이 크면 SA가 못 찾는 구간

### 정보 누출

Planted field는 대각 성분의 부호를 통해 target 정보를 누출:
- Q[(i,i)] < 0이면 target_i = 1일 가능성 높음
- 단, clause penalty가 대각에 기여하는 양이 O(alpha)이므로
- field << O(alpha)이면 노이즈에 묻힘

---

## 5. 다른 생성기와의 비교

| 특성 | Zero-Exp | Hard Mode | Wishart | Quiet Planting |
|------|----------|-----------|---------|----------------|
| SA 난이도 | 쉬움 | 쉬움 | 어려움 (alpha~0.7) | field=0.5: N<300 쉬움, N>500 어려움 |
| 구별 불가능성 | E[q_ij]=0 | 구조 노출 | low-rank 노출 | alpha < 3.86 보장 |
| QUBO 크기 | n(n-1)/2 | sparse | n^2 | n(1+alpha) |
| Ground state 유일성 | 유일 (확률적) | 유일 | Z2 대칭 | 축퇴 (field 필요) |
| 난이도 제어 | density | noise_ratio | alpha | alpha + field |

---

## 6. 사용법

```bash
# 기본 생성 (field 없음 — 축퇴 있음)
python3 qubo_quiet_planted.py 10110 4.2

# planted field로 축퇴 해소
python3 -c "
from qubo_quiet_planted import create_qubo_quiet_planted
Q, clauses, info = create_qubo_quiet_planted('10110', alpha=4.2, field_strength=0.5)
print(f'Target energy: {info[\"target_energy\"]:.4f}')
print(f'QUBO size: {info[\"total_vars\"]}')
"

# alpha sweep 실험
python3 test_quiet_planted.py 50 10

# 4-way 비교
python3 test_quiet_planted.py --compare
```

---

## 7. 참고문헌

- Krzakala, F. & Zdeborova, L. (2009). "Hiding Quiet Solutions in Random Constraint Satisfaction Problems." *Physical Review Letters*, 102(23), 238701.
- Rosenberg, I. G. (1975). "Reduction of bivalent maximization to the quadratic case." *Cahiers du Centre d'Etudes de Recherche Operationnelle*, 17, 71-74.
