# Quiet Planting QUBO 생성기

## 개요

**Krzakala & Zdeborova (2009)**의 quiet planting 기법을 기반으로, planted random 3-SAT 문제를 Rosenberg reduction으로 QUBO로 변환하는 방법론. Clause density alpha = m/n으로 난이도를 제어하며, QUBO 크기는 n(1+alpha) (보조변수 포함).

## 이론적 배경

### Quiet Planting이란?

일반적인 planted 문제(목표 해가 알려진 문제)는 생성 과정에서 목표 해에 대한 정보가 노출될 수 있다. **Quiet planting**은 이를 방지하는 기법:

> Planted random 3-SAT에서, 모든 clause는 목표 해를 만족하도록 생성되지만, 생성된 인스턴스는 순수 랜덤 3-SAT과 **통계적으로 구별 불가능**하다.

### 3-SAT 상전이

| alpha 범위 | 의미 | SAT 상태 |
|:----------:|:----:|:--------:|
| < 3.86 | Condensation threshold 이하 | SAT, quiet planting 보장 |
| 3.86 ~ 4.27 | Condensation 이상, SAT threshold 이하 | SAT, 어려움 |
| ~ 4.27 | **SAT/UNSAT threshold** | 경계 |
| > 4.27 | 랜덤은 UNSAT, planted는 여전히 SAT | 어려움 |

alpha < 3.86에서 planted 인스턴스는 랜덤 인스턴스와 구별 불가능 (다항 시간 알고리즘 없음).

### Rosenberg Reduction (Cubic → Quadratic)

3-SAT clause는 3개 변수의 곱(z1*z2*z3)을 포함하여 cubic. 이를 quadratic으로 변환:

**보조변수 도입**: y = z1 * z2

**Rosenberg linearization** (P=1):
```
g(x, y) = y*z3 + z1*z2 - 2*z1*y - 2*z2*y + 3*y
```

여기서 zi는 리터럴에 따라 xi 또는 (1-xi).

각 clause마다 1개의 보조변수가 추가되므로, QUBO 크기 = n + m (n: 원래 변수, m: clause 수).

### SAT 해 축퇴 문제와 Planted Field

**문제**: 모든 SAT 해가 모든 clause를 만족하므로, clause 위반 penalty = 0. 따라서 모든 SAT 해가 동일한 QUBO 에너지를 가짐 → 목표 해를 특정할 수 없음.

**해결**: **Planted field** (작은 선형 편향)를 추가:
```python
if target[i] == '1':
    Q[(i,i)] -= eps_i   # x_i=1 선호
else:
    Q[(i,i)] += eps_i   # x_i=0 선호
```

field_strength 권장값: 0.1 ~ 1.0.

## 구현 방식

### 전체 파이프라인

```
1. Planted 3-SAT 생성
   - m = alpha * n개의 clause 생성
   - 각 clause: 3개 변수 무작위 선택 → target 할당과 다른 7가지 패턴 중 하나를 violating assignment로 선택
   - 결과: target이 모든 clause를 만족

2. Rosenberg Reduction
   - 각 clause k에 보조변수 y_k = z_a * z_b 도입
   - 5개 항으로 분해: y*z_c, z_a*z_b, -2*z_a*y, -2*z_b*y, 3*y
   - 각 항을 x, (1-x) 대입으로 QUBO 형태로 변환

3. Clause 가중치 (선택)
   - clause_weight_range=(lo, hi)로 clause별 랜덤 가중치 부여
   - 주의: SAT 해 간 축퇴를 깨지 않음 (모든 SAT 해에서 penalty=0)

4. Planted Field (선택)
   - field_strength > 0이면 target 방향 선형 편향 추가
   - 이것이 실제로 축퇴를 깸

5. 검증
   - target이 모든 clause를 만족하는지 확인
   - target + 최적 보조변수에서 QUBO 에너지 계산
```

### 핵심 파라미터

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `alpha` | 4.2 | Clause density m/n. 높을수록 어려움 |
| `clause_weight_range` | None | Clause별 가중치 범위 (축퇴 깨지 않음) |
| `field_strength` | 0.0 | Planted field 강도 (축퇴 깸, 권장 0.1~1.0) |
| `seed` | None | 재현성을 위한 난수 시드 |

### 주요 함수

| 함수 | 설명 |
|------|------|
| `create_planted_3sat(target, alpha, seed)` | Planted random 3-SAT 생성 |
| `verify_3sat_solution(clauses, assignment)` | 할당이 모든 clause를 만족하는지 검증 |
| `clause_to_qubo(clause, aux_idx)` | 단일 clause → Rosenberg reduction → QUBO |
| `compute_auxiliary_values(clauses, assignment, n)` | 주어진 할당의 최적 보조변수 값 계산 |
| `create_qubo_quiet_planted(target, alpha, seed, ...)` | **메인 진입점** |
| `extract_original_solution(sample, n)` | SA 결과에서 원래 n개 변수만 추출 |

### 반환값

```python
Q, clauses, info = create_qubo_quiet_planted(target, alpha=4.2)
# Q: QUBO 딕셔너리 {(i,j): weight}
# clauses: 3-SAT clause 리스트
# info: {'n', 'm', 'total_vars', 'alpha', 'target_energy', ...}
```

## SA 난이도 특성

### N 스케일링 (alpha=4.2, field_strength=0.5)

```
N=100:  ~100% 성공
N=300:  ~70%  성공
N=500:  ~20%  성공
N=750:  ~0%   성공
```

field_strength가 감소할수록 SA가 더 어려워지며, field=0이면 축퇴로 인해 target 복원 불가.

### QUBO 크기 영향

QUBO 변수 수 = n + m = n(1 + alpha). alpha=4.2일 때 QUBO는 원래 문제의 5.2배 크기. 이 확장된 변수 공간이 SA의 탐색을 어렵게 만든다.

## 파일 구성

| 파일 | 역할 |
|------|------|
| `qubo_quiet_planted.py` | 생성기 |
| `test_quiet_planted.py` | SA 실험 (alpha sweep, N scaling, 4-way comparison) |

## 참고 문헌

### 핵심 논문

1. **Krzakala, F. & Zdeborova, L.** "Hiding quiet solutions in random constraint satisfaction problems." *Physical Review Letters*, 102(23), 238701, 2009.
   - Quiet planting 개념 제시
   - Condensation threshold (alpha < 3.86) 이하에서 planted 인스턴스가 랜덤과 구별 불가능
   - [arXiv:0901.2130](https://arxiv.org/abs/0901.2130)

### Rosenberg Reduction

2. **Rosenberg, I. G.** "Reduction of bivalent maximization to the quadratic case." *Cahiers du Centre d'Etudes de Recherche Operationnelle*, 17, 71-74, 1975.
   - Cubic 이상의 의사불 함수를 quadratic으로 변환
   - 보조변수 y = z1*z2 도입, penalty term으로 y의 정확성 강제

3. **Boros, E. & Hammer, P. L.** "Pseudo-Boolean optimization." *Discrete Applied Mathematics*, 123(1-3), 155-225, 2002. — Rosenberg reduction의 일반화.

### 3-SAT 상전이

4. **Mezard, M., Parisi, G., & Zecchina, R.** "Analytic and algorithmic solution of random satisfiability problems." *Science*, 297(5582), 812-815, 2002.
   - Random 3-SAT의 condensation/satisfiability threshold 분석

5. **Krzakala, F., Montanari, A., Ricci-Tersenghi, F., Semerjian, G., & Zdeborova, L.** "Gibbs states and the set of solutions of random constraint satisfaction problems." *Proceedings of the National Academy of Sciences*, 104(25), 10318-10323, 2007.

### QUBO 변환

6. 표준 QUBO 변환: 3-SAT → penalty function → Rosenberg reduction → quadratic QUBO.

## 사용법

```bash
# 기본 생성 (target="10110", alpha=4.2)
python3 quiet_planting/qubo_quiet_planted.py 10110 4.2

# 시드 지정
python3 quiet_planting/qubo_quiet_planted.py 10110 4.2 42

# 길이로 랜덤 목표 생성
python3 quiet_planting/qubo_quiet_planted.py 50 3.86

# Alpha sweep 실험
python3 quiet_planting/test_quiet_planted.py 50 10

# N 스케일링 실험
python3 quiet_planting/test_quiet_planted.py --scaling 4.2

# 3-way 비교 (Quiet vs Wishart vs ZeroExp)
python3 quiet_planting/test_quiet_planted.py --compare
```
