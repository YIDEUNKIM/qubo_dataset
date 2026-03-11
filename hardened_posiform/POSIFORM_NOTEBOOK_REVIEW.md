# posiform.ipynb 검토 피드백

## 검토 대상
- `hardened_posiform/posiform.ipynb`
- 비교 논문: Hahn, Pelofske, Djidjev (2023) "Posiform Planting", Pelofske, Hahn, Djidjev (2024) "Increasing the Hardness of Posiform Planting Using Random QUBOs"
- 비교 구현: `posiform/qubo_posiform.py`, `hardened_posiform/qubo_posiform_hardened.py`

---

## 1. 수학적으로 올바른 부분

### posiform 항 전개 (정확)

4개 조건과 행렬 업데이트 모두 올바름:

| Wrong tuple | 조건 | Indicator | 행렬 변경 |
|---|---|---|---|
| (0,0) | `opt[i]!=0 or opt[j]!=0` | `a(1-xi)(1-xj)` | `Q[i,i]-=a, Q[j,j]-=a, Q[i,j]+=a` |
| (0,1) | `opt[i]!=0 or opt[j]!=1` | `a(1-xi)xj` | `Q[j,j]+=a, Q[i,j]-=a` |
| (1,0) | `opt[i]!=1 or opt[j]!=0` | `axi(1-xj)` | `Q[i,i]+=a, Q[i,j]-=a` |
| (1,1) | `opt[i]!=1 or opt[j]!=1` | `axixj` | `Q[i,j]+=a` |

- 각 indicator는 non-negative이고, target (opt[i], opt[j])에서 정확히 0
- 조건 `opt[i] != a or opt[j] != b`는 `(opt[i], opt[j]) != (a, b)`과 동치 → 3개 wrong tuple만 정확히 선택
- ground state 보존은 수학적으로 보장됨

---

## 2. 논문과의 차이점

### 2.1 한 쌍당 3개 wrong tuple 전부 추가 (논문: 1개만) — **해결 완료**

**논문 (Hahn 2023, Section 2.3)**:
- 랜덤 pair (i,j) 선택 → 3개 wrong tuple 중 **1개만 랜덤 선택** → 2-SAT clause로 추가

**노트북 (수정 전)**:
- 매 pair마다 **3개 wrong tuple 전부** 동시 추가

수학적으로 유효하지만, 논문의 clause-by-clause 구성과 다름. 3개를 동시에 추가하면 pair (i,j)의 net 행렬 변경이 target에 따라 결정적으로 고정됨:

| Target | net Q[i,i] | net Q[j,j] | net Q[i,j] |
|---|---|---|---|
| (0,0) | +a | +a | -a |
| (0,1) | 0 | -a | +a |
| (1,0) | -a | 0 | +a |
| (1,1) | 0 | 0 | -a |

**해결**: `posiform_planting` 함수를 논문 방식으로 수정.
- `all_tuples = [(0,0), (0,1), (1,0), (1,1)]`에서 `target_tuple`을 제외한 3개 중 1개만 랜덤 선택
- `int(opt[i])` 변환으로 float 비교 위험도 함께 해소

```python
all_tuples = [(0, 0), (0, 1), (1, 0), (1, 1)]
target_tuple = (int(opt[i]), int(opt[j]))
wrong_tuples = [t for t in all_tuples if t != target_tuple]
wi, wj = random.choice(wrong_tuples)
```

### 2.2 2-SAT uniqueness 검증 없음 (가장 중요한 차이) — **해결 완료**

**논문**: clause 추가 후 **2-SAT solver (MiniSat/Tarjan SCC)로 유일성 확인**. Target이 유일한 해가 될 때까지 반복. posiform P(x)가 x\*에서만 0이 됨을 보장.

**노트북 (수정 전)**: uniqueness 검증 없이 t번 반복만 수행.

**문제**: 어떤 변수가 한 번도 pair에 선택되지 않으면, 그 변수를 flip해도 posiform 값 불변.
- n=500, t=1000일 때 변수 i 미선택 확률 ≈ e^(-4) ≈ 1.8%
- 약 9개 변수가 미포함 가능
- 해당 변수의 subgraph QUBO가 축퇴(degenerate)이면 **ground state 유일성 깨짐**

**해결**: 논문과 동일하게 **MiniSat (`python-sat`)** 으로 2-SAT uniqueness 검증 도입. `t` 파라미터는 제거하고, 내부적으로 `max_clauses = 10 * n`을 상한으로 사용 (2-SAT phase transition은 O(n)).

1. **MiniSat uniqueness check** (`posiform_planting` 수정):
   - `t` 파라미터 제거 → MiniSat이 유일성 확보 시 자동 종료
   - clause 추가 시 CNF clause도 함께 구성
   - 주기적으로 MiniSat에 blocking clause(target 차단)를 넣어 다른 해 존재 여부 확인
   - 유일하면 조기 종료, 아니면 clause 계속 추가 (상한 10*n)
   - 기존 커버리지 보충 로직은 제거 (MiniSat이 더 강한 보장)

```python
from pysat.solvers import Minisat22

def check_uniqueness():
    with Minisat22() as solver:
        for clause in clauses_cnf:
            solver.add_clause(clause)
        # target 차단: "적어도 하나의 변수가 target과 달라야 한다"
        blocking = [-(i + 1) if int(opt[i]) == 1 else (i + 1) for i in range(n)]
        solver.add_clause(blocking)
        return not solver.solve()  # UNSAT → 유일
```

2. **축퇴도 검증** (`find_opt_brute_force` + 검증 셀 수정):
   - `find_opt_brute_force`에 `num_degenerate` 반환 추가 (동일 에너지 상태 개수)
   - 검증 셀에서 `deg > 1`이면 Warning 출력

```python
opt, opt_val, deg = find_opt_brute_force(qubo, False)
if deg > 1:
    print(f"[Warning] Degenerate ground state: {deg} solutions")
```

### 2.3 고정 계수 a vs. 랜덤 양수 계수

**논문 (Hahn 2023, Section 2.3)**: "the coefficient b_{zz'} > 0 of the posiform is actually freely choosable (as long as it is positive)" → 각 clause마다 **독립적인 랜덤 양수 계수** b > 0 부여.

**노트북**: 모든 term에 **동일한 a = 0.1** 사용.

**영향**: 계수가 동일하면 QUBO의 에너지 landscape가 균일해져 solver에게 추가 정보를 줄 수 있음. 논문처럼 랜덤 계수를 사용하면 landscape가 더 복잡해짐.

### 2.4 posiform scaling factor (α) 부재

**Pelofske 2024의 핵심**: `Q_final = Σ R_i + α × P` — α는 별도의 전역 스케일링 파라미터 (실험값: 0.1, 0.01). α가 작을수록 문제가 더 어려워짐.

**노트북**: posiform이 random QUBO와 같은 행렬에 직접 더해짐 (α=1에 해당). per-term 계수 a=0.1이 다른 역할을 함.

**수정 방향**: posiform QUBO를 별도로 생성한 뒤, `Q_final = Q_random + α * Q_posiform`으로 결합. `qubo_posiform_hardened.py`의 방식 참고.

### 2.5 연속 vs. 이산 계수 random QUBO

**논문 (Pelofske 2024, Section 2.1)**: discrete coefficient {-1, +1} (lin2) 또는 {-1, -0.9, ..., 0.9, 1} (lin20).

**노트북**: `np.random.uniform(-1, 1)` 연속 균일분포.

**영향**: 이산 계수(특히 lin2)가 SA에 대해 더 어려운 문제를 생성함 (논문 실험 결과).

### 2.6 Full matrix vs. upper-triangular — **해결 완료**

**노트북 (수정 전)**: full n×n matrix 사용. `x@mat@x`에서 (i,j)와 (j,i) 모두 에너지에 기여. random QUBO에서 mat[i][j]와 mat[j][i]가 독립적으로 생성됨.

**논문/기존 코드**: upper-triangular dict `{(i,j): weight}` (i≤j). 각 pair의 coupling이 단일 값.

**해결**: 상삼각 행렬 형식으로 전환.
- `gen_random_qubo`: `np.triu()`로 하삼각을 0으로 처리
- `posiform_planting`: off-diagonal 업데이트 시 `lo, hi = min(i,j), max(i,j)`로 항상 상삼각 위치에 접근

**참고**: `x @ mat @ x`에서 full matrix는 coupling = `Q[i][j] + Q[j][i]` (두 번 기여), upper-triangular는 coupling = `Q[i][j]` (한 번)이므로, 같은 값이 들어있어도 에너지가 달라진다. D-Wave `neal`/`dimod`의 표준 입력도 upper-triangular `{(i,j): weight}` (i ≤ j) 형식이므로 상삼각이 맞다.

---

## 3. 잠재적 버그/위험

### 3.1 float 비교
`opt[i] != 0`, `opt[i] != 1`로 float를 비교. 현재 코드에서는 `np.array(bits)`가 정확한 0/1 정수를 반환하므로 문제없지만, 부동소수점 연산이 개입되면 조건 판정이 깨질 수 있음.

**수정 방향**: `int(opt[i])` 변환 또는 `opt[i] < 0.5` / `opt[i] > 0.5` 비교 사용.

### 3.2 상수항 누락
`(1-xi)(1-xj)` 전개 시 상수 `+a`가 행렬에 반영되지 않음. argmin에는 영향 없지만, 보고되는 `opt_val = opt@mat@opt`은 실제 posiform 에너지와 다름.

### 3.3 검증 루프의 한계
Cell 5에서 `np.array_equal(opt_init, opt)`로 최적해 변경만 확인하지, 축퇴도(degeneracy)는 검사하지 않음. 여러 개의 ground state가 있어도 enumeration 순서에 따라 같은 해가 반환되면 통과함.

**수정 방향**: brute force에서 ground state 개수 (degeneracy count)도 함께 확인.

---

## 4. 비교 요약표

| 항목 | 논문 (Pelofske 2024) | 노트북 (posiform.ipynb) | 기존 구현 (qubo_posiform_hardened.py) |
|---|---|---|---|
| Wrong tuple/pair | 1개 랜덤 | 3개 전부 | 1개 랜덤 (via qubo_posiform.py) |
| Uniqueness 검증 | MiniSat 2-SAT solver | 없음 | Tarjan SCC + uniqueness check |
| Posiform 계수 | 랜덤 양수 (논문: 1 고정) | 고정 a | 랜덤 양수 (coeff_range) |
| Scaling α | 별도 파라미터 (0.1, 0.01) | 없음 (α=1) | 별도 파라미터 (posiform_scale) |
| Random QUBO 계수 | 이산 (lin2/lin20) | 연속 uniform | 이산 (lin2/lin20) |
| 행렬 형식 | upper-triangular | upper-triangular (수정 완료) | upper-triangular dict |
| Subgraph 분할 | Kernighan-Lin bisection | 랜덤 크기 순차 분할 | 균등 순차 분할 |
| GS 유일성 보장 | 증명 (Section 2.2) | 미보장 | 보장 (posiform uniqueness) |

---

## 5. 결론

노트북의 posiform planting은 **posiform의 non-negativity와 target에서의 zero 속성은 유지**하므로 ground state 보존은 수학적으로 맞다. 그러나 **uniqueness 보장이 없고, 논문의 핵심 구조 (2-SAT + uniqueness check + α scaling + 이산 계수)가 빠져 있어** Pelofske 2024 논문의 정확한 재현이라고 보기 어렵다.

논문을 정확히 재현한 구현은 `hardened_posiform/qubo_posiform_hardened.py`이며, 노트북은 posiform planting의 기본 아이디어를 프로토타이핑한 수준으로 이해해야 한다.

---

## 6. 수정 방법 제안

### 6.1 [필수] Wrong tuple 선택: 3개 전부 → 1개 랜덤

**현재 (노트북)**:
```python
# 매 pair마다 3개 wrong tuple 전부 추가
if opt[i] != 0 or opt[j] != 0:  # (0,0) 배제
    ...
if opt[i] != 0 or opt[j] != 1:  # (0,1) 배제
    ...
if opt[i] != 1 or opt[j] != 0:  # (1,0) 배제
    ...
if opt[i] != 1 or opt[j] != 1:  # (1,1) 배제
    ...
```

**수정안** (현재 적용됨):
```python
def posiform_planting(mat, opt, a):
    n = mat.shape[0]
    all_tuples = [(0, 0), (0, 1), (1, 0), (1, 1)]
    clauses_cnf = []
    max_clauses = 10 * n  # t 대신 자동 상한

    for step in range(max_clauses):
        i, j = random.sample(range(n), 2)
        target_tuple = (int(opt[i]), int(opt[j]))
        wrong_tuples = [t for t in all_tuples if t != target_tuple]
        wi, wj = random.choice(wrong_tuples)

        # posiform 항 추가 + CNF clause 추가
        # ... (MiniSat uniqueness check로 자동 종료)
```

**효과**: 논문과 동일한 clause-by-clause 구성. 에너지 landscape의 다양성 증가.

### 6.2 [필수] Uniqueness 검증 추가

**방법 A — 간단한 변수 커버리지 확인** (최소한의 수정):
```python
def posiform_planting(mat, opt, a, t):
    n = mat.shape[0]
    covered = set()

    for _ in range(t):
        i, j = random.sample(range(n), 2)
        covered.add(i)
        covered.add(j)
        # ... (posiform 항 추가)

    # 미커버 변수 보충
    uncovered = set(range(n)) - covered
    for i in uncovered:
        j = random.choice([x for x in range(n) if x != i])
        # posiform 항 추가 (위와 동일)
```

**방법 B — 2-SAT uniqueness 검증 도입** (논문 재현):
```python
# posiform/qubo_posiform.py의 함수 활용
from posiform.qubo_posiform import create_planted_2sat, posiform_to_qubo

def gen_posiform_qubo_v2(n, min_sub_graph_size, max_sub_graph_size, alpha):
    mat, opt = gen_concatenated_random_qubo(n, min_sub_graph_size, max_sub_graph_size)
    target = ''.join(str(int(b)) for b in opt)

    # 유일성 보장된 posiform QUBO 생성
    clauses, is_unique = create_planted_2sat(target)
    Q_posiform, _ = posiform_to_qubo(n, clauses, coeff_range=(1.0, 1.0))

    # Q_final = Q_random + alpha * Q_posiform
    for (i, j), w in Q_posiform.items():
        mat[i][j] += alpha * w  # upper-triangular → full matrix 변환 주의

    return mat, opt, opt @ mat @ opt, is_unique
```

### 6.3 [권장] posiform scaling factor (α) 분리

**현재**: posiform과 random QUBO가 같은 행렬에 혼합 (α=1).

**수정안**:
```python
def gen_posiform_qubo(n, min_sub_graph_size, max_sub_graph_size, alpha):
    # 1. Random QUBO 생성
    mat_random, opt = gen_concatenated_random_qubo(n, min_sub_graph_size, max_sub_graph_size)

    # 2. Posiform QUBO를 별도 행렬로 생성
    mat_posiform = np.zeros((n, n))
    posiform_planting(mat_posiform, opt, a=1.0)  # MiniSat이 자동 종료

    # 3. 결합: Q_final = Q_random + α * Q_posiform
    mat = mat_random + alpha * mat_posiform

    return mat, opt, opt @ mat @ opt
```

**효과**: α를 독립적으로 조절하여 난이도 튜닝 가능 (α=0.01이 가장 어려움).

### 6.4 [권장] 랜덤 양수 계수 도입

**현재**: 모든 posiform term에 고정 계수 `a`.

**수정안**:
```python
def posiform_planting(mat, opt, coeff_range):
    n = mat.shape[0]
    lo, hi = coeff_range
    max_clauses = 10 * n
    for step in range(max_clauses):
        i, j = random.sample(range(n), 2)
        a = random.uniform(lo, hi)  # 매 clause마다 랜덤 계수
        # ... (MiniSat uniqueness check로 자동 종료)
```

**효과**: 에너지 landscape가 더 복잡해져 solver에게 패턴을 노출하지 않음.

### 6.5 [권장] 이산 계수 random QUBO

**현재**:
```python
def gen_random_qubo(n):
    return np.random.uniform(-coeff_boundary, coeff_boundary, (n, n))
```

**수정안**:
```python
COEFF_LIN2 = [-1, 1]
COEFF_LIN20 = [round(-1 + 0.1 * i, 1) for i in range(21)]

def gen_random_qubo(n, coeff_type='lin2'):
    coeffs = COEFF_LIN2 if coeff_type == 'lin2' else COEFF_LIN20
    return np.array([[random.choice(coeffs) for _ in range(n)] for _ in range(n)])
```

**효과**: 논문과 동일한 이산 계수. lin2가 SA에 대해 더 어려운 문제를 생성.

### 6.6 [선택] float 비교 안전화

**현재**:
```python
if opt[i] != 0 or opt[j] != 0:
```

**수정안**:
```python
opt = opt.astype(int)  # 함수 시작 시 정수 변환
# 이후 opt[i] != 0 비교는 정수 비교로 안전
```

### 6.7 [선택] 검증 루프 강화

**현재**: 최적해 변경만 확인.

**수정안**:
```python
def find_opt_brute_force(mat, debug=False):
    n = mat.shape[0]
    best_x = None
    best_val = float('inf')
    num_degenerate = 0  # 축퇴도 추가

    for bits in product([0, 1], repeat=n):
        x = np.array(bits)
        cur_val = x @ mat @ x

        if cur_val < best_val - 1e-12:
            best_val = cur_val
            best_x = x
            num_degenerate = 1
        elif abs(cur_val - best_val) < 1e-12:
            num_degenerate += 1

    if debug:
        print("opt_x:", best_x, "\nopt_val:", best_val, "\ndegeneracy:", num_degenerate)

    return best_x, best_val, num_degenerate

# 검증 루프에서
opt, opt_val, deg = find_opt_brute_force(qubo, False)
if not np.array_equal(opt_init, opt):
    print("[Error] Optimum changed!")
if deg > 1:
    print(f"[Warning] Degenerate ground state: {deg} solutions")
```

### 6.8 [선택] Upper-triangular 행렬 형식 전환 — **해결 완료**

**수정 전**: full n×n matrix → `x @ mat @ x`에서 대칭 pair (i,j)와 (j,i) 모두 기여.

**해결**:

`gen_random_qubo` — `np.triu()`로 상삼각만 유지:
```python
def gen_random_qubo(n):
    random_qubo = np.random.uniform(-coeff_boundary, coeff_boundary, (n, n))
    return np.triu(random_qubo)  # [수정] 상삼각만 유지, 하삼각은 0
```

`posiform_planting` — off-diagonal 접근 시 `mat[lo][hi]` 사용:
```python
# [수정] 상삼각 형식 유지: off-diagonal은 항상 mat[작은][큰]에 접근
lo, hi = min(i, j), max(i, j)

if wi == 0 and wj == 0:   # (1-xi)(1-xj) * a
    mat[i][i] -= a
    mat[j][j] -= a
    mat[lo][hi] += a      # [수정] mat[i][j] → mat[lo][hi]
elif wi == 0 and wj == 1: # (1-xi)xj * a
    mat[j][j] += a
    mat[lo][hi] -= a      # [수정] mat[i][j] → mat[lo][hi]
elif wi == 1 and wj == 0: # xi(1-xj) * a
    mat[i][i] += a
    mat[lo][hi] -= a      # [수정] mat[i][j] → mat[lo][hi]
else:                     # xixj * a
    mat[lo][hi] += a      # [수정] mat[i][j] → mat[lo][hi]
```

---

## 7. 수정 우선순위

| 우선순위 | 항목 | 난이도 | 영향도 |
|---|---|---|---|
| **P0** | 6.2 Uniqueness 검증 | 중 | 높음 — GS 유일성 보장 |
| **P0** | 6.1 Wrong tuple 1개 선택 | 낮 | 높음 — 논문 방법론 일치 |
| **P1** | 6.3 α scaling 분리 | 낮 | 중 — 난이도 튜닝 |
| **P1** | 6.4 랜덤 계수 | 낮 | 중 — landscape 다양성 |
| **P2** | 6.5 이산 계수 | 낮 | 중 — 논문 재현 |
| **P2** | 6.6 float 비교 | 낮 | 낮 — 안전성 |
| **P2** | 6.7 검증 강화 | 낮 | 낮 — 디버깅 |
| ~~**P3**~~ | ~~6.8 행렬 형식~~ | ~~중~~ | ~~낮~~ — **해결 완료** |
