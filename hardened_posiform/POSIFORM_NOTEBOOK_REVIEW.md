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

### 2.3 고정 계수 a — **수정 불필요**

**논문 (Hahn 2023, Section 2.3)**: "the coefficient b_{zz'} > 0 of the posiform is actually freely choosable (as long as it is positive)" → 양수이기만 하면 고정값도 유효. 실험에서는 {1, 2} 집합에서 선택.

**논문 (Pelofske 2024, Section 2.2)**: "All posiform coefficients corresponding to a 2-SAT clause are chosen as **1**." → 명시적으로 1 고정.

**노트북**: 모든 term에 **동일한 a = 0.1** 사용. 논문과 모순 없음 (양수 고정값).

**참고**: `qubo_posiform_hardened.py`도 기본값 `posiform_coeff_range=(1.0, 1.0)` (고정 계수).

### 2.4 posiform scaling factor (α) 부재 — **해결 완료**

**Pelofske 2024의 핵심**: `Q_final = Σ R_i + α × P` — α는 별도의 전역 스케일링 파라미터 (실험값: 0.1, 0.01). α가 작을수록 문제가 더 어려워짐.

**노트북 (수정 전)**: posiform이 random QUBO와 같은 행렬에 직접 더해짐 (α=1에 해당). per-term 계수 a=0.1이 다른 역할을 함.

```python
# 수정 전: mat(random QUBO)를 직접 수정 → α=1 고정
def posiform_planting(mat, opt, a):    # mat를 직접 변경
    ...
    mat[i][j] += a                     # random QUBO에 바로 더함
    ...
    # 반환값 없음 (void)

def gen_posiform_qubo(n, min_sub_graph_size, max_sub_graph_size, a):
    mat, opt = gen_concatenated_random_qubo(n, min_sub_graph_size, max_sub_graph_size)
    posiform_planting(mat, opt, a)     # mat를 직접 수정 → α=1 고정
    return mat, opt, opt@mat@opt
```

**해결**: posiform QUBO를 별도 행렬로 생성한 뒤, `Q_final = Q_random + α * Q_posiform`으로 결합.

```python
# [수정 2.4] posiform을 별도 행렬에 쌓아서 반환
def posiform_planting(opt, n, a):              # [수정 2.4] mat 제거, n 추가, 반환값 있음
    mat_posiform = np.zeros((n, n))            # [수정 2.4] 별도의 posiform 전용 행렬
    ...
    # [수정 2.4] mat → mat_posiform: posiform 전용 행렬에 기록
    mat_posiform[lo][hi] += a
    ...
    return mat_posiform                        # [수정 2.4] 별도 행렬 반환

# [수정 2.4] posiform_scale(α) 파라미터 추가
def gen_posiform_qubo(n, min_sub_graph_size, max_sub_graph_size, a, posiform_scale=1.0):
    mat_random, opt = gen_concatenated_random_qubo(n, min_sub_graph_size, max_sub_graph_size)
    mat_posiform = posiform_planting(opt, n, a)                   # [수정 2.4]
    mat = mat_random + posiform_scale * mat_posiform              # [수정 2.4]
    return mat, opt, opt @ mat @ opt
```

### 2.5 연속 vs. 이산 계수 random QUBO — **해결 완료**

**논문 (Pelofske 2024, Section 2.1)**: discrete coefficient {-1, +1} (lin2) 또는 {-1, -0.9, ..., 0.9, 1} (lin20).

**노트북 (수정 전)**: `np.random.uniform(-1, 1)` 연속 균일분포.

```python
# 수정 전: 연속 균일분포
coeff_boundary = 1

def gen_random_qubo(n):
    random_qubo = np.random.uniform(-coeff_boundary, coeff_boundary, (n, n))
    return np.triu(random_qubo)
```

**해결**: `coeff_type` 파라미터 추가. 이산 계수 집합 `COEFF_LIN2`, `COEFF_LIN20` 정의.

```python
# [수정 2.5] 이산 계수 집합 정의
COEFF_LIN2 = [-1, 1]                                              # [수정 2.5]
COEFF_LIN20 = [round(-1 + 0.1 * i, 1) for i in range(21)]        # [수정 2.5]

# [수정 2.5] coeff_type 파라미터 추가
def gen_random_qubo(n, coeff_type='lin2'):                         # [수정 2.5]
    coeffs = COEFF_LIN2 if coeff_type == 'lin2' else COEFF_LIN20   # [수정 2.5]
    random_qubo = np.array([[random.choice(coeffs) for _ in range(n)] for _ in range(n)])  # [수정 2.5]
    return np.triu(random_qubo)

# 호출 체인에 coeff_type 전파
def gen_concatenated_random_qubo(n, ..., coeff_type='lin2', ...):  # [수정 2.5]
    ...
    cur_mat = gen_random_qubo(cur_n, coeff_type)                   # [수정 2.5]
    ...

def gen_posiform_qubo(n, ..., coeff_type='lin2'):                  # [수정 2.5]
    mat_random, opt = gen_concatenated_random_qubo(..., coeff_type) # [수정 2.5]
    ...
```

### 2.6 Full matrix vs. upper-triangular — **해결 완료**

**노트북 (수정 전)**: full n×n matrix 사용. `x@mat@x`에서 (i,j)와 (j,i) 모두 에너지에 기여. random QUBO에서 mat[i][j]와 mat[j][i]가 독립적으로 생성됨.

**논문/기존 코드**: upper-triangular dict `{(i,j): weight}` (i≤j). 각 pair의 coupling이 단일 값.

```python
# 수정 전: full matrix
def gen_random_qubo(n):
    random_qubo = np.random.uniform(-coeff_boundary, coeff_boundary, (n, n))
    return random_qubo                 # 하삼각도 포함

def posiform_planting(mat, opt, a):
    ...
    mat[i][j] += a                     # i > j일 때 하삼각에 기록됨
```

**해결**: 상삼각 행렬 형식으로 전환.

```python
# [수정] gen_random_qubo — np.triu()로 상삼각만 유지
def gen_random_qubo(n):
    random_qubo = np.random.uniform(-coeff_boundary, coeff_boundary, (n, n))
    return np.triu(random_qubo)        # [수정] 상삼각만 유지, 하삼각은 0

# [수정] posiform_planting — off-diagonal 접근 시 mat[lo][hi] 사용
def posiform_planting(...):
    ...
    lo, hi = min(i, j), max(i, j)      # [수정] 상삼각 형식 유지

    if wi == 0 and wj == 0:            # (1-xi)(1-xj) * a
        mat_posiform[i][i] -= a
        mat_posiform[j][j] -= a
        mat_posiform[lo][hi] += a      # [수정] mat[i][j] → mat[lo][hi]
    elif wi == 0 and wj == 1:          # (1-xi)xj * a
        mat_posiform[j][j] += a
        mat_posiform[lo][hi] -= a      # [수정] mat[i][j] → mat[lo][hi]
    elif wi == 1 and wj == 0:          # xi(1-xj) * a
        mat_posiform[i][i] += a
        mat_posiform[lo][hi] -= a      # [수정] mat[i][j] → mat[lo][hi]
    else:                              # xixj * a
        mat_posiform[lo][hi] += a      # [수정] mat[i][j] → mat[lo][hi]
```

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
| Wrong tuple/pair | 1개 랜덤 | 1개 랜덤 (수정 완료) | 1개 랜덤 (via qubo_posiform.py) |
| Uniqueness 검증 | MiniSat 2-SAT solver | MiniSat (수정 완료) | Tarjan SCC + uniqueness check |
| Posiform 계수 | 자유 선택 (양수, 고정도 유효) | 고정 a | coeff_range (기본 1.0 고정) |
| Scaling α | 별도 파라미터 (0.1, 0.01) | 별도 파라미터 (posiform_scale, 수정 완료) | 별도 파라미터 (posiform_scale) |
| Random QUBO 계수 | 이산 (lin2/lin20) | 이산 (lin2/lin20, 수정 완료) | 이산 (lin2/lin20) |
| 행렬 형식 | upper-triangular | upper-triangular (수정 완료) | upper-triangular dict |
| Subgraph 분할 | Kernighan-Lin bisection | 균등 분할 (크기 차이 최대 1, 수정 완료) | 균등 순차 분할 |
| GS 유일성 보장 | 증명 (Section 2.2) | 보장 (MiniSat, 수정 완료) | 보장 (posiform uniqueness) |

---

## 5. 결론

노트북의 posiform planting은 **posiform의 non-negativity와 target에서의 zero 속성은 유지**하므로 ground state 보존은 수학적으로 맞다.

주요 수정 사항 (해결 완료):
- **2.1** Wrong tuple 1개 랜덤 선택 (논문 방식)
- **2.2** MiniSat 2-SAT uniqueness 검증 (GS 유일성 보장)
- **2.4** posiform scaling factor α 분리 (`Q_final = Q_random + α × Q_posiform`)
- **2.5** 이산 계수 random QUBO (lin2/lin20)
- **2.6** Upper-triangular 행렬 형식
- Subgraph 균등 분할 (`max_sub_graph_size`만 사용, 크기 차이 최대 1)

남은 차이점:
- **Subgraph 분할 알고리즘**: 논문은 Kernighan-Lin bisection (hardware graph용), 노트북은 순차 균등 분할. complete graph에서는 분할 알고리즘 차이가 성능에 영향 미미.

정밀 재현 구현은 `hardened_posiform/qubo_posiform_hardened.py`이며, 노트북은 핵심 구조를 대부분 반영한 학습/프로토타이핑 용도로 사용할 수 있다.

---

## 6. 수정 방법 제안

### 6.1 [필수] Wrong tuple 선택: 3개 전부 → 1개 랜덤 — **해결 완료**

**수정 전 (노트북)**:
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

### 6.2 [필수] Uniqueness 검증 추가 — **해결 완료** (방법 B 적용)

**방법 A — 간단한 변수 커버리지 확인** (미적용):
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

**방법 B — 2-SAT uniqueness 검증 도입** (적용됨):
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

### 6.3 [권장] posiform scaling factor (α) 분리 — **해결 완료**

**수정 전**: posiform과 random QUBO가 같은 행렬에 혼합 (α=1).

**적용된 수정**:
```python
# posiform_planting: mat 제거, 별도 행렬 생성 후 반환
def posiform_planting(opt, n, a):              # [수정 2.4] mat 제거, n 추가, 반환값 있음
    mat_posiform = np.zeros((n, n))            # [수정 2.4] 별도의 posiform 전용 행렬
    ...
    return mat_posiform                        # [수정 2.4] 별도 행렬 반환

# gen_posiform_qubo: α로 결합
def gen_posiform_qubo(n, min_sub_graph_size, max_sub_graph_size, a, posiform_scale=1.0):
    mat_random, opt = gen_concatenated_random_qubo(n, min_sub_graph_size, max_sub_graph_size)
    mat_posiform = posiform_planting(opt, n, a)                   # [수정 2.4]
    mat = mat_random + posiform_scale * mat_posiform              # [수정 2.4]
    return mat, opt, opt @ mat @ opt
```

**효과**: α를 독립적으로 조절하여 난이도 튜닝 가능 (α=0.01이 가장 어려움).

### 6.4 ~~[권장] 랜덤 양수 계수 도입~~ — **수정 불필요**

논문은 "freely choosable (as long as it is positive)"로 고정 계수도 유효. `qubo_posiform_hardened.py`도 기본값 `(1.0, 1.0)` 고정. 현재 노트북의 고정 `a = 0.1`은 논문과 모순 없음.

### 6.5 [권장] 이산 계수 random QUBO — **해결 완료**

**수정 전**:
```python
coeff_boundary = 1
def gen_random_qubo(n):
    return np.random.uniform(-coeff_boundary, coeff_boundary, (n, n))
```

**적용된 수정**:
```python
COEFF_LIN2 = [-1, 1]                                              # [수정 2.5]
COEFF_LIN20 = [round(-1 + 0.1 * i, 1) for i in range(21)]        # [수정 2.5]

def gen_random_qubo(n, coeff_type='lin2'):                         # [수정 2.5]
    coeffs = COEFF_LIN2 if coeff_type == 'lin2' else COEFF_LIN20
    random_qubo = np.array([[random.choice(coeffs) for _ in range(n)] for _ in range(n)])
    return np.triu(random_qubo)
```

**효과**: 논문과 동일한 이산 계수. lin2가 SA에 대해 더 어려운 문제를 생성.

### 6.6 [선택] float 비교 안전화 — **해결 완료** (2.1에서 함께 해결)

**수정 전**: `opt[i] != 0` float 비교.

**해결**: `int(opt[i])` 변환 사용. `target_tuple = (int(opt[i]), int(opt[j]))` 형태로 정수 비교.

### 6.7 [선택] 검증 루프 강화 — **해결 완료** (2.2에서 함께 해결)

**수정 전**: 최적해 변경만 확인.

**해결**: `find_opt_brute_force`에 `num_degenerate` 반환 추가. 검증 셀에서 `deg > 1`이면 Warning 출력.

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
| ~~**P0**~~ | ~~6.2 Uniqueness 검증~~ | ~~중~~ | ~~높음~~ — **해결 완료** |
| ~~**P0**~~ | ~~6.1 Wrong tuple 1개 선택~~ | ~~낮~~ | ~~높음~~ — **해결 완료** |
| ~~**P1**~~ | ~~6.3 α scaling 분리~~ | ~~낮~~ | ~~중~~ — **해결 완료** |
| ~~**P1**~~ | ~~6.4 랜덤 계수~~ | ~~낮~~ | ~~중~~ — **수정 불필요** |
| ~~**P2**~~ | ~~6.5 이산 계수~~ | ~~낮~~ | ~~중~~ — **해결 완료** |
| ~~**P2**~~ | ~~6.6 float 비교~~ | ~~낮~~ | ~~낮~~ — **해결 완료** (2.1) |
| ~~**P2**~~ | ~~6.7 검증 강화~~ | ~~낮~~ | ~~낮~~ — **해결 완료** (2.2) |
| ~~**P3**~~ | ~~6.8 행렬 형식~~ | ~~중~~ | ~~낮~~ — **해결 완료** |
