"""
Posiform Planting QUBO 생성기 (2-SAT 기반)

Hahn, Pelofske, Djidjev (2023) 기반: planted 2-SAT → posiform → QUBO.
보조변수 없이 직접 n-bit QUBO를 생성.

핵심 아이디어:
  1. Target x*로부터 2-SAT clause를 반복 생성 (랜덤 변수쌍 → wrong tuple 배제)
  2. 2-SAT이 x*만을 유일한 해로 가질 때까지 clause 추가 (Tarjan SCC + uniqueness)
  3. 각 clause에 랜덤 양수 계수 부여 → posiform P(x) 구성
  4. P(x)를 x̄=1-x 대입으로 표준 QUBO Q 행렬로 변환
  5. P(x*) = 0이므로 x*가 유일한 ground state

장점:
  - 보조변수 없음 (QUBO 크기 = n)
  - x*가 수학적으로 유일한 ground state임이 보장
  - 2-SAT은 다항 시간에 풀리지만, QUBO로 변환 후에는 NP-hard
"""

import random
import itertools
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from qubo_utils import calculate_energy, save_qubo_edgelist, print_q_matrix, print_qubo_formula


def tarjan_scc(adj, n_nodes):
    """
    반복적(iterative) Tarjan SCC 알고리즘.

    Python 재귀 제한 회피를 위해 명시적 스택 사용.

    Args:
        adj: 인접 리스트 (list of lists)
        n_nodes: 노드 수

    Returns:
        comp: 각 노드의 SCC 번호 (역위상 순서: 작은 번호가 위상적으로 뒤)
    """
    index_counter = [0]
    stack = []
    on_stack = [False] * n_nodes
    index = [-1] * n_nodes
    lowlink = [-1] * n_nodes
    comp = [-1] * n_nodes
    comp_counter = [0]

    def strongconnect(v):
        # 반복적 구현: (node, neighbor_index) 스택
        call_stack = [(v, 0)]
        index[v] = lowlink[v] = index_counter[0]
        index_counter[0] += 1
        stack.append(v)
        on_stack[v] = True

        while call_stack:
            node, ni = call_stack[-1]

            if ni < len(adj[node]):
                call_stack[-1] = (node, ni + 1)
                w = adj[node][ni]
                if index[w] == -1:
                    # 미방문: 재귀 진입
                    index[w] = lowlink[w] = index_counter[0]
                    index_counter[0] += 1
                    stack.append(w)
                    on_stack[w] = True
                    call_stack.append((w, 0))
                elif on_stack[w]:
                    lowlink[node] = min(lowlink[node], index[w])
            else:
                # 모든 이웃 처리 완료
                if lowlink[node] == index[node]:
                    # SCC 루트: 스택에서 팝
                    c = comp_counter[0]
                    comp_counter[0] += 1
                    while True:
                        w = stack.pop()
                        on_stack[w] = False
                        comp[w] = c
                        if w == node:
                            break

                call_stack.pop()
                if call_stack:
                    parent = call_stack[-1][0]
                    lowlink[parent] = min(lowlink[parent], lowlink[node])

    for v in range(n_nodes):
        if index[v] == -1:
            strongconnect(v)

    return comp


def build_implication_graph(n, clauses):
    """
    2-SAT 함축 그래프 구축.

    리터럴 인코딩: 2*i = x_i, 2*i+1 = ¬x_i
    각 clause (l₁ ∨ l₂): ¬l₁ → l₂, ¬l₂ → l₁

    Args:
        n: 변수 수
        clauses: [(lit1, lit2), ...] — lit = (var_idx, is_positive)

    Returns:
        adj: 인접 리스트 (2*n 노드)
    """
    num_nodes = 2 * n
    adj = [[] for _ in range(num_nodes)]

    for (v1, p1), (v2, p2) in clauses:
        # 리터럴 → 그래프 노드
        l1 = 2 * v1 + (0 if p1 else 1)      # l₁
        l2 = 2 * v2 + (0 if p2 else 1)      # l₂
        neg_l1 = l1 ^ 1                      # ¬l₁
        neg_l2 = l2 ^ 1                      # ¬l₂

        # ¬l₁ → l₂
        adj[neg_l1].append(l2)
        # ¬l₂ → l₁
        adj[neg_l2].append(l1)

    return adj


def solve_2sat(n, clauses):
    """
    Tarjan SCC 기반 2-SAT 풀기.

    Args:
        n: 변수 수
        clauses: [(lit1, lit2), ...] — lit = (var_idx, is_positive)

    Returns:
        solution: [bool, ...] 길이 n (None if UNSAT)
    """
    adj = build_implication_graph(n, clauses)
    comp = tarjan_scc(adj, 2 * n)

    # UNSAT 판정: x_i와 ¬x_i가 같은 SCC
    for i in range(n):
        if comp[2 * i] == comp[2 * i + 1]:
            return None

    # 해 구성: Tarjan 역위상 순서 규약
    # comp 번호가 작을수록 위상적으로 뒤 → comp[2*i] > comp[2*i+1]이면 x_i = True
    solution = [comp[2 * i] > comp[2 * i + 1] for i in range(n)]
    return solution


def check_2sat_uniqueness(n, clauses, solution):
    """
    2-SAT 해의 유일성 검사.

    각 변수를 현재 값의 반대로 강제하는 unit clause 추가 후 satisfiability 확인.
    모든 변수에서 UNSAT이면 유일한 해.

    Args:
        n: 변수 수
        clauses: 기존 clause 리스트
        solution: [bool, ...] 현재 해

    Returns:
        is_unique: True if solution is the only satisfying assignment
        flippable: 반대로 강제해도 SAT인 변수 인덱스 목록
    """
    flippable = []
    for i in range(n):
        # x_i를 반대로 강제: (x_i = !solution[i]) → unit clause (¬x_i ∨ ¬x_i) 또는 (x_i ∨ x_i)
        forced_val = not solution[i]
        # unit clause: 변수 i가 forced_val이어야 함 → (i, forced_val) ∨ (i, forced_val)
        test_clauses = clauses + [((i, forced_val), (i, forced_val))]
        result = solve_2sat(n, test_clauses)
        if result is not None:
            flippable.append(i)
    return len(flippable) == 0, flippable


def create_planted_2sat(target, max_clauses_factor=10, seed=None):
    """
    Target이 유일한 해인 planted 2-SAT 생성.

    2단계 접근:
      1. 랜덤 clause 추가 (대부분의 변수 고정)
      2. Targeted 단계: 아직 flippable한 변수에 대해 집중 clause 추가

    Args:
        target: 목표 비트스트링 (예: "10110")
        max_clauses_factor: 최대 clause 수 = factor * n
        seed: 난수 시드

    Returns:
        clauses: [(lit1, lit2), ...] — lit = (var_idx, is_positive)
        is_unique: 유일한 해 여부
    """
    n = len(target)
    if n < 2:
        raise ValueError(f"변수 수가 2 미만입니다: n={n}")

    rng = random.Random(seed)
    target_bits = [int(b) for b in target]
    target_sol = [b == 1 for b in target_bits]
    max_clauses = max_clauses_factor * n

    clauses = []
    added = set()  # 중복 clause 방지

    def add_clause_for_pair(i, j, wrong=None):
        """변수 쌍 (i,j)에 대해 clause 추가. wrong이 None이면 랜덤 선택."""
        ti, tj = target_bits[i], target_bits[j]
        if wrong is None:
            wrong_tuples = [(vi, vj) for vi in (0, 1) for vj in (0, 1)
                            if (vi, vj) != (ti, tj)]
            wrong = rng.choice(wrong_tuples)

        lit1 = (i, wrong[0] == 0)
        lit2 = (j, wrong[1] == 0)
        clause_key = (lit1, lit2)
        if clause_key in added:
            return False
        added.add(clause_key)
        clauses.append((lit1, lit2))
        return True

    # --- Phase 1: 랜덤 clause 추가 ---
    check_interval = max(1, n // 4)

    for iteration in range(max_clauses):
        i, j = sorted(rng.sample(range(n), 2))
        add_clause_for_pair(i, j)

        if (len(clauses) % check_interval == 0):
            is_unique, flippable = check_2sat_uniqueness(n, clauses, target_sol)
            if is_unique:
                return clauses, True

    # --- Phase 2: Targeted clause 추가 ---
    # 아직 flippable한 변수에 대해 집중적으로 clause 추가
    is_unique, flippable = check_2sat_uniqueness(n, clauses, target_sol)
    if is_unique:
        return clauses, True

    max_targeted = 3 * n  # targeted 단계 최대 추가 clause
    for _ in range(max_targeted):
        if not flippable:
            break

        # flippable 변수 하나 선택
        fi = rng.choice(flippable)
        # fi와 다른 변수를 쌍으로 clause 추가
        others = [j for j in range(n) if j != fi]
        partner = rng.choice(others)

        i, j = min(fi, partner), max(fi, partner)
        # 가능한 모든 wrong tuple 시도
        ti, tj = target_bits[i], target_bits[j]
        wrong_tuples = [(vi, vj) for vi in (0, 1) for vj in (0, 1)
                        if (vi, vj) != (ti, tj)]
        rng.shuffle(wrong_tuples)
        for wrong in wrong_tuples:
            if add_clause_for_pair(i, j, wrong):
                break

        # 주기적 재확인
        if len(clauses) % max(1, len(flippable)) == 0:
            is_unique, flippable = check_2sat_uniqueness(n, clauses, target_sol)
            if is_unique:
                return clauses, True

    # 최종 확인
    is_unique, _ = check_2sat_uniqueness(n, clauses, target_sol)
    return clauses, is_unique


def posiform_to_qubo(n, clauses, coeff_range=(1.0, 3.0), seed=None):
    """
    2-SAT clause들을 posiform으로 구성 후 QUBO로 변환.

    각 clause는 하나의 wrong tuple (vi, vj)를 배제:
      Exclude (0,0): b*(1-x_i)(1-x_j) → Q[i,i]-=b, Q[j,j]-=b, Q[i,j]+=b, const+=b
      Exclude (0,1): b*(1-x_i)*x_j → Q[j,j]+=b, Q[i,j]-=b
      Exclude (1,0): b*x_i*(1-x_j) → Q[i,i]+=b, Q[i,j]-=b
      Exclude (1,1): b*x_i*x_j → Q[i,j]+=b

    Args:
        n: 변수 수
        clauses: [(lit1, lit2), ...]
        coeff_range: 랜덤 양수 계수 범위 (lo, hi)
        seed: 난수 시드

    Returns:
        Q: QUBO 딕셔너리 {(i,j): weight} (i <= j)
        constant: 상수항 (에너지 오프셋)
    """
    rng = random.Random(seed)
    Q = {}
    constant = 0.0

    for (v1, p1), (v2, p2) in clauses:
        i, j = min(v1, v2), max(v1, v2)
        # 계수
        lo, hi = coeff_range
        b = rng.uniform(lo, hi)

        # clause가 배제하는 할당 결정
        # lit = (var, is_positive): is_positive=True → x_i, is_positive=False → ¬x_i
        # clause: (x_i if p1 else ¬x_i) ∨ (x_j if p2 else ¬x_j)
        # 위반(배제) 조건: 두 리터럴 모두 False
        #   x_i가 p1이면 True → 위반 시 x_i = (not p1)
        #   x_j가 p2이면 True → 위반 시 x_j = (not p2)
        wrong_i = 0 if p1 else 1  # p1=True(x_i) → 위반 시 x_i=0, p1=False(¬x_i) → 위반 시 x_i=1
        wrong_j = 0 if p2 else 1

        # 변수 인덱스 정렬 반영
        if v1 <= v2:
            wi, wj = wrong_i, wrong_j
        else:
            wi, wj = wrong_j, wrong_i

        # posiform 항 → QUBO 변환
        if wi == 0 and wj == 0:
            # b*(1-x_i)*(1-x_j) = b - b*x_i - b*x_j + b*x_i*x_j
            constant += b
            Q[(i, i)] = Q.get((i, i), 0) - b
            Q[(j, j)] = Q.get((j, j), 0) - b
            Q[(i, j)] = Q.get((i, j), 0) + b
        elif wi == 0 and wj == 1:
            # b*(1-x_i)*x_j = b*x_j - b*x_i*x_j
            Q[(j, j)] = Q.get((j, j), 0) + b
            Q[(i, j)] = Q.get((i, j), 0) - b
        elif wi == 1 and wj == 0:
            # b*x_i*(1-x_j) = b*x_i - b*x_i*x_j
            Q[(i, i)] = Q.get((i, i), 0) + b
            Q[(i, j)] = Q.get((i, j), 0) - b
        else:
            # b*x_i*x_j
            Q[(i, j)] = Q.get((i, j), 0) + b

    # 0에 가까운 항 제거
    Q = {k: v for k, v in Q.items() if abs(v) > 1e-15}

    return Q, constant


def create_qubo_posiform(target, coeff_range=(1.0, 3.0), max_clauses_factor=10, seed=None):
    """
    Posiform Planting QUBO 생성 (메인 진입점).

    planted 2-SAT → posiform → QUBO.
    보조변수 없이 직접 n-bit QUBO 생성.

    Args:
        target: 목표 비트스트링 (예: "10110")
        coeff_range: posiform 계수 범위 (lo, hi)
        max_clauses_factor: 최대 clause 수 = factor * n
        seed: 난수 시드

    Returns:
        Q: QUBO 딕셔너리 {(i,j): weight} (i <= j)
        info: 메타정보 딕셔너리
    """
    n = len(target)

    # 2-SAT 시드와 계수 시드 분리
    if seed is not None:
        sat_seed = seed
        coeff_seed = seed + 1000
    else:
        sat_seed = None
        coeff_seed = None

    # planted 2-SAT 생성
    clauses, is_unique = create_planted_2sat(target, max_clauses_factor=max_clauses_factor,
                                              seed=sat_seed)

    # posiform → QUBO 변환
    Q, constant = posiform_to_qubo(n, clauses, coeff_range=coeff_range, seed=coeff_seed)

    # target 에너지 계산
    target_energy = calculate_energy(target, Q)

    info = {
        'n': n,
        'num_clauses': len(clauses),
        'is_unique': is_unique,
        'coeff_range': coeff_range,
        'constant_offset': constant,
        'target_energy': target_energy,
        'clauses': clauses,
    }

    return Q, info


def verify_ground_state(Q, target, num_random_samples=10000):
    """
    target이 ground state인지 검증.

    1. Single-flip 이웃 검사 (완전)
    2. 랜덤 샘플 검사

    Returns:
        (is_local_min, is_likely_global, stats_dict)
    """
    n = len(target)
    target_energy = calculate_energy(target, Q)

    # 1. Single-flip 이웃 검사
    min_flip_delta = float('inf')
    for i in range(n):
        flipped = list(target)
        flipped[i] = '0' if flipped[i] == '1' else '1'
        flipped_str = ''.join(flipped)
        flipped_energy = calculate_energy(flipped_str, Q)
        delta = flipped_energy - target_energy
        min_flip_delta = min(min_flip_delta, delta)

    is_local_min = min_flip_delta > -1e-10

    # 2. 랜덤 샘플 검사
    lower_count = 0
    min_random_energy = float('inf')

    for _ in range(num_random_samples):
        rand_state = ''.join(str(random.randint(0, 1)) for _ in range(n))
        rand_energy = calculate_energy(rand_state, Q)
        if rand_energy < target_energy - 1e-10:
            lower_count += 1
        min_random_energy = min(min_random_energy, rand_energy)

    is_likely_global = (lower_count == 0)

    stats = {
        'target_energy': target_energy,
        'min_flip_delta': min_flip_delta,
        'is_local_min': is_local_min,
        'lower_count': lower_count,
        'num_random_samples': num_random_samples,
        'min_random_energy': min_random_energy,
        'is_likely_global': is_likely_global,
    }

    return is_local_min, is_likely_global, stats


def verify_brute_force(Q, target, n):
    """n <= 20일 때 brute force로 ground state 완전 검증."""
    if n > 20:
        return None

    target_energy = calculate_energy(target, Q)
    best_energy = float('inf')
    best_state = None
    num_degenerate = 0

    for bits in itertools.product([0, 1], repeat=n):
        state = ''.join(map(str, bits))
        energy = calculate_energy(state, Q)
        if energy < best_energy - 1e-10:
            best_energy = energy
            best_state = state
            num_degenerate = 1
        elif abs(energy - best_energy) < 1e-10:
            num_degenerate += 1

    return {
        'target_energy': target_energy,
        'best_energy': best_energy,
        'best_state': best_state,
        'is_ground_state': abs(best_energy - target_energy) < 1e-10,
        'energy_gap': target_energy - best_energy if not abs(best_energy - target_energy) < 1e-10 else 0.0,
        'num_degenerate': num_degenerate,
    }


if __name__ == "__main__":
    target = "10110"
    coeff_range = (1.0, 3.0)
    seed = None

    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if set(arg) <= {'0', '1'} and len(arg) >= 2:
            # 0과 1로만 구성된 문자열 → 이진 타겟으로 해석
            target = arg
        elif arg.isdigit():
            # 순수 숫자 (2~9 포함) → N 길이로 해석
            length = int(arg)
            if length < 2:
                length = 2
            random.seed(42)
            target = ''.join(str(random.randint(0, 1)) for _ in range(length))
            print(f"[설정] 길이 {length}의 랜덤 목표 해 생성")
        else:
            target = arg

    if len(sys.argv) > 2:
        seed = int(sys.argv[2])

    n = len(target)
    print("=" * 60)
    print("Posiform Planting QUBO 생성기 (Planted 2-SAT → Posiform)")
    print("=" * 60)
    print(f"Target: {target} (N={n})")
    print(f"계수 범위: {coeff_range}")
    if seed is not None:
        print(f"Seed: {seed}")

    # QUBO 생성
    Q, info = create_qubo_posiform(target, coeff_range=coeff_range, seed=seed)

    print(f"\n[2-SAT 정보]")
    print(f"  Clause 수: {info['num_clauses']}")
    print(f"  유일한 해: {'OK' if info['is_unique'] else 'FAIL (clause 부족)'}")
    print(f"  QUBO 변수 수: {n} (보조변수 없음)")
    print(f"  QUBO 비영 항 수: {len(Q)}")

    # 행렬 출력 (작을 때만)
    if n <= 12:
        print_q_matrix(Q, n)
        print_qubo_formula(Q)

    # 에너지 검증
    target_energy = info['target_energy']
    print(f"\n[에너지 검증]")
    print(f"  Target 에너지: {target_energy:.6f}")

    # Brute force 검증
    if n <= 20:
        print(f"\n[Brute Force 검증] (N={n})")
        bf = verify_brute_force(Q, target, n)
        print(f"  Target 에너지: {bf['target_energy']:.6f}")
        print(f"  최소 에너지:   {bf['best_energy']:.6f}")
        print(f"  최소 상태:     {bf['best_state']}")
        print(f"  축퇴도:        {bf['num_degenerate']}")
        if bf['is_ground_state']:
            print("  Ground State 검증 성공!")
        else:
            print(f"  Ground State 검증 실패 (gap: {bf['energy_gap']:.6f})")
    else:
        print(f"\n[통계적 검증] (N={n})")
        is_local, is_global, stats = verify_ground_state(Q, target)
        print(f"  Target 에너지: {stats['target_energy']:.6f}")
        print(f"  최소 flip delta: {stats['min_flip_delta']:.6f}")
        print(f"  Local minimum: {'OK' if is_local else 'FAIL'}")
        print(f"  랜덤 {stats['num_random_samples']}개 중 더 낮은 에너지: {stats['lower_count']}개")
        print(f"  Global minimum (추정): {'OK' if is_global else 'FAIL'}")

    # 결과 저장
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(output_dir, exist_ok=True)
    filename_target = f"{target[:5]}_{n}"
    output_file = os.path.join(output_dir, f"posiform_{filename_target}.txt")
    save_qubo_edgelist(Q, output_file, target)
    print(f"\n[저장 완료] {output_file}")
