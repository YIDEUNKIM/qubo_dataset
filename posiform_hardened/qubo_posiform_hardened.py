"""
Hardened Posiform Planting QUBO 생성기

Pelofske, Hahn, Djidjev (2024) 기반:
"Increasing the Hardness of Posiform Planting Using Random QUBOs
 for Programmable Quantum Annealer Benchmarking"

핵심 아이디어:
  1. n개 변수를 k개의 disjoint subgraph로 분할
  2. 각 subgraph에 discrete coefficient random QUBO 생성
  3. 각 random QUBO의 ground state를 정확히 계산 (brute force)
  4. 모든 subproblem 최적해를 concatenate → target bitstring
  5. 이 target으로 전체 그래프에 posiform planted QUBO 생성
  6. Q_final = Σ R_i + α × P  (α = posiform scaling coefficient)

Ground State 유일성 증명 (논문 Section 2.2):
  sub_i(X*) minimizes R_i (by construction)
  X* uniquely minimizes P (posiform guarantee)
  → Q(X*) = Σ R_i(X*) + αP(X*) < Σ R_i(X̂) + αP(X̂) for any X̂ ≠ X*

난이도 조절 파라미터:
  - posiform_scale (α): 작을수록 어려움 (논문: 0.01이 가장 어려움)
  - coeff_type: 'lin2' ({-1,+1}) vs 'lin20' ({-1,-0.9,...,0.9,1})
  - subgraph_size: 각 subgraph의 변수 수
"""

import random
import itertools
import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from qubo_utils import calculate_energy, save_qubo_edgelist, print_q_matrix, print_qubo_formula
from posiform.qubo_posiform import create_qubo_posiform


# ============================================================
# Discrete Coefficient Sets (논문 Section 2.1)
# ============================================================
COEFF_LIN2 = [-1, 1]
COEFF_LIN20 = [round(-1 + 0.1 * i, 1) for i in range(21)]  # -1.0, -0.9, ..., 0.9, 1.0


def partition_variables(n, max_subgraph_size):
    """
    n개 변수를 max_subgraph_size 이하의 disjoint 그룹으로 분할.

    논문: Kernighan-Lin recursive bisection (hardware graph용).
    여기서는 complete graph이므로 순차적 동일 크기 분할.

    Args:
        n: 총 변수 수
        max_subgraph_size: 각 subgraph 최대 크기

    Returns:
        partitions: [[var_indices], ...] — disjoint 변수 그룹
    """
    k = max(1, -(-n // max_subgraph_size))  # ceil division
    partitions = []
    for i in range(k):
        start = i * n // k
        end = (i + 1) * n // k
        partitions.append(list(range(start, end)))
    return partitions


def generate_random_qubo(variables, coeff_type='lin2', seed=None):
    """
    Disjoint subgraph에 대한 random discrete coefficient QUBO 생성.

    논문 Step 2: "Choose linear and quadratic random coefficients for all
    of the edges and nodes within the disjoint induced subgraphs."

    Args:
        variables: 변수 인덱스 리스트
        coeff_type: 'lin2' or 'lin20'
        seed: 난수 시드

    Returns:
        Q: QUBO 딕셔너리 {(i,j): weight} (i <= j)
    """
    rng = random.Random(seed)
    coeffs = COEFF_LIN2 if coeff_type == 'lin2' else COEFF_LIN20

    Q = {}
    # Linear terms
    for v in variables:
        Q[(v, v)] = rng.choice(coeffs)

    # Quadratic terms (complete graph within subgraph)
    for idx_a in range(len(variables)):
        for idx_b in range(idx_a + 1, len(variables)):
            i, j = variables[idx_a], variables[idx_b]
            Q[(i, j)] = rng.choice(coeffs)

    return Q


def solve_qubo_brute_force(Q, variables):
    """
    Brute force로 subproblem QUBO의 ground state 계산.

    논문 Step 3: "Use an exact classical solver to compute an optimal
    variable assignment for each of the random subproblems."
    (논문은 CPLEX, 여기서는 brute force — 인덱스 사전매핑으로 최적화)

    Args:
        Q: QUBO 딕셔너리
        variables: 변수 인덱스 리스트

    Returns:
        best_assignment: {var_idx: 0 or 1}
        best_energy: float
        num_degenerate: ground state 축퇴도
    """
    n = len(variables)
    if n > 23:
        raise ValueError(f"Brute force 불가: subgraph size {n} > 23 (2^23 = 8M)")

    # 사전 매핑: variable index → bit position shift amount
    var_to_shift = {var: (n - 1 - idx) for idx, var in enumerate(variables)}

    # Q 항을 (shift_i, shift_j, weight) 튜플로 사전 변환
    linear_terms = []  # (shift_i, weight)
    quad_terms = []    # (shift_i, shift_j, weight)
    for (i, j), w in Q.items():
        if i == j:
            linear_terms.append((var_to_shift[i], w))
        else:
            quad_terms.append((var_to_shift[i], var_to_shift[j], w))

    best_energy = float('inf')
    best_bits = 0
    num_degenerate = 0

    for bits in range(1 << n):
        energy = 0.0
        for shift, w in linear_terms:
            if (bits >> shift) & 1:
                energy += w
        for si, sj, w in quad_terms:
            if ((bits >> si) & 1) and ((bits >> sj) & 1):
                energy += w

        if energy < best_energy - 1e-12:
            best_energy = energy
            best_bits = bits
            num_degenerate = 1
        elif abs(energy - best_energy) < 1e-12:
            num_degenerate += 1

    best_assignment = {}
    for idx, var in enumerate(variables):
        best_assignment[var] = (best_bits >> (n - 1 - idx)) & 1

    return best_assignment, best_energy, num_degenerate


def create_qubo_hardened_posiform(n, max_subgraph_size=15, coeff_type='lin2',
                                  posiform_scale=0.1,
                                  posiform_coeff_range=(1.0, 1.0),
                                  seed=None):
    """
    Hardened Posiform Planting QUBO 생성 (메인 진입점).

    논문 Algorithm (Section 2.1):
      1. Partition variables into disjoint subgraphs
      2. Generate random discrete-coefficient QUBO on each subgraph
      3. Solve each subproblem exactly (brute force)
      4. Concatenate solutions → target bitstring
      5. Generate posiform planted QUBO with this target
      6. Q_final = Σ R_i + α × P_posiform

    Args:
        n: 총 변수 수
        max_subgraph_size: 각 random QUBO subgraph 최대 크기 (≤23)
        coeff_type: 'lin2' ({-1,+1}) or 'lin20' ({-1,-0.9,...,0.9,1})
        posiform_scale: α — posiform QUBO 스케일링 (작을수록 어려움)
        posiform_coeff_range: posiform clause 계수 범위 (논문: (1.0, 1.0))
        seed: 난수 시드

    Returns:
        Q_final: 결합된 QUBO 딕셔너리
        info: 메타정보 딕셔너리
    """
    rng = random.Random(seed)

    # Step 1: 변수 분할
    partitions = partition_variables(n, max_subgraph_size)
    k = len(partitions)

    # Step 2 & 3: 각 partition에 random QUBO 생성 + ground state 계산
    random_qubos = []
    target_bits = [0] * n
    total_random_energy = 0.0
    total_degenerate = 1

    for variables in partitions:
        part_seed = rng.randint(0, 10**9)

        R = generate_random_qubo(variables, coeff_type=coeff_type, seed=part_seed)
        assignment, energy, deg = solve_qubo_brute_force(R, variables)

        for var, val in assignment.items():
            target_bits[var] = val

        random_qubos.append(R)
        total_random_energy += energy
        total_degenerate *= deg

    # Step 4: Target bitstring
    target = ''.join(map(str, target_bits))

    # Step 5: Posiform planted QUBO 생성
    posiform_seed = rng.randint(0, 10**9)
    Q_posiform, posiform_info = create_qubo_posiform(
        target,
        coeff_range=posiform_coeff_range,
        seed=posiform_seed
    )

    if not posiform_info['is_unique']:
        # 재시도 (다른 시드)
        for retry in range(5):
            posiform_seed = rng.randint(0, 10**9)
            Q_posiform, posiform_info = create_qubo_posiform(
                target,
                coeff_range=posiform_coeff_range,
                max_clauses_factor=20,
                seed=posiform_seed
            )
            if posiform_info['is_unique']:
                break

    # Step 6: 결합 — Q_final = Σ R_i + α × P
    Q_final = {}

    for R in random_qubos:
        for key, val in R.items():
            Q_final[key] = Q_final.get(key, 0) + val

    for key, val in Q_posiform.items():
        Q_final[key] = Q_final.get(key, 0) + posiform_scale * val

    Q_final = {k: v for k, v in Q_final.items() if abs(v) > 1e-15}

    # 에너지 계산
    target_energy = calculate_energy(target, Q_final)

    info = {
        'n': n,
        'target': target,
        'num_partitions': k,
        'partition_sizes': [len(p) for p in partitions],
        'coeff_type': coeff_type,
        'posiform_scale': posiform_scale,
        'posiform_coeff_range': posiform_coeff_range,
        'posiform_is_unique': posiform_info['is_unique'],
        'posiform_num_clauses': posiform_info['num_clauses'],
        'random_total_energy': total_random_energy,
        'random_total_degenerate': total_degenerate,
        'posiform_energy_at_target': posiform_info['target_energy'],
        'target_energy': target_energy,
        'num_qubo_terms': len(Q_final),
    }

    return Q_final, info


def verify_hardened_ground_state(Q, target, n):
    """
    Hardened QUBO의 ground state 검증.

    n ≤ 20: brute force 완전 검증
    n > 20: single-flip + random sampling
    """
    target_energy = calculate_energy(target, Q)

    if n <= 20:
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
            'num_degenerate': num_degenerate,
        }

    # n > 20: statistical verification
    min_flip_delta = float('inf')
    for i in range(n):
        flipped = list(target)
        flipped[i] = '0' if flipped[i] == '1' else '1'
        delta = calculate_energy(''.join(flipped), Q) - target_energy
        min_flip_delta = min(min_flip_delta, delta)

    lower_count = 0
    num_samples = 10000
    for _ in range(num_samples):
        rand_state = ''.join(str(random.randint(0, 1)) for _ in range(n))
        if calculate_energy(rand_state, Q) < target_energy - 1e-10:
            lower_count += 1

    return {
        'target_energy': target_energy,
        'min_flip_delta': min_flip_delta,
        'is_local_min': min_flip_delta > -1e-10,
        'lower_count': lower_count,
        'is_likely_global': lower_count == 0,
    }


if __name__ == "__main__":
    n = 30
    max_subgraph_size = 15
    coeff_type = 'lin2'
    posiform_scale = 0.1
    seed = 42

    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg.isdigit():
            n = int(arg)

    if len(sys.argv) > 2:
        coeff_type = sys.argv[2]  # lin2 or lin20

    if len(sys.argv) > 3:
        posiform_scale = float(sys.argv[3])

    if len(sys.argv) > 4:
        seed = int(sys.argv[4])

    print("=" * 70)
    print("Hardened Posiform Planting QUBO 생성기")
    print("(Pelofske, Hahn, Djidjev 2024)")
    print("=" * 70)
    print(f"N={n}, subgraph_size={max_subgraph_size}, coeff_type={coeff_type}")
    print(f"posiform_scale(α)={posiform_scale}, seed={seed}")

    t0 = time.time()
    Q, info = create_qubo_hardened_posiform(
        n, max_subgraph_size=max_subgraph_size,
        coeff_type=coeff_type, posiform_scale=posiform_scale,
        seed=seed
    )
    gen_time = time.time() - t0

    target = info['target']
    print(f"\n[생성 결과]")
    print(f"  Target: {target}")
    print(f"  Partitions: {info['num_partitions']}개 (sizes={info['partition_sizes']})")
    print(f"  Posiform unique: {'OK' if info['posiform_is_unique'] else 'FAIL'}")
    print(f"  Posiform clauses: {info['posiform_num_clauses']}")
    print(f"  Random subproblem 총 축퇴도: {info['random_total_degenerate']}")
    print(f"  QUBO 비영 항 수: {info['num_qubo_terms']}")
    print(f"  생성 시간: {gen_time:.2f}s")

    print(f"\n[에너지 분석]")
    print(f"  Σ R_i(x*) = {info['random_total_energy']:.6f}")
    print(f"  α·P(x*)  = {posiform_scale * info['posiform_energy_at_target']:.6f}")
    print(f"  Q(x*)     = {info['target_energy']:.6f}")

    if n <= 12:
        print_q_matrix(Q, n)

    # Ground state 검증
    print(f"\n[Ground State 검증]")
    result = verify_hardened_ground_state(Q, target, n)
    if n <= 20:
        print(f"  Best energy:  {result['best_energy']:.6f}")
        print(f"  Best state:   {result['best_state']}")
        print(f"  축퇴도:       {result['num_degenerate']}")
        print(f"  Ground state: {'OK' if result['is_ground_state'] else 'FAIL'}")
    else:
        print(f"  Target energy:   {result['target_energy']:.6f}")
        print(f"  Min flip delta:  {result['min_flip_delta']:.6f}")
        print(f"  Local minimum:   {'OK' if result['is_local_min'] else 'FAIL'}")
        print(f"  Global (추정):   {'OK' if result['is_likely_global'] else 'FAIL'}")

    # 저장
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(output_dir, exist_ok=True)
    fname = f"hardened_{target[:5]}_{n}_{coeff_type}_a{posiform_scale}.txt"
    output_file = os.path.join(output_dir, fname)
    save_qubo_edgelist(Q, output_file, target)
    print(f"\n[저장 완료] {output_file}")
