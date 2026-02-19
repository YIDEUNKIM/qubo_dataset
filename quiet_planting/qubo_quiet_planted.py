"""
Quiet Planting QUBO 생성기 (Planted 3-SAT → QUBO)

Krzakala & Zdeborova (2009) 기반: planted random 3-SAT를 Rosenberg reduction으로
QUBO로 변환. clause density alpha = m/n으로 난이도 제어.

핵심 아이디어:
  1. target이 반드시 만족하는 3-SAT clause만 생성 (quiet planting)
  2. 각 clause의 위반 페널티 z1*z2*z3 (cubic)를 Rosenberg reduction으로 quadratic화
  3. 보조변수 y = z1*z2를 강제하는 penalty 추가 → QUBO 크기 = n + m

상전이:
  - alpha < 3.86: condensation threshold 이하, quiet planting 보장 (구별 불가능)
  - alpha ≈ 4.27: satisfiability threshold (SAT/UNSAT 경계)
  - alpha > 4.27: 랜덤 3-SAT는 UNSAT이지만 planted는 여전히 SAT
"""

import random
import itertools
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from qubo_utils import calculate_energy, save_qubo_edgelist, print_q_matrix, print_qubo_formula


def create_planted_3sat(target, alpha=4.2, seed=None):
    """
    Planted random 3-SAT 생성.

    target이 반드시 만족하는 clause만 생성. 각 clause에 대해:
    - 3개 변수 무작위 선택
    - target 할당과 다른 7가지 violating pattern 중 하나를 균일 샘플

    Args:
        target: 목표 비트스트링 (예: "10110")
        alpha: clause density m/n (난이도 파라미터)
        seed: 난수 시드

    Returns:
        clauses: [(vars_tuple, violating_tuple), ...]
            vars_tuple: (a, b, c) — 변수 인덱스 (정렬됨, a < b < c)
            violating_tuple: (va, vb, vc) — 이 clause를 위반하는 할당
    """
    n = len(target)
    if n < 3:
        raise ValueError(f"변수 수가 3 미만입니다: n={n}")

    m = int(alpha * n)
    rng = random.Random(seed)
    target_bits = [int(b) for b in target]

    all_patterns = [(i, j, k) for i in (0, 1) for j in (0, 1) for k in (0, 1)]
    clauses = []

    for _ in range(m):
        # 3개 변수 무작위 선택 (정렬)
        vars_triple = tuple(sorted(rng.sample(range(n), 3)))
        a, b, c = vars_triple

        # target의 이 변수들에 대한 할당
        target_abc = (target_bits[a], target_bits[b], target_bits[c])

        # violating assignment: target과 다른 7가지 중 균일 샘플
        valid_patterns = [p for p in all_patterns if p != target_abc]
        violating = rng.choice(valid_patterns)

        clauses.append((vars_triple, violating))

    return clauses


def verify_3sat_solution(clauses, assignment):
    """
    할당이 모든 clause를 만족하는지 검증.

    Args:
        clauses: create_planted_3sat의 반환값
        assignment: 비트스트링 (예: "10110")

    Returns:
        (all_satisfied, num_satisfied, total_clauses)
    """
    bits = [int(b) for b in assignment]
    satisfied = 0

    for (a, b, c), (va, vb, vc) in clauses:
        # clause는 할당이 violating assignment와 같을 때만 위반됨
        if (bits[a], bits[b], bits[c]) != (va, vb, vc):
            satisfied += 1

    return satisfied == len(clauses), satisfied, len(clauses)


def clause_to_qubo(clause, aux_idx):
    """
    단일 clause를 Rosenberg reduction으로 QUBO 항으로 변환.

    Rosenberg reduction (P=1):
      g(x,y) = y*z3 + z1*z2 - 2*z1*y - 2*z2*y + 3*y
      여기서 zi = xi (vi=1) 또는 (1-xi) (vi=0)

    Args:
        clause: ((a, b, c), (va, vb, vc))
        aux_idx: 보조변수 y의 인덱스 (>= n)

    Returns:
        Q_terms: {(i, j): weight} — QUBO 기여분 (i <= j)
        constant: 상수항 (에너지 오프셋)
    """
    (a, b, c), (va, vb, vc) = clause
    y = aux_idx
    Q = {}
    constant = 0.0

    def add_term(i, j, val):
        """Q[(min(i,j), max(i,j))] += val"""
        key = (min(i, j), max(i, j))
        Q[key] = Q.get(key, 0) + val

    def add_linear(i, val):
        """Q[(i, i)] += val (선형항)"""
        Q[(i, i)] = Q.get((i, i), 0) + val

    # --- Term 1: y * z_c ---
    # z_c = x_c (vc=1) 또는 (1-x_c) (vc=0)
    if vc == 1:
        # y * x_c
        add_term(c, y, 1.0)
    else:
        # y * (1 - x_c) = y - y*x_c
        add_linear(y, 1.0)
        add_term(c, y, -1.0)

    # --- Term 2: z_a * z_b ---
    if va == 1 and vb == 1:
        # x_a * x_b
        add_term(a, b, 1.0)
    elif va == 1 and vb == 0:
        # x_a * (1 - x_b) = x_a - x_a*x_b
        add_linear(a, 1.0)
        add_term(a, b, -1.0)
    elif va == 0 and vb == 1:
        # (1 - x_a) * x_b = x_b - x_a*x_b
        add_linear(b, 1.0)
        add_term(a, b, -1.0)
    else:
        # (1-x_a) * (1-x_b) = 1 - x_a - x_b + x_a*x_b
        constant += 1.0
        add_linear(a, -1.0)
        add_linear(b, -1.0)
        add_term(a, b, 1.0)

    # --- Term 3: -2 * z_a * y ---
    if va == 1:
        # -2 * x_a * y
        add_term(a, y, -2.0)
    else:
        # -2 * (1-x_a) * y = -2y + 2*x_a*y
        add_linear(y, -2.0)
        add_term(a, y, 2.0)

    # --- Term 4: -2 * z_b * y ---
    if vb == 1:
        # -2 * x_b * y
        add_term(b, y, -2.0)
    else:
        # -2 * (1-x_b) * y = -2y + 2*x_b*y
        add_linear(y, -2.0)
        add_term(b, y, 2.0)

    # --- Term 5: 3 * y ---
    add_linear(y, 3.0)

    return Q, constant


def compute_auxiliary_values(clauses, assignment, n):
    """
    주어진 할당에 대해 최적의 보조변수 값 계산.

    각 clause k의 보조변수 y_k = z_a * z_b (z_i는 clause의 처음 2개 변수에 대응).

    Args:
        clauses: create_planted_3sat의 반환값
        assignment: 원래 n개 변수의 비트스트링
        n: 원래 변수 수

    Returns:
        aux_str: 보조변수 비트스트링 (길이 m)
    """
    bits = [int(b) for b in assignment]
    aux_bits = []

    for (a, b, _c), (va, vb, _vc) in clauses:
        z_a = bits[a] if va == 1 else (1 - bits[a])
        z_b = bits[b] if vb == 1 else (1 - bits[b])
        aux_bits.append(str(z_a * z_b))

    return ''.join(aux_bits)


def create_qubo_quiet_planted(target, alpha=4.2, seed=None,
                               clause_weight_range=None, field_strength=0.0):
    """
    Quiet Planting QUBO 생성 (메인 진입점).

    planted 3-SAT → Rosenberg reduction → QUBO.

    Args:
        target: 목표 비트스트링 (예: "10110")
        alpha: clause density m/n
        seed: 난수 시드
        clause_weight_range: clause별 랜덤 가중치 범위 (예: (1.0, 3.0)).
            None이면 모든 clause에 동일 가중치 1.0 적용.
            주의: clause 가중치만으로는 SAT 해 간 축퇴를 깰 수 없음
            (모든 SAT 해가 모든 clause를 만족 → 0 * w = 0)
        field_strength: planted field 강도. 각 변수에 target 방향으로
            미세한 선형 편향 추가. 이것이 실제로 축퇴를 깸.
            0이면 field 없음. 권장: 0.1 ~ 1.0

    Returns:
        Q: QUBO 딕셔너리 {(i,j): weight} (i <= j)
        clauses: 생성된 3-SAT clause 리스트
        info: 메타정보 딕셔너리
    """
    n = len(target)
    rng = random.Random(seed)
    clauses = create_planted_3sat(target, alpha=alpha, seed=seed)
    m = len(clauses)

    Q = {}
    total_constant = 0.0

    for k, clause in enumerate(clauses):
        aux_idx = n + k
        q_terms, const = clause_to_qubo(clause, aux_idx)

        # clause별 랜덤 가중치 적용
        if clause_weight_range is not None:
            w_lo, w_hi = clause_weight_range
            w_k = rng.uniform(w_lo, w_hi)
        else:
            w_k = 1.0

        for key, val in q_terms.items():
            Q[key] = Q.get(key, 0) + val * w_k
        total_constant += const * w_k

    # planted field: target 방향으로 미세한 선형 편향
    if field_strength > 0:
        for i in range(n):
            eps_i = rng.uniform(field_strength * 0.5, field_strength * 1.5)
            if target[i] == '1':
                # x_i = 1을 선호: -eps * x_i (x_i=1일 때 에너지 감소)
                Q[(i, i)] = Q.get((i, i), 0) - eps_i
            else:
                # x_i = 0을 선호: +eps * x_i (x_i=0일 때 에너지 감소)
                Q[(i, i)] = Q.get((i, i), 0) + eps_i

    # 0에 가까운 항 제거
    Q = {k: v for k, v in Q.items() if abs(v) > 1e-15}

    # 검증: target이 모든 clause를 만족하는지
    all_sat, num_sat, total = verify_3sat_solution(clauses, target)
    assert all_sat, f"target이 {total - num_sat}개 clause를 위반!"

    # target + 최적 보조변수에서의 에너지 계산
    aux_str = compute_auxiliary_values(clauses, target, n)
    full_target = target + aux_str
    target_energy = calculate_energy(full_target, Q)

    info = {
        'n': n,
        'm': m,
        'total_vars': n + m,
        'alpha': alpha,
        'constant_offset': total_constant,
        'target_energy': target_energy,
        'clauses': clauses,
        'clause_weight_range': clause_weight_range,
        'field_strength': field_strength,
    }

    return Q, clauses, info


def extract_original_solution(sample, n):
    """
    SA 결과에서 원래 n개 변수만 추출.

    Args:
        sample: SA sample dict {var_idx: 0 or 1}
        n: 원래 변수 수

    Returns:
        비트스트링 (길이 n)
    """
    return ''.join(str(sample.get(k, 0)) for k in range(n))


def verify_ground_state(Q, target, clauses, n, num_random_samples=10000):
    """
    target + 최적 aux가 ground state인지 검증.

    1. Single-flip 이웃 검사 (원래 변수 + 보조변수 모두)
    2. 랜덤 샘플 검사

    Returns:
        (is_local_min, is_likely_global, stats_dict)
    """
    m = len(clauses)
    total_n = n + m
    aux_str = compute_auxiliary_values(clauses, target, n)
    full_target = target + aux_str
    target_energy = calculate_energy(full_target, Q)

    # 1. Single-flip 이웃 검사
    min_flip_delta = float('inf')
    for i in range(total_n):
        flipped = list(full_target)
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
        # 원래 변수는 랜덤, 보조변수는 최적값 사용
        rand_orig = ''.join(str(random.randint(0, 1)) for _ in range(n))
        rand_aux = compute_auxiliary_values(clauses, rand_orig, n)
        rand_state = rand_orig + rand_aux
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
        'full_target': full_target,
        'total_vars': total_n,
    }

    return is_local_min, is_likely_global, stats


def verify_brute_force(Q, target, clauses, n):
    """
    Brute force로 ground state 완전 검증 (n+m <= 20일 때만).

    Returns:
        결과 딕셔너리 또는 None (크기 초과 시)
    """
    m = len(clauses)
    total_n = n + m

    if total_n > 20:
        return None

    aux_str = compute_auxiliary_values(clauses, target, n)
    full_target = target + aux_str
    target_energy = calculate_energy(full_target, Q)

    best_energy = float('inf')
    best_state = None
    num_degenerate = 0

    for bits in itertools.product([0, 1], repeat=total_n):
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
        'best_original': best_state[:n] if best_state else None,
        'is_ground_state': abs(best_energy - target_energy) < 1e-10,
        'energy_gap': target_energy - best_energy if not abs(best_energy - target_energy) < 1e-10 else 0.0,
        'num_degenerate': num_degenerate,
        'total_vars': total_n,
    }


if __name__ == "__main__":
    target = "10110"
    alpha = 4.2
    seed = None

    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if set(arg) <= {'0', '1'} and len(arg) >= 3:
            target = arg
        elif arg.isdigit():
            length = int(arg)
            if length < 3:
                length = 3
            random.seed(42)
            target = ''.join(str(random.randint(0, 1)) for _ in range(length))
            print(f"[설정] 길이 {length}의 랜덤 목표 해 생성")
        else:
            target = arg

    if len(sys.argv) > 2:
        alpha = float(sys.argv[2])

    if len(sys.argv) > 3:
        seed = int(sys.argv[3])

    n = len(target)
    m = int(alpha * n)

    print("=" * 60)
    print("Quiet Planting QUBO 생성기 (Planted 3-SAT → Rosenberg)")
    print("=" * 60)
    print(f"Target: {target} (N={n})")
    print(f"Alpha: {alpha} (M={m}, QUBO 크기={n + m})")
    if seed is not None:
        print(f"Seed: {seed}")

    # QUBO 생성
    Q, clauses, info = create_qubo_quiet_planted(target, alpha=alpha, seed=seed)

    print(f"\n[3-SAT 정보]")
    print(f"  Clause 수: {info['m']}")
    print(f"  전체 변수 수: {info['total_vars']} (원래 {n} + 보조 {info['m']})")
    print(f"  QUBO 비영 항 수: {len(Q)}")

    # 3-SAT 검증
    all_sat, num_sat, total = verify_3sat_solution(clauses, target)
    print(f"\n[3-SAT 검증]")
    print(f"  Target 만족: {num_sat}/{total} clauses ({'OK' if all_sat else 'FAIL'})")

    # 행렬 출력 (작을 때만)
    if info['total_vars'] <= 12:
        print_q_matrix(Q, info['total_vars'])
        print_qubo_formula(Q)

    # 에너지 검증
    aux_str = compute_auxiliary_values(clauses, target, n)
    full_target = target + aux_str
    target_energy = calculate_energy(full_target, Q)
    print(f"\n[에너지 검증]")
    print(f"  Target 에너지: {target_energy:.6f}")
    print(f"  보조변수 값: {aux_str[:20]}{'...' if len(aux_str) > 20 else ''}")

    # Brute force 검증
    total_n = info['total_vars']
    if total_n <= 20:
        print(f"\n[Brute Force 검증] (전체 변수 {total_n}개)")
        bf = verify_brute_force(Q, target, clauses, n)
        print(f"  Target 에너지: {bf['target_energy']:.6f}")
        print(f"  최소 에너지:   {bf['best_energy']:.6f}")
        print(f"  최소 상태:     {bf['best_state']}")
        print(f"  원래 변수:     {bf['best_original']}")
        print(f"  축퇴도:        {bf['num_degenerate']}")
        if bf['is_ground_state']:
            print("  Ground State 검증 성공!")
        else:
            print(f"  Ground State 검증 실패 (gap: {bf['energy_gap']:.6f})")
    else:
        print(f"\n[통계적 검증] (전체 변수 {total_n}개)")
        is_local, is_global, stats = verify_ground_state(Q, target, clauses, n)
        print(f"  Target 에너지: {stats['target_energy']:.6f}")
        print(f"  최소 flip delta: {stats['min_flip_delta']:.6f}")
        print(f"  Local minimum: {'OK' if is_local else 'FAIL'}")
        print(f"  랜덤 {stats['num_random_samples']}개 중 더 낮은 에너지: {stats['lower_count']}개")
        print(f"  Global minimum (추정): {'OK' if is_global else 'FAIL'}")

    # 결과 저장
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(output_dir, exist_ok=True)
    filename_target = f"{target[:5]}_{n}"
    output_file = os.path.join(output_dir, f"quiet_{filename_target}_a{alpha}.txt")
    save_qubo_edgelist(Q, output_file, full_target)
    print(f"\n[저장 완료] {output_file}")
