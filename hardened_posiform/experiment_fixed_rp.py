"""
R/P 고정 α 비교 실험

핵심: 인스턴스당 R_i들과 P를 딱 한 번만 생성하고,
α = {0, 0.001, 0.01, 0.1}만 바꿔서 Q = Σ R_i + α × P 를 결합.
동일한 landscape 위에서 α의 순수한 영향만 관찰한다.
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hardened_posiform.qubo_posiform_hardened import (
    partition_variables, generate_random_qubo, solve_qubo_brute_force
)
from posiform.qubo_posiform import create_qubo_posiform
from qubo_utils import calculate_energy


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


# ============================================================
# Step 1: R과 P를 분리 생성
# ============================================================
def generate_components(n, max_subgraph_size=15, coeff_type='lin2', seed=None):
    """
    R_i 리스트와 P를 각각 반환. 결합(α 적용)은 하지 않는다.

    Returns:
        R_sum: dict — Σ R_i (모든 subgraph random QUBO의 합)
        Q_posiform: dict — posiform QUBO (스케일링 전, 즉 P)
        target: str — 각 R_i 최적해를 concatenate한 비트스트링
    """
    rng = random.Random(seed)

    # 변수 분할
    partitions = partition_variables(n, max_subgraph_size)

    # 각 partition에 random QUBO 생성 + brute-force 풀이
    R_sum = {}
    target_bits = [0] * n

    for variables in partitions:
        part_seed = rng.randint(0, 10**9)
        R = generate_random_qubo(variables, coeff_type=coeff_type, seed=part_seed)
        assignment, energy, deg = solve_qubo_brute_force(R, variables)

        for var, val in assignment.items():
            target_bits[var] = val

        # R_sum에 누적
        for key, val in R.items():
            R_sum[key] = R_sum.get(key, 0) + val

    target = ''.join(map(str, target_bits))

    # Posiform 생성
    posiform_seed = rng.randint(0, 10**9)
    Q_posiform, posiform_info = create_qubo_posiform(
        target, coeff_range=(1.0, 1.0), seed=posiform_seed
    )

    return R_sum, Q_posiform, target


# ============================================================
# Step 2: α를 적용해서 결합
# ============================================================
def combine(R_sum, Q_posiform, alpha):
    """Q = Σ R_i + α × P"""
    Q = dict(R_sum)  # 복사
    for key, val in Q_posiform.items():
        Q[key] = Q.get(key, 0) + alpha * val
    return Q


# ============================================================
# Step 3: 실험 본체
# ============================================================
def run_experiment(n_bits=500, num_instances=10, num_reads=50):
    alpha_list = [0, 0.001, 0.01, 0.1]
    sweep_list = [100, 500, 1000, 5000]
    coeff_types = ['lin2', 'lin20']
    sampler = neal.SimulatedAnnealingSampler()

    print("=" * 100)
    print("R/P 고정 α 비교 실험")
    print(f"N={n_bits}, instances={num_instances}, reads={num_reads}")
    print(f"α: {alpha_list}")
    print(f"sweeps: {sweep_list}")
    print("=" * 100)

    for ct in coeff_types:
        print(f"\n{'━' * 90}")
        print(f"  coeff_type = {ct}")
        print(f"{'━' * 90}")

        # ── 인스턴스별로 R, P 한 번만 생성 ──
        components = []
        t0 = time.time()
        for run in range(num_instances):
            R_sum, Q_posiform, target = generate_components(
                n_bits, coeff_type=ct, seed=run * 53
            )
            components.append((R_sum, Q_posiform, target))
        gen_time = time.time() - t0
        print(f"  R, P 생성 완료: {num_instances}개 인스턴스 ({gen_time:.1f}s)")

        # ── sanity check: α=0.1에서 target이 local min인지 ──
        print(f"\n  [Sanity Check] α=0.1에서 target이 local minimum인지:")
        for idx, (R_sum, Q_posiform, target) in enumerate(components):
            Q = combine(R_sum, Q_posiform, 0.1)
            target_energy = calculate_energy(target, Q)
            min_delta = float('inf')
            for i in range(n_bits):
                flipped = list(target)
                flipped[i] = '0' if flipped[i] == '1' else '1'
                delta = calculate_energy(''.join(flipped), Q) - target_energy
                min_delta = min(min_delta, delta)
            ok = min_delta > -1e-10
            if not ok:
                print(f"    inst {idx}: FAIL (min_flip_delta={min_delta:.6f})")
        print(f"    전부 OK") if all(True for _ in components) else None

        # ── 각 α × sweep 조합으로 SA 실행 ──
        # results[alpha][sweeps] = {'target_found': int, 'total': int, 'hammings': list}
        results = {a: {s: {'target_found': 0, 'total': 0, 'hammings': []}
                       for s in sweep_list}
                   for a in alpha_list}

        for alpha in alpha_list:
            print(f"\n  α = {alpha}")

            for idx, (R_sum, Q_posiform, target) in enumerate(components):
                Q = combine(R_sum, Q_posiform, alpha)

                for sweeps in sweep_list:
                    ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=sweeps)

                    for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                        found = ''.join(str(sample[k]) for k in range(n_bits))
                        h = hamming_distance(target, found)
                        r = results[alpha][sweeps]
                        r['total'] += 1
                        r['hammings'].append(h)
                        if found == target:
                            r['target_found'] += 1

                # 진행 표시
                print(f"    inst {idx} done")

        # ── 결과 출력 ──
        print(f"\n  {'─' * 80}")
        print(f"  [결과] Target Sampling Rate (%)")
        print(f"  {'─' * 80}")
        header = f"  {'α':<8} |"
        for s in sweep_list:
            header += f" sw={s:<5}"
        print(header)
        print(f"  {'-' * (10 + 8 * len(sweep_list))}")

        for alpha in alpha_list:
            row = f"  {alpha:<8} |"
            for s in sweep_list:
                r = results[alpha][s]
                rate = 100.0 * r['target_found'] / r['total'] if r['total'] > 0 else 0
                row += f" {rate:>5.1f}% "
            print(row)

        print(f"\n  [결과] Hamming Distance (median / avg)")
        print(f"  {'─' * 80}")
        header = f"  {'α':<8} |"
        for s in sweep_list:
            header += f"  sw={s:<10}"
        print(header)
        print(f"  {'-' * (10 + 14 * len(sweep_list))}")

        for alpha in alpha_list:
            row = f"  {alpha:<8} |"
            for s in sweep_list:
                h = np.array(results[alpha][s]['hammings'])
                row += f"  {np.median(h):>3.0f}/{h.mean():>5.1f}   "
            print(row)


if __name__ == "__main__":
    n_bits = 500
    num_instances = 10

    if len(sys.argv) > 1:
        n_bits = int(sys.argv[1])
    if len(sys.argv) > 2:
        num_instances = int(sys.argv[2])

    run_experiment(n_bits=n_bits, num_instances=num_instances)
