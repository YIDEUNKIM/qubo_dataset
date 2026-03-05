"""
R/P 고정 α 비교 실험 — 에너지 분석

experiment_fixed_rp.py와 동일한 R/P 생성 후,
SA 해의 에너지가 target 에너지와 일치하는지 분석.

핵심 측정:
  - GS 에너지 도달률: SA 해의 에너지 == target 에너지인 비율
  - Target 일치률: SA 해의 비트스트링 == target인 비율
  - GS인데 target 아닌 수: 축퇴 지표
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


def generate_components(n, max_subgraph_size=15, coeff_type='lin2', seed=None):
    """R_sum과 Q_posiform을 분리 반환."""
    rng = random.Random(seed)
    partitions = partition_variables(n, max_subgraph_size)
    R_sum = {}
    target_bits = [0] * n

    for variables in partitions:
        part_seed = rng.randint(0, 10**9)
        R = generate_random_qubo(variables, coeff_type=coeff_type, seed=part_seed)
        assignment, energy, deg = solve_qubo_brute_force(R, variables)
        for var, val in assignment.items():
            target_bits[var] = val
        for key, val in R.items():
            R_sum[key] = R_sum.get(key, 0) + val

    target = ''.join(map(str, target_bits))
    posiform_seed = rng.randint(0, 10**9)
    Q_posiform, _ = create_qubo_posiform(target, coeff_range=(1.0, 1.0), seed=posiform_seed)
    return R_sum, Q_posiform, target


def combine(R_sum, Q_posiform, alpha):
    """Q = Σ R_i + α × P"""
    Q = dict(R_sum)
    for key, val in Q_posiform.items():
        Q[key] = Q.get(key, 0) + alpha * val
    return Q


def run_energy_analysis(n_bits=500, num_instances=100, num_reads=50):
    alpha_list = [0, 0.001, 0.01, 0.1]
    sweep_list = [1000, 5000]
    coeff_types = ['lin2', 'lin20']
    sampler = neal.SimulatedAnnealingSampler()

    print("=" * 100)
    print("R/P 고정 α 비교 — 에너지 분석")
    print(f"N={n_bits}, instances={num_instances}, reads={num_reads}")
    print(f"α: {alpha_list}, sweeps: {sweep_list}")
    print("=" * 100)

    for ct in coeff_types:
        print(f"\n{'━' * 90}")
        print(f"  coeff_type = {ct}")
        print(f"{'━' * 90}")

        # R, P 생성
        components = []
        t0 = time.time()
        for run in range(num_instances):
            R_sum, Q_posiform, target = generate_components(
                n_bits, coeff_type=ct, seed=run * 53
            )
            components.append((R_sum, Q_posiform, target))
        print(f"  R, P 생성 완료: {num_instances}개 ({time.time()-t0:.1f}s)")

        for sweeps in sweep_list:
            print(f"\n  ── sweep = {sweeps} ──")
            print(f"  {'α':<8} | {'GS에너지%':>9} | {'Target%':>9} | "
                  f"{'GS≠Target':>10} | {'Ham med':>8} | {'Ham avg':>8}")
            print(f"  {'-' * 70}")

            for alpha in alpha_list:
                gs_energy_match = 0
                target_match = 0
                total = 0
                hammings = []

                t0 = time.time()
                for idx, (R_sum, Q_posiform, target) in enumerate(components):
                    Q = combine(R_sum, Q_posiform, alpha)
                    target_energy = calculate_energy(target, Q)

                    ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=sweeps)

                    for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                        total += 1
                        found = ''.join(str(sample[k]) for k in range(n_bits))

                        if abs(energy - target_energy) < 1e-6:
                            gs_energy_match += 1
                        if found == target:
                            target_match += 1

                        hammings.append(hamming_distance(target, found))

                elapsed = time.time() - t0
                h = np.array(hammings)
                gs_pct = 100.0 * gs_energy_match / total
                tgt_pct = 100.0 * target_match / total
                gs_not_tgt = gs_energy_match - target_match

                print(f"  {alpha:<8} | {gs_pct:>8.1f}% | {tgt_pct:>8.1f}% | "
                      f"{gs_not_tgt:>10} | {np.median(h):>8.0f} | {h.mean():>8.1f}"
                      f"  ({elapsed:.0f}s)")


if __name__ == "__main__":
    n_bits = 500
    num_instances = 100

    if len(sys.argv) > 1:
        n_bits = int(sys.argv[1])
    if len(sys.argv) > 2:
        num_instances = int(sys.argv[2])

    run_energy_analysis(n_bits=n_bits, num_instances=num_instances)
