"""
α=0 vs α=0.01 vs α=0.1 비교 실험

Hardened Posiform에서 posiform_scale(α)를 0으로 설정하면
Q_final = Σ R_i (순수 random subproblem QUBO, posiform 기여 없음).

비교 조건:
  - N=500, subgraph_size=15
  - coeff_type: lin2, lin20
  - α: 0, 0.01, 0.1
  - SA sweep 수 변화에 따른 성공률 전이 (S-curve)
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hardened_posiform.qubo_posiform_hardened import create_qubo_hardened_posiform
from qubo_utils import calculate_energy


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def run_alpha_comparison(n_bits=500, num_instances=10, reads_per_sweep=50):
    """
    α=0 vs α=0.01 vs α=0.1 sweep 전이 비교.
    """
    sweep_values = [1, 5, 10, 50, 100, 500, 1000, 5000]
    configs = [
        ('lin2', 0.0),
        ('lin2', 0.01),
        ('lin2', 0.1),
        ('lin20', 0.0),
        ('lin20', 0.01),
        ('lin20', 0.1),
    ]

    print("=" * 120)
    print(f"α=0 vs α=0.01 vs α=0.1 비교 실험 (Sweep 전이)")
    print(f"N={n_bits}, num_instances={num_instances}, reads_per_sweep={reads_per_sweep}")
    print(f"Sweep values: {sweep_values}")
    print("=" * 120)

    sampler = neal.SimulatedAnnealingSampler()
    all_results = []

    for ct, alpha in configs:
        print(f"\n{'─' * 100}")
        print(f"  {ct}, α={alpha}")
        print(f"{'─' * 100}")

        # QUBO 인스턴스 미리 생성
        t0 = time.time()
        instances = []
        for run in range(num_instances):
            Q, info = create_qubo_hardened_posiform(
                n_bits, max_subgraph_size=15, coeff_type=ct,
                posiform_scale=alpha, seed=run * 53
            )
            # α=0일 때도 posiform_is_unique 체크 (생성 자체는 됨)
            if alpha == 0.0 or info['posiform_is_unique']:
                instances.append((Q, info))
        gen_time = time.time() - t0
        actual = len(instances)
        print(f"  인스턴스 생성: {actual}/{num_instances} (생성 시간: {gen_time:.1f}s)")

        # 축퇴도 통계 (α=0일 때 특히 중요)
        degeneracies = [info['random_total_degenerate'] for _, info in instances]
        print(f"  Random subproblem 평균 축퇴도: {np.mean(degeneracies):.1f} "
              f"(min={min(degeneracies)}, max={max(degeneracies)})")

        for sweeps in sweep_values:
            total_samples = 0
            ground_state_found = 0
            total_hamming = 0

            t0 = time.time()
            for Q, info in instances:
                target = info['target']

                ss = sampler.sample_qubo(Q, num_reads=reads_per_sweep,
                                         num_sweeps=sweeps)

                for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                    total_samples += 1
                    found = ''.join(str(sample[k]) for k in range(n_bits))
                    if found == target:
                        ground_state_found += 1
                    total_hamming += hamming_distance(target, found)

            solve_time = time.time() - t0
            rate = 100.0 * ground_state_found / total_samples if total_samples > 0 else 0
            avg_h = total_hamming / total_samples if total_samples > 0 else float('nan')

            print(f"  sweeps={sweeps:<6} | GS rate: {ground_state_found:>4}/{total_samples} "
                  f"({rate:>6.2f}%) | avg Hamming: {avg_h:>6.1f} | time: {solve_time:.1f}s")

            all_results.append({
                'coeff_type': ct, 'scale': alpha, 'sweeps': sweeps,
                'gs_rate': rate, 'avg_hamming': avg_h,
                'gs_found': ground_state_found, 'total_samples': total_samples,
            })

    # ==========================================
    # 요약 테이블 1: Ground-State Sampling Rate
    # ==========================================
    print("\n" + "=" * 120)
    print("요약 1: Ground-State Sampling Rate (%)")
    print("=" * 120)
    header = f"{'Config':<16} |"
    for s in sweep_values:
        header += f" {s:>7}"
    print(header)
    print("-" * (18 + 8 * len(sweep_values)))

    for ct, alpha in configs:
        label = f"{ct},α={alpha}"
        row = f"{label:<16} |"
        for s in sweep_values:
            matching = [r for r in all_results
                        if r['coeff_type'] == ct and r['scale'] == alpha
                        and r['sweeps'] == s]
            if matching:
                row += f" {matching[0]['gs_rate']:>6.1f}%"
            else:
                row += f" {'N/A':>7}"
        print(row)

    # ==========================================
    # 요약 테이블 2: Average Hamming Distance
    # ==========================================
    print("\n" + "=" * 120)
    print("요약 2: Average Hamming Distance")
    print("=" * 120)
    header = f"{'Config':<16} |"
    for s in sweep_values:
        header += f" {s:>7}"
    print(header)
    print("-" * (18 + 8 * len(sweep_values)))

    for ct, alpha in configs:
        label = f"{ct},α={alpha}"
        row = f"{label:<16} |"
        for s in sweep_values:
            matching = [r for r in all_results
                        if r['coeff_type'] == ct and r['scale'] == alpha
                        and r['sweeps'] == s]
            if matching:
                row += f" {matching[0]['avg_hamming']:>7.1f}"
            else:
                row += f" {'N/A':>7}"
        print(row)

    # ==========================================
    # 요약 테이블 3: α별 비교 (lin2 / lin20 각각)
    # ==========================================
    print("\n" + "=" * 120)
    print("요약 3: α값 비교 — 난이도 순서")
    print("=" * 120)

    for ct in ['lin2', 'lin20']:
        print(f"\n  [{ct}]")
        print(f"  {'α':<6} | ", end="")
        for s in sweep_values:
            print(f" sw={s:<5}", end="")
        print()
        print(f"  " + "-" * (8 + 8 * len(sweep_values)))

        for alpha in [0.0, 0.01, 0.1]:
            print(f"  {alpha:<6} | ", end="")
            for s in sweep_values:
                matching = [r for r in all_results
                            if r['coeff_type'] == ct and r['scale'] == alpha
                            and r['sweeps'] == s]
                if matching:
                    print(f" {matching[0]['gs_rate']:>5.1f}%", end=" ")
                else:
                    print(f" {'N/A':>6}", end=" ")
            print()

    return all_results


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    n_bits = 500
    num_instances = 10

    if len(sys.argv) > 1:
        n_bits = int(sys.argv[1])
    if len(sys.argv) > 2:
        num_instances = int(sys.argv[2])

    run_alpha_comparison(n_bits=n_bits, num_instances=num_instances, reads_per_sweep=50)
