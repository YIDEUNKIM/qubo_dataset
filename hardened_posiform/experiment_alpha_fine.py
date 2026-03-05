"""
α=0.001 fine-grained 실험

α=0 (축퇴 지배) → α=0.01 (genuine hardness) 사이의 전환점 탐색.
α=0.001에서 posiform signal이 축퇴를 깨기에 충분한지 확인.
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


def run_experiment(n_bits=500, num_instances=10, reads_per_sweep=50):
    sweep_values = [1, 5, 10, 50, 100, 500, 1000, 5000]
    alpha_values = [0.0, 0.001, 0.01, 0.1]
    coeff_types = ['lin2', 'lin20']

    configs = [(ct, a) for ct in coeff_types for a in alpha_values]

    print("=" * 120)
    print("α=0.001 Fine-Grained 실험")
    print(f"N={n_bits}, num_instances={num_instances}, reads_per_sweep={reads_per_sweep}")
    print(f"α values: {alpha_values}")
    print("=" * 120)

    sampler = neal.SimulatedAnnealingSampler()
    all_results = []

    for ct, alpha in configs:
        print(f"\n{'━' * 100}")
        print(f"  {ct}, α={alpha}")
        print(f"{'━' * 100}")

        t0 = time.time()
        instances = []
        for run in range(num_instances):
            Q, info = create_qubo_hardened_posiform(
                n_bits, max_subgraph_size=15, coeff_type=ct,
                posiform_scale=alpha, seed=run * 53
            )
            instances.append((Q, info))
        gen_time = time.time() - t0
        print(f"  인스턴스 생성: {len(instances)}개 ({gen_time:.1f}s)")

        degeneracies = [info['random_total_degenerate'] for _, info in instances]
        print(f"  Random subproblem 평균 축퇴도: {np.mean(degeneracies):.1f}")

        # local min 검사
        print(f"\n  [Target GS 여부 검사]")
        for idx, (Q, info) in enumerate(instances):
            target = info['target']
            target_energy = calculate_energy(target, Q)
            min_flip_delta = float('inf')
            for i in range(n_bits):
                flipped = list(target)
                flipped[i] = '0' if flipped[i] == '1' else '1'
                delta = calculate_energy(''.join(flipped), Q) - target_energy
                min_flip_delta = min(min_flip_delta, delta)
            is_local_min = min_flip_delta > -1e-10
            print(f"    inst {idx}: target_E={target_energy:>10.4f} | "
                  f"min_flip_delta={min_flip_delta:>+10.6f} | "
                  f"local_min={'O' if is_local_min else 'X'}")

        # sweep 전이
        for sweeps in sweep_values:
            total_samples = 0
            target_found = 0
            hamming_list = []

            t0 = time.time()
            for Q, info in instances:
                target = info['target']
                ss = sampler.sample_qubo(Q, num_reads=reads_per_sweep,
                                         num_sweeps=sweeps)
                for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                    total_samples += 1
                    found = ''.join(str(sample[k]) for k in range(n_bits))
                    if found == target:
                        target_found += 1
                    hamming_list.append(hamming_distance(target, found))

            solve_time = time.time() - t0
            rate = 100.0 * target_found / total_samples
            h_arr = np.array(hamming_list)

            print(f"  sw={sweeps:<6} | target rate: {target_found:>4}/{total_samples} "
                  f"({rate:>6.2f}%) | Hamming avg={h_arr.mean():>6.1f} "
                  f"med={np.median(h_arr):>5.0f} "
                  f"min={h_arr.min():>3} max={h_arr.max():>3} | {solve_time:.1f}s")

            all_results.append({
                'coeff_type': ct, 'alpha': alpha, 'sweeps': sweeps,
                'target_rate': rate, 'avg_hamming': float(h_arr.mean()),
                'median_hamming': float(np.median(h_arr)),
            })

    # 에너지 분석 (sweep=1000 고정)
    print("\n" + "=" * 120)
    print("에너지 분석 (sweep=1000)")
    print("=" * 120)

    energy_summary = []
    for ct, alpha in configs:
        print(f"\n  [{ct}, α={alpha}]")

        gs_energy_match = 0
        target_match = 0
        total_samples = 0
        hamming_all = []

        for run in range(num_instances):
            Q, info = create_qubo_hardened_posiform(
                n_bits, max_subgraph_size=15, coeff_type=ct,
                posiform_scale=alpha, seed=run * 53
            )
            target = info['target']
            target_energy = info['target_energy']

            ss = sampler.sample_qubo(Q, num_reads=reads_per_sweep, num_sweeps=1000)
            for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                total_samples += 1
                found = ''.join(str(sample[k]) for k in range(n_bits))
                if abs(energy - target_energy) < 1e-10:
                    gs_energy_match += 1
                if found == target:
                    target_match += 1
                hamming_all.append(hamming_distance(target, found))

        gs_pct = 100.0 * gs_energy_match / total_samples
        tgt_pct = 100.0 * target_match / total_samples
        h_arr = np.array(hamming_all)

        print(f"  GS 에너지 일치: {gs_energy_match}/{total_samples} ({gs_pct:.1f}%)")
        print(f"  Target 일치:    {target_match}/{total_samples} ({tgt_pct:.1f}%)")
        print(f"  GS찾았지만 target 아닌 수: {gs_energy_match - target_match}")
        print(f"  Hamming: mean={h_arr.mean():.1f}, median={np.median(h_arr):.0f}")

        energy_summary.append({
            'config': f"{ct},α={alpha}",
            'gs_pct': gs_pct, 'tgt_pct': tgt_pct,
            'hamming_med': float(np.median(h_arr)),
            'gs_not_target': gs_energy_match - target_match,
        })

    # 요약 테이블
    print("\n" + "=" * 120)
    print("요약 1: Target Sampling Rate (%)")
    print("=" * 120)
    for ct in coeff_types:
        print(f"\n  [{ct}]")
        header = f"  {'α':<8} |"
        for s in sweep_values:
            header += f" sw={s:<5}"
        print(header)
        print(f"  " + "-" * (10 + 8 * len(sweep_values)))
        for alpha in alpha_values:
            row = f"  {alpha:<8} |"
            for s in sweep_values:
                matching = [r for r in all_results
                            if r['coeff_type'] == ct and r['alpha'] == alpha
                            and r['sweeps'] == s]
                if matching:
                    row += f" {matching[0]['target_rate']:>5.1f}% "
                else:
                    row += f" {'N/A':>6} "
            print(row)

    print("\n" + "=" * 120)
    print("요약 2: Median Hamming Distance (sw=5000)")
    print("=" * 120)
    for ct in coeff_types:
        print(f"\n  [{ct}]")
        for alpha in alpha_values:
            matching = [r for r in all_results
                        if r['coeff_type'] == ct and r['alpha'] == alpha
                        and r['sweeps'] == 5000]
            if matching:
                print(f"  α={alpha:<8} | Hamming med={matching[0]['median_hamming']:>6.0f} "
                      f"avg={matching[0]['avg_hamming']:>6.1f}")

    print("\n" + "=" * 120)
    print("요약 3: 에너지 분석 (sweep=1000)")
    print("=" * 120)
    print(f"  {'Config':<16} | {'GS에너지%':>8} | {'Target%':>8} | "
          f"{'Ham med':>7} | {'GS찾았지만 target아닌':>20}")
    print(f"  " + "-" * 80)
    for s in energy_summary:
        print(f"  {s['config']:<16} | {s['gs_pct']:>7.1f}% | {s['tgt_pct']:>7.1f}% | "
              f"{s['hamming_med']:>7.0f} | {s['gs_not_target']:>20}")


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    n_bits = 500
    num_instances = 10

    if len(sys.argv) > 1:
        n_bits = int(sys.argv[1])
    if len(sys.argv) > 2:
        num_instances = int(sys.argv[2])

    run_experiment(n_bits=n_bits, num_instances=num_instances, reads_per_sweep=50)
