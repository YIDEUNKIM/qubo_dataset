"""
α=0에서 SA가 ground state 에너지는 찾는지 확인.

가설: SA가 target을 못 찾는 이유는 ground state를 못 찾아서가 아니라,
      축퇴된 ground state 중 target이 아닌 다른 해를 찾기 때문.

검증: SA가 찾은 해의 에너지 vs target 에너지 비교.
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


def run_energy_analysis(n_bits=500, num_instances=10, num_reads=50, num_sweeps=1000):
    """
    α=0, 0.01, 0.1에서 SA가 찾은 해의 에너지를 target 에너지와 비교.
    """
    configs = [
        ('lin2', 0.0),
        ('lin2', 0.01),
        ('lin2', 0.1),
        ('lin20', 0.0),
        ('lin20', 0.01),
        ('lin20', 0.1),
    ]

    print("=" * 120)
    print(f"SA 에너지 분석: ground state를 찾지만 target이 아닌 해를 찾는 것인지 검증")
    print(f"N={n_bits}, num_instances={num_instances}, num_reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 120)

    sampler = neal.SimulatedAnnealingSampler()

    for ct, alpha in configs:
        print(f"\n{'─' * 100}")
        print(f"  {ct}, α={alpha}")
        print(f"{'─' * 100}")

        exact_gs_count = 0       # 에너지가 target과 정확히 같은 해
        target_found_count = 0   # target 자체를 찾은 횟수
        total_samples = 0
        energy_diffs = []        # (SA 에너지 - target 에너지) 리스트
        hamming_when_exact = []  # 에너지가 같을 때 hamming distance

        t0 = time.time()
        for run in range(num_instances):
            Q, info = create_qubo_hardened_posiform(
                n_bits, max_subgraph_size=15, coeff_type=ct,
                posiform_scale=alpha, seed=run * 53
            )
            target = info['target']
            target_energy = info['target_energy']
            degeneracy = info['random_total_degenerate']

            ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)

            instance_exact = 0
            instance_target = 0

            for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                total_samples += 1
                found = ''.join(str(sample[k]) for k in range(n_bits))
                diff = energy - target_energy

                energy_diffs.append(diff)

                if abs(diff) < 1e-10:
                    exact_gs_count += 1
                    instance_exact += 1
                    hd = hamming_distance(target, found)
                    hamming_when_exact.append(hd)

                if found == target:
                    target_found_count += 1
                    instance_target += 1

            print(f"  inst {run}: 축퇴도={degeneracy:>20,} | "
                  f"GS에너지 일치: {instance_exact}/{num_reads} | "
                  f"target 일치: {instance_target}/{num_reads} | "
                  f"target_E={target_energy:.4f}")

        elapsed = time.time() - t0
        energy_diffs = np.array(energy_diffs)

        print(f"\n  *** 종합 ({ct}, α={alpha}) ***")
        print(f"  전체 샘플: {total_samples}")
        print(f"  Ground state 에너지 일치: {exact_gs_count}/{total_samples} "
              f"({100.0*exact_gs_count/total_samples:.1f}%)")
        print(f"  Target bitstring 일치:    {target_found_count}/{total_samples} "
              f"({100.0*target_found_count/total_samples:.1f}%)")

        if exact_gs_count > 0 and exact_gs_count > target_found_count:
            print(f"  → GS 에너지 찾았지만 target 아닌 경우: "
                  f"{exact_gs_count - target_found_count}개")
            print(f"  → GS 에너지 일치 시 평균 Hamming: {np.mean(hamming_when_exact):.1f}")
            # hamming 분포
            h_arr = np.array(hamming_when_exact)
            print(f"     Hamming 분포: min={h_arr.min()}, max={h_arr.max()}, "
                  f"median={np.median(h_arr):.0f}")
            print(f"     Hamming=0 (target): {np.sum(h_arr == 0)}개")

        print(f"  에너지 차이 (SA - target): "
              f"mean={energy_diffs.mean():.4f}, min={energy_diffs.min():.4f}, "
              f"median={np.median(energy_diffs):.4f}")
        near_zero = np.sum(np.abs(energy_diffs) < 1e-6)
        print(f"  |차이| < 1e-6인 샘플: {near_zero}/{total_samples}")
        print(f"  시간: {elapsed:.1f}s")

    # 최종 요약
    print("\n" + "=" * 120)
    print("최종 요약")
    print("=" * 120)
    print(f"{'Config':<16} | {'GS에너지일치':<16} | {'Target일치':<16} | {'차이(GS찾았지만 target아님)'}")
    print("-" * 80)


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    n_bits = 500
    num_instances = 10

    if len(sys.argv) > 1:
        n_bits = int(sys.argv[1])
    if len(sys.argv) > 2:
        num_instances = int(sys.argv[2])

    run_energy_analysis(n_bits=n_bits, num_instances=num_instances,
                        num_reads=50, num_sweeps=1000)
