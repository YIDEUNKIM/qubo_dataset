"""
Zero-Expectation QUBO SA 실험

N을 점진적으로 늘려가며 SA 성공률, 에너지 정확도, 해밍거리를 측정.
Wishart 실험(test_wishart.py)과 동일한 프레임워크로 비교 가능.
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from zero_expectation.qubo_zero_expectation import create_qubo_precise, DefaultZeroExpectationModel
from qubo_utils import calculate_energy


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def classify_result(target, found_solution, target_energy, found_energy):
    """결과 분류: EXACT / SYM_MATCH / ENERGY_MATCH / FAIL"""
    exact_match = (found_solution == target)
    inverse_target = ''.join('0' if b == '1' else '1' for b in target)
    sym_match = (found_solution == inverse_target)
    energy_match = abs(found_energy - target_energy) < 1e-4

    if exact_match:
        return "EXACT"
    elif sym_match:
        return "SYM_MATCH"
    elif energy_match:
        return "ENERGY_MATCH"
    else:
        return "FAIL"


def run_scaling_experiment(sizes=None, num_runs=10, num_reads=200, num_sweeps=1000):
    """
    N 스케일링 실험: Zero-Expectation QUBO에 대해 N 증가 시 SA 성능 변화.
    """
    if sizes is None:
        sizes = [10, 20, 50, 100, 200, 300, 500]

    print("=" * 90)
    print("Zero-Expectation QUBO — SA Scaling 실험")
    print(f"Sizes={sizes}, num_runs={num_runs}, num_reads={num_reads}, num_sweeps={num_sweeps}")
    print("=" * 90)

    sampler = neal.SimulatedAnnealingSampler()
    model = DefaultZeroExpectationModel()
    results = []

    for n in sizes:
        print(f"\n--- N = {n} ---")
        print(f"{'Run':<4} | {'Result':<12} | {'Target E':<12} | {'Found E':<12} | {'E Ratio':<10} | {'Hamming':<8} | {'Time':<6}")
        print("-" * 75)

        success_count = 0
        energy_ratios = []
        hamming_dists = []
        solve_times = []

        for run in range(num_runs):
            target = ''.join(str(random.randint(0, 1)) for _ in range(n))

            # QUBO 생성 (Zero-Expectation)
            Q = create_qubo_precise(target, density=1.0, model=model)
            target_energy = calculate_energy(target, Q)

            # SA 풀기
            t0 = time.time()
            sampleset = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
            solve_time = time.time() - t0

            best_sample = sampleset.first.sample
            best_energy = sampleset.first.energy
            found = ''.join(str(best_sample[k]) for k in range(n))

            # 분류
            result = classify_result(target, found, target_energy, best_energy)
            if result in ("EXACT", "SYM_MATCH"):
                success_count += 1

            # 메트릭
            ratio = best_energy / target_energy if abs(target_energy) > 1e-10 else float('nan')
            hdist = hamming_distance(target, found)
            energy_ratios.append(ratio)
            hamming_dists.append(hdist)
            solve_times.append(solve_time)

            print(f"{run+1:<4} | {result:<12} | {target_energy:<12.2f} | {best_energy:<12.2f} | {ratio:<10.4f} | {hdist:<8} | {solve_time:<6.2f}")

        # N별 요약
        success_rate = 100.0 * success_count / num_runs
        avg_ratio = np.nanmean(energy_ratios)
        avg_hamming = np.mean(hamming_dists)
        avg_time = np.mean(solve_times)

        print(f"\n  [N={n} 요약]")
        print(f"    성공률: {success_count}/{num_runs} ({success_rate:.1f}%)")
        print(f"    평균 에너지비: {avg_ratio:.4f}")
        print(f"    평균 해밍거리: {avg_hamming:.1f} (N의 {100*avg_hamming/n:.1f}%)")
        print(f"    평균 시간: {avg_time:.2f}s")

        results.append({
            'n': n,
            'success_rate': success_rate,
            'success_count': success_count,
            'avg_energy_ratio': avg_ratio,
            'avg_hamming': avg_hamming,
            'hamming_ratio': avg_hamming / n,
            'avg_time': avg_time,
        })

    # 전체 요약 테이블
    print("\n" + "=" * 90)
    print("전체 요약 (Zero-Expectation QUBO — N Scaling)")
    print("=" * 90)
    print(f"{'N':<8} | {'Success%':<10} | {'Count':<8} | {'Avg E Ratio':<12} | {'Avg Hamming':<12} | {'Hamming/N':<10} | {'Avg Time':<10}")
    print("-" * 85)
    for r in results:
        print(f"{r['n']:<8} | {r['success_rate']:<10.1f} | "
              f"{r['success_count']}/{num_runs:<5} | {r['avg_energy_ratio']:<12.4f} | "
              f"{r['avg_hamming']:<12.1f} | {r['hamming_ratio']:<10.2f} | {r['avg_time']:<10.2f}")

    return results


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    sizes = [10, 20, 50, 100, 200, 300, 500]
    num_runs = 10
    num_reads = 200
    num_sweeps = 1000

    if len(sys.argv) > 1:
        # 커스텀 크기 지정: python3 test_zero_expectation.py 10,20,50,100
        sizes = [int(x) for x in sys.argv[1].split(',')]
    if len(sys.argv) > 2:
        num_runs = int(sys.argv[2])

    run_scaling_experiment(
        sizes=sizes,
        num_runs=num_runs,
        num_reads=num_reads,
        num_sweeps=num_sweeps,
    )
