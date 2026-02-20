"""
Quiet Planting SA 실험 프레임워크

실험 1: Alpha sweep — clause density에 따른 SA 성공률 변화
실험 2: Scaling — 고정 alpha에서 N 증가에 따른 성공률 변화
실험 3: 3-way 비교 — Quiet vs Wishart vs ZeroExp
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from quiet_planting.qubo_quiet_planted import (
    create_qubo_quiet_planted,
    compute_auxiliary_values,
    extract_original_solution,
)
from wishart.qubo_wishart import create_qubo_wishart
from qubo_utils import calculate_energy
from zero_expectation.qubo_zero_expectation import create_qubo_precise


def classify_result(target, found_solution, target_energy, found_energy):
    """결과 분류: EXACT / ENERGY_MATCH / FAIL"""
    if found_solution == target:
        return "EXACT"
    if abs(found_energy - target_energy) < 1e-4:
        return "ENERGY_MATCH"
    return "FAIL"


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def run_alpha_sweep(n_bits=50, alpha_values=None, num_runs=10,
                    num_reads=200, num_sweeps=None):
    """
    Alpha sweep 실험: clause density에 따른 SA 성공률 측정.

    3-SAT 상전이 관련 alpha 값:
      - alpha < 3.86: condensation 이하, quiet planting 보장
      - alpha ~ 4.27: SAT/UNSAT 경계
      - alpha > 4.27: planted는 여전히 SAT, 하지만 어려움
    """
    if alpha_values is None:
        alpha_values = [2.0, 3.0, 3.5, 3.86, 4.0, 4.2, 4.5, 5.0]

    print("=" * 95)
    print(f"Quiet Planting — Alpha Sweep 실험")
    print(f"N={n_bits}, num_runs={num_runs}, num_reads={num_reads}")
    print("=" * 95)

    sampler = neal.SimulatedAnnealingSampler()
    results_summary = []

    for alpha in alpha_values:
        m = int(alpha * n_bits)
        total_vars = n_bits + m
        sweeps = num_sweeps or max(1000, 10 * total_vars)

        print(f"\n--- alpha={alpha} (M={m}, QUBO 크기={total_vars}, sweeps={sweeps}) ---")
        print(f"{'Run':<4} | {'Result':<12} | {'Target E':<12} | {'Found E':<12} | "
              f"{'Ratio':<8} | {'Hamming':<8} | {'Time':<6}")
        print("-" * 80)

        counts = {"EXACT": 0, "ENERGY_MATCH": 0, "FAIL": 0}
        energy_ratios = []
        hamming_dists = []
        solve_times = []

        for run in range(num_runs):
            target = ''.join(str(random.randint(0, 1)) for _ in range(n_bits))

            Q, clauses, info = create_qubo_quiet_planted(target, alpha=alpha)
            target_energy = info['target_energy']

            t0 = time.time()
            sampleset = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=sweeps)
            solve_time = time.time() - t0

            best_sample = sampleset.first.sample
            best_energy = sampleset.first.energy
            found_orig = extract_original_solution(best_sample, n_bits)

            result = classify_result(target, found_orig, target_energy, best_energy)
            counts[result] += 1

            ratio = best_energy / target_energy if abs(target_energy) > 1e-10 else float('nan')
            hdist = hamming_distance(target, found_orig)
            energy_ratios.append(ratio)
            hamming_dists.append(hdist)
            solve_times.append(solve_time)

            print(f"{run+1:<4} | {result:<12} | {target_energy:<12.2f} | {best_energy:<12.2f} | "
                  f"{ratio:<8.4f} | {hdist:<8} | {solve_time:<6.2f}")

        success = counts["EXACT"]
        success_rate = 100.0 * success / num_runs
        avg_ratio = np.nanmean(energy_ratios)
        avg_hamming = np.mean(hamming_dists)
        avg_time = np.mean(solve_times)

        print(f"\n  [요약] alpha={alpha}")
        print(f"    성공률: {success}/{num_runs} ({success_rate:.1f}%)")
        print(f"    EXACT={counts['EXACT']}, ENERGY={counts['ENERGY_MATCH']}, FAIL={counts['FAIL']}")
        print(f"    평균 에너지비: {avg_ratio:.4f}, 평균 해밍거리: {avg_hamming:.1f}, "
              f"평균 시간: {avg_time:.2f}s")

        results_summary.append({
            'alpha': alpha,
            'success_rate': success_rate,
            'counts': dict(counts),
            'avg_energy_ratio': avg_ratio,
            'avg_hamming': avg_hamming,
            'avg_time': avg_time,
        })

    # 전체 요약
    print("\n" + "=" * 95)
    print("전체 요약 (Alpha Sweep)")
    print("=" * 95)
    print(f"{'Alpha':<8} | {'Success%':<10} | {'EXACT':<6} | {'ENERGY':<6} | "
          f"{'FAIL':<6} | {'Avg E Ratio':<12} | {'Avg Hamming':<12}")
    print("-" * 80)
    for r in results_summary:
        c = r['counts']
        print(f"{r['alpha']:<8} | {r['success_rate']:<10.1f} | {c.get('EXACT',0):<6} | "
              f"{c.get('ENERGY_MATCH',0):<6} | {c.get('FAIL',0):<6} | "
              f"{r['avg_energy_ratio']:<12.4f} | {r['avg_hamming']:<12.1f}")

    return results_summary


def run_scaling_experiment(sizes=None, alpha=4.2, num_runs=10,
                           num_reads=200, num_sweeps=None):
    """
    스케일링 실험: 고정 alpha에서 N 증가에 따른 성공률 변화.
    """
    if sizes is None:
        sizes = [10, 20, 50, 100]

    print("=" * 90)
    print(f"Quiet Planting — Scaling 실험 (alpha={alpha})")
    print(f"Sizes={sizes}, num_runs={num_runs}, num_reads={num_reads}")
    print("=" * 90)

    sampler = neal.SimulatedAnnealingSampler()
    results = []

    for n in sizes:
        m = int(alpha * n)
        total_vars = n + m
        sweeps = num_sweeps or max(1000, 10 * total_vars)

        print(f"\n--- N={n} (QUBO 크기={total_vars}, sweeps={sweeps}) ---")
        success_count = 0
        energy_ratios = []
        hamming_dists = []

        for run in range(num_runs):
            target = ''.join(str(random.randint(0, 1)) for _ in range(n))
            Q, clauses, info = create_qubo_quiet_planted(target, alpha=alpha)
            target_energy = info['target_energy']

            sampleset = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=sweeps)
            best_sample = sampleset.first.sample
            best_energy = sampleset.first.energy
            found_orig = extract_original_solution(best_sample, n)

            result = classify_result(target, found_orig, target_energy, best_energy)
            if result == "EXACT":
                success_count += 1

            ratio = best_energy / target_energy if abs(target_energy) > 1e-10 else float('nan')
            energy_ratios.append(ratio)
            hamming_dists.append(hamming_distance(target, found_orig))

            print(f"  Run {run+1}: {result} (E_ratio={ratio:.4f}, Hamming={hamming_dists[-1]})")

        success_rate = 100.0 * success_count / num_runs
        avg_ratio = np.nanmean(energy_ratios)
        avg_hamming = np.mean(hamming_dists)

        print(f"  [N={n}] 성공률: {success_count}/{num_runs} ({success_rate:.1f}%), "
              f"평균 E비: {avg_ratio:.4f}, 평균 해밍: {avg_hamming:.1f}")

        results.append({
            'n': n,
            'success_rate': success_rate,
            'success_count': success_count,
            'avg_energy_ratio': avg_ratio,
            'avg_hamming': avg_hamming,
        })

    # 요약
    print("\n" + "=" * 70)
    print(f"Scaling 요약 (alpha={alpha})")
    print("=" * 70)
    print(f"{'N':<8} | {'QUBO Size':<10} | {'Success%':<10} | {'Count':<8} | "
          f"{'Avg E Ratio':<12} | {'Avg Hamming':<12}")
    print("-" * 70)
    for r in results:
        qubo_size = r['n'] + int(alpha * r['n'])
        print(f"{r['n']:<8} | {qubo_size:<10} | {r['success_rate']:<10.1f} | "
              f"{r['success_count']}/{num_runs:<5} | {r['avg_energy_ratio']:<12.4f} | "
              f"{r['avg_hamming']:<12.1f}")

    return results


def run_comparison_experiment(n_bits=50, num_runs=10, quiet_alpha=4.2,
                               wishart_alpha=0.7,
                               num_reads=200, num_sweeps=1000):
    """
    3-way 비교: Quiet Planting vs Wishart vs ZeroExp.
    """
    print("=" * 100)
    print(f"3-Way 비교 실험")
    print(f"N={n_bits}, num_runs={num_runs}")
    print(f"  Quiet: alpha={quiet_alpha}")
    print(f"  Wishart: alpha={wishart_alpha}")
    print(f"  ZeroExp: density=1.0")
    print("=" * 100)

    sampler = neal.SimulatedAnnealingSampler()

    quiet_m = int(quiet_alpha * n_bits)
    quiet_total = n_bits + quiet_m
    quiet_sweeps = max(1000, 10 * quiet_total)

    methods = ['Quiet', 'Wishart', 'ZeroExp']
    counts = {m: {"EXACT": 0, "ENERGY_MATCH": 0, "FAIL": 0} for m in methods}
    hammings = {k: [] for k in counts}

    print(f"\n{'Run':<4} | {'Quiet':<12} | {'Q_Ham':<6} | {'Wishart':<12} | {'W_Ham':<6} | "
          f"{'ZeroExp':<12} | {'Z_Ham':<6}")
    print("-" * 80)

    for run in range(num_runs):
        target = ''.join(str(random.randint(0, 1)) for _ in range(n_bits))

        # --- Quiet Planting ---
        Q_q, clauses_q, info_q = create_qubo_quiet_planted(target, alpha=quiet_alpha)
        te_q = info_q['target_energy']
        ss_q = sampler.sample_qubo(Q_q, num_reads=num_reads, num_sweeps=quiet_sweeps)
        found_q = extract_original_solution(ss_q.first.sample, n_bits)
        result_q = classify_result(target, found_q, te_q, ss_q.first.energy)
        hdist_q = hamming_distance(target, found_q)

        # --- Wishart ---
        Q_w = create_qubo_wishart(target, alpha=wishart_alpha)
        te_w = calculate_energy(target, Q_w)
        ss_w = sampler.sample_qubo(Q_w, num_reads=num_reads, num_sweeps=num_sweeps)
        found_w = ''.join(str(ss_w.first.sample[k]) for k in range(n_bits))
        result_w = classify_result(target, found_w, te_w, ss_w.first.energy)
        inverse_w = ''.join('0' if b == '1' else '1' for b in target)
        if found_w == inverse_w:
            result_w = "SYM_MATCH"
        hdist_w = hamming_distance(target, found_w)

        # --- Zero Expectation ---
        Q_z = create_qubo_precise(target, density=1.0)
        te_z = calculate_energy(target, Q_z)
        ss_z = sampler.sample_qubo(Q_z, num_reads=num_reads, num_sweeps=num_sweeps)
        found_z = ''.join(str(ss_z.first.sample[k]) for k in range(n_bits))
        result_z = classify_result(target, found_z, te_z, ss_z.first.energy)
        hdist_z = hamming_distance(target, found_z)

        # 집계
        for name, result, hdist in [
            ('Quiet', result_q, hdist_q),
            ('Wishart', result_w, hdist_w),
            ('ZeroExp', result_z, hdist_z),
        ]:
            r_key = "EXACT" if result in ("EXACT", "SYM_MATCH") else result
            counts[name][r_key] += 1
            hammings[name].append(hdist)

        print(f"{run+1:<4} | {result_q:<12} | {hdist_q:<6} | {result_w:<12} | {hdist_w:<6} | "
              f"{result_z:<12} | {hdist_z:<6}")

    # 요약
    print("\n" + "=" * 80)
    print("비교 요약")
    print("=" * 80)
    print(f"{'Method':<12} | {'Success%':<10} | {'EXACT':<6} | {'ENERGY':<6} | "
          f"{'FAIL':<6} | {'Avg Hamming':<12}")
    print("-" * 65)

    for name in methods:
        c = counts[name]
        success = c.get('EXACT', 0)
        success_rate = 100.0 * success / num_runs
        avg_h = np.mean(hammings[name])
        print(f"{name:<12} | {success_rate:<10.1f} | {success:<6} | "
              f"{c.get('ENERGY_MATCH',0):<6} | {c.get('FAIL',0):<6} | {avg_h:<12.1f}")


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    n_bits = 50
    num_runs = 10
    mode = "alpha"
    alpha_arg = 4.2

    if len(sys.argv) > 1:
        if sys.argv[1] == "--scaling":
            mode = "scaling"
            if len(sys.argv) > 2:
                alpha_arg = float(sys.argv[2])
            else:
                alpha_arg = 4.2
        elif sys.argv[1] == "--compare":
            mode = "compare"
        else:
            n_bits = int(sys.argv[1])

    if len(sys.argv) > 2 and mode == "alpha":
        num_runs = int(sys.argv[2])

    if mode == "alpha":
        run_alpha_sweep(
            n_bits=n_bits,
            alpha_values=[2.0, 3.0, 3.5, 3.86, 4.0, 4.2, 4.5, 5.0],
            num_runs=num_runs,
            num_reads=200,
        )

    elif mode == "scaling":
        run_scaling_experiment(
            sizes=[10, 20, 50, 100],
            alpha=alpha_arg,
            num_runs=10,
            num_reads=200,
        )

    elif mode == "compare":
        run_comparison_experiment(
            n_bits=n_bits,
            num_runs=num_runs,
        )
