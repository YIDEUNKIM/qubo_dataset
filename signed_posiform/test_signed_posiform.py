"""
Signed Posiform SA 실험 프레임워크

실험 1: Scaling — N 증가에 따른 SA 성공률 변화
실험 2: negative_ratio sweep — 음수 비율에 따른 난이도 변화
실험 3: 비교 — Signed Posiform vs Plain Posiform

Note: Signed Posiform은 gap 계산에 O(2^n) 필요 → 생성 가능 범위 n ≤ 25
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from signed_posiform.qubo_signed_posiform import create_qubo_signed_posiform
from posiform.qubo_posiform import create_qubo_posiform
from qubo_utils import calculate_energy


def classify_result(target, found_solution, target_energy, found_energy):
    """결과 분류: EXACT / ENERGY_MATCH / FAIL"""
    if found_solution == target:
        return "EXACT"
    if abs(found_energy - target_energy) < 1e-4:
        return "ENERGY_MATCH"
    return "FAIL"


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def run_scaling_experiment(sizes=None, num_runs=10,
                           coeff_range=(1.0, 3.0), negative_ratio=0.3,
                           num_reads=200, num_sweeps=None):
    """
    스케일링 실험: N 증가에 따른 SA 성공률 변화.
    Signed Posiform은 생성에 O(2^n) 필요 → 실용 범위 n ≤ 25.
    """
    if sizes is None:
        sizes = [10, 15, 20]

    print("=" * 90)
    print(f"Signed Posiform — Scaling 실험")
    print(f"Sizes={sizes}, num_runs={num_runs}, num_reads={num_reads}")
    print(f"계수 범위: {coeff_range}, 음수 비율: {negative_ratio:.0%}")
    print("=" * 90)

    sampler = neal.SimulatedAnnealingSampler()
    results = []

    for n in sizes:
        sweeps = num_sweeps or max(1000, 10 * n)

        print(f"\n--- N={n} (sweeps={sweeps}) ---")
        success_count = 0
        energy_ratios = []
        hamming_dists = []
        gen_times = []
        neg_ratios_actual = []

        for run in range(num_runs):
            target = ''.join(str(random.randint(0, 1)) for _ in range(n))

            t0 = time.time()
            Q, info = create_qubo_signed_posiform(
                target, coeff_range=coeff_range, negative_ratio=negative_ratio,
                seed=run, verbose=False)
            gen_time = time.time() - t0
            gen_times.append(gen_time)

            if not info['is_unique']:
                print(f"  Run {run+1}: SKIP (유일성 실패)")
                continue

            target_energy = info['target_energy'] + info['constant_offset']
            q_total = info['q_pos_count'] + info['q_neg_count']
            neg_pct = info['q_neg_count'] / q_total * 100 if q_total > 0 else 0
            neg_ratios_actual.append(neg_pct)

            sampleset = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=sweeps)
            best_sample = sampleset.first.sample
            best_energy = sampleset.first.energy
            found = ''.join(str(best_sample[k]) for k in range(n))

            result = classify_result(target, found, target_energy, best_energy)
            if result == "EXACT":
                success_count += 1

            ratio = best_energy / target_energy if abs(target_energy) > 1e-10 else float('nan')
            energy_ratios.append(ratio)
            hamming_dists.append(hamming_distance(target, found))

            print(f"  Run {run+1}: {result} (Hamming={hamming_dists[-1]}, "
                  f"Q_neg={neg_pct:.0f}%, gen={gen_time:.2f}s)")

        success_rate = 100.0 * success_count / num_runs
        avg_ratio = np.nanmean(energy_ratios) if energy_ratios else float('nan')
        avg_hamming = np.mean(hamming_dists) if hamming_dists else float('nan')
        avg_gen = np.mean(gen_times) if gen_times else float('nan')
        avg_neg = np.mean(neg_ratios_actual) if neg_ratios_actual else 0

        print(f"  [N={n}] 성공률: {success_count}/{num_runs} ({success_rate:.1f}%), "
              f"평균 해밍: {avg_hamming:.1f}, Q_neg: {avg_neg:.1f}%, "
              f"생성시간: {avg_gen:.2f}s")

        results.append({
            'n': n,
            'success_rate': success_rate,
            'success_count': success_count,
            'avg_energy_ratio': avg_ratio,
            'avg_hamming': avg_hamming,
            'avg_gen_time': avg_gen,
            'avg_neg_pct': avg_neg,
        })

    # 요약
    print("\n" + "=" * 90)
    print(f"Scaling 요약 (Signed Posiform, neg_ratio={negative_ratio:.0%})")
    print("=" * 90)
    print(f"{'N':<6} | {'Success%':<10} | {'Count':<8} | "
          f"{'Avg Hamming':<12} | {'Q Neg%':<8} | {'Gen Time':<10}")
    print("-" * 70)
    for r in results:
        print(f"{r['n']:<6} | {r['success_rate']:<10.1f} | "
              f"{r['success_count']}/{num_runs:<5} | {r['avg_hamming']:<12.1f} | "
              f"{r['avg_neg_pct']:<8.1f} | {r['avg_gen_time']:<10.2f}")

    return results


def run_neg_ratio_sweep(n_bits=15, neg_ratios=None, num_runs=10,
                        coeff_range=(1.0, 3.0),
                        num_reads=200, num_sweeps=None):
    """
    음수 비율 sweep 실험: negative_ratio에 따른 난이도 변화.
    """
    if neg_ratios is None:
        neg_ratios = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7]

    sweeps = num_sweeps or max(1000, 10 * n_bits)

    print("=" * 95)
    print(f"Signed Posiform — 음수 비율 Sweep 실험")
    print(f"N={n_bits}, num_runs={num_runs}, num_reads={num_reads}, sweeps={sweeps}")
    print("=" * 95)

    sampler = neal.SimulatedAnnealingSampler()
    results_summary = []

    for neg_ratio in neg_ratios:
        print(f"\n--- negative_ratio={neg_ratio:.1f} ---")

        counts = {"EXACT": 0, "ENERGY_MATCH": 0, "FAIL": 0}
        hamming_dists = []
        q_neg_pcts = []
        gaps = []

        for run in range(num_runs):
            target = ''.join(str(random.randint(0, 1)) for _ in range(n_bits))

            Q, info = create_qubo_signed_posiform(
                target, coeff_range=coeff_range, negative_ratio=neg_ratio,
                seed=run, verbose=False)

            target_energy = info['target_energy'] + info['constant_offset']
            gaps.append(info['energy_gap'])

            q_total = info['q_pos_count'] + info['q_neg_count']
            q_neg_pcts.append(info['q_neg_count'] / q_total * 100 if q_total > 0 else 0)

            sampleset = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=sweeps)
            best_sample = sampleset.first.sample
            best_energy = sampleset.first.energy
            found = ''.join(str(best_sample[k]) for k in range(n_bits))

            result = classify_result(target, found, target_energy, best_energy)
            counts[result] += 1
            hamming_dists.append(hamming_distance(target, found))

            print(f"  Run {run+1}: {result} (Hamming={hamming_dists[-1]}, "
                  f"Q_neg={q_neg_pcts[-1]:.0f}%, gap={gaps[-1]:.4f})")

        success = counts["EXACT"]
        success_rate = 100.0 * success / num_runs
        avg_hamming = np.mean(hamming_dists)
        avg_neg = np.mean(q_neg_pcts)
        avg_gap = np.mean(gaps)

        print(f"  [요약] neg_ratio={neg_ratio:.1f}: "
              f"성공률={success_rate:.1f}%, Q_neg={avg_neg:.1f}%, gap={avg_gap:.4f}")

        results_summary.append({
            'neg_ratio': neg_ratio,
            'success_rate': success_rate,
            'counts': dict(counts),
            'avg_hamming': avg_hamming,
            'avg_q_neg_pct': avg_neg,
            'avg_gap': avg_gap,
        })

    # 전체 요약
    print("\n" + "=" * 95)
    print("전체 요약 (음수 비율 Sweep)")
    print("=" * 95)
    print(f"{'Neg Ratio':<10} | {'Success%':<10} | {'EXACT':<6} | {'ENERGY':<6} | "
          f"{'FAIL':<6} | {'Q Neg%':<8} | {'Avg Gap':<10} | {'Avg Hamming':<12}")
    print("-" * 85)
    for r in results_summary:
        c = r['counts']
        print(f"{r['neg_ratio']:<10.1f} | {r['success_rate']:<10.1f} | "
              f"{c.get('EXACT',0):<6} | {c.get('ENERGY_MATCH',0):<6} | "
              f"{c.get('FAIL',0):<6} | {r['avg_q_neg_pct']:<8.1f} | "
              f"{r['avg_gap']:<10.4f} | {r['avg_hamming']:<12.1f}")

    return results_summary


def run_comparison_experiment(n_bits=20, num_runs=10,
                               coeff_range=(1.0, 3.0), negative_ratio=0.3,
                               num_reads=200, num_sweeps=None):
    """
    비교 실험: Signed Posiform vs Plain Posiform.
    동일 target, 동일 SA 설정으로 비교.
    """
    sweeps = num_sweeps or max(1000, 10 * n_bits)

    print("=" * 90)
    print(f"Signed Posiform vs Plain Posiform 비교 실험")
    print(f"N={n_bits}, num_runs={num_runs}, num_reads={num_reads}, sweeps={sweeps}")
    print(f"Signed: coeff_range={coeff_range}, neg_ratio={negative_ratio}")
    print(f"Plain:  coeff_range={coeff_range}")
    print("=" * 90)

    sampler = neal.SimulatedAnnealingSampler()

    methods = ['SignedPosiform', 'Posiform']
    counts = {m: {"EXACT": 0, "ENERGY_MATCH": 0, "FAIL": 0} for m in methods}
    hammings = {m: [] for m in methods}

    print(f"\n{'Run':<4} | {'SignedPosiform':<14} {'S_Ham':<6} {'S_Neg%':<7} | "
          f"{'Posiform':<14} {'P_Ham':<6}")
    print("-" * 70)

    for run in range(num_runs):
        target = ''.join(str(random.randint(0, 1)) for _ in range(n_bits))

        # --- Signed Posiform ---
        Q_s, info_s = create_qubo_signed_posiform(
            target, coeff_range=coeff_range, negative_ratio=negative_ratio,
            seed=run, verbose=False)
        te_s = info_s['target_energy'] + info_s['constant_offset']
        q_total_s = info_s['q_pos_count'] + info_s['q_neg_count']
        neg_pct_s = info_s['q_neg_count'] / q_total_s * 100 if q_total_s > 0 else 0

        ss_s = sampler.sample_qubo(Q_s, num_reads=num_reads, num_sweeps=sweeps)
        found_s = ''.join(str(ss_s.first.sample[k]) for k in range(n_bits))
        result_s = classify_result(target, found_s, te_s, ss_s.first.energy)
        hdist_s = hamming_distance(target, found_s)

        # --- Plain Posiform ---
        Q_p, info_p = create_qubo_posiform(target, coeff_range=coeff_range, seed=run)
        te_p = info_p['target_energy']

        ss_p = sampler.sample_qubo(Q_p, num_reads=num_reads, num_sweeps=sweeps)
        found_p = ''.join(str(ss_p.first.sample[k]) for k in range(n_bits))
        result_p = classify_result(target, found_p, te_p, ss_p.first.energy)
        hdist_p = hamming_distance(target, found_p)

        # 집계
        for name, result, hdist in [
            ('SignedPosiform', result_s, hdist_s),
            ('Posiform', result_p, hdist_p),
        ]:
            counts[name][result] += 1
            hammings[name].append(hdist)

        print(f"{run+1:<4} | {result_s:<14} {hdist_s:<6} {neg_pct_s:<7.0f} | "
              f"{result_p:<14} {hdist_p:<6}")

    # 요약
    print("\n" + "=" * 80)
    print("비교 요약")
    print("=" * 80)
    print(f"{'Method':<16} | {'Success%':<10} | {'EXACT':<6} | {'ENERGY':<6} | "
          f"{'FAIL':<6} | {'Avg Hamming':<12}")
    print("-" * 70)

    for name in methods:
        c = counts[name]
        success = c.get('EXACT', 0)
        success_rate = 100.0 * success / num_runs
        avg_h = np.mean(hammings[name])
        print(f"{name:<16} | {success_rate:<10.1f} | {success:<6} | "
              f"{c.get('ENERGY_MATCH',0):<6} | {c.get('FAIL',0):<6} | {avg_h:<12.1f}")


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    num_runs = 10
    mode = "scaling"

    if len(sys.argv) > 1:
        if sys.argv[1] == "--scaling":
            mode = "scaling"
        elif sys.argv[1] == "--neg-sweep":
            mode = "neg_sweep"
        elif sys.argv[1] == "--compare":
            mode = "compare"

    if len(sys.argv) > 2:
        num_runs = int(sys.argv[2])

    if mode == "scaling":
        run_scaling_experiment(
            sizes=[10, 15, 20],
            num_runs=num_runs,
            num_reads=200,
        )

    elif mode == "neg_sweep":
        run_neg_ratio_sweep(
            n_bits=15,
            num_runs=num_runs,
            num_reads=200,
        )

    elif mode == "compare":
        run_comparison_experiment(
            n_bits=20,
            num_runs=num_runs,
            num_reads=200,
        )
