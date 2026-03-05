"""
Truth Table Concat QUBO SA 실험 프레임워크

실험 1: h-Scaling — 블록 수 증가에 따른 SA 성공률 (k=7 고정)
실험 2: 8-way 비교 — 동일 변수 수에서 기존 방법론과 SA 성공률 비교
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from truthtable_concat.qubo_truthtable_concat import create_qubo_concat
from qubo_utils import calculate_energy


def run_sa_on_concat(Q, info, num_reads=50, num_sweeps=1000):
    """Concat QUBO에 SA를 실행하고 원래 변수(k*h)만 추출하여 target과 비교"""
    n_original = info['n_original']
    target = info['target']

    sampler = neal.SimulatedAnnealingSampler()
    ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)

    gs_found = 0
    total_hamming = 0

    for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
        # 원래 k*h개 변수만 추출 (블록별 원래 변수)
        found_bits = extract_original_bits(sample, info)
        if found_bits == target:
            gs_found += 1
        total_hamming += sum(1 for a, b in zip(target, found_bits) if a != b)

    return gs_found, num_reads, total_hamming / num_reads


def extract_original_bits(sample, info):
    """SA 결과에서 원래 변수(k*h)만 추출"""
    k = info['n_block']
    blocks = info['blocks']
    bits = []
    for block in blocks:
        offset = block['offset']
        for j in range(k):
            bits.append(str(sample.get(offset + j, 0)))
    return ''.join(bits)


# ─────────────────────────────────────────────────
#  실험 1: h-Scaling
# ─────────────────────────────────────────────────

def run_scaling(k=7, h_values=None, num_runs=10, num_reads=100, num_sweeps=5000,
                mode='approx'):
    """
    블록 수(h) 증가에 따른 SA 성공률.
    k=7 고정, h=1,2,5,10,20으로 총 변수 수 7~140 확장.
    """
    if h_values is None:
        h_values = [1, 2, 5, 10, 20]

    print("=" * 90)
    print(f"실험 1: h-Scaling (k={k}, mode={mode})")
    print(f"h_values={h_values}, runs={num_runs}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 90)

    all_results = []

    for h in h_values:
        gs_total = 0
        samples_total = 0
        hamming_total = 0
        qubo_sizes = []

        t0 = time.time()
        for run in range(num_runs):
            target = ''.join([str(random.randint(0, 1)) for _ in range(k)])
            Q, info = create_qubo_concat(
                target, h, gap=2.0, noise_scale=1.0,
                mode=mode, seed=run * 1000, verbose=False)
            qubo_sizes.append(info['n_total'])

            gs, reads, avg_h = run_sa_on_concat(
                Q, info, num_reads=num_reads, num_sweeps=num_sweeps)
            gs_total += gs
            samples_total += reads
            hamming_total += avg_h

        elapsed = time.time() - t0
        rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
        avg_hamming = hamming_total / num_runs
        n_vars = k * h

        print(f"  h={h:<3} (N={n_vars:<4}) | QUBO={np.mean(qubo_sizes):>5.0f} | "
              f"GS rate: {gs_total:>4}/{samples_total} ({rate:>6.2f}%) | "
              f"avg Hamming: {avg_hamming:>5.2f} | time: {elapsed:.1f}s")

        all_results.append({
            'h': h, 'n_vars': n_vars, 'gs_rate': rate,
            'avg_hamming': avg_hamming, 'avg_qubo_size': np.mean(qubo_sizes),
        })

    # 요약
    print(f"\n{'='*60}")
    print(f"h-Scaling 요약 (k={k})")
    print(f"{'='*60}")
    print(f"{'h':>4} | {'N':>5} | {'GS Rate':>10} | {'Avg Hamming':>12}")
    print("-" * 40)
    for r in all_results:
        bar = '#' * int(r['gs_rate'] / 2)
        print(f"{r['h']:>4} | {r['n_vars']:>5} | {r['gs_rate']:>9.2f}% | "
              f"{r['avg_hamming']:>12.2f} | {bar}")

    return all_results


# ─────────────────────────────────────────────────
#  실험 2: 8-way 비교
# ─────────────────────────────────────────────────

def run_comparison(n_bits=35, num_runs=10, num_reads=100, num_sweeps=1000):
    """
    동일 변수 수(N≈35)에서 Concat vs Hardened 비교.
    """
    k = 7
    h = n_bits // k

    print("=" * 90)
    print(f"실험 2: Hardened vs Concat 비교 (N≈{n_bits}, Concat: k={k},h={h})")
    print(f"runs={num_runs}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 90)

    sampler = neal.SimulatedAnnealingSampler()
    results = {}

    # --- Concat Approx ---
    gs_total, samples_total = 0, 0
    t0 = time.time()
    for run in range(num_runs):
        target = ''.join([str(random.randint(0, 1)) for _ in range(k)])
        Q, info = create_qubo_concat(
            target, h, gap=2.0, mode='approx', seed=run * 1000, verbose=False)
        gs, reads, _ = run_sa_on_concat(Q, info, num_reads, num_sweeps)
        gs_total += gs
        samples_total += reads
    rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
    print(f"  Concat-Approx    | {rate:>6.2f}% ({time.time()-t0:.1f}s)")
    results['Concat-Approx'] = rate

    # --- Concat Exact ---
    gs_total, samples_total = 0, 0
    t0 = time.time()
    for run in range(num_runs):
        target = ''.join([str(random.randint(0, 1)) for _ in range(k)])
        Q, info = create_qubo_concat(
            target, h, gap=2.0, mode='exact', seed=run * 1000, verbose=False)
        gs, reads, _ = run_sa_on_concat(Q, info, num_reads, num_sweeps)
        gs_total += gs
        samples_total += reads
    rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
    print(f"  Concat-Exact     | {rate:>6.2f}% ({time.time()-t0:.1f}s)")
    results['Concat-Exact'] = rate

    # --- Hardened Posiform ---
    n = n_bits
    try:
        from hardened_posiform.qubo_posiform_hardened import create_qubo_hardened_posiform
        gs_total, samples_total = 0, 0
        t0 = time.time()
        for run in range(num_runs):
            Q, hp_info = create_qubo_hardened_posiform(
                n, coeff_type='lin2', posiform_scale=0.1, seed=run)
            target_hp = hp_info['target']
            ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
            best = ss.first
            found = ''.join(str(best.sample[k_]) for k_ in range(n))
            if found == target_hp:
                gs_total += 1
            samples_total += 1
        rate = 100.0 * gs_total / samples_total
        print(f"  Hardened(lin2)   | {rate:>6.2f}% ({time.time()-t0:.1f}s)")
        results['Hardened'] = rate
    except Exception as e:
        print(f"  Hardened         | skip: {e}")

    # 요약
    print(f"\n{'='*50}")
    print(f"Hardened vs Concat 비교 요약 (N≈{n_bits})")
    print(f"{'='*50}")
    for method, rate in sorted(results.items(), key=lambda x: x[1]):
        bar = '#' * int(rate / 2)
        print(f"  {method:<18} | {rate:>6.2f}% | {bar}")

    return results


# ─────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────

def main():
    random.seed(42)
    np.random.seed(42)

    if len(sys.argv) < 2:
        print("사용법:")
        print("  python truthtable_concat/test_truthtable_concat.py --scaling [num_runs]")
        print("  python truthtable_concat/test_truthtable_concat.py --compare [num_runs]")
        print()
        print("  --scaling: h 증가에 따른 SA 성공률 (k=7 고정, h=1,2,5,10,20)")
        print("  --compare: 8-way 비교 (N≈35, Concat vs 기존 방법론)")
        return

    mode = sys.argv[1]
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 10

    qubo_mode = 'exact' if '--exact' in sys.argv else 'approx'

    if mode == '--scaling':
        run_scaling(k=7, num_runs=count, num_sweeps=1000, mode=qubo_mode)
    elif mode == '--compare':
        run_comparison(n_bits=35, num_runs=count, num_sweeps=1000)
    else:
        print(f"알 수 없는 모드: {mode}")


if __name__ == '__main__':
    main()
