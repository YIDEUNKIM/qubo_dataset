"""
Truth Table QUBO SA 실험 프레임워크

실험 1: 에너지 갭 프리셋 — gap 크기별 SA 성공률 (gap ↓ = 난이도 ↑)
실험 2: 다중 계곡 프리셋 — local minima 수 vs SA 성공률
실험 3: N-Scaling — n 증가에 따른 QUBO 크기 / SA 성능 변화
실험 4: 7-way 비교 — 기존 방법론과 SA 성공률 비교
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from truthtable.qubo_truthtable import (
    create_qubo_truthtable, preset_energy_gap, preset_multi_valley,
    compute_aux_values
)
from qubo_utils import calculate_energy


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def run_sa_on_truthtable(Q, info, num_reads=50, num_sweeps=1000):
    """Truth Table QUBO에 SA를 실행하고 결과를 반환"""
    n = info['n_original']
    total = info['n_total']
    target = info['target']
    aux_info = info['aux_info']

    sampler = neal.SimulatedAnnealingSampler()
    ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)

    gs_found = 0
    total_hamming = 0

    for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
        # 원래 n개 변수만 추출
        found_bits = ''.join(str(sample[k]) for k in range(n))
        if found_bits == target:
            gs_found += 1
        total_hamming += hamming_distance(target, found_bits)

    return gs_found, num_reads, total_hamming / num_reads


# ─────────────────────────────────────────────────
#  실험 1: 에너지 갭 프리셋 — gap sweep
# ─────────────────────────────────────────────────

def run_gap_sweep(n_bits=5, num_instances=20, num_reads=100, num_sweeps=5000):
    """
    에너지 갭 크기별 SA 성공률.
    gap이 작을수록 ground state 찾기 어려움 (양자 어닐링 minimum gap과 직결).
    """
    gap_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

    print("=" * 90)
    print(f"실험 1: Energy Gap Sweep")
    print(f"n={n_bits}, instances={num_instances}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 90)

    all_results = []

    for gap in gap_values:
        gs_total = 0
        samples_total = 0
        hamming_total = 0
        qubo_sizes = []

        t0 = time.time()
        for run in range(num_instances):
            target = ''.join([str(random.randint(0, 1)) for _ in range(n_bits)])
            tt = preset_energy_gap(n_bits, target, gap=gap, noise_scale=1.0, seed=run)
            Q, info = create_qubo_truthtable(tt, verbose=False)
            qubo_sizes.append(info['n_total'])

            gs, reads, avg_h = run_sa_on_truthtable(
                Q, info, num_reads=num_reads, num_sweeps=num_sweeps)
            gs_total += gs
            samples_total += reads
            hamming_total += avg_h

        elapsed = time.time() - t0
        rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
        avg_hamming = hamming_total / num_instances
        avg_size = np.mean(qubo_sizes)

        print(f"  gap={gap:<5} | GS rate: {gs_total:>4}/{samples_total} "
              f"({rate:>6.2f}%) | avg Hamming: {avg_hamming:>5.2f} | "
              f"QUBO size: {avg_size:.0f} | time: {elapsed:.1f}s")

        all_results.append({
            'gap': gap, 'gs_rate': rate, 'avg_hamming': avg_hamming,
            'avg_qubo_size': avg_size,
        })

    # 요약
    print(f"\n{'='*60}")
    print("Gap Sweep 요약")
    print(f"{'='*60}")
    print(f"{'Gap':>8} | {'GS Rate':>10} | {'Avg Hamming':>12}")
    print("-" * 36)
    for r in all_results:
        print(f"{r['gap']:>8.1f} | {r['gs_rate']:>9.2f}% | {r['avg_hamming']:>12.2f}")

    return all_results


# ─────────────────────────────────────────────────
#  실험 2: 다중 계곡 프리셋 — local minima 수 sweep
# ─────────────────────────────────────────────────

def run_valley_sweep(n_bits=5, num_instances=20, num_reads=100, num_sweeps=5000):
    """
    Local minima 수 증가 → SA가 local minimum에 갇힐 확률 증가.
    Quantum tunneling 능력 벤치마크.
    """
    valley_counts = [1, 2, 3, 4]

    print("=" * 90)
    print(f"실험 2: Multi-Valley Sweep (local minima 수)")
    print(f"n={n_bits}, instances={num_instances}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 90)

    all_results = []

    for n_valleys in valley_counts:
        gs_total = 0
        samples_total = 0
        hamming_total = 0

        t0 = time.time()
        for run in range(num_instances):
            rng = np.random.default_rng(run * 100)
            # 랜덤 targets 생성 (서로 Hamming 거리가 먼)
            targets = []
            for _ in range(n_valleys):
                for attempt in range(100):
                    t = ''.join([str(rng.integers(0, 2)) for _ in range(n_bits)])
                    # 기존 targets와 충분히 먼지 확인
                    if all(hamming_distance(t, prev) >= n_bits // 3 for prev in targets):
                        targets.append(t)
                        break
                else:
                    targets.append(t)

            tt = preset_multi_valley(n_bits, targets, gap=0.5,
                                     barrier_height=5.0, seed=run)
            Q, info = create_qubo_truthtable(tt, verbose=False)

            gs, reads, avg_h = run_sa_on_truthtable(
                Q, info, num_reads=num_reads, num_sweeps=num_sweeps)
            gs_total += gs
            samples_total += reads
            hamming_total += avg_h

        elapsed = time.time() - t0
        rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
        avg_hamming = hamming_total / num_instances

        print(f"  valleys={n_valleys} | GS rate: {gs_total:>4}/{samples_total} "
              f"({rate:>6.2f}%) | avg Hamming: {avg_hamming:>5.2f} | time: {elapsed:.1f}s")

        all_results.append({
            'n_valleys': n_valleys, 'gs_rate': rate, 'avg_hamming': avg_hamming,
        })

    return all_results


# ─────────────────────────────────────────────────
#  실험 3: N-Scaling
# ─────────────────────────────────────────────────

def run_scaling(sizes=None, num_runs=10, num_reads=100, num_sweeps=5000):
    """
    n 증가에 따른 QUBO 크기 폭발 + SA 성능 변화.
    """
    if sizes is None:
        sizes = [3, 4, 5, 6, 7]

    print("=" * 90)
    print(f"실험 3: N-Scaling (gap=2.0)")
    print(f"Sizes={sizes}, runs={num_runs}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 90)

    all_results = []

    for n in sizes:
        gs_total = 0
        samples_total = 0
        qubo_sizes = []
        aux_counts = []

        t0 = time.time()
        for run in range(num_runs):
            target = ''.join([str(random.randint(0, 1)) for _ in range(n)])
            tt = preset_energy_gap(n, target, gap=2.0, noise_scale=1.0, seed=run)
            Q, info = create_qubo_truthtable(tt, verbose=False)
            qubo_sizes.append(info['n_total'])
            aux_counts.append(info['n_aux'])

            gs, reads, _ = run_sa_on_truthtable(
                Q, info, num_reads=num_reads, num_sweeps=num_sweeps)
            gs_total += gs
            samples_total += reads

        elapsed = time.time() - t0
        rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
        avg_size = np.mean(qubo_sizes)
        avg_aux = np.mean(aux_counts)

        print(f"  n={n:<3} | QUBO={avg_size:>7.0f} (aux={avg_aux:>6.0f}) | "
              f"GS rate: {gs_total:>4}/{samples_total} ({rate:>6.2f}%) | time: {elapsed:.1f}s")

        all_results.append({
            'n': n, 'qubo_size': avg_size, 'aux_count': avg_aux,
            'gs_rate': rate,
        })

    return all_results


# ─────────────────────────────────────────────────
#  실험 4: 7-way 비교
# ─────────────────────────────────────────────────

def run_comparison(n_bits=5, num_runs=10, num_reads=100, num_sweeps=5000):
    """
    Truth Table (gap, valley) vs 기존 6개 방법론 SA 성공률 비교.
    n이 작아야 truth table 방법이 동작하므로 n=8 기준.
    """
    print("=" * 90)
    print(f"실험 4: 7-way 비교 (n={n_bits})")
    print(f"runs={num_runs}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 90)

    sampler = neal.SimulatedAnnealingSampler()
    results = {}

    # --- Truth Table: Energy Gap ---
    gs_total, samples_total = 0, 0
    t0 = time.time()
    for run in range(num_runs):
        target = ''.join([str(random.randint(0, 1)) for _ in range(n_bits)])
        tt = preset_energy_gap(n_bits, target, gap=1.0, seed=run)
        Q, info = create_qubo_truthtable(tt, verbose=False)
        gs, reads, _ = run_sa_on_truthtable(Q, info, num_reads, num_sweeps)
        gs_total += gs
        samples_total += reads
    rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
    print(f"  TruthTable-Gap   | {rate:>6.2f}% ({time.time()-t0:.1f}s)")
    results['TT-Gap'] = rate

    # --- Truth Table: Multi-Valley ---
    gs_total, samples_total = 0, 0
    t0 = time.time()
    for run in range(num_runs):
        rng = np.random.default_rng(run)
        t1 = ''.join([str(rng.integers(0, 2)) for _ in range(n_bits)])
        t2 = ''.join([str(1 - int(c)) for c in t1])  # 반대 bitstring
        tt = preset_multi_valley(n_bits, [t1, t2], gap=0.5, seed=run)
        Q, info = create_qubo_truthtable(tt, verbose=False)
        gs, reads, _ = run_sa_on_truthtable(Q, info, num_reads, num_sweeps)
        gs_total += gs
        samples_total += reads
    rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
    print(f"  TruthTable-Valley| {rate:>6.2f}% ({time.time()-t0:.1f}s)")
    results['TT-Valley'] = rate

    # --- 기존 방법론들 (n=8에서 실행 가능한 것만) ---
    try:
        from zero_expectation.qubo_zero_expectation import create_qubo_precise
        gs_total, samples_total = 0, 0
        t0 = time.time()
        for run in range(num_runs):
            target = ''.join([str(random.randint(0, 1)) for _ in range(n_bits)])
            Q = create_qubo_precise(target)
            ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
            best = ss.first
            found = ''.join(str(best.sample[k]) for k in range(n_bits))
            if found == target:
                gs_total += 1
            samples_total += 1
        rate = 100.0 * gs_total / samples_total
        print(f"  ZeroExpectation  | {rate:>6.2f}% ({time.time()-t0:.1f}s)")
        results['ZeroExp'] = rate
    except Exception as e:
        print(f"  ZeroExpectation  | 실패: {e}")

    try:
        from wishart.qubo_wishart import create_qubo_wishart
        gs_total, samples_total = 0, 0
        t0 = time.time()
        for run in range(num_runs):
            target = ''.join([str(random.randint(0, 1)) for _ in range(n_bits)])
            Q = create_qubo_wishart(target, alpha=0.7, seed=run)
            ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
            best = ss.first
            found = ''.join(str(best.sample[k]) for k in range(n_bits))
            if found == target:
                gs_total += 1
            samples_total += 1
        rate = 100.0 * gs_total / samples_total
        print(f"  Wishart(α=0.7)   | {rate:>6.2f}% ({time.time()-t0:.1f}s)")
        results['Wishart'] = rate
    except Exception as e:
        print(f"  Wishart          | 실패: {e}")

    try:
        from posiform.qubo_posiform import create_qubo_posiform
        gs_total, samples_total = 0, 0
        t0 = time.time()
        for run in range(num_runs):
            target = ''.join([str(random.randint(0, 1)) for _ in range(n_bits)])
            Q, info = create_qubo_posiform(target, seed=run)
            ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
            best = ss.first
            found = ''.join(str(best.sample[k]) for k in range(n_bits))
            if found == target:
                gs_total += 1
            samples_total += 1
        rate = 100.0 * gs_total / samples_total
        print(f"  Posiform         | {rate:>6.2f}% ({time.time()-t0:.1f}s)")
        results['Posiform'] = rate
    except Exception as e:
        print(f"  Posiform         | 실패: {e}")

    # 요약
    print(f"\n{'='*50}")
    print("7-way 비교 요약")
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
        print("  python test_truthtable.py --gap [num_instances]")
        print("  python test_truthtable.py --valley [num_instances]")
        print("  python test_truthtable.py --scaling [num_runs]")
        print("  python test_truthtable.py --compare [num_runs]")
        return

    mode = sys.argv[1]
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 10

    if mode == '--gap':
        run_gap_sweep(n_bits=5, num_instances=count)
    elif mode == '--valley':
        run_valley_sweep(n_bits=5, num_instances=count)
    elif mode == '--scaling':
        run_scaling(num_runs=count)
    elif mode == '--compare':
        run_comparison(n_bits=5, num_runs=count)
    else:
        print(f"알 수 없는 모드: {mode}")


if __name__ == '__main__':
    main()
