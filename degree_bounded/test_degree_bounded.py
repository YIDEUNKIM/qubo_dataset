"""
Degree-Bounded Möbius QUBO SA 실험 프레임워크

실험 1: Ground State 검증 (brute force, n ≤ 20)
실험 2: d-sweep — degree별 SA 성공률
실험 3: n-scaling — n 증가에 따른 SA 성공률 및 생성 시간
실험 4: 비교 — 기존 방법론과 SA 성공률 비교
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from degree_bounded.qubo_degree_bounded import (
    create_qubo_degree_bounded, verify_ground_state, compute_aux_values
)
from qubo_utils import calculate_energy


def run_sa(Q, info, num_reads=100, num_sweeps=1000):
    """SA 실행 후 target 일치 여부 평가 (원래 변수만 비교)"""
    n = info['n']
    target = info['target']
    aux_info = info['aux_info']

    sampler = neal.SimulatedAnnealingSampler()
    ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)

    gs_found = 0
    total_hamming = 0

    for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
        found = ''.join(str(sample.get(i, 0)) for i in range(n))
        if found == target:
            gs_found += 1
        total_hamming += sum(1 for a, b in zip(target, found) if a != b)

    return gs_found, num_reads, total_hamming / num_reads


# ─────────────────────────────────────────────────
#  실험 1: Ground State 검증
# ─────────────────────────────────────────────────

def run_verification(num_runs=20):
    """n ≤ 15에서 brute force로 ground state 보장 검증"""
    print("=" * 70)
    print("실험 1: Ground State 검증 (brute force)")
    print("=" * 70)

    configs = [
        (8, 2), (8, 3), (10, 2), (10, 3), (10, 4),
        (12, 2), (12, 3), (15, 2), (15, 3),
    ]

    for n, d in configs:
        ok_count = 0
        unique_count = 0
        for run in range(num_runs):
            rng = random.Random(42 + run)
            target = ''.join(str(rng.randint(0, 1)) for _ in range(n))
            Q, info = create_qubo_degree_bounded(
                target, d=d, seed=run * 100, verbose=False)
            result = verify_ground_state(
                Q, target, n, info['aux_info'], info['offset'])
            if result['is_ground_state']:
                ok_count += 1
            if result.get('is_unique', False):
                unique_count += 1

        print(f"  n={n:<3} d={d} | GS: {ok_count}/{num_runs} | "
              f"Unique: {unique_count}/{num_runs} | "
              f"aux={info['n_aux']}")

    print()


# ─────────────────────────────────────────────────
#  실험 2: d-sweep
# ─────────────────────────────────────────────────

def run_d_sweep(n=30, num_runs=10, num_reads=100, num_sweeps=1000):
    """degree별 SA 성공률"""
    d_values = [2, 3, 4, 5]

    print("=" * 70)
    print(f"실험 2: d-sweep (n={n})")
    print(f"runs={num_runs}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 70)

    results = []

    for d in d_values:
        gs_total = 0
        samples_total = 0
        hamming_total = 0.0
        qubo_sizes = []
        gen_times = []

        for run in range(num_runs):
            rng = random.Random(42 + run)
            target = ''.join(str(rng.randint(0, 1)) for _ in range(n))

            t0 = time.time()
            Q, info = create_qubo_degree_bounded(
                target, d=d, density=0.5, seed=run * 100, verbose=False)
            gen_times.append(time.time() - t0)
            qubo_sizes.append(info['n_total'])

            gs, reads, avg_h = run_sa(
                Q, info, num_reads=num_reads, num_sweeps=num_sweeps)
            gs_total += gs
            samples_total += reads
            hamming_total += avg_h * reads

        rate = 100.0 * gs_total / samples_total
        avg_hamming = hamming_total / samples_total
        avg_qubo = np.mean(qubo_sizes)
        avg_gen = np.mean(gen_times)

        print(f"  d={d} | QUBO={avg_qubo:>6.0f} (aux={avg_qubo-n:>5.0f}) | "
              f"GS: {gs_total:>4}/{samples_total} ({rate:>6.2f}%) | "
              f"Hamming: {avg_hamming:>5.2f} | gen: {avg_gen:.3f}s")

        results.append({
            'd': d, 'gs_rate': rate, 'avg_qubo': avg_qubo,
            'avg_hamming': avg_hamming, 'avg_gen_time': avg_gen,
        })

    return results


# ─────────────────────────────────────────────────
#  실험 3: n-scaling
# ─────────────────────────────────────────────────

def run_scaling(d=2, n_values=None, num_runs=10, num_reads=100, num_sweeps=1000):
    """n 증가에 따른 SA 성공률 및 생성 시간"""
    if n_values is None:
        if d == 2:
            n_values = [20, 50, 100, 200, 500]
        else:
            n_values = [20, 30, 50, 75, 100]

    print("=" * 70)
    print(f"실험 3: n-scaling (d={d})")
    print(f"n_values={n_values}, runs={num_runs}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 70)

    results = []

    for n in n_values:
        gs_total = 0
        samples_total = 0
        hamming_total = 0.0
        qubo_sizes = []
        gen_times = []

        dens = min(1.0, max(0.05, 1000 / max(1, n**d // 6)))

        for run in range(num_runs):
            rng = random.Random(42 + run)
            target = ''.join(str(rng.randint(0, 1)) for _ in range(n))

            t0 = time.time()
            Q, info = create_qubo_degree_bounded(
                target, d=d, density=dens, seed=run * 100, verbose=False)
            gen_times.append(time.time() - t0)
            qubo_sizes.append(info['n_total'])

            gs, reads, avg_h = run_sa(
                Q, info, num_reads=num_reads, num_sweeps=num_sweeps)
            gs_total += gs
            samples_total += reads
            hamming_total += avg_h * reads

        rate = 100.0 * gs_total / samples_total
        avg_hamming = hamming_total / samples_total
        avg_qubo = np.mean(qubo_sizes)
        avg_gen = np.mean(gen_times)

        print(f"  n={n:<4} (dens={dens:.2f}) | QUBO={avg_qubo:>7.0f} | "
              f"GS: {gs_total:>4}/{samples_total} ({rate:>6.2f}%) | "
              f"Hamming: {avg_hamming:>5.2f} | gen: {avg_gen:.3f}s")

        results.append({
            'n': n, 'gs_rate': rate, 'avg_qubo': avg_qubo,
            'avg_hamming': avg_hamming, 'avg_gen_time': avg_gen,
            'density': dens,
        })

    return results


# ─────────────────────────────────────────────────
#  실험 4: 기존 방법론 비교
# ─────────────────────────────────────────────────

def run_comparison(n=50, num_runs=10, num_reads=100, num_sweeps=1000):
    """DegreeBounded vs Hardened vs Concat 비교"""
    print("=" * 70)
    print(f"실험 4: 방법론 비교 (N≈{n})")
    print(f"runs={num_runs}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 70)

    sampler = neal.SimulatedAnnealingSampler()
    results = []

    # --- Degree-Bounded d=2 ---
    gs_total, samples_total = 0, 0
    t0 = time.time()
    for run in range(num_runs):
        rng = random.Random(42 + run)
        target = ''.join(str(rng.randint(0, 1)) for _ in range(n))
        Q, info = create_qubo_degree_bounded(
            target, d=2, seed=run * 100, verbose=False)
        gs, reads, _ = run_sa(Q, info, num_reads, num_sweeps)
        gs_total += gs
        samples_total += reads
    rate = 100.0 * gs_total / samples_total
    elapsed = time.time() - t0
    print(f"  DegreeBounded(d=2)     | {rate:>6.2f}% | {gs_total:>4}/{samples_total} | {elapsed:.1f}s")
    results.append({'name': 'DegreeBounded(d=2)', 'gs_rate': rate,
                    'gs_count': gs_total, 'total': samples_total, 'time': elapsed})

    # --- Degree-Bounded d=3 ---
    gs_total, samples_total = 0, 0
    t0 = time.time()
    for run in range(num_runs):
        rng = random.Random(42 + run)
        target = ''.join(str(rng.randint(0, 1)) for _ in range(n))
        Q, info = create_qubo_degree_bounded(
            target, d=3, density=0.3, seed=run * 100, verbose=False)
        gs, reads, _ = run_sa(Q, info, num_reads, num_sweeps)
        gs_total += gs
        samples_total += reads
    rate = 100.0 * gs_total / samples_total
    elapsed = time.time() - t0
    print(f"  DegreeBounded(d=3)     | {rate:>6.2f}% | {gs_total:>4}/{samples_total} | {elapsed:.1f}s")
    results.append({'name': 'DegreeBounded(d=3)', 'gs_rate': rate,
                    'gs_count': gs_total, 'total': samples_total, 'time': elapsed})

    # --- Degree-Bounded d=2 + posiform ---
    gs_total, samples_total = 0, 0
    t0 = time.time()
    for run in range(num_runs):
        rng = random.Random(42 + run)
        target = ''.join(str(rng.randint(0, 1)) for _ in range(n))
        Q, info = create_qubo_degree_bounded(
            target, d=2, posiform_scale=0.1, seed=run * 100, verbose=False)
        gs, reads, _ = run_sa(Q, info, num_reads, num_sweeps)
        gs_total += gs
        samples_total += reads
    rate = 100.0 * gs_total / samples_total
    elapsed = time.time() - t0
    print(f"  DegreeBounded(d=2,α=0.1)| {rate:>6.2f}% | {gs_total:>4}/{samples_total} | {elapsed:.1f}s")
    results.append({'name': 'DegreeBounded(d=2,α=0.1)', 'gs_rate': rate,
                    'gs_count': gs_total, 'total': samples_total, 'time': elapsed})

    # --- Hardened Posiform ---
    try:
        from hardened_posiform.qubo_posiform_hardened import create_qubo_hardened_posiform

        for alpha in [0.01, 0.1]:
            gs_total, samples_total = 0, 0
            t0 = time.time()
            for run in range(num_runs):
                Q, hp_info = create_qubo_hardened_posiform(
                    n, coeff_type='lin2', posiform_scale=alpha, seed=run * 53)
                target_hp = hp_info['target']
                ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
                for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                    found = ''.join(str(sample[j]) for j in range(n))
                    if found == target_hp:
                        gs_total += 1
                    samples_total += 1
            rate = 100.0 * gs_total / samples_total
            elapsed = time.time() - t0
            print(f"  Hardened(α={alpha})       | {rate:>6.2f}% | {gs_total:>4}/{samples_total} | {elapsed:.1f}s")
            results.append({'name': f'Hardened(α={alpha})', 'gs_rate': rate,
                            'gs_count': gs_total, 'total': samples_total, 'time': elapsed})
    except Exception as e:
        print(f"  Hardened | skip: {e}")

    # 요약
    results.sort(key=lambda x: x['gs_rate'])
    print(f"\n{'='*60}")
    print(f"비교 요약 (N={n})")
    print(f"{'='*60}")
    max_rate = max(r['gs_rate'] for r in results) if results else 1
    for r in results:
        bar = '#' * int(r['gs_rate'] / max(max_rate, 1) * 30) if max_rate > 0 else ''
        print(f"  {r['name']:<28} | {r['gs_rate']:>6.2f}% | {bar}")

    return results


# ─────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────

def main():
    random.seed(42)
    np.random.seed(42)

    if len(sys.argv) < 2:
        print("사용법:")
        print("  python degree_bounded/test_degree_bounded.py --verify [num_runs]")
        print("  python degree_bounded/test_degree_bounded.py --d-sweep [num_runs]")
        print("  python degree_bounded/test_degree_bounded.py --scaling [num_runs]")
        print("  python degree_bounded/test_degree_bounded.py --compare [num_runs]")
        print()
        print("  --verify:   Ground State brute force 검증 (n≤15)")
        print("  --d-sweep:  degree별 SA 성공률 (n=30)")
        print("  --scaling:  n 증가에 따른 SA 성공률")
        print("  --compare:  기존 방법론 비교 (N=50)")
        return

    mode = sys.argv[1]
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 10

    if mode == '--verify':
        run_verification(num_runs=count)
    elif mode == '--d-sweep':
        run_d_sweep(num_runs=count, num_sweeps=1000)
    elif mode == '--scaling':
        run_scaling(d=2, num_runs=count, num_sweeps=1000)
        print()
        run_scaling(d=3, num_runs=count, num_sweeps=1000)
    elif mode == '--compare':
        run_comparison(n=50, num_runs=count, num_sweeps=1000)
    else:
        print(f"알 수 없는 모드: {mode}")


if __name__ == '__main__':
    main()
