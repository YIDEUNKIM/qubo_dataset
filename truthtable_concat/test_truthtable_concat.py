"""
Truth Table Concat QUBO SA 실험 프레임워크

실험 1: GS 검증 — random landscape brute force ground state 검증
실험 2: h-Scaling — 블록 수 증가에 따른 SA 성공률 (k=7 고정, 6 configs)
실험 3: 14-way 비교 — Concat-Approx/Greedy/Random × α + Hardened × 3α
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from truthtable_concat.qubo_truthtable_concat import (
    create_qubo_concat, verify_ground_state
)
from qubo_utils import calculate_energy

# Lazy imports
try:
    from hardened_posiform.qubo_posiform_hardened import \
        create_qubo_hardened_posiform
except ImportError:
    create_qubo_hardened_posiform = None


def run_sa_on_concat(Q, info, num_reads=100, num_sweeps=1000):
    """Concat QUBO에 SA를 실행하고 원래 변수(k*h)만 추출하여 target과 비교"""
    target = info['target']

    sampler = neal.SimulatedAnnealingSampler()
    ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)

    gs_found = 0
    total_hamming = 0

    for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
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
#  실험 1: GS 검증 (random landscape)
# ─────────────────────────────────────────────────

def run_verification(num_runs=20):
    """Brute force로 random landscape GS 보장 검증 (n <= 21)"""
    print("=" * 70)
    print("실험 1: Ground State 검증 (random landscape, brute force)")
    print("=" * 70)

    configs = [
        (7, 7, 0.01),
        (7, 7, 0.1),
        (14, 7, 0.01),
        (14, 7, 0.1),
        (21, 7, 0.01),
        (21, 7, 0.1),
    ]

    for n, k, alpha in configs:
        h = n // k
        ok_count = 0
        unique_count = 0
        gap_sum = 0.0
        lm_sum = 0

        for run in range(num_runs):
            target_placeholder = '0' * k
            Q, info = create_qubo_concat(
                target_placeholder, h,
                landscape='random', posiform_scale=alpha,
                seed=run * 100, verbose=False)

            result = verify_ground_state(
                Q, info['target'], info['n_total'])
            if result['is_ground_state']:
                ok_count += 1
            if result['is_unique']:
                unique_count += 1
            gap_sum += result['energy_gap']
            if info.get('total_local_minima'):
                lm_sum += info['total_local_minima']

        avg_gap = gap_sum / num_runs
        avg_lm = lm_sum / num_runs
        print(f"  n={n:<3} k={k} alpha={alpha:<5} | "
              f"GS: {ok_count}/{num_runs} | "
              f"Unique: {unique_count}/{num_runs} | "
              f"gap: {avg_gap:.4f} | LM: {avg_lm:.1f}")

    print()


# ─────────────────────────────────────────────────
#  실험 2: h-Scaling (6 configs)
# ─────────────────────────────────────────────────

def run_scaling(k=7, h_values=None, num_runs=10, num_reads=100,
                num_sweeps=1000):
    """
    블록 수(h) 증가에 따른 SA 성공률.
    k=7 고정, 6개 config:
      Approx, Approx(α=0.01), Random(α=0.01), Random(α=0.1),
      Greedy, Greedy(α=0.01)
    """
    if h_values is None:
        h_values = [1, 2, 5, 10, 20]

    # (name, mode, alpha, landscape)
    configs = [
        ('Approx',         'approx', None, 'planted'),
        ('Approx(α=0.01)', 'approx', 0.01, 'planted'),
        ('Random(α=0.01)', 'approx', 0.01, 'random'),
        ('Random(α=0.1)',  'approx', 0.1,  'random'),
        ('Greedy',         'exact',  None,  'planted'),
        ('Greedy(α=0.01)', 'exact',  0.01,  'planted'),
    ]

    print("=" * 100)
    print(f"실험 2: h-Scaling (k={k}, 6 configs)")
    print(f"h_values={h_values}, runs={num_runs}, reads={num_reads}, "
          f"sweeps={num_sweeps}")
    print("=" * 100)

    all_results = {}

    for cfg_name, mode, alpha, landscape in configs:
        print(f"\n--- {cfg_name} (mode={mode}, α={alpha}, "
              f"landscape={landscape}) ---")
        cfg_results = []

        for h in h_values:
            gs_total = 0
            samples_total = 0
            hamming_total = 0
            qubo_sizes = []

            t0 = time.time()
            for run in range(num_runs):
                rng = random.Random(42 + run)
                target = ''.join(
                    [str(rng.randint(0, 1)) for _ in range(k)])

                kwargs = dict(
                    mode=mode,
                    posiform_scale=alpha,
                    landscape=landscape,
                    seed=run * 1000,
                    verbose=False,
                )
                if mode == 'exact' and landscape == 'planted':
                    kwargs['reduce_strategy'] = 'greedy'

                Q, info = create_qubo_concat(target, h, **kwargs)
                qubo_sizes.append(info['n_total'])

                gs, reads, avg_h = run_sa_on_concat(
                    Q, info, num_reads=num_reads,
                    num_sweeps=num_sweeps)
                gs_total += gs
                samples_total += reads
                hamming_total += avg_h

            elapsed = time.time() - t0
            rate = (100.0 * gs_total / samples_total
                    if samples_total > 0 else 0)
            avg_hamming = hamming_total / num_runs
            n_vars = k * h

            print(f"  h={h:<3} (N={n_vars:<4}) | "
                  f"QUBO={np.mean(qubo_sizes):>5.0f} | "
                  f"GS rate: {gs_total:>4}/{samples_total} "
                  f"({rate:>6.2f}%) | "
                  f"avg Hamming: {avg_hamming:>5.2f} | "
                  f"time: {elapsed:.1f}s")

            cfg_results.append({
                'h': h, 'n_vars': n_vars, 'gs_rate': rate,
                'avg_hamming': avg_hamming,
                'avg_qubo_size': np.mean(qubo_sizes),
            })

        all_results[cfg_name] = cfg_results

    # 요약
    print(f"\n{'='*100}")
    print(f"h-Scaling 요약 (k={k})")
    print(f"{'='*100}")
    header = f"{'h':>4} | {'N':>5}"
    for cfg_name in all_results:
        header += f" | {cfg_name:>18}"
    print(header)
    print("-" * len(header))

    for hi, h in enumerate(h_values):
        row = f"{h:>4} | {k*h:>5}"
        for cfg_name in all_results:
            r = all_results[cfg_name][hi]
            row += f" | {r['gs_rate']:>17.2f}%"
        print(row)

    return all_results


# ─────────────────────────────────────────────────
#  실험 3: 14-way 비교
# ─────────────────────────────────────────────────

def run_comparison(num_runs=10, num_reads=100, num_sweeps=1000):
    """
    14-way 비교 (N=35, k=7, h=5):
      Concat-Approx × 4α + Concat-Greedy × 4α + Random × 3α + Hardened × 3α
    """
    k = 7
    h = 5
    n_bits = k * h  # 35

    alphas = [None, 0.001, 0.01, 0.1]

    print("=" * 90)
    print(f"실험 3: 14-way 비교 (N={n_bits}, Concat: k={k},h={h})")
    print(f"runs={num_runs}, reads={num_reads}, sweeps={num_sweeps}")
    print(f"총 샘플/방법: {num_runs * num_reads}")
    print("=" * 90)

    sampler = neal.SimulatedAnnealingSampler()
    results = []

    # ─── Concat configs (Approx × 4 + Greedy × 4) ───
    concat_configs = []
    for alpha in alphas:
        alpha_str = f"(α={alpha})" if alpha is not None else ""
        concat_configs.append(
            ('approx', alpha, 'planted',
             f"Concat-Approx{alpha_str}"))
    for alpha in alphas:
        alpha_str = f"(α={alpha})" if alpha is not None else ""
        concat_configs.append(
            ('exact', alpha, 'planted',
             f"Concat-Greedy{alpha_str}"))

    for mode, alpha, landscape, name in concat_configs:
        gs_total = 0
        samples_total = 0
        hamming_total = 0.0
        qubo_sizes = []

        t0 = time.time()
        for run in range(num_runs):
            rng = random.Random(42 + run)
            target = ''.join(
                [str(rng.randint(0, 1)) for _ in range(k)])

            kwargs = dict(
                mode=mode,
                posiform_scale=alpha,
                landscape=landscape,
                seed=run * 1000,
                verbose=False,
            )
            if mode == 'exact':
                kwargs['reduce_strategy'] = 'greedy'

            Q, info = create_qubo_concat(target, h, **kwargs)
            qubo_sizes.append(info['n_total'])

            gs, reads, avg_h = run_sa_on_concat(
                Q, info, num_reads=num_reads,
                num_sweeps=num_sweeps)
            gs_total += gs
            samples_total += reads
            hamming_total += avg_h * num_reads

        elapsed = time.time() - t0
        rate = (100.0 * gs_total / samples_total
                if samples_total > 0 else 0)
        avg_hamming = (hamming_total / samples_total
                       if samples_total > 0 else 0)

        print(f"  {name:<28} | {rate:>6.2f}% | "
              f"{gs_total:>4}/{samples_total} | "
              f"QUBO={np.mean(qubo_sizes):>5.0f} | {elapsed:.1f}s")

        results.append({
            'name': name, 'gs_rate': rate, 'gs_count': gs_total,
            'total_samples': samples_total,
            'avg_hamming': avg_hamming,
            'avg_qubo_size': np.mean(qubo_sizes), 'time': elapsed,
        })

    # ─── Random landscape (3 alphas) ───
    random_alphas = [0.001, 0.01, 0.1]

    for alpha in random_alphas:
        name = f"Random(α={alpha})"
        gs_total = 0
        samples_total = 0
        hamming_total = 0.0
        qubo_sizes = []

        t0 = time.time()
        for run in range(num_runs):
            rng = random.Random(42 + run)
            target = ''.join(
                [str(rng.randint(0, 1)) for _ in range(k)])

            Q, info = create_qubo_concat(
                target, h,
                landscape='random', posiform_scale=alpha,
                seed=run * 1000, verbose=False)
            qubo_sizes.append(info['n_total'])

            gs, reads, avg_h = run_sa_on_concat(
                Q, info, num_reads=num_reads,
                num_sweeps=num_sweeps)
            gs_total += gs
            samples_total += reads
            hamming_total += avg_h * num_reads

        elapsed = time.time() - t0
        rate = (100.0 * gs_total / samples_total
                if samples_total > 0 else 0)
        avg_hamming = (hamming_total / samples_total
                       if samples_total > 0 else 0)

        print(f"  {name:<28} | {rate:>6.2f}% | "
              f"{gs_total:>4}/{samples_total} | "
              f"QUBO={np.mean(qubo_sizes):>5.0f} | {elapsed:.1f}s")

        results.append({
            'name': name, 'gs_rate': rate, 'gs_count': gs_total,
            'total_samples': samples_total,
            'avg_hamming': avg_hamming,
            'avg_qubo_size': np.mean(qubo_sizes), 'time': elapsed,
        })

    # ─── Hardened Posiform (3 alphas) ───
    hardened_alphas = [0.001, 0.01, 0.1]

    if create_qubo_hardened_posiform is not None:
        for alpha in hardened_alphas:
            name = f"Hardened(α={alpha})"
            gs_total = 0
            samples_total = 0
            hamming_total = 0.0
            qubo_sizes = []

            t0 = time.time()
            for run in range(num_runs):
                Q, hp_info = create_qubo_hardened_posiform(
                    n_bits, coeff_type='lin2',
                    posiform_scale=alpha, seed=run * 53)
                target_hp = hp_info['target']
                qubo_sizes.append(n_bits)

                ss = sampler.sample_qubo(
                    Q, num_reads=num_reads,
                    num_sweeps=num_sweeps)

                for sample, energy, _ in ss.data(
                        ['sample', 'energy', 'num_occurrences']):
                    found = ''.join(
                        str(sample[j]) for j in range(n_bits))
                    if found == target_hp:
                        gs_total += 1
                    hamming_total += sum(
                        1 for a, b in zip(target_hp, found)
                        if a != b)
                    samples_total += 1

            elapsed = time.time() - t0
            rate = (100.0 * gs_total / samples_total
                    if samples_total > 0 else 0)
            avg_hamming = (hamming_total / samples_total
                           if samples_total > 0 else 0)

            print(f"  {name:<28} | {rate:>6.2f}% | "
                  f"{gs_total:>4}/{samples_total} | "
                  f"QUBO={np.mean(qubo_sizes):>5.0f} | "
                  f"{elapsed:.1f}s")

            results.append({
                'name': name, 'gs_rate': rate,
                'gs_count': gs_total,
                'total_samples': samples_total,
                'avg_hamming': avg_hamming,
                'avg_qubo_size': np.mean(qubo_sizes),
                'time': elapsed,
            })
    else:
        print("  [Hardened Posiform] import 실패 — 건너뜀")

    # ─── 요약 테이블 (성공률 오름차순) ───
    results.sort(key=lambda x: x['gs_rate'])

    print(f"\n{'='*85}")
    print(f"14-way 비교 요약 (N={n_bits}, k={k}, h={h})")
    print(f"{'='*85}")
    print(f"{'Method':<28} | {'GS Rate':>8} | {'Count':>10} | "
          f"{'QUBO':>5} | {'Time':>6} |")
    print(f"{'-'*85}")

    max_rate = max(r['gs_rate'] for r in results) if results else 1
    for r in results:
        bar_len = (int(r['gs_rate'] / max(max_rate, 1) * 30)
                   if max_rate > 0 else 0)
        bar = '#' * bar_len
        print(f"{r['name']:<28} | {r['gs_rate']:>7.2f}% | "
              f"{r['gs_count']:>4}/{r['total_samples']:<4} | "
              f"{r['avg_qubo_size']:>5.0f} | {r['time']:>5.1f}s | "
              f"{bar}")

    return results


# ─────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────

def main():
    random.seed(42)
    np.random.seed(42)

    if len(sys.argv) < 2:
        print("사용법:")
        print("  python truthtable_concat/test_truthtable_concat.py "
              "--verify [num_runs]")
        print("  python truthtable_concat/test_truthtable_concat.py "
              "--scaling [num_runs]")
        print("  python truthtable_concat/test_truthtable_concat.py "
              "--compare [num_runs]")
        print()
        print("  --verify:  GS brute force 검증 "
              "(random landscape, n<=21)")
        print("  --scaling: h 증가에 따른 SA 성공률 "
              "(k=7, 6 configs)")
        print("  --compare: 14-way 비교 "
              "(N=35, Approx/Greedy/Random/Hardened)")
        return

    mode = sys.argv[1]
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 10

    if mode == '--verify':
        run_verification(num_runs=count)
    elif mode == '--scaling':
        run_scaling(k=7, num_runs=count, num_sweeps=1000)
    elif mode == '--compare':
        run_comparison(num_runs=count, num_sweeps=1000)
    else:
        print(f"알 수 없는 모드: {mode}")


if __name__ == '__main__':
    main()
