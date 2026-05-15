"""
Truth Table QUBO SA 실험 프레임워크

실험 1: 다중 계곡 프리셋 — local minima 수 vs SA 성공률
실험 2: N-Scaling — n 증가에 따른 QUBO 크기 / SA 성능 변화
실험 3: 7-way 비교 — 기존 방법론과 SA 성공률 비교
실험 4: 차수축소 전략 비교 — original/cache/greedy
실험 5: Greedy 확장 스케일링 — exact greedy vs approx
실험 6: Sweep 전이 — S-curve
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from etc.truthtable.qubo_truthtable import (
    create_qubo_truthtable, create_qubo_approx,
    create_qubo_approx_optimized,
    preset_random_landscape, preset_multi_valley
)


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
#  실험 1: 다중 계곡 프리셋 — local minima 수 sweep
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
    print(f"실험 2: N-Scaling (random landscape)")
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
            tt = preset_random_landscape(n, target, seed=run)
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
    Truth Table (random, valley) vs 기존 방법론 SA 성공률 비교.
    n이 작아야 truth table 방법이 동작하므로 n=5 기준.
    """
    print("=" * 90)
    print(f"실험 3: 7-way 비교 (n={n_bits})")
    print(f"runs={num_runs}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 90)

    sampler = neal.SimulatedAnnealingSampler()
    results = {}

    # --- Approx: Random Landscape ---
    gs_total, samples_total = 0, 0
    t0 = time.time()
    for run in range(num_runs):
        target = ''.join([str(random.randint(0, 1)) for _ in range(n_bits)])
        tt = preset_random_landscape(n_bits, target, seed=run)
        Q, info = create_qubo_approx(tt, verbose=False)
        gs, reads, _ = run_sa_on_truthtable(Q, info, num_reads, num_sweeps)
        gs_total += gs
        samples_total += reads
    rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
    print(f"  Approx-Random    | {rate:>6.2f}% ({time.time()-t0:.1f}s)")
    results['Approx-Random'] = rate

    # --- Approx: Multi-Valley ---
    gs_total, samples_total = 0, 0
    t0 = time.time()
    for run in range(num_runs):
        rng = np.random.default_rng(run)
        t1 = ''.join([str(rng.integers(0, 2)) for _ in range(n_bits)])
        t2 = ''.join([str(1 - int(c)) for c in t1])
        tt = preset_multi_valley(n_bits, [t1, t2], gap=0.5, seed=run)
        Q, info = create_qubo_approx(tt, verbose=False)
        gs, reads, _ = run_sa_on_truthtable(Q, info, num_reads, num_sweeps)
        gs_total += gs
        samples_total += reads
    rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
    print(f"  Approx-Valley    | {rate:>6.2f}% ({time.time()-t0:.1f}s)")
    results['Approx-Valley'] = rate

    # --- Exact (Rosenberg) — n≤5만 실용적 ---
    gs_total, samples_total = 0, 0
    t0 = time.time()
    for run in range(num_runs):
        target = ''.join([str(random.randint(0, 1)) for _ in range(n_bits)])
        tt = preset_random_landscape(n_bits, target, seed=run)
        Q, info = create_qubo_truthtable(tt, verbose=False)
        gs, reads, _ = run_sa_on_truthtable(Q, info, num_reads, num_sweeps)
        gs_total += gs
        samples_total += reads
    rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
    print(f"  Exact-Rosenberg  | {rate:>6.2f}% ({time.time()-t0:.1f}s)")
    results['Exact'] = rate

    # --- 기존 방법론들 (n=8에서 실행 가능한 것만) ---
    try:
        from etc.zero_expectation.qubo_zero_expectation import create_qubo_precise
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
        from etc.posiform.qubo_posiform import create_qubo_posiform
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
#  실험 4: 차수축소 전략 비교
# ─────────────────────────────────────────────────

def run_strategy_comparison(sizes=None, num_runs=5, num_reads=100, num_sweeps=5000):
    """
    3가지 Rosenberg 차수축소 전략 비교:
      original — 매번 새 보조변수
      cache    — 동일 쌍 재활용
      greedy   — 빈도 기반 쌍 선택 + 재활용

    각 전략별 보조변수 수, QUBO 크기, SA 성공률, 생성 시간 측정.
    """
    if sizes is None:
        sizes = [3, 4, 5, 6, 7, 8, 9]

    strategies = ['original', 'cache', 'greedy']

    print("=" * 100)
    print(f"실험 4: Rosenberg 차수축소 전략 비교")
    print(f"Sizes={sizes}, runs={num_runs}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 100)

    all_results = []

    for n in sizes:
        print(f"\n--- n = {n} ---")
        row = {'n': n}

        for strategy in strategies:
            aux_counts = []
            qubo_sizes = []
            gs_total = 0
            samples_total = 0
            gen_times = []

            for run in range(num_runs):
                target = ''.join([str(random.randint(0, 1)) for _ in range(n)])
                tt = preset_random_landscape(n, target, seed=run)

                t0 = time.time()
                Q, info = create_qubo_truthtable(
                    tt, verbose=False, reduce_strategy=strategy)
                gen_time = time.time() - t0
                gen_times.append(gen_time)

                aux_counts.append(info['n_aux'])
                qubo_sizes.append(info['n_total'])

                gs, reads, _ = run_sa_on_truthtable(
                    Q, info, num_reads=num_reads, num_sweeps=num_sweeps)
                gs_total += gs
                samples_total += reads

            rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
            avg_aux = np.mean(aux_counts)
            avg_size = np.mean(qubo_sizes)
            avg_gen = np.mean(gen_times)

            print(f"  {strategy:<8} | aux={avg_aux:>7.0f} | QUBO={avg_size:>7.0f} | "
                  f"GS rate: {rate:>6.2f}% | gen: {avg_gen:.3f}s")

            row[f'{strategy}_aux'] = avg_aux
            row[f'{strategy}_qubo'] = avg_size
            row[f'{strategy}_rate'] = rate
            row[f'{strategy}_gen'] = avg_gen

        all_results.append(row)

    # 요약 표
    print(f"\n{'='*100}")
    print("전략 비교 요약 — 보조변수 수")
    print(f"{'='*100}")
    print(f"{'n':>3} | {'original':>10} | {'cache':>10} | {'greedy':>10} | {'greedy 절감률':>14}")
    print("-" * 60)
    for r in all_results:
        orig = r.get('original_aux', 0)
        greedy = r.get('greedy_aux', 0)
        saving = 100.0 * (1 - greedy / orig) if orig > 0 else 0
        print(f"{r['n']:>3} | {orig:>10.0f} | {r.get('cache_aux', 0):>10.0f} | "
              f"{greedy:>10.0f} | {saving:>13.1f}%")

    print(f"\n전략 비교 요약 — SA 성공률")
    print(f"{'='*100}")
    print(f"{'n':>3} | {'original':>10} | {'cache':>10} | {'greedy':>10}")
    print("-" * 45)
    for r in all_results:
        print(f"{r['n']:>3} | {r.get('original_rate', 0):>9.2f}% | "
              f"{r.get('cache_rate', 0):>9.2f}% | {r.get('greedy_rate', 0):>9.2f}%")

    return all_results


# ─────────────────────────────────────────────────
#  실험 5: Greedy 확장 스케일링 (exact greedy vs approx)
# ─────────────────────────────────────────────────

def run_greedy_scaling(sizes=None, num_runs=100, num_reads=100, num_sweeps=5000):
    """
    Greedy 전략으로 정확 모드의 실용 한계 확장 측정.
    같은 n에서 approx 모드도 실행하여 1:1 비교.
    """
    if sizes is None:
        sizes = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

    print("=" * 100, flush=True)
    print(f"실험 5: Greedy 확장 스케일링 (exact greedy vs approx)", flush=True)
    print(f"Sizes={sizes}, runs={num_runs}, reads={num_reads}, sweeps={num_sweeps}", flush=True)
    print("=" * 100, flush=True)

    all_results = []

    for n in sizes:
        print(f"\n--- n = {n} ---", flush=True)
        row = {'n': n}

        for mode in ['greedy', 'approx']:
            qubo_sizes = []
            aux_counts = []
            gs_total = 0
            samples_total = 0
            gen_times = []

            for run in range(num_runs):
                target = ''.join([str(random.randint(0, 1)) for _ in range(n)])
                tt = preset_random_landscape(n, target, seed=run)

                t0 = time.time()
                if mode == 'greedy':
                    Q, info = create_qubo_truthtable(
                        tt, verbose=False, reduce_strategy='greedy')
                else:
                    Q, info = create_qubo_approx_optimized(tt, verbose=False)
                gen_time = time.time() - t0
                gen_times.append(gen_time)

                qubo_sizes.append(info['n_total'])
                aux_counts.append(info['n_aux'])

                gs, reads, _ = run_sa_on_truthtable(
                    Q, info, num_reads=num_reads, num_sweeps=num_sweeps)
                gs_total += gs
                samples_total += reads

            rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0
            avg_size = np.mean(qubo_sizes)
            avg_aux = np.mean(aux_counts)
            avg_gen = np.mean(gen_times)

            print(f"  {mode:<8} | QUBO={avg_size:>7.0f} (aux={avg_aux:>6.0f}) | "
                  f"GS rate: {rate:>6.2f}% | gen: {avg_gen:.3f}s", flush=True)

            row[f'{mode}_qubo'] = avg_size
            row[f'{mode}_aux'] = avg_aux
            row[f'{mode}_rate'] = rate
            row[f'{mode}_gen'] = avg_gen

        all_results.append(row)

    # 요약 표
    print(f"\n{'='*100}")
    print("Greedy 확장 스케일링 요약")
    print(f"{'='*100}")
    print(f"{'n':>3} | {'Greedy QUBO':>12} | {'Greedy Rate':>12} | "
          f"{'Approx QUBO':>12} | {'Approx Rate':>12} | {'Greedy Gen':>10} | {'Approx Gen':>10}")
    print("-" * 85)
    for r in all_results:
        print(f"{r['n']:>3} | {r.get('greedy_qubo', 0):>12.0f} | "
              f"{r.get('greedy_rate', 0):>11.2f}% | "
              f"{r.get('approx_qubo', 0):>12.0f} | "
              f"{r.get('approx_rate', 0):>11.2f}% | "
              f"{r.get('greedy_gen', 0):>9.3f}s | "
              f"{r.get('approx_gen', 0):>9.3f}s")

    return all_results


# ─────────────────────────────────────────────────
#  실험 6: Sweep 전이 (S-curve)
# ─────────────────────────────────────────────────

def run_sweep_transition(n=8, num_instances=100, num_reads=100):
    """
    Greedy로 축소된 QUBO의 SA hardness 특성 분석.
    sweep 수를 변화시키며 SA 성공률의 S-curve 전이 측정.
    """
    sweep_values = [50, 100, 200, 500, 1000, 2000, 5000, 10000]

    print("=" * 100)
    print(f"실험 6: Sweep 전이 (S-curve)")
    print(f"n={n}, instances={num_instances}, reads={num_reads}")
    print(f"sweeps={sweep_values}")
    print("=" * 100)

    # 인스턴스 사전 생성 (모든 sweep 값에 동일 인스턴스 사용)
    print(f"\n[1] QUBO 인스턴스 {num_instances}개 사전 생성 (greedy)...")
    instances = []
    qubo_sizes = []
    t0 = time.time()
    for run in range(num_instances):
        target = ''.join([str(random.randint(0, 1)) for _ in range(n)])
        tt = preset_random_landscape(n, target, seed=run)
        Q, info = create_qubo_truthtable(tt, verbose=False, reduce_strategy='greedy')
        instances.append((Q, info))
        qubo_sizes.append(info['n_total'])
    gen_time = time.time() - t0
    print(f"  생성 완료: avg QUBO size={np.mean(qubo_sizes):.0f}, time={gen_time:.1f}s")

    # sweep별 SA 실행
    print(f"\n[2] Sweep별 SA 실행...")
    all_results = []

    for sweeps in sweep_values:
        gs_total = 0
        samples_total = 0

        t0 = time.time()
        for Q, info in instances:
            gs, reads, _ = run_sa_on_truthtable(
                Q, info, num_reads=num_reads, num_sweeps=sweeps)
            gs_total += gs
            samples_total += reads

        elapsed = time.time() - t0
        rate = 100.0 * gs_total / samples_total if samples_total > 0 else 0

        print(f"  sweeps={sweeps:>6} | GS rate: {gs_total:>5}/{samples_total} "
              f"({rate:>6.2f}%) | time: {elapsed:.1f}s")

        all_results.append({
            'sweeps': sweeps, 'gs_rate': rate,
            'gs_found': gs_total, 'samples': samples_total,
        })

    # 요약
    print(f"\n{'='*60}")
    print(f"Sweep 전이 요약 (n={n}, greedy)")
    print(f"{'='*60}")
    print(f"{'Sweeps':>8} | {'GS Rate':>10}")
    print("-" * 22)
    for r in all_results:
        bar = '#' * int(r['gs_rate'] / 2)
        print(f"{r['sweeps']:>8} | {r['gs_rate']:>9.2f}% | {bar}")

    return all_results


# ─────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────

def main():
    random.seed(42)
    np.random.seed(42)

    if len(sys.argv) < 2:
        print("사용법:")
        print("  python test_truthtable.py --valley [num_instances]")
        print("  python test_truthtable.py --scaling [num_runs]")
        print("  python test_truthtable.py --compare [num_runs]")
        print("  python test_truthtable.py --strategy [num_runs] [sizes]")
        print("  python test_truthtable.py --greedy-scaling [num_runs] [sizes]")
        print("  python test_truthtable.py --sweep [num_instances]")
        print()
        print("  --strategy: 3가지 Rosenberg 전략 비교 (original/cache/greedy)")
        print("    예: --strategy 5          (n=3..9, 5 runs)")
        print("    예: --strategy 5 3,4,5,6  (n=3..6, 5 runs)")
        print("  --greedy-scaling: greedy 확장 스케일링 (exact greedy vs approx)")
        print("    예: --greedy-scaling 100              (n=3..12, 100 runs)")
        print("    예: --greedy-scaling 100 3,4,5,6,7,8  (n=3..8, 100 runs)")
        print("  --sweep: sweep 전이 S-curve (n=8, greedy)")
        print("    예: --sweep 100  (100 instances)")
        return

    mode = sys.argv[1]
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 10

    if mode == '--valley':
        run_valley_sweep(n_bits=5, num_instances=count)
    elif mode == '--scaling':
        run_scaling(num_runs=count)
    elif mode == '--compare':
        run_comparison(n_bits=5, num_runs=count)
    elif mode == '--strategy':
        sizes = None
        if len(sys.argv) > 3:
            sizes = [int(x) for x in sys.argv[3].split(',')]
        run_strategy_comparison(sizes=sizes, num_runs=count)
    elif mode == '--greedy-scaling':
        sizes = None
        if len(sys.argv) > 3:
            sizes = [int(x) for x in sys.argv[3].split(',')]
        run_greedy_scaling(sizes=sizes, num_runs=count)
    elif mode == '--sweep':
        run_sweep_transition(n=8, num_instances=count)
    else:
        print(f"알 수 없는 모드: {mode}")


if __name__ == '__main__':
    main()
