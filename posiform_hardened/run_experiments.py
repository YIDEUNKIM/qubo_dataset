"""
Hardened Posiform 실험 스크립트: sweep 전이 + N-scaling + 비교
"""
import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from posiform_hardened.qubo_posiform_hardened import create_qubo_hardened_posiform
from posiform.qubo_posiform import create_qubo_posiform
from qubo_utils import calculate_energy


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def run_sweep_transition():
    """실험 1: Sweep 전이 (논문 Fig 7-8 재현)"""
    N = 500
    num_inst = 8
    reads = 100
    sweep_vals = [1, 5, 10, 50, 100, 500, 1000, 5000]

    configs = [
        ('lin2', 0.01),
        ('lin2', 0.1),
        ('lin20', 0.01),
        ('lin20', 0.1),
    ]

    sampler = neal.SimulatedAnnealingSampler()

    print("=" * 110)
    print(f"실험 1: Sweep 전이 (N={N}, instances={num_inst}, reads/sweep={reads})")
    print(f"Sweep values: {sweep_vals}")
    print("=" * 110)

    summary = {}

    for ct, alpha in configs:
        label = f"{ct},a={alpha}"
        print(f"\n--- {label} ---")
        summary[label] = {}

        t0 = time.time()
        insts = []
        for i in range(num_inst):
            Q, info = create_qubo_hardened_posiform(
                N, max_subgraph_size=15, coeff_type=ct,
                posiform_scale=alpha, seed=i * 53
            )
            if info['posiform_is_unique']:
                insts.append((Q, info))
        print(f"  생성: {len(insts)}/{num_inst} ({time.time()-t0:.1f}s)")

        for sweeps in sweep_vals:
            gs = 0
            total = 0
            total_h = 0

            t0 = time.time()
            for Q, info in insts:
                target = info['target']
                ss = sampler.sample_qubo(Q, num_reads=reads, num_sweeps=sweeps)
                for s in ss.samples():
                    found = ''.join(str(s[k]) for k in range(N))
                    total += 1
                    if found == target:
                        gs += 1
                    total_h += hamming(target, found)
            dt = time.time() - t0

            rate = 100.0 * gs / total if total else 0
            avg_h = total_h / total if total else 0
            summary[label][sweeps] = rate

            print(f"  sweeps={sweeps:<6} | GS: {gs:>4}/{total} ({rate:>6.2f}%) | "
                  f"H={avg_h:>6.1f} | {dt:.1f}s")

    # 요약 테이블
    print("\n" + "=" * 110)
    print("Sweep 전이 요약 (Ground-State Sampling Rate %)")
    print("=" * 110)
    header = f"{'Config':<14} |"
    for s in sweep_vals:
        header += f" {s:>7}"
    print(header)
    print("-" * (16 + 8 * len(sweep_vals)))

    for label in summary:
        row = f"{label:<14} |"
        for s in sweep_vals:
            row += f" {summary[label].get(s, 0):>6.1f}%"
        print(row)


def run_scaling_and_compare():
    """실험 2+3: N-Scaling + Hardened vs Plain 비교"""
    sizes = [50, 100, 200, 300, 500]
    num_runs = 15
    reads = 100
    sweeps = 500
    sampler = neal.SimulatedAnnealingSampler()

    configs = [
        ('lin2', 0.01, 'Hard(lin2,0.01)'),
        ('lin2', 0.1, 'Hard(lin2,0.1)'),
        ('lin20', 0.01, 'Hard(lin20,0.01)'),
        ('lin20', 0.1, 'Hard(lin20,0.1)'),
        (None, None, 'Plain Posiform'),
    ]

    print("\n" + "=" * 100)
    print(f"실험 2+3: N-Scaling + Hardened vs Plain")
    print(f"Sizes={sizes}, runs={num_runs}, reads={reads}, sweeps={sweeps}")
    print("=" * 100)

    summary = {}

    for ct, alpha, label in configs:
        print(f"\n--- {label} ---")
        summary[label] = {}

        for n in sizes:
            gs = 0
            total = 0
            total_h = 0

            t0 = time.time()
            for run in range(num_runs):
                if ct is not None:
                    Q, info = create_qubo_hardened_posiform(
                        n, max_subgraph_size=15, coeff_type=ct,
                        posiform_scale=alpha, seed=run * 31 + n
                    )
                    target = info['target']
                    ok = info['posiform_is_unique']
                else:
                    rng = random.Random(run * 31 + n)
                    target = ''.join(str(rng.randint(0, 1)) for _ in range(n))
                    Q, pinfo = create_qubo_posiform(
                        target, coeff_range=(1.0, 3.0), seed=run * 31 + n
                    )
                    ok = pinfo['is_unique']

                if not ok:
                    continue

                ss = sampler.sample_qubo(Q, num_reads=reads, num_sweeps=sweeps)
                found = ''.join(str(ss.first.sample[k]) for k in range(n))
                total += 1
                if found == target:
                    gs += 1
                total_h += hamming(target, found)

            dt = time.time() - t0
            rate = 100.0 * gs / total if total else 0
            avg_h = total_h / total if total else 0
            summary[label][n] = rate

            print(f"  N={n:<4} | {gs:>2}/{total} ({rate:>5.1f}%) | "
                  f"H={avg_h:>5.1f} | {dt:.1f}s")

    # 요약 테이블
    print("\n" + "=" * 100)
    print("N-Scaling 요약 (SA Best-Sample Success %)")
    print("=" * 100)
    header = f"{'Method':<20} |"
    for n in sizes:
        header += f" N={n:<5}"
    print(header)
    print("-" * (22 + 8 * len(sizes)))

    for _, _, label in configs:
        row = f"{label:<20} |"
        for n in sizes:
            row += f" {summary[label].get(n, 0):>5.1f}%"
        print(row)


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    mode = sys.argv[1] if len(sys.argv) > 1 else "all"

    if mode in ("sweep", "all"):
        run_sweep_transition()

    if mode in ("scaling", "all"):
        run_scaling_and_compare()

    print("\n모든 실험 완료.")
