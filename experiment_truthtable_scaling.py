"""
Truth Table 정확 모드: n=3~16 생성 시간 + SA 시간 분석
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from truthtable.qubo_truthtable import (
    create_qubo_truthtable, preset_random_landscape, compute_aux_values
)
from qubo_utils import calculate_energy


def run_scaling_timing(sizes=None):
    if sizes is None:
        sizes = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    num_reads = 50
    num_sweeps = 1000
    sampler = neal.SimulatedAnnealingSampler()

    print("=" * 110)
    print(f"Truth Table 정확 모드: N-Scaling 시간 분석")
    print(f"reads={num_reads}, sweeps={num_sweeps}, 1 instance per n")
    print("=" * 110)
    print(f"{'n':>4} | {'보조변수':>10} | {'QUBO크기':>10} | {'생성시간':>10} | {'SA시간':>10} | {'총시간':>10} | {'SA성공률':>10}")
    print("-" * 110)

    for n in sizes:
        random.seed(42)
        np.random.seed(42)

        target = ''.join([str(random.randint(0, 1)) for _ in range(n)])
        tt = preset_random_landscape(n, target, seed=42)

        # 생성 시간 측정
        t_gen_start = time.time()
        try:
            Q, info = create_qubo_truthtable(tt, verbose=False)
        except MemoryError:
            print(f"{n:>4} | {'MemoryError':>10} |")
            break
        t_gen = time.time() - t_gen_start

        n_aux = info['n_aux']
        n_total = info['n_total']

        # SA 시간 측정
        t_sa_start = time.time()
        ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
        t_sa = time.time() - t_sa_start

        # 성공률 계산
        gs_found = 0
        for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
            found = ''.join(str(sample[k]) for k in range(n))
            if found == target:
                gs_found += 1

        rate = 100.0 * gs_found / num_reads
        t_total = t_gen + t_sa

        # 시간 포맷팅
        def fmt_time(t):
            if t < 1:
                return f"{t:.2f}s"
            elif t < 60:
                return f"{t:.1f}s"
            elif t < 3600:
                return f"{t/60:.1f}m"
            else:
                return f"{t/3600:.1f}h"

        print(f"{n:>4} | {n_aux:>10,} | {n_total:>10,} | {fmt_time(t_gen):>10} | "
              f"{fmt_time(t_sa):>10} | {fmt_time(t_total):>10} | {rate:>9.1f}%")

        # n=16 이상에서 생성만 10분 넘으면 중단
        if t_gen > 600:
            print("  (생성 시간 10분 초과 → 중단)")
            break


if __name__ == "__main__":
    run_scaling_timing()
