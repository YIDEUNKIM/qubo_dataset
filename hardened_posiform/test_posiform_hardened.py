"""
Hardened Posiform SA 실험 프레임워크

Pelofske et al. (2024) 논문 기반 실험:
  실험 1: Sweep 전이 — SA sweep 수 변화에 따른 성공률 S-curve (N=1000)
  실험 2: N-Scaling — 문제 크기 증가에 따른 SA 성공률 감소
  실험 3: Hardened vs Plain — hardened가 plain보다 어려운지 비교
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hardened_posiform.qubo_posiform_hardened import create_qubo_hardened_posiform
from posiform.qubo_posiform import create_qubo_posiform
from qubo_utils import calculate_energy


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def run_sweep_transition(n_bits=1000, num_instances=10, reads_per_sweep=50):
    """
    실험 1: SA sweep 수 변화에 따른 성공률 전이.

    논문 Fig 7-8 재현: sweep 수 증가에 따른 ground-state sampling rate.
    핵심: α=0.01이 α=0.1보다 더 많은 sweep을 요구.

    각 인스턴스 × 각 sweep 값 × reads_per_sweep개 SA sample.
    성공률 = (ground state를 찾은 sample 수) / (전체 sample 수).
    """
    sweep_values = [1, 5, 10, 50, 100, 500, 1000, 5000, 10000]
    configs = [
        ('lin2', 0.01),
        ('lin2', 0.1),
        ('lin20', 0.01),
        ('lin20', 0.1),
    ]

    print("=" * 110)
    print(f"실험 1: SA Sweep Count 전이 실험 (논문 Fig 7-8 재현)")
    print(f"N={n_bits}, num_instances={num_instances}, reads_per_sweep={reads_per_sweep}")
    print(f"Sweep values: {sweep_values}")
    print("=" * 110)

    sampler = neal.SimulatedAnnealingSampler()
    all_results = []

    for ct, alpha in configs:
        print(f"\n{'─' * 90}")
        print(f"  {ct}, α={alpha}")
        print(f"{'─' * 90}")

        # QUBO 인스턴스 미리 생성
        t0 = time.time()
        instances = []
        for run in range(num_instances):
            Q, info = create_qubo_hardened_posiform(
                n_bits, max_subgraph_size=15, coeff_type=ct,
                posiform_scale=alpha, seed=run * 53
            )
            if info['posiform_is_unique']:
                instances.append((Q, info))
        gen_time = time.time() - t0
        actual = len(instances)
        print(f"  인스턴스 생성: {actual}/{num_instances} (생성 시간: {gen_time:.1f}s)")

        for sweeps in sweep_values:
            total_samples = 0
            ground_state_found = 0
            total_hamming = 0

            t0 = time.time()
            for Q, info in instances:
                target = info['target']
                te = info['target_energy']

                ss = sampler.sample_qubo(Q, num_reads=reads_per_sweep,
                                         num_sweeps=sweeps)

                for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                    total_samples += 1
                    found = ''.join(str(sample[k]) for k in range(n_bits))
                    if found == target:
                        ground_state_found += 1
                    total_hamming += hamming_distance(target, found)

            solve_time = time.time() - t0
            rate = 100.0 * ground_state_found / total_samples if total_samples > 0 else 0
            avg_h = total_hamming / total_samples if total_samples > 0 else float('nan')

            print(f"  sweeps={sweeps:<6} | GS rate: {ground_state_found:>4}/{total_samples} "
                  f"({rate:>6.2f}%) | avg Hamming: {avg_h:>6.1f} | time: {solve_time:.1f}s")

            all_results.append({
                'coeff_type': ct, 'scale': alpha, 'sweeps': sweeps,
                'gs_rate': rate, 'avg_hamming': avg_h,
                'gs_found': ground_state_found, 'total_samples': total_samples,
            })

    # 요약 테이블
    print("\n" + "=" * 110)
    print("Sweep 전이 요약 (Ground-State Sampling Rate %)")
    print("=" * 110)
    header = f"{'Config':<14} |"
    for s in sweep_values:
        header += f" {s:>7}"
    print(header)
    print("-" * (16 + 8 * len(sweep_values)))

    for ct, alpha in configs:
        label = f"{ct},α={alpha}"
        row = f"{label:<14} |"
        for s in sweep_values:
            matching = [r for r in all_results
                        if r['coeff_type'] == ct and r['scale'] == alpha
                        and r['sweeps'] == s]
            if matching:
                row += f" {matching[0]['gs_rate']:>6.1f}%"
            else:
                row += f" {'N/A':>7}"
        print(row)

    return all_results


def run_scaling_experiment(sizes=None, num_runs=30, num_reads=50, num_sweeps=500):
    """
    실험 2: N-Scaling.

    고정 sweeps에서 N 증가 → SA 성공률 감소 관찰.
    sweeps를 일부러 적게 설정하여 차이를 극대화.
    """
    if sizes is None:
        sizes = [100, 200, 500, 1000, 2000]

    configs = [
        ('lin2', 0.01),
        ('lin2', 0.1),
        ('lin20', 0.01),
    ]

    print("=" * 100)
    print(f"실험 2: N-Scaling")
    print(f"Sizes={sizes}, num_runs={num_runs}, num_reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 100)

    sampler = neal.SimulatedAnnealingSampler()
    all_results = []

    for ct, alpha in configs:
        print(f"\n{'─' * 80}")
        print(f"  {ct}, α={alpha}")
        print(f"{'─' * 80}")

        for n in sizes:
            total_samples = 0
            gs_found = 0
            total_hamming = 0

            t0 = time.time()
            for run in range(num_runs):
                Q, info = create_qubo_hardened_posiform(
                    n, max_subgraph_size=15, coeff_type=ct,
                    posiform_scale=alpha, seed=run * 31 + n
                )
                if not info['posiform_is_unique']:
                    continue

                target = info['target']
                te = info['target_energy']
                ss = sampler.sample_qubo(Q, num_reads=num_reads,
                                         num_sweeps=num_sweeps)

                best = ss.first
                found = ''.join(str(best.sample[k]) for k in range(n))
                total_samples += 1
                if found == target:
                    gs_found += 1
                total_hamming += hamming_distance(target, found)

            total_time = time.time() - t0
            rate = 100.0 * gs_found / total_samples if total_samples > 0 else 0
            avg_h = total_hamming / total_samples if total_samples > 0 else float('nan')

            print(f"  N={n:<5} | 성공: {gs_found:>2}/{total_samples} ({rate:>5.1f}%) | "
                  f"Hamming: {avg_h:>6.1f} | time: {total_time:.1f}s")

            all_results.append({
                'coeff_type': ct, 'scale': alpha, 'n': n,
                'success_rate': rate, 'avg_hamming': avg_h,
            })

    # 요약
    print("\n" + "=" * 100)
    print("N-Scaling 요약")
    print("=" * 100)
    print(f"{'Config':<14} |", end="")
    for n in sizes:
        print(f" N={n:<5}", end="")
    print("  (success %)")
    print("-" * (16 + 8 * len(sizes)))

    for ct, alpha in configs:
        label = f"{ct},α={alpha}"
        print(f"{label:<14} |", end="")
        for n in sizes:
            matching = [r for r in all_results
                        if r['coeff_type'] == ct and r['scale'] == alpha
                        and r['n'] == n]
            if matching:
                print(f" {matching[0]['success_rate']:>5.1f}%", end=" ")
            else:
                print(f" {'N/A':>6}", end=" ")
        print()

    return all_results


def run_hardened_vs_plain(n_bits=1000, num_runs=30, num_reads=50, num_sweeps=500):
    """
    실험 3: Hardened vs Plain Posiform 비교.

    동일 N, 동일 SA 자원에서 hardened가 얼마나 더 어려운지.
    """
    print("=" * 100)
    print(f"실험 3: Hardened vs Plain Posiform 비교")
    print(f"N={n_bits}, num_runs={num_runs}, num_reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 100)

    sampler = neal.SimulatedAnnealingSampler()

    configs = [
        ('Hard(lin2,α=0.01)', 'lin2', 0.01),
        ('Hard(lin2,α=0.1)',  'lin2', 0.1),
        ('Hard(lin20,α=0.01)', 'lin20', 0.01),
        ('Hard(lin20,α=0.1)', 'lin20', 0.1),
        ('Plain Posiform',    None, None),
    ]

    results = {name: {'gs': 0, 'total': 0, 'hammings': []} for name, _, _ in configs}

    for run in range(num_runs):
        first_target = None

        for name, ct, alpha in configs:
            if ct is not None:
                Q, info = create_qubo_hardened_posiform(
                    n_bits, max_subgraph_size=15, coeff_type=ct,
                    posiform_scale=alpha, seed=run * 77
                )
                target = info['target']
                te = info['target_energy']
                if first_target is None:
                    first_target = target
            else:
                target = first_target
                Q, info = create_qubo_posiform(target, coeff_range=(1.0, 3.0),
                                               seed=run * 77)
                te = info['target_energy']

            ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
            found = ''.join(str(ss.first.sample[k]) for k in range(n_bits))

            hdist = hamming_distance(target, found)
            results[name]['total'] += 1
            if found == target:
                results[name]['gs'] += 1
            results[name]['hammings'].append(hdist)

        if (run + 1) % 5 == 0:
            print(f"  Run {run+1}/{num_runs} 완료...")

    # 요약
    print("\n" + "=" * 100)
    print("Hardened vs Plain 요약")
    print("=" * 100)
    print(f"{'Method':<22} | {'Success%':<10} | {'Count':<10} | {'Avg Hamming':<12}")
    print("-" * 60)

    for name, _, _ in configs:
        r = results[name]
        rate = 100.0 * r['gs'] / r['total'] if r['total'] > 0 else 0
        avg_h = np.mean(r['hammings']) if r['hammings'] else float('nan')
        print(f"{name:<22} | {rate:<10.1f} | {r['gs']:>3}/{r['total']:<5} | {avg_h:<12.1f}")

    return results


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    mode = "sweep"
    n_bits = 1000
    num_runs = 10

    if len(sys.argv) > 1:
        if sys.argv[1] == "--sweep":
            mode = "sweep"
        elif sys.argv[1] == "--scaling":
            mode = "scaling"
        elif sys.argv[1] == "--compare":
            mode = "compare"
        else:
            n_bits = int(sys.argv[1])

    if len(sys.argv) > 2:
        num_runs = int(sys.argv[2])

    if mode == "sweep":
        run_sweep_transition(n_bits=n_bits, num_instances=num_runs, reads_per_sweep=50)
    elif mode == "scaling":
        run_scaling_experiment(num_runs=num_runs)
    elif mode == "compare":
        run_hardened_vs_plain(n_bits=n_bits, num_runs=num_runs)
