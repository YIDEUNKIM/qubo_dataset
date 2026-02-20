"""
Wishart Planted Ensemble SA 실험 프레임워크

실험 1: Alpha sweep — easy-hard-easy 상전이 프로파일 확인
실험 2: Scaling — hard alpha에서 N 증가 시 성공률 지수적 감소 확인
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from wishart.qubo_wishart import create_qubo_wishart, verify_ground_state
from qubo_utils import calculate_energy


def classify_result(target, found_solution, target_energy, found_energy):
    """결과 분류: EXACT / SYM_MATCH / ENERGY_MATCH / FAIL"""
    exact_match = (found_solution == target)
    inverse_target = ''.join('0' if b == '1' else '1' for b in target)
    sym_match = (found_solution == inverse_target)
    energy_match = abs(found_energy - target_energy) < 1e-4

    if exact_match:
        return "EXACT"
    elif sym_match:
        return "SYM_MATCH"
    elif energy_match:
        return "ENERGY_MATCH"
    else:
        return "FAIL"


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def run_wishart_experiment(num_runs=10, n_bits=100, alpha_values=None,
                           num_reads=200, num_sweeps=1000):
    """
    Alpha sweep 실험: 각 alpha 값에서 SA 성공률 측정.

    Args:
        num_runs: alpha당 실험 횟수
        n_bits: 문제 크기
        alpha_values: 테스트할 alpha 값 리스트
        num_reads: SA sample 수
        num_sweeps: SA sweep 수
    """
    if alpha_values is None:
        alpha_values = [0.3, 0.5, 0.7, 0.8, 1.0, 1.5]

    print("=" * 90)
    print(f"Wishart Planted Ensemble — Alpha Sweep 실험")
    print(f"N={n_bits}, num_runs={num_runs}, num_reads={num_reads}, num_sweeps={num_sweeps}")
    print("=" * 90)

    sampler = neal.SimulatedAnnealingSampler()
    results_summary = []

    for alpha in alpha_values:
        print(f"\n--- alpha = {alpha} (M={int(alpha * n_bits)}) ---")
        print(f"{'Run':<4} | {'Result':<12} | {'Target E':<12} | {'Found E':<12} | {'Ratio':<8} | {'Hamming':<8} | {'Time':<6}")
        print("-" * 75)

        counts = {"EXACT": 0, "SYM_MATCH": 0, "ENERGY_MATCH": 0, "FAIL": 0}
        energy_ratios = []
        hamming_dists = []
        solve_times = []

        for run in range(num_runs):
            # 랜덤 타겟 생성
            target = ''.join(str(random.randint(0, 1)) for _ in range(n_bits))

            # QUBO 생성
            Q = create_qubo_wishart(target, alpha=alpha)
            target_energy = calculate_energy(target, Q)

            # SA 풀기
            t0 = time.time()
            sampleset = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
            solve_time = time.time() - t0

            best_sample = sampleset.first.sample
            best_energy = sampleset.first.energy
            found = ''.join(str(best_sample[k]) for k in range(n_bits))

            # 분류
            result = classify_result(target, found, target_energy, best_energy)
            counts[result] += 1

            # 메트릭
            ratio = best_energy / target_energy if abs(target_energy) > 1e-10 else float('nan')
            hdist = hamming_distance(target, found)
            energy_ratios.append(ratio)
            hamming_dists.append(hdist)
            solve_times.append(solve_time)

            print(f"{run+1:<4} | {result:<12} | {target_energy:<12.2f} | {best_energy:<12.2f} | {ratio:<8.4f} | {hdist:<8} | {solve_time:<6.2f}")

        # alpha 요약
        success = counts["EXACT"] + counts["SYM_MATCH"]
        success_rate = 100.0 * success / num_runs
        avg_ratio = np.nanmean(energy_ratios)
        avg_hamming = np.mean(hamming_dists)
        avg_time = np.mean(solve_times)

        print(f"\n  [요약] alpha={alpha}")
        print(f"    성공률: {success}/{num_runs} ({success_rate:.1f}%)")
        print(f"    EXACT={counts['EXACT']}, SYM={counts['SYM_MATCH']}, "
              f"ENERGY={counts['ENERGY_MATCH']}, FAIL={counts['FAIL']}")
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

    # 전체 요약 테이블
    print("\n" + "=" * 90)
    print("전체 요약 (Alpha Sweep)")
    print("=" * 90)
    print(f"{'Alpha':<8} | {'Success%':<10} | {'EXACT':<6} | {'SYM':<6} | {'ENERGY':<6} | {'FAIL':<6} | {'Avg E Ratio':<12} | {'Avg Hamming':<12}")
    print("-" * 90)
    for r in results_summary:
        c = r['counts']
        print(f"{r['alpha']:<8} | {r['success_rate']:<10.1f} | {c.get('EXACT',0):<6} | "
              f"{c.get('SYM_MATCH',0):<6} | {c.get('ENERGY_MATCH',0):<6} | {c.get('FAIL',0):<6} | "
              f"{r['avg_energy_ratio']:<12.4f} | {r['avg_hamming']:<12.1f}")

    return results_summary


def run_scaling_experiment(sizes=None, alpha=0.7, num_runs=10,
                           num_reads=200, num_sweeps=1000):
    """
    스케일링 실험: 고정 alpha에서 N 증가에 따른 성공률 변화.

    Args:
        sizes: 테스트할 N 리스트
        alpha: 고정 alpha 값
        num_runs: 크기당 실험 횟수
    """
    if sizes is None:
        sizes = [50, 100, 150, 200, 300, 500]

    print("=" * 90)
    print(f"Wishart Planted Ensemble — Scaling 실험 (alpha={alpha})")
    print(f"Sizes={sizes}, num_runs={num_runs}, num_reads={num_reads}, num_sweeps={num_sweeps}")
    print("=" * 90)

    sampler = neal.SimulatedAnnealingSampler()
    results = []

    for n in sizes:
        print(f"\n--- N = {n} ---")
        success_count = 0
        energy_ratios = []
        hamming_dists = []

        for run in range(num_runs):
            target = ''.join(str(random.randint(0, 1)) for _ in range(n))
            Q = create_qubo_wishart(target, alpha=alpha)
            target_energy = calculate_energy(target, Q)

            sampleset = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
            best_sample = sampleset.first.sample
            best_energy = sampleset.first.energy
            found = ''.join(str(best_sample[k]) for k in range(n))

            result = classify_result(target, found, target_energy, best_energy)
            if result in ("EXACT", "SYM_MATCH"):
                success_count += 1

            ratio = best_energy / target_energy if abs(target_energy) > 1e-10 else float('nan')
            energy_ratios.append(ratio)
            hamming_dists.append(hamming_distance(target, found))

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

    # 요약 테이블
    print("\n" + "=" * 70)
    print(f"Scaling 요약 (alpha={alpha})")
    print("=" * 70)
    print(f"{'N':<8} | {'Success%':<10} | {'Count':<8} | {'Avg E Ratio':<12} | {'Avg Hamming':<12}")
    print("-" * 60)
    for r in results:
        print(f"{r['n']:<8} | {r['success_rate']:<10.1f} | "
              f"{r['success_count']}/{num_runs:<5} | {r['avg_energy_ratio']:<12.4f} | "
              f"{r['avg_hamming']:<12.1f}")

    return results


def hardness_metrics(Q, target, sampler, num_reads=1000, num_sweeps=1000):
    """
    난이도 메트릭 계산.

    Returns:
        dict with: tts_99 (Time-to-Solution 99%), energy_ratio, spectral_gap_estimate
    """
    n = len(target)
    target_energy = calculate_energy(target, Q)

    # 다수 SA run으로 성공 확률 추정
    t0 = time.time()
    sampleset = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
    wall_time = time.time() - t0
    time_per_read = wall_time / num_reads

    # 성공 횟수 세기
    success_count = 0
    energies = []
    for sample, energy in zip(sampleset.samples(), sampleset.record.energy):
        found = ''.join(str(sample[k]) for k in range(n))
        inverse = ''.join('0' if b == '1' else '1' for b in target)
        if found == target or found == inverse:
            success_count += 1
        energies.append(energy)

    p_success = success_count / num_reads if num_reads > 0 else 0.0

    # TTS (Time-to-Solution) at 99% confidence
    if p_success > 0 and p_success < 1.0:
        tts_99 = time_per_read * np.log(1 - 0.99) / np.log(1 - p_success)
    elif p_success >= 1.0:
        tts_99 = time_per_read
    else:
        tts_99 = float('inf')

    # 에너지 비율 (found / target)
    best_energy = min(energies)
    energy_ratio = best_energy / target_energy if abs(target_energy) > 1e-10 else float('nan')

    # 스펙트럴 갭 추정 (에너지 분포의 표준편차로 간접 추정)
    spectral_gap_estimate = np.std(energies)

    return {
        'p_success': p_success,
        'tts_99': tts_99,
        'energy_ratio': energy_ratio,
        'best_energy': best_energy,
        'target_energy': target_energy,
        'energy_std': spectral_gap_estimate,
        'wall_time': wall_time,
    }


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    n_bits = 100
    num_runs = 10
    mode = "alpha"  # alpha / scaling / metrics

    if len(sys.argv) > 1:
        if sys.argv[1] == "--scaling":
            mode = "scaling"
            if len(sys.argv) > 2:
                alpha_arg = float(sys.argv[2])
            else:
                alpha_arg = 0.7
        elif sys.argv[1] == "--metrics":
            mode = "metrics"
        else:
            n_bits = int(sys.argv[1])

    if len(sys.argv) > 2 and mode == "alpha":
        num_runs = int(sys.argv[2])

    if mode == "alpha":
        run_wishart_experiment(
            num_runs=num_runs,
            n_bits=n_bits,
            alpha_values=[0.3, 0.5, 0.7, 0.8, 1.0, 1.5],
            num_reads=200,
            num_sweeps=1000,
        )

    elif mode == "scaling":
        run_scaling_experiment(
            sizes=[50, 100, 150, 200, 300, 500],
            alpha=alpha_arg,
            num_runs=10,
            num_reads=200,
            num_sweeps=1000,
        )

    elif mode == "metrics":
        print("=== Hardness Metrics 측정 ===")
        sampler = neal.SimulatedAnnealingSampler()
        target = ''.join(str(random.randint(0, 1)) for _ in range(n_bits))
        Q = create_qubo_wishart(target, alpha=0.7)
        metrics = hardness_metrics(Q, target, sampler, num_reads=1000, num_sweeps=1000)
        print(f"  N={n_bits}, alpha=0.7")
        print(f"  P(success): {metrics['p_success']:.4f}")
        print(f"  TTS(99%): {metrics['tts_99']:.2f}s")
        print(f"  Energy ratio: {metrics['energy_ratio']:.6f}")
        print(f"  Energy std: {metrics['energy_std']:.4f}")
        print(f"  Wall time: {metrics['wall_time']:.2f}s")
