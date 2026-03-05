"""
α < 0 실험: Posiform이 anti-signal로 작용하는지 검증

가설:
  α > 0: posiform이 target 방향으로 gradient signal 제공 (SA 성공률 증가)
  α = 0: posiform 제거 → signal 없음 → 축퇴에 의한 실패
  α < 0: posiform이 target **반대 방향**으로 signal 제공 (SA가 target에서 멀어짐)

  이 세 조건을 비교하면 posiform이 단순한 축퇴 제거 장치가 아니라
  방향성 있는 signal임을 직접 증명할 수 있다.

측정 항목:
  1. Target sampling rate (sweep 전이)
  2. GS 에너지 도달률 vs target 일치율 (에너지 분석)
  3. SA 해의 target 대비 Hamming distance 분포
  4. 에너지 분해: Σ R_i 성분 vs αP 성분
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hardened_posiform.qubo_posiform_hardened import create_qubo_hardened_posiform
from qubo_utils import calculate_energy


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def bitstring_complement(s):
    return ''.join('1' if c == '0' else '0' for c in s)


def run_negative_alpha_experiment(n_bits=500, num_instances=10, reads_per_sweep=50):
    """
    α = -0.1, -0.01, 0, 0.01, 0.1 sweep 전이 + 에너지 분석 + Hamming 분포.
    """
    sweep_values = [1, 5, 10, 50, 100, 500, 1000, 5000]
    alpha_values = [-0.1, -0.01, 0.0, 0.01, 0.1]
    coeff_types = ['lin2', 'lin20']

    configs = [(ct, a) for ct in coeff_types for a in alpha_values]

    print("=" * 120)
    print("α < 0 실험: Posiform Anti-Signal 검증")
    print(f"N={n_bits}, num_instances={num_instances}, reads_per_sweep={reads_per_sweep}")
    print(f"α values: {alpha_values}")
    print(f"Sweep values: {sweep_values}")
    print("=" * 120)

    sampler = neal.SimulatedAnnealingSampler()
    all_results = []

    for ct, alpha in configs:
        print(f"\n{'━' * 100}")
        print(f"  {ct}, α={alpha}")
        print(f"{'━' * 100}")

        # QUBO 인스턴스 생성
        t0 = time.time()
        instances = []
        for run in range(num_instances):
            Q, info = create_qubo_hardened_posiform(
                n_bits, max_subgraph_size=15, coeff_type=ct,
                posiform_scale=alpha, seed=run * 53
            )
            instances.append((Q, info))
        gen_time = time.time() - t0
        print(f"  인스턴스 생성: {len(instances)}개 ({gen_time:.1f}s)")

        # 축퇴도 통계
        degeneracies = [info['random_total_degenerate'] for _, info in instances]
        print(f"  Random subproblem 평균 축퇴도: {np.mean(degeneracies):.1f}")

        # target이 ground state인지 검사 (α<0이면 아닐 수 있음)
        print(f"\n  [Target GS 여부 검사 — single-flip test]")
        for idx, (Q, info) in enumerate(instances):
            target = info['target']
            target_energy = calculate_energy(target, Q)
            min_flip_delta = float('inf')
            for i in range(n_bits):
                flipped = list(target)
                flipped[i] = '0' if flipped[i] == '1' else '1'
                delta = calculate_energy(''.join(flipped), Q) - target_energy
                min_flip_delta = min(min_flip_delta, delta)
            is_local_min = min_flip_delta > -1e-10
            print(f"    inst {idx}: target_E={target_energy:>10.4f} | "
                  f"min_flip_delta={min_flip_delta:>+8.4f} | "
                  f"local_min={'O' if is_local_min else 'X'}")

        # Sweep 전이
        for sweeps in sweep_values:
            total_samples = 0
            target_found = 0
            total_hamming = 0
            hamming_list = []

            t0 = time.time()
            for Q, info in instances:
                target = info['target']

                ss = sampler.sample_qubo(Q, num_reads=reads_per_sweep,
                                         num_sweeps=sweeps)

                for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                    total_samples += 1
                    found = ''.join(str(sample[k]) for k in range(n_bits))
                    if found == target:
                        target_found += 1
                    hd = hamming_distance(target, found)
                    total_hamming += hd
                    hamming_list.append(hd)

            solve_time = time.time() - t0
            rate = 100.0 * target_found / total_samples if total_samples > 0 else 0
            avg_h = total_hamming / total_samples if total_samples > 0 else float('nan')
            h_arr = np.array(hamming_list)

            print(f"  sw={sweeps:<6} | target rate: {target_found:>4}/{total_samples} "
                  f"({rate:>6.2f}%) | Hamming avg={avg_h:>6.1f} "
                  f"med={np.median(h_arr):>5.0f} "
                  f"min={h_arr.min():>3} max={h_arr.max():>3} | {solve_time:.1f}s")

            all_results.append({
                'coeff_type': ct, 'alpha': alpha, 'sweeps': sweeps,
                'target_rate': rate, 'avg_hamming': avg_h,
                'median_hamming': float(np.median(h_arr)),
                'target_found': target_found, 'total_samples': total_samples,
            })

    # ==========================================
    # 요약 테이블 1: Target Sampling Rate
    # ==========================================
    print("\n" + "=" * 120)
    print("요약 1: Target Sampling Rate (%)")
    print("=" * 120)

    for ct in coeff_types:
        print(f"\n  [{ct}]")
        header = f"  {'α':<8} |"
        for s in sweep_values:
            header += f" sw={s:<5}"
        print(header)
        print(f"  " + "-" * (10 + 8 * len(sweep_values)))

        for alpha in alpha_values:
            row = f"  {alpha:<8} |"
            for s in sweep_values:
                matching = [r for r in all_results
                            if r['coeff_type'] == ct and r['alpha'] == alpha
                            and r['sweeps'] == s]
                if matching:
                    row += f" {matching[0]['target_rate']:>5.1f}% "
                else:
                    row += f" {'N/A':>6} "
            print(row)

    # ==========================================
    # 요약 테이블 2: Median Hamming Distance
    # ==========================================
    print("\n" + "=" * 120)
    print("요약 2: Median Hamming Distance (SA 해 vs target)")
    print("=" * 120)

    for ct in coeff_types:
        print(f"\n  [{ct}]")
        header = f"  {'α':<8} |"
        for s in sweep_values:
            header += f" sw={s:<5}"
        print(header)
        print(f"  " + "-" * (10 + 8 * len(sweep_values)))

        for alpha in alpha_values:
            row = f"  {alpha:<8} |"
            for s in sweep_values:
                matching = [r for r in all_results
                            if r['coeff_type'] == ct and r['alpha'] == alpha
                            and r['sweeps'] == s]
                if matching:
                    row += f" {matching[0]['median_hamming']:>6.0f} "
                else:
                    row += f" {'N/A':>6} "
            print(row)

    # ==========================================
    # 요약 테이블 3: 대칭성 비교 (α vs -α)
    # ==========================================
    print("\n" + "=" * 120)
    print("요약 3: α vs -α 대칭성 비교 (sw=5000)")
    print("=" * 120)
    print(f"  N={n_bits} (random Hamming ≈ {n_bits/2:.0f})")
    print()
    print(f"  {'Config':<12} | {'target rate':<12} | {'Hamming med':<12} | 해석")
    print(f"  " + "-" * 70)

    for ct in coeff_types:
        for alpha in alpha_values:
            matching = [r for r in all_results
                        if r['coeff_type'] == ct and r['alpha'] == alpha
                        and r['sweeps'] == 5000]
            if matching:
                r = matching[0]
                # 해석
                if alpha > 0:
                    interp = "posiform signal → target 방향"
                elif alpha == 0:
                    interp = "signal 없음 → 축퇴"
                else:
                    interp = "anti-signal → target 반대 방향?"
                print(f"  {ct},α={alpha:<5} | {r['target_rate']:>6.1f}%     "
                      f"| {r['median_hamming']:>6.0f}       | {interp}")
        print()

    return all_results


def run_energy_decomposition(n_bits=500, num_instances=10, num_reads=50, num_sweeps=1000):
    """
    α < 0에서 에너지 분해 분석.

    SA가 찾은 해에서:
      - Σ R_i 성분 (random QUBO 에너지)
      - αP 성분 (posiform 에너지)
    을 분리하여, anti-signal이 어떤 방향으로 작용하는지 확인.
    """
    alpha_values = [-0.1, -0.01, 0.0, 0.01, 0.1]
    coeff_types = ['lin2', 'lin20']
    configs = [(ct, a) for ct in coeff_types for a in alpha_values]

    print("\n" + "=" * 120)
    print("에너지 분해 분석: Σ R_i 성분 vs αP 성분")
    print(f"N={n_bits}, instances={num_instances}, reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 120)

    sampler = neal.SimulatedAnnealingSampler()
    summary = []

    for ct, alpha in configs:
        print(f"\n{'─' * 100}")
        print(f"  {ct}, α={alpha}")
        print(f"{'─' * 100}")

        gs_energy_match = 0
        target_match = 0
        total_samples = 0
        hamming_all = []
        energy_diffs = []

        # 에너지 분해용
        random_energy_gaps = []  # Σ R_i(SA해) - Σ R_i(target)
        posiform_energy_gaps = []  # P(SA해) - P(target)

        for run in range(num_instances):
            Q, info = create_qubo_hardened_posiform(
                n_bits, max_subgraph_size=15, coeff_type=ct,
                posiform_scale=alpha, seed=run * 53
            )
            target = info['target']
            target_energy = info['target_energy']

            # Σ R_i와 P를 분리해서 다시 생성 (에너지 분해용)
            Q_no_posiform, info_no_p = create_qubo_hardened_posiform(
                n_bits, max_subgraph_size=15, coeff_type=ct,
                posiform_scale=0.0, seed=run * 53
            )
            # P = (Q - Q_no_posiform) / α if α ≠ 0
            # 직접 posiform 에너지 계산: P(x) = (E_Q(x) - E_R(x)) / α

            target_R_energy = calculate_energy(target, Q_no_posiform)

            ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)

            for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                total_samples += 1
                found = ''.join(str(sample[k]) for k in range(n_bits))
                diff = energy - target_energy
                energy_diffs.append(diff)
                hd = hamming_distance(target, found)
                hamming_all.append(hd)

                if abs(diff) < 1e-10:
                    gs_energy_match += 1
                if found == target:
                    target_match += 1

                # 에너지 분해
                found_R_energy = calculate_energy(found, Q_no_posiform)
                random_energy_gaps.append(found_R_energy - target_R_energy)
                if abs(alpha) > 1e-15:
                    found_P_energy = (energy - found_R_energy) / alpha
                    target_P_energy = (target_energy - target_R_energy) / alpha
                    posiform_energy_gaps.append(found_P_energy - target_P_energy)

            print(f"  inst {run}: target_E={target_energy:>10.4f} | "
                  f"Σ R(target)={target_R_energy:>10.4f}")

        h_arr = np.array(hamming_all)
        e_arr = np.array(energy_diffs)
        r_arr = np.array(random_energy_gaps)

        gs_pct = 100.0 * gs_energy_match / total_samples
        tgt_pct = 100.0 * target_match / total_samples

        print(f"\n  *** 종합 ({ct}, α={alpha}) ***")
        print(f"  GS 에너지 일치:  {gs_energy_match}/{total_samples} ({gs_pct:.1f}%)")
        print(f"  Target 일치:     {target_match}/{total_samples} ({tgt_pct:.1f}%)")
        print(f"  Hamming: mean={h_arr.mean():.1f}, median={np.median(h_arr):.0f}, "
              f"min={h_arr.min()}, max={h_arr.max()}")
        print(f"  에너지 차이 (Q(SA)-Q(target)): mean={e_arr.mean():.4f}, "
              f"min={e_arr.min():.4f}, median={np.median(e_arr):.4f}")
        print(f"  Σ R_i gap (R(SA)-R(target)):   mean={r_arr.mean():.4f}, "
              f"min={r_arr.min():.4f}, median={np.median(r_arr):.4f}")
        if len(posiform_energy_gaps) > 0:
            p_arr = np.array(posiform_energy_gaps)
            print(f"  P gap (P(SA)-P(target)):       mean={p_arr.mean():.4f}, "
                  f"min={p_arr.min():.4f}, median={np.median(p_arr):.4f}")

        summary.append({
            'config': f"{ct},α={alpha}",
            'gs_pct': gs_pct, 'tgt_pct': tgt_pct,
            'hamming_med': float(np.median(h_arr)),
            'R_gap_mean': float(r_arr.mean()),
            'P_gap_mean': float(np.mean(posiform_energy_gaps)) if posiform_energy_gaps else None,
        })

    # 최종 요약
    print("\n" + "=" * 120)
    print("최종 요약: 에너지 분해")
    print("=" * 120)
    print(f"  {'Config':<16} | {'GS에너지%':>8} | {'Target%':>8} | "
          f"{'Ham med':>7} | {'R gap':>10} | {'P gap':>10}")
    print(f"  " + "-" * 80)
    for s in summary:
        p_str = f"{s['P_gap_mean']:>10.4f}" if s['P_gap_mean'] is not None else "     N/A  "
        print(f"  {s['config']:<16} | {s['gs_pct']:>7.1f}% | {s['tgt_pct']:>7.1f}% | "
              f"{s['hamming_med']:>7.0f} | {s['R_gap_mean']:>10.4f} | {p_str}")

    return summary


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    n_bits = 500
    num_instances = 10

    if len(sys.argv) > 1:
        n_bits = int(sys.argv[1])
    if len(sys.argv) > 2:
        num_instances = int(sys.argv[2])

    # 실험 1: Sweep 전이 + Hamming 분포
    run_negative_alpha_experiment(
        n_bits=n_bits, num_instances=num_instances, reads_per_sweep=50
    )

    # 실험 2: 에너지 분해 분석
    run_energy_decomposition(
        n_bits=n_bits, num_instances=num_instances,
        num_reads=50, num_sweeps=1000
    )
