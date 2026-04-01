"""
Δ_R vs α 난이도 실험: lin2 vs lin20 비교

예상:
  α=0: lin2(Δ_R=1.0)가 쉬움, lin20(Δ_R=0.1)이 어려움
  α>0: lin2(ρ 큼)가 어려움, lin20(ρ 작음)이 쉬움
  → α 어딘가에서 역전이 발생
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hardened_posiform.qubo_posiform_hardened import (
    create_qubo_hardened_posiform,
    partition_variables,
    generate_random_qubo,
    solve_qubo_brute_force,
)
from qubo_utils import calculate_energy


def compute_delta_r_and_degeneracy(n, max_sub, coeff_type, seed):
    """각 subblock의 GS갭과 축퇴도를 계산."""
    rng = random.Random(seed)
    partitions = partition_variables(n, max_sub)

    delta_r = float('inf')
    total_deg = 1
    block_gaps = []

    for variables in partitions:
        part_seed = rng.randint(0, 10**9)
        R = generate_random_qubo(variables, coeff_type=coeff_type, seed=part_seed)

        # brute force로 모든 에너지 계산
        var_list = sorted(variables)
        k = len(var_list)
        var_to_shift = {var: (k - 1 - idx) for idx, var in enumerate(var_list)}

        linear_terms = []
        quad_terms = []
        for (i, j), w in R.items():
            if i == j:
                linear_terms.append((var_to_shift[i], w))
            else:
                quad_terms.append((var_to_shift[i], var_to_shift[j], w))

        energies = []
        for bits in range(1 << k):
            energy = 0.0
            for shift, w in linear_terms:
                if (bits >> shift) & 1:
                    energy += w
            for si, sj, w in quad_terms:
                if ((bits >> si) & 1) and ((bits >> sj) & 1):
                    energy += w
            energies.append(energy)

        energies.sort()
        gs_e = energies[0]
        deg = sum(1 for e in energies if abs(e - gs_e) < 1e-12)
        total_deg *= deg

        # first excited state 갭
        gap = float('inf')
        for e in energies:
            if e > gs_e + 1e-12:
                gap = e - gs_e
                break
        block_gaps.append(gap)
        delta_r = min(delta_r, gap)

    return delta_r, total_deg, block_gaps


def run_experiment():
    """메인 실험: lin2 vs lin20, α별 Energy GSP 비교."""
    n_bits = 200
    max_sub = 10
    num_instances = 100
    num_reads = 50
    num_sweeps = 1000
    alphas = [0, 0.001, 0.005, 0.01, 0.03, 0.05, 0.1]

    sampler = neal.SimulatedAnnealingSampler()

    print("=" * 100)
    print(f"Δ_R vs α 난이도 실험: lin2 vs lin20")
    print(f"n={n_bits}, subgraph_size={max_sub}, instances={num_instances}, "
          f"reads={num_reads}, sweeps={num_sweeps}")
    print("=" * 100)

    # ── Part 1: Δ_R, 축퇴도 통계 ──
    print("\n[Part 1] Δ_R 및 축퇴도 통계 (10 인스턴스 샘플)")
    print("-" * 70)
    for ct in ['lin2', 'lin20']:
        dr_list = []
        deg_list = []
        for seed in range(10):
            dr, deg, _ = compute_delta_r_and_degeneracy(n_bits, max_sub, ct, seed * 53)
            dr_list.append(dr)
            deg_list.append(deg)
        print(f"  {ct:5s}: Δ_R = {np.mean(dr_list):.4f} (±{np.std(dr_list):.4f}), "
              f"축퇴도 중앙값 = {np.median(deg_list):.0f}, "
              f"범위 = [{min(deg_list):.0f}, {max(deg_list):.0f}]")

    # ── Part 2: α별 Energy GSP ──
    print(f"\n[Part 2] α별 Energy GSP (에너지 기준 성공률)")
    print("-" * 100)

    # 같은 R/P 구조에서 α만 변경 (Fixed R/P 설계)
    results = {ct: {} for ct in ['lin2', 'lin20']}

    for ct in ['lin2', 'lin20']:
        print(f"\n  === {ct} ===")

        # 인스턴스 생성 (α=0.1로 생성, 나중에 α만 변경)
        instances = []
        t0 = time.time()
        for seed in range(num_instances):
            Q, info = create_qubo_hardened_posiform(
                n_bits, max_subgraph_size=max_sub, coeff_type=ct,
                posiform_scale=0.1, seed=seed * 53
            )
            if info['posiform_is_unique']:
                # R과 P를 분리 저장
                # R 재구성
                rng = random.Random(seed * 53)
                partitions = partition_variables(n_bits, max_sub)
                R_combined = {}
                for variables in partitions:
                    part_seed = rng.randint(0, 10**9)
                    R = generate_random_qubo(variables, coeff_type=ct, seed=part_seed)
                    for key, val in R.items():
                        R_combined[key] = R_combined.get(key, 0) + val

                # P = (Q_final - R) / 0.1
                P = {}
                all_keys = set(Q.keys()) | set(R_combined.keys())
                for key in all_keys:
                    q_val = Q.get(key, 0)
                    r_val = R_combined.get(key, 0)
                    p_val = (q_val - r_val) / 0.1
                    if abs(p_val) > 1e-12:
                        P[key] = p_val

                instances.append({
                    'R': R_combined,
                    'P': P,
                    'target': info['target'],
                    'target_energy_R': info['random_total_energy'],
                })
        gen_time = time.time() - t0
        print(f"  인스턴스: {len(instances)}/{num_instances} ({gen_time:.1f}s)")

        # 각 α에 대해 Q = R + α·P 재구성하여 SA 실행
        for alpha in alphas:
            total_samples = 0
            energy_success = 0
            target_success = 0

            t0 = time.time()
            for inst in instances:
                # Q_final = R + α·P
                Q_final = {}
                for key, val in inst['R'].items():
                    Q_final[key] = val
                for key, val in inst['P'].items():
                    Q_final[key] = Q_final.get(key, 0) + alpha * val
                Q_final = {k: v for k, v in Q_final.items() if abs(v) > 1e-15}

                # target 에너지 계산
                target_e = calculate_energy(inst['target'], Q_final)

                ss = sampler.sample_qubo(Q_final, num_reads=num_reads,
                                         num_sweeps=num_sweeps)

                for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                    total_samples += 1
                    # 에너지 기준 성공
                    if abs(energy - target_e) < 1e-4:
                        energy_success += 1
                    # target 일치 성공
                    found = ''.join(str(sample[k]) for k in range(n_bits))
                    if found == inst['target']:
                        target_success += 1

            solve_time = time.time() - t0
            e_rate = 100.0 * energy_success / total_samples if total_samples else 0
            t_rate = 100.0 * target_success / total_samples if total_samples else 0

            results[ct][alpha] = (e_rate, t_rate)
            print(f"  α={alpha:<6} | Energy GSP: {e_rate:>6.1f}% | "
                  f"Target GSP: {t_rate:>6.1f}% | time: {solve_time:.1f}s")

    # ── Part 3: 비교 요약 ──
    print(f"\n{'=' * 100}")
    print(f"[Part 3] 비교 요약 — lin2 vs lin20")
    print(f"{'=' * 100}")
    print(f"{'α':>8} | {'lin2 E-GSP':>12} {'lin2 T-GSP':>12} | "
          f"{'lin20 E-GSP':>12} {'lin20 T-GSP':>12} | {'E-GSP 우위':>10}")
    print("-" * 90)
    for alpha in alphas:
        e1, t1 = results['lin2'].get(alpha, (0, 0))
        e2, t2 = results['lin20'].get(alpha, (0, 0))
        winner = "lin2" if e1 > e2 else ("lin20" if e2 > e1 else "동일")
        print(f"{alpha:>8} | {e1:>10.1f}% {t1:>10.1f}% | "
              f"{e2:>10.1f}% {t2:>10.1f}% | {winner:>10}")

    # 예측 검증
    print(f"\n{'=' * 100}")
    print("[Part 4] 예측 검증")
    print("=" * 100)

    e_lin2_a0 = results['lin2'].get(0, (0,0))[0]
    e_lin20_a0 = results['lin20'].get(0, (0,0))[0]
    print(f"  예측 1) α=0에서 lin2가 더 쉬움 (Δ_R 크고 축퇴도 높으므로)")
    print(f"    → lin2 E-GSP={e_lin2_a0:.1f}%, lin20 E-GSP={e_lin20_a0:.1f}%  "
          f"{'✓ 맞음' if e_lin2_a0 > e_lin20_a0 else '✗ 틀림'}")

    # α=0.1에서 비교
    e_lin2_high = results['lin2'].get(0.1, (0,0))[0]
    e_lin20_high = results['lin20'].get(0.1, (0,0))[0]
    print(f"  예측 2) 큰 α에서 lin20이 더 쉬움 (ρ가 10배 작으므로)")
    print(f"    → lin2 E-GSP={e_lin2_high:.1f}%, lin20 E-GSP={e_lin20_high:.1f}%  "
          f"{'✓ 맞음' if e_lin20_high >= e_lin2_high else '✗ 틀림 (둘 다 높아 차이 미미할 수 있음)'}")

    # 역전 지점 탐색
    print(f"\n  역전 분석:")
    for alpha in alphas:
        e1 = results['lin2'].get(alpha, (0,0))[0]
        e2 = results['lin20'].get(alpha, (0,0))[0]
        diff = e1 - e2
        bar = "█" * int(abs(diff) / 2) if abs(diff) > 1 else "·"
        direction = f"lin2 +{diff:.1f}%" if diff > 0 else f"lin20 +{-diff:.1f}%"
        print(f"    α={alpha:<6}: {direction:>16}  {bar}")


if __name__ == '__main__':
    run_experiment()
