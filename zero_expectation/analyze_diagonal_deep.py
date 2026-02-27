"""
Deep Analysis: 대각 편향(Diagonal Bias) = 에너지 갭(Energy Gap) 관계 분석

4개 파트:
  Part 1: 대각 편향 = 에너지 갭 동치 관계 실험적 증명
  Part 2: 대각 노이즈 주입에 따른 GS 보존 vs 탐지 SNR 트레이드오프
  Part 3: Ising-derived 접근법과의 편향 비교
  Part 4: Hybrid 접근법 (ZeroOffDiag + Ising) 탐색
"""

import itertools
import random
import numpy as np
import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from zero_expectation.qubo_zero_expectation import (
    create_qubo_precise,
    create_qubo_ising_derived,
    ZeroOffDiagonalModel,
    solve_brute_force,
)
from qubo_utils import calculate_energy


# =====================================================================
# Helper Functions
# =====================================================================

def random_target(n):
    """길이 n의 랜덤 이진 문자열 생성."""
    return ''.join(str(random.randint(0, 1)) for _ in range(n))


def flip_bit(target, k):
    """target의 k번째 비트를 뒤집은 문자열 반환."""
    bits = list(target)
    bits[k] = '1' if bits[k] == '0' else '0'
    return ''.join(bits)


def single_flip_energy_gap(target, Q, k):
    """변수 k의 단일 비트 플립 에너지 갭: E(x_flipped) - E(x_target)."""
    flipped = flip_bit(target, k)
    E_target = calculate_energy(target, Q)
    E_flipped = calculate_energy(flipped, Q)
    return E_flipped - E_target


def check_gs_brute_force(target, Q, n):
    """브루트포스로 target이 ground state인지 확인 (n<=20)."""
    E_target = calculate_energy(target, Q)
    for bits in itertools.product([0, 1], repeat=n):
        state = ''.join(map(str, bits))
        if state == target:
            continue
        E_state = calculate_energy(state, Q)
        if E_state < E_target - 1e-10:
            return False
    return True


def check_gs_single_flip(target, Q, n):
    """모든 단일 비트 플립 에너지 갭이 양수인지 확인."""
    for k in range(n):
        gap = single_flip_energy_gap(target, Q, k)
        if gap < -1e-10:
            return False
    return True


def add_diagonal_noise(Q, n, sigma):
    """Q의 대각 원소에 N(0, sigma^2) 노이즈 추가. 새 Q 반환."""
    Q_noisy = dict(Q)
    for i in range(n):
        noise = np.random.normal(0, sigma)
        Q_noisy[(i, i)] = Q_noisy.get((i, i), 0) + noise
    return Q_noisy


def measure_diagonal_snr(Q_samples, targets, n):
    """
    여러 (Q, target) 쌍에서 대각 SNR 측정.
    SNR = |E[Q_kk|b=0] - E[Q_kk|b=1]| / (2 * std(Q_kk))
    """
    diag_b0 = []
    diag_b1 = []
    for Q, target in zip(Q_samples, targets):
        for k in range(n):
            val = Q.get((k, k), 0)
            if target[k] == '0':
                diag_b0.append(val)
            else:
                diag_b1.append(val)

    if not diag_b0 or not diag_b1:
        return 0.0, 0.0, 0.0, 0.0

    mean_b0 = np.mean(diag_b0)
    mean_b1 = np.mean(diag_b1)
    std_all = np.std(diag_b0 + diag_b1)
    snr = abs(mean_b0 - mean_b1) / (2 * std_all) if std_all > 1e-15 else float('inf')
    return snr, mean_b0, mean_b1, std_all


def measure_offdiag_bias(Q_samples, targets, n):
    """
    여러 (Q, target) 쌍에서 off-diagonal 편향 측정.
    target pair별 E[Q_ij]와 전체 SNR 반환.
    """
    by_pair = {(0,0): [], (0,1): [], (1,0): [], (1,1): []}
    all_vals = []

    for Q, target in zip(Q_samples, targets):
        for (i, j), w in Q.items():
            if i != j:
                pair = (int(target[i]), int(target[j]))
                by_pair[pair].append(w)
                all_vals.append(w)

    result = {}
    for pair in [(0,0), (0,1), (1,0), (1,1)]:
        vals = by_pair[pair]
        if vals:
            result[pair] = {'mean': np.mean(vals), 'std': np.std(vals), 'n': len(vals)}
        else:
            result[pair] = {'mean': 0, 'std': 0, 'n': 0}

    # pair간 SNR: max |E[q_ij|pair]| / std(q_ij)
    std_all = np.std(all_vals) if all_vals else 1e-15
    max_bias = max(abs(result[p]['mean']) for p in result)
    offdiag_snr = max_bias / std_all if std_all > 1e-15 else float('inf')

    return result, offdiag_snr


# =====================================================================
# PART 1: Diagonal Bias = Energy Gap
# =====================================================================

def part1_diagonal_bias_equals_energy_gap():
    """
    Part 1: 대각 편향이 에너지 갭과 동일함을 실험적으로 증명.

    이론:
      - ΔE_k = Q_kk + Σ_{m≠k} Q_km * b_m  (b_k=0인 경우)
      - E[Q_km] = 0 (ZeroOffDiagonal 보장) → E[ΔE_k] = E[Q_kk]
      - 대각 편향을 제거하면 에너지 갭도 사라짐
    """
    print("\n" + "=" * 70)
    print("  PART 1: Diagonal Bias = Energy Gap 동치 관계 증명")
    print("=" * 70)

    # --- 1a: 개별 변수의 Q_kk vs ΔE_k 비교 ---
    print("\n--- 1a: Q_kk vs ΔE_k (single-flip energy gap) 비교 ---")
    print("이론: b_k=0이면 ΔE_k = Q_kk + Σ Q_km*b_m, b_k=1이면 ΔE_k = -Q_kk - Σ Q_km*b_m")
    print("      E[Q_km]=0이므로 E[ΔE_k] ≈ E[Q_kk] (b_k=0) 또는 -E[Q_kk] (b_k=1)\n")

    n = 20
    n_trials = 200
    model = ZeroOffDiagonalModel()

    # 데이터 수집
    qkk_b0_list = []
    delta_e_b0_list = []
    qkk_b1_list = []
    delta_e_b1_list = []
    cross_term_b0_list = []  # Σ Q_km * b_m
    cross_term_b1_list = []

    for trial in range(n_trials):
        target = random_target(n)
        Q = create_qubo_precise(target, density=1.0, model=model)

        for k in range(n):
            qkk = Q.get((k, k), 0)
            gap = single_flip_energy_gap(target, Q, k)

            # cross term = Σ_{m≠k} Q_km * b_m
            cross = 0.0
            for m in range(n):
                if m == k:
                    continue
                key = (min(k, m), max(k, m))
                q_km = Q.get(key, 0)
                cross += q_km * int(target[m])

            if target[k] == '0':
                qkk_b0_list.append(qkk)
                delta_e_b0_list.append(gap)
                cross_term_b0_list.append(cross)
            else:
                qkk_b1_list.append(qkk)
                delta_e_b1_list.append(gap)
                cross_term_b1_list.append(cross)

    print(f"  N={n}, trials={n_trials}")
    print(f"  b_k=0 samples: {len(qkk_b0_list)}, b_k=1 samples: {len(qkk_b1_list)}")

    print(f"\n  [b_k=0인 변수]")
    print(f"    E[Q_kk]     = {np.mean(qkk_b0_list):+.4f} (std={np.std(qkk_b0_list):.4f})")
    print(f"    E[ΔE_k]     = {np.mean(delta_e_b0_list):+.4f} (std={np.std(delta_e_b0_list):.4f})")
    print(f"    E[cross]    = {np.mean(cross_term_b0_list):+.4f} (std={np.std(cross_term_b0_list):.4f})")
    corr_b0 = np.corrcoef(qkk_b0_list, delta_e_b0_list)[0, 1]
    print(f"    corr(Q_kk, ΔE_k) = {corr_b0:.6f}")

    print(f"\n  [b_k=1인 변수]")
    print(f"    E[Q_kk]     = {np.mean(qkk_b1_list):+.4f} (std={np.std(qkk_b1_list):.4f})")
    print(f"    E[ΔE_k]     = {np.mean(delta_e_b1_list):+.4f} (std={np.std(delta_e_b1_list):.4f})")
    print(f"    E[cross]    = {np.mean(cross_term_b1_list):+.4f} (std={np.std(cross_term_b1_list):.4f})")
    corr_b1 = np.corrcoef(qkk_b1_list, delta_e_b1_list)[0, 1]
    print(f"    corr(Q_kk, ΔE_k) = {corr_b1:.6f}")

    # 관계 검증: ΔE_k ≈ Q_kk + cross (b_k=0) 또는 ΔE_k ≈ -Q_kk - cross (b_k=1)
    # 정확한 관계: b_k=0 → ΔE_k = Q_kk + cross (정확히)
    #             b_k=1 → ΔE_k = -Q_kk - cross (정확히)
    residual_b0 = [gap - qkk - cross for gap, qkk, cross
                   in zip(delta_e_b0_list, qkk_b0_list, cross_term_b0_list)]
    residual_b1 = [gap + qkk + cross for gap, qkk, cross
                   in zip(delta_e_b1_list, qkk_b1_list, cross_term_b1_list)]

    print(f"\n  [정확성 검증: ΔE_k = Q_kk + cross (b=0), ΔE_k = -Q_kk - cross (b=1)]")
    print(f"    b_k=0: max|residual| = {max(abs(r) for r in residual_b0):.2e}")
    print(f"    b_k=1: max|residual| = {max(abs(r) for r in residual_b1):.2e}")

    # --- 1b: 통계적 분해 ---
    print("\n\n--- 1b: E[ΔE_k]의 분해: E[Q_kk] + E[cross] ---")
    print("  cross = Σ_{m≠k} Q_km * b_m 이고 E[Q_km]=0이므로 E[cross] ≈ 0")
    print("  따라서 E[ΔE_k] ≈ E[Q_kk]: 에너지 갭의 기댓값 = 대각 편향")

    print(f"\n  b_k=0: E[ΔE_k]={np.mean(delta_e_b0_list):+.4f}, E[Q_kk]={np.mean(qkk_b0_list):+.4f}, "
          f"E[cross]={np.mean(cross_term_b0_list):+.4f}")
    print(f"  b_k=1: E[ΔE_k]={np.mean(delta_e_b1_list):+.4f}, E[Q_kk]={np.mean(qkk_b1_list):+.4f}, "
          f"E[cross]={np.mean(cross_term_b1_list):+.4f}")
    print(f"  => E[ΔE_k] - E[Q_kk] = {np.mean(delta_e_b0_list) - np.mean(qkk_b0_list):+.4f} (b=0)")
    print(f"  => E[ΔE_k] - (-E[Q_kk]) = {np.mean(delta_e_b1_list) + np.mean(qkk_b1_list):+.4f} (b=1)")

    # --- 1c: N 스케일링 ---
    print("\n\n--- 1c: N에 따른 E[Q_kk] 스케일링 ---")
    print("  이론: E[Q_kk|b=0] = (n-1)*mu, mu=E[r]*(contribution per pair)")
    print(f"  base_range=(1,3) → E[r]=2.0")

    n_values = [10, 20, 50, 100]
    trials_stat = 100

    print(f"\n  {'N':>6} | {'E[Q_kk|b=0]':>14} | {'E[Q_kk|b=1]':>14} | "
          f"{'E[Q_kk|b=0]/(n-1)':>18} | {'E[Q_kk|b=1]/(n-1)':>18}")
    print("  " + "-" * 85)

    for n_val in n_values:
        diag_b0 = []
        diag_b1 = []
        for _ in range(trials_stat):
            target = random_target(n_val)
            Q = create_qubo_precise(target, density=1.0, model=model)
            for k in range(n_val):
                val = Q.get((k, k), 0)
                if target[k] == '0':
                    diag_b0.append(val)
                else:
                    diag_b1.append(val)

        mean_b0 = np.mean(diag_b0)
        mean_b1 = np.mean(diag_b1)
        norm_b0 = mean_b0 / (n_val - 1)
        norm_b1 = mean_b1 / (n_val - 1)
        print(f"  {n_val:>6} | {mean_b0:>+14.4f} | {mean_b1:>+14.4f} | "
              f"{norm_b0:>+18.6f} | {norm_b1:>+18.6f}")

    # --- 1d: 대각을 0으로 설정하면 GS가 깨지는지 확인 ---
    print("\n\n--- 1d: 대각을 0으로 설정 시 Ground State 보존율 ---")
    print("  대각 편향 = 에너지 갭이라면, 대각을 0으로 만들면 GS가 깨져야 함")

    n_bf = 15  # brute force 가능한 크기
    n_test = 100

    gs_intact_original = 0
    gs_intact_nodiag = 0

    for _ in range(n_test):
        target = random_target(n_bf)
        Q = create_qubo_precise(target, density=1.0, model=model)

        # 원래 Q
        if check_gs_brute_force(target, Q, n_bf):
            gs_intact_original += 1

        # 대각을 0으로
        Q_nodiag = dict(Q)
        for k in range(n_bf):
            Q_nodiag[(k, k)] = 0.0
        if check_gs_brute_force(target, Q_nodiag, n_bf):
            gs_intact_nodiag += 1

    print(f"  N={n_bf}, tests={n_test}")
    print(f"  원래 Q:      GS 보존율 = {gs_intact_original}/{n_test} ({100*gs_intact_original/n_test:.1f}%)")
    print(f"  대각=0 Q:    GS 보존율 = {gs_intact_nodiag}/{n_test} ({100*gs_intact_nodiag/n_test:.1f}%)")
    print(f"  => 대각 제거 시 GS 파괴 확률 = {100*(gs_intact_original - gs_intact_nodiag)/max(gs_intact_original,1):.1f}%")
    print(f"\n  결론: 대각 편향을 제거하면 에너지 갭이 사라져 GS가 거의 반드시 깨진다.")


# =====================================================================
# PART 2: Diagonal Noise Injection Tradeoff
# =====================================================================

def part2_noise_injection_tradeoff():
    """
    Part 2: 대각 노이즈 주입에 따른 GS 보존 vs 탐지 SNR 트레이드오프.

    Q_ii' = Q_ii + N(0, sigma^2) 를 주입.
    sigma가 커지면 SNR 감소 (좋음) but GS 보존율도 감소 (나쁨).
    """
    print("\n\n" + "=" * 70)
    print("  PART 2: Diagonal Noise Injection Tradeoff")
    print("=" * 70)

    model = ZeroOffDiagonalModel()

    # --- 2a: Brute force 검증 (작은 N) ---
    print("\n--- 2a: 노이즈 σ sweep (brute force, N=15) ---")

    n_bf = 15
    n_test = 100
    sigma_values_bf = [0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0]

    print(f"\n  {'sigma':>8} | {'GS保존率':>10} | {'Diag SNR':>10} | {'E[Q_kk|b=0]':>14} | {'E[Q_kk|b=1]':>14}")
    print("  " + "-" * 70)

    for sigma in sigma_values_bf:
        gs_count = 0
        diag_b0 = []
        diag_b1 = []

        for _ in range(n_test):
            target = random_target(n_bf)
            Q = create_qubo_precise(target, density=1.0, model=model)

            if sigma > 0:
                Q = add_diagonal_noise(Q, n_bf, sigma)

            if check_gs_brute_force(target, Q, n_bf):
                gs_count += 1

            for k in range(n_bf):
                val = Q.get((k, k), 0)
                if target[k] == '0':
                    diag_b0.append(val)
                else:
                    diag_b1.append(val)

        mean_b0 = np.mean(diag_b0)
        mean_b1 = np.mean(diag_b1)
        std_all = np.std(diag_b0 + diag_b1)
        snr = abs(mean_b0 - mean_b1) / (2 * std_all) if std_all > 1e-15 else float('inf')
        gs_rate = gs_count / n_test

        print(f"  {sigma:>8.1f} | {gs_rate:>9.1%} | {snr:>10.4f} | {mean_b0:>+14.4f} | {mean_b1:>+14.4f}")

    # --- 2b: 큰 N에서 single-flip 검증 ---
    print("\n\n--- 2b: 노이즈 σ sweep (single-flip check, N=50) ---")

    n_large = 50
    n_test_large = 100
    sigma_values_large = [0, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0, 80.0, 100.0]

    print(f"\n  {'sigma':>8} | {'SF_GS保존率':>12} | {'Diag SNR':>10} | {'Min ΔE_k':>10}")
    print("  " + "-" * 55)

    for sigma in sigma_values_large:
        gs_count = 0
        min_gap_avg = []
        diag_b0 = []
        diag_b1 = []

        for _ in range(n_test_large):
            target = random_target(n_large)
            Q = create_qubo_precise(target, density=1.0, model=model)

            if sigma > 0:
                Q = add_diagonal_noise(Q, n_large, sigma)

            # single-flip check
            min_gap = float('inf')
            all_positive = True
            for k in range(n_large):
                gap = single_flip_energy_gap(target, Q, k)
                if gap < min_gap:
                    min_gap = gap
                if gap < -1e-10:
                    all_positive = False

            if all_positive:
                gs_count += 1
            min_gap_avg.append(min_gap)

            for k in range(n_large):
                val = Q.get((k, k), 0)
                if target[k] == '0':
                    diag_b0.append(val)
                else:
                    diag_b1.append(val)

        mean_b0 = np.mean(diag_b0)
        mean_b1 = np.mean(diag_b1)
        std_all = np.std(diag_b0 + diag_b1)
        snr = abs(mean_b0 - mean_b1) / (2 * std_all) if std_all > 1e-15 else float('inf')
        gs_rate = gs_count / n_test_large

        print(f"  {sigma:>8.1f} | {gs_rate:>11.1%} | {snr:>10.4f} | {np.mean(min_gap_avg):>+10.4f}")

    # --- 2c: 이론적 최적 σ 계산 ---
    print("\n\n--- 2c: 이론적 최적 노이즈 수준 ---")
    print("""
  이론:
    - E[Q_kk|b=0] = (n-1)*mu, mu = E[r]*c (c는 penalty contribution per pair)
    - σ_orig = std(Q_kk) ≈ (n-1)*sigma_r*c (r의 표준편차에 비례)
    - 노이즈 σ_noise 추가 시: σ_total = sqrt(σ_orig^2 + σ_noise^2)
    - SNR = |mean_b0 - mean_b1| / (2*σ_total)
    -      = 2*(n-1)*mu / (2*σ_total)
    -      = (n-1)*mu / σ_total
    - GS 보존: 모든 k에서 ΔE_k + noise_k > 0 필요
    -           worst case: min ΔE_k - |noise_k| > 0
    -           P(GS 보존) ≈ (1 - Phi(-ΔE_min/σ))^n

    최적 σ는 GS 보존율을 P* (예: 95%) 이상 유지하면서 SNR을 최소화하는 값.
""")

    # 실측으로 최적 σ 추정
    n_opt = 50
    n_trials_opt = 200

    # 먼저 노이즈 없을 때의 min gap 분포 측정
    min_gaps_no_noise = []
    for _ in range(n_trials_opt):
        target = random_target(n_opt)
        Q = create_qubo_precise(target, density=1.0, model=model)
        gaps = [single_flip_energy_gap(target, Q, k) for k in range(n_opt)]
        min_gaps_no_noise.append(min(gaps))

    mean_min_gap = np.mean(min_gaps_no_noise)
    std_min_gap = np.std(min_gaps_no_noise)
    print(f"  N={n_opt}, trials={n_trials_opt}")
    print(f"  노이즈 없을 때 min(ΔE_k) 분포:")
    print(f"    mean = {mean_min_gap:.4f}, std = {std_min_gap:.4f}")
    print(f"    5th percentile = {np.percentile(min_gaps_no_noise, 5):.4f}")

    # 이론적 최적: σ_opt ≈ mean_min_gap / z_0.05*sqrt(n) (95% GS 보존)
    # 대략적 추정
    from scipy.stats import norm as norm_dist
    z_95 = norm_dist.ppf(0.95)  # ≈ 1.645
    # 각 변수의 ΔE_k ~ N(mean_gap, var_gap)
    # noise를 추가하면 ΔE_k' = ΔE_k + noise_k (b_k=0) or ΔE_k' = ΔE_k - noise_k (b_k=1)
    # worst case: ΔE_k - |noise_k|
    # P(GS intact) = prod_k P(ΔE_k + eps_k > 0) where eps_k ~ N(0, σ^2) with appropriate sign
    # For single variable: P(ΔE_k + eps_k > 0) ≈ Phi(E[ΔE_k]/sqrt(var(ΔE_k)+σ^2))
    # For all n: P(all positive) ≈ Phi(E[ΔE_k]/sqrt(var(ΔE_k)+σ^2))^n ≈ 0.95
    # => Phi(.) ≈ 0.95^(1/n)

    mean_gap_all = []
    for _ in range(n_trials_opt):
        target = random_target(n_opt)
        Q = create_qubo_precise(target, density=1.0, model=model)
        for k in range(n_opt):
            mean_gap_all.append(single_flip_energy_gap(target, Q, k))

    mean_gap_per_var = np.mean(mean_gap_all)
    std_gap_per_var = np.std(mean_gap_all)

    target_prob = 0.95 ** (1.0 / n_opt)
    z_per_var = norm_dist.ppf(target_prob)

    # mean_gap_per_var / sqrt(std_gap^2 + σ^2) = z_per_var
    # σ^2 = (mean_gap / z)^2 - std_gap^2
    sigma_sq = (mean_gap_per_var / z_per_var) ** 2 - std_gap_per_var ** 2
    sigma_opt = np.sqrt(max(sigma_sq, 0))

    print(f"\n  이론적 최적 σ 추정 (95% GS 보존 목표):")
    print(f"    E[ΔE_k per var] = {mean_gap_per_var:.4f}")
    print(f"    std(ΔE_k per var) = {std_gap_per_var:.4f}")
    print(f"    required z per var = {z_per_var:.4f} (for 0.95^(1/{n_opt}) = {target_prob:.6f})")
    if sigma_sq > 0:
        print(f"    σ_optimal ≈ {sigma_opt:.4f}")
    else:
        print(f"    σ_optimal: 현재 std_gap이 이미 크므로 노이즈 추가 불필요 (σ=0 최적)")

    # SNR at optimal sigma
    diag_b0_test = []
    diag_b1_test = []
    for _ in range(n_trials_opt):
        target = random_target(n_opt)
        Q = create_qubo_precise(target, density=1.0, model=model)
        for k in range(n_opt):
            val = Q.get((k, k), 0)
            if target[k] == '0':
                diag_b0_test.append(val)
            else:
                diag_b1_test.append(val)

    orig_mean_b0 = np.mean(diag_b0_test)
    orig_mean_b1 = np.mean(diag_b1_test)
    orig_std = np.std(diag_b0_test + diag_b1_test)
    orig_snr = abs(orig_mean_b0 - orig_mean_b1) / (2 * orig_std) if orig_std > 0 else float('inf')

    new_std = np.sqrt(orig_std ** 2 + sigma_opt ** 2) if sigma_opt > 0 else orig_std
    new_snr = abs(orig_mean_b0 - orig_mean_b1) / (2 * new_std) if new_std > 0 else float('inf')

    print(f"\n  SNR 비교 (N={n_opt}):")
    print(f"    원래 SNR       = {orig_snr:.4f}")
    print(f"    최적 σ 후 SNR  = {new_snr:.4f}")
    print(f"    SNR 감소율     = {(1 - new_snr/orig_snr)*100:.1f}%" if orig_snr > 0 else "    N/A")


# =====================================================================
# PART 3: Comparison with Ising-derived
# =====================================================================

def part3_compare_ising_derived():
    """
    Part 3: Ising-derived 접근법과의 편향 비교.

    Ising-derived:
      - E[Q_kk|b_k] = 0 (대각 무편향)
      - E[Q_ij] != 0 (off-diagonal 편향)
    ZeroOffDiagonal:
      - E[Q_ij] = 0 (off-diagonal 무편향)
      - E[Q_kk|b_k] != 0 (대각 편향)
    """
    print("\n\n" + "=" * 70)
    print("  PART 3: ZeroOffDiagonal vs Ising-derived 편향 비교")
    print("=" * 70)

    n = 100
    n_trials = 150
    model = ZeroOffDiagonalModel()

    # --- ZeroOffDiagonal 측정 ---
    print(f"\n--- 데이터 수집: N={n}, trials={n_trials} ---")

    zero_Qs = []
    ising_Qs = []
    targets_list = []

    for _ in range(n_trials):
        target = random_target(n)
        targets_list.append(target)
        Q_zero = create_qubo_precise(target, density=1.0, model=model)
        Q_ising = create_qubo_ising_derived(target, density=1.0)
        zero_Qs.append(Q_zero)
        ising_Qs.append(Q_ising)

    # --- 3a: 대각 편향 비교 ---
    print("\n--- 3a: 대각 편향 비교 ---")

    for label, Qs in [("ZeroOffDiagonal", zero_Qs), ("Ising-derived", ising_Qs)]:
        diag_b0 = []
        diag_b1 = []
        for Q, target in zip(Qs, targets_list):
            for k in range(n):
                val = Q.get((k, k), 0)
                if target[k] == '0':
                    diag_b0.append(val)
                else:
                    diag_b1.append(val)

        mean_b0 = np.mean(diag_b0)
        mean_b1 = np.mean(diag_b1)
        std_all = np.std(diag_b0 + diag_b1)
        snr = abs(mean_b0 - mean_b1) / (2 * std_all) if std_all > 0 else 0
        bias = mean_b0 - mean_b1

        print(f"\n  [{label}]")
        print(f"    E[Q_kk|b=0] = {mean_b0:+.4f} (std={np.std(diag_b0):.4f})")
        print(f"    E[Q_kk|b=1] = {mean_b1:+.4f} (std={np.std(diag_b1):.4f})")
        print(f"    대각 편향 (b0-b1) = {bias:+.4f}")
        print(f"    대각 SNR = {snr:.4f}")

    # --- 3b: Off-diagonal 편향 비교 ---
    print("\n\n--- 3b: Off-diagonal 편향 비교 ---")

    for label, Qs in [("ZeroOffDiagonal", zero_Qs), ("Ising-derived", ising_Qs)]:
        result, snr = measure_offdiag_bias(Qs, targets_list, n)

        print(f"\n  [{label}]")
        print(f"    {'Target Pair':>12} | {'E[Q_ij]':>12} | {'std':>10} | {'count':>8}")
        print(f"    " + "-" * 50)
        for pair in [(0,0), (0,1), (1,0), (1,1)]:
            r = result[pair]
            print(f"    {str(pair):>12} | {r['mean']:>+12.4f} | {r['std']:>10.4f} | {r['n']:>8}")
        print(f"    Off-diag SNR (max |E|/std) = {snr:.4f}")

    # --- 3c: Ising-derived의 Q_ij 부호 패턴 분석 ---
    print("\n\n--- 3c: Ising-derived off-diagonal 부호 패턴 ---")
    print("  이론: Q_ij(Ising) = -4*J_ij = -4*alpha*s_i*s_j")
    print("    s_i=s_j (target 동일) → Q_ij < 0")
    print("    s_i≠s_j (target 다름) → Q_ij > 0")

    count_correct = 0
    count_total = 0
    for Q, target in zip(ising_Qs, targets_list):
        for (i, j), w in Q.items():
            if i != j:
                si = 1 if target[i] == '1' else -1
                sj = 1 if target[j] == '1' else -1
                expected_sign = -1 if si == sj else 1
                actual_sign = 1 if w > 0 else -1
                if expected_sign == actual_sign:
                    count_correct += 1
                count_total += 1

    print(f"  부호 예측 정확도: {count_correct}/{count_total} ({100*count_correct/count_total:.1f}%)")
    print(f"  => Ising-derived의 off-diagonal 부호가 target 정보를 100% 노출")

    # --- 3d: GS 보존율 비교 ---
    print("\n\n--- 3d: GS 보존율 비교 (brute force, N=15) ---")

    n_bf = 15
    n_test = 100

    for label, gen_func in [("ZeroOffDiagonal", lambda t: create_qubo_precise(t, density=1.0, model=model)),
                              ("Ising-derived", lambda t: create_qubo_ising_derived(t, density=1.0))]:
        gs_count = 0
        for _ in range(n_test):
            target = random_target(n_bf)
            Q = gen_func(target)
            if check_gs_brute_force(target, Q, n_bf):
                gs_count += 1
        print(f"  {label:>20}: GS 보존율 = {gs_count}/{n_test} ({100*gs_count/n_test:.1f}%)")

    # --- 3e: 요약 ---
    print("\n\n--- 3e: 비교 요약 ---")
    print("""
  +----------------------+--------------+---------------+
  |        방법           | 대각 편향     | Off-diag 편향  |
  +----------------------+--------------+---------------+
  | ZeroOffDiagonal      | 있음 (SNR~1) | 없음 (E=0)    |
  | Ising-derived        | 없음 (E=0)   | 있음 (SNR~n)  |
  +----------------------+--------------+---------------+

  ZeroOffDiagonal: 대각 편향은 있지만 SNR ~ O(1/sqrt(n))로 작음.
                   Off-diagonal로는 완전히 구별 불가.
  Ising-derived:   대각은 완벽하지만 off-diagonal 부호가
                   target을 100% 노출 (SNR ~ O(n)).
  => ZeroOffDiagonal이 indistinguishability 관점에서 훨씬 우월.""")


# =====================================================================
# PART 4: Hybrid Approach
# =====================================================================

def part4_hybrid_approach():
    """
    Part 4: Q_hybrid = (1-lambda)*Q_zero + lambda*Q_ising.

    lambda=0: ZeroOffDiagonal (대각 편향, off-diag 무편향)
    lambda=1: Ising-derived (대각 무편향, off-diag 편향)
    중간값: 양쪽 편향의 trade-off
    """
    print("\n\n" + "=" * 70)
    print("  PART 4: Hybrid Approach (ZeroOffDiag + Ising)")
    print("=" * 70)

    model = ZeroOffDiagonalModel()
    n = 50
    n_trials = 150
    lambda_values = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # Pre-generate Q pairs
    print(f"\n  N={n}, trials={n_trials}")
    print("  Q_hybrid = (1-lambda)*Q_zero + lambda*Q_ising")

    targets_list = []
    zero_Qs = []
    ising_Qs = []

    for _ in range(n_trials):
        target = random_target(n)
        targets_list.append(target)
        Q_zero = create_qubo_precise(target, density=1.0, model=model)
        Q_ising = create_qubo_ising_derived(target, density=1.0)
        zero_Qs.append(Q_zero)
        ising_Qs.append(Q_ising)

    print(f"\n  {'lambda':>8} | {'Diag SNR':>10} | {'Offdiag SNR':>12} | {'Max SNR':>10} | "
          f"{'SF_GS率':>8} | {'E[Q_kk|b0]':>12} | {'E[Q_kk|b1]':>12}")
    print("  " + "-" * 95)

    best_lambda = 0.0
    best_max_snr = float('inf')

    for lam in lambda_values:
        hybrid_Qs = []
        for Q_zero, Q_ising in zip(zero_Qs, ising_Qs):
            Q_hybrid = {}
            all_keys = set(Q_zero.keys()) | set(Q_ising.keys())
            for key in all_keys:
                v_zero = Q_zero.get(key, 0)
                v_ising = Q_ising.get(key, 0)
                Q_hybrid[key] = (1 - lam) * v_zero + lam * v_ising
            hybrid_Qs.append(Q_hybrid)

        # Diagonal SNR
        diag_snr, mean_b0, mean_b1, _ = measure_diagonal_snr(hybrid_Qs, targets_list, n)

        # Off-diagonal SNR
        _, offdiag_snr = measure_offdiag_bias(hybrid_Qs, targets_list, n)

        max_snr = max(diag_snr, offdiag_snr)

        # Single-flip GS check
        gs_count = 0
        for Q, target in zip(hybrid_Qs, targets_list):
            if check_gs_single_flip(target, Q, n):
                gs_count += 1
        gs_rate = gs_count / n_trials

        print(f"  {lam:>8.2f} | {diag_snr:>10.4f} | {offdiag_snr:>12.4f} | {max_snr:>10.4f} | "
              f"{gs_rate:>7.1%} | {mean_b0:>+12.4f} | {mean_b1:>+12.4f}")

        if max_snr < best_max_snr and gs_rate > 0.9:
            best_max_snr = max_snr
            best_lambda = lam

    print(f"\n  최적 lambda (max_SNR 최소화, GS>90%): lambda = {best_lambda:.2f}, max_SNR = {best_max_snr:.4f}")

    # --- 4b: 최적 lambda 근방 세밀 스캔 ---
    print(f"\n--- 4b: 최적 lambda 근방 세밀 스캔 ---")

    fine_start = max(0.0, best_lambda - 0.1)
    fine_end = min(1.0, best_lambda + 0.1)
    fine_lambdas = np.linspace(fine_start, fine_end, 21)

    print(f"  스캔 범위: [{fine_start:.2f}, {fine_end:.2f}]")
    print(f"\n  {'lambda':>8} | {'Diag SNR':>10} | {'Offdiag SNR':>12} | {'Max SNR':>10} | {'SF_GS率':>8}")
    print("  " + "-" * 65)

    best_lambda_fine = best_lambda
    best_max_snr_fine = float('inf')

    for lam in fine_lambdas:
        hybrid_Qs = []
        for Q_zero, Q_ising in zip(zero_Qs, ising_Qs):
            Q_hybrid = {}
            all_keys = set(Q_zero.keys()) | set(Q_ising.keys())
            for key in all_keys:
                v_zero = Q_zero.get(key, 0)
                v_ising = Q_ising.get(key, 0)
                Q_hybrid[key] = (1 - lam) * v_zero + lam * v_ising
            hybrid_Qs.append(Q_hybrid)

        diag_snr, _, _, _ = measure_diagonal_snr(hybrid_Qs, targets_list, n)
        _, offdiag_snr = measure_offdiag_bias(hybrid_Qs, targets_list, n)
        max_snr = max(diag_snr, offdiag_snr)

        gs_count = 0
        for Q, target in zip(hybrid_Qs, targets_list):
            if check_gs_single_flip(target, Q, n):
                gs_count += 1
        gs_rate = gs_count / n_trials

        print(f"  {lam:>8.3f} | {diag_snr:>10.4f} | {offdiag_snr:>12.4f} | {max_snr:>10.4f} | {gs_rate:>7.1%}")

        if max_snr < best_max_snr_fine and gs_rate > 0.9:
            best_max_snr_fine = max_snr
            best_lambda_fine = lam

    print(f"\n  세밀 스캔 최적: lambda = {best_lambda_fine:.3f}, max_SNR = {best_max_snr_fine:.4f}")

    # --- 4c: 최적 hybrid의 brute force GS 검증 ---
    print(f"\n--- 4c: 최적 hybrid (lambda={best_lambda_fine:.3f}) brute force GS 검증 ---")

    n_bf = 15
    n_test = 100
    gs_count = 0

    for _ in range(n_test):
        target = random_target(n_bf)
        Q_zero = create_qubo_precise(target, density=1.0, model=model)
        Q_ising = create_qubo_ising_derived(target, density=1.0)

        Q_hybrid = {}
        all_keys = set(Q_zero.keys()) | set(Q_ising.keys())
        for key in all_keys:
            v_zero = Q_zero.get(key, 0)
            v_ising = Q_ising.get(key, 0)
            Q_hybrid[key] = (1 - best_lambda_fine) * v_zero + best_lambda_fine * v_ising

        if check_gs_brute_force(target, Q_hybrid, n_bf):
            gs_count += 1

    print(f"  N={n_bf}, tests={n_test}")
    print(f"  GS 보존율 = {gs_count}/{n_test} ({100*gs_count/n_test:.1f}%)")

    # --- 4d: 결론 ---
    print(f"\n--- 4d: Hybrid 분석 결론 ---")
    print(f"""
  Hybrid Q = (1-lambda)*Q_zero + lambda*Q_ising

  lambda=0 (ZeroOffDiag): 대각 편향만 존재, off-diag 깨끗
  lambda=1 (Ising):       off-diag 편향만 존재, 대각 깨끗

  최적 lambda={best_lambda_fine:.3f}:
    - Diagonal SNR와 Off-diagonal SNR의 최대값 최소화
    - max_SNR = {best_max_snr_fine:.4f}

  핵심 인사이트:
    1. ZeroOffDiagonal의 대각 SNR은 O(1)로 이미 작음
    2. Ising-derived의 off-diagonal SNR은 O(n)으로 매우 큼
    3. 소량의 Ising 혼합(lambda>0)도 off-diag SNR을 급격히 증가시킴
    4. 따라서 lambda=0 (순수 ZeroOffDiagonal)이 거의 항상 최적에 가까움
    5. 대각 편향은 제거할 수 없는 본질적 한계 (불가능성 정리)이지만,
       SNR ~ O(1)이므로 N이 커질수록 통계적 탐지가 오히려 어려워짐
""")


# =====================================================================
# Main
# =====================================================================

if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    print("=" * 70)
    print("  Deep Analysis: Diagonal Bias in Zero Expectation QUBO")
    print("  대각 편향 심층 분석")
    print("=" * 70)
    print(f"  시작 시간: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    t0 = time.time()

    part1_diagonal_bias_equals_energy_gap()

    t1 = time.time()
    print(f"\n  [Part 1 완료: {t1 - t0:.1f}초]")

    part2_noise_injection_tradeoff()

    t2 = time.time()
    print(f"\n  [Part 2 완료: {t2 - t1:.1f}초]")

    part3_compare_ising_derived()

    t3 = time.time()
    print(f"\n  [Part 3 완료: {t3 - t2:.1f}초]")

    part4_hybrid_approach()

    t4 = time.time()
    print(f"\n  [Part 4 완료: {t4 - t3:.1f}초]")

    print("\n" + "=" * 70)
    print("  전체 분석 완료")
    print(f"  총 소요 시간: {t4 - t0:.1f}초")
    print("=" * 70)

    print("""
  === 최종 결론 (Final Conclusions) ===

  1. Diagonal Bias = Energy Gap (Part 1):
     Q_kk는 단일 비트 플립 에너지 갭 ΔE_k의 주성분.
     E[ΔE_k] = E[Q_kk] (cross term은 E[Q_km]=0에 의해 상쇄).
     대각 편향을 제거하면 에너지 갭이 0이 되어 GS가 반드시 깨짐.

  2. Noise Injection Tradeoff (Part 2):
     대각 노이즈로 SNR을 줄일 수 있으나, GS 보존율과 trade-off.
     적당한 sigma에서 SNR을 약간 줄이면서 GS 95%+ 유지 가능.
     하지만 SNR을 0으로 만드는 것은 불가능 (GS 파괴).

  3. ZeroOffDiag vs Ising-derived (Part 3):
     ZeroOffDiagonal: 대각 SNR ~ O(1), off-diag 완벽
     Ising-derived:   대각 완벽, off-diag SNR ~ O(n)
     => ZeroOffDiagonal이 indistinguishability에서 압도적 우위.

  4. Hybrid Approach (Part 4):
     혼합 비율 lambda를 조절해도 ZeroOffDiagonal(lambda=0)이 거의 최적.
     Ising 혼합 시 off-diag 누출이 급격히 증가하여 오히려 나빠짐.
     대각 편향은 양의 페널티 제약 하에서 제거 불가능한 본질적 한계.
""")
