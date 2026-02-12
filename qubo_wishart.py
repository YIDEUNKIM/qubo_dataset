"""
Wishart Planted Ensemble QUBO 생성기

Hamze et al. (2020) 기반: alpha = M/N 파라미터로 난이도 제어.
1차 상전이 근처 (alpha ~ 0.5-0.8)에서 SA가 실패하는 어려운 QUBO 생성.

핵심 아이디어:
  1. 목표 비트스트링 t → spin 벡터 (+1/-1)
  2. t에 직교하는 M개의 Gaussian 벡터 W 생성 (M = alpha * N)
  3. J = -(1/N) W W^T (대각 제거)
  4. W^T t = 0이므로 t가 ground state임이 수학적으로 보장됨
"""

import random
import itertools
import sys
import os
import numpy as np

from qubo_zero_expectation import calculate_energy, save_qubo_edgelist, print_q_matrix, print_qubo_formula


def create_wishart_ising(n, alpha, target, seed=None):
    """
    Wishart Planted Ensemble Ising 모델 생성.

    Args:
        n: 큐빗 수
        alpha: M/N 비율 (난이도 파라미터, 0.5~0.8이 가장 어려움)
        target: 목표 비트스트링 (예: "10110")
        seed: 난수 시드 (재현성)

    Returns:
        J_dict: Ising coupling {(i,j): J_ij} (i < j)
        h_dict: 외부 자기장 {i: h_i} (Wishart에서는 모두 0)
    """
    rng = np.random.default_rng(seed)
    M = int(alpha * n)

    # spin 벡터: +1/-1
    t = np.array([1.0 if b == '1' else -1.0 for b in target])
    t_norm_sq = np.dot(t, t)  # = n

    # M개의 N차원 Gaussian 벡터 생성
    W = rng.standard_normal((n, M))

    # 각 열을 t에 직교 투영: w_mu → w_mu - (w_mu · t / ||t||^2) * t
    for mu in range(M):
        proj = np.dot(W[:, mu], t) / t_norm_sq
        W[:, mu] -= proj * t

    # J = -(1/N) W W^T
    J_matrix = -(1.0 / n) * (W @ W.T)

    # 대각 제거 + 상삼각만 추출
    J_dict = {}
    h_dict = {i: 0.0 for i in range(n)}

    for i in range(n):
        for j in range(i + 1, n):
            J_val = J_matrix[i, j]
            if abs(J_val) > 1e-15:
                J_dict[(i, j)] = J_val

    return J_dict, h_dict


def ising_to_qubo(J_dict, h_dict, n):
    """
    Ising → QUBO 변환.

    H_Ising = -sum_{i<j} J_ij s_i s_j - sum_i h_i s_i
    s_i = 2x_i - 1

    Returns:
        Q: QUBO 딕셔너리 {(i,j): weight}
    """
    Q = {}

    # J coupling 변환
    for (i, j), J_ij in J_dict.items():
        # -J_ij s_i s_j = -J_ij (4x_ix_j - 2x_i - 2x_j + 1)
        #               = -4J_ij x_ix_j + 2J_ij x_i + 2J_ij x_j - J_ij
        Q[(i, j)] = Q.get((i, j), 0) + (-4 * J_ij)
        Q[(i, i)] = Q.get((i, i), 0) + (2 * J_ij)
        Q[(j, j)] = Q.get((j, j), 0) + (2 * J_ij)

    # h 외부 자기장 변환
    for i, h_i in h_dict.items():
        if abs(h_i) > 1e-15:
            # -h_i s_i = -h_i (2x_i - 1) = -2h_i x_i + h_i
            Q[(i, i)] = Q.get((i, i), 0) + (-2 * h_i)

    return Q


def create_qubo_wishart(target, alpha=0.7, seed=None):
    """
    Wishart Planted Ensemble QUBO 생성 (메인 진입점).

    Args:
        target: 목표 비트스트링 (예: "10110")
        alpha: M/N 비율 (0.5~0.8 = 어려움, 0.3 = 쉬움, 1.5 = 쉬움)
        seed: 난수 시드

    Returns:
        Q: QUBO 딕셔너리 {(i,j): weight}
    """
    n = len(target)
    J_dict, h_dict = create_wishart_ising(n, alpha, target, seed=seed)
    Q = ising_to_qubo(J_dict, h_dict, n)
    return Q


def verify_ground_state(Q, target, num_random_samples=10000):
    """
    target이 ground state인지 검증.

    1. Single-flip 이웃 검사 (완전): 모든 비트 하나 뒤집어서 에너지 증가 확인
    2. 랜덤 샘플 검사: 무작위 상태들과 비교

    Returns:
        (is_local_min, is_likely_global, stats_dict)
    """
    n = len(target)
    target_energy = calculate_energy(target, Q)

    # 1. Single-flip 이웃 검사
    min_flip_delta = float('inf')
    for i in range(n):
        flipped = list(target)
        flipped[i] = '0' if flipped[i] == '1' else '1'
        flipped_str = ''.join(flipped)
        flipped_energy = calculate_energy(flipped_str, Q)
        delta = flipped_energy - target_energy
        min_flip_delta = min(min_flip_delta, delta)

    is_local_min = min_flip_delta > -1e-10

    # 2. 랜덤 샘플 검사
    lower_count = 0
    min_random_energy = float('inf')

    for _ in range(num_random_samples):
        rand_state = ''.join(str(random.randint(0, 1)) for _ in range(n))
        rand_energy = calculate_energy(rand_state, Q)
        if rand_energy < target_energy - 1e-10:
            lower_count += 1
        min_random_energy = min(min_random_energy, rand_energy)

    is_likely_global = (lower_count == 0)

    stats = {
        'target_energy': target_energy,
        'min_flip_delta': min_flip_delta,
        'is_local_min': is_local_min,
        'lower_count': lower_count,
        'num_random_samples': num_random_samples,
        'min_random_energy': min_random_energy,
        'is_likely_global': is_likely_global,
    }

    return is_local_min, is_likely_global, stats


def verify_brute_force(Q, target, n):
    """n <= 20일 때 brute force로 ground state 완전 검증."""
    if n > 20:
        return None

    target_energy = calculate_energy(target, Q)
    best_energy = float('inf')
    best_state = None

    for bits in itertools.product([0, 1], repeat=n):
        state = ''.join(map(str, bits))
        energy = calculate_energy(state, Q)
        if energy < best_energy:
            best_energy = energy
            best_state = state

    return {
        'target_energy': target_energy,
        'best_energy': best_energy,
        'best_state': best_state,
        'is_ground_state': abs(best_energy - target_energy) < 1e-10,
        'energy_gap': best_energy - target_energy if best_state != target else 0.0,
    }


if __name__ == "__main__":
    target = "10110"
    alpha = 0.7
    seed = None

    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if set(arg) <= {'0', '1'}:
            # 0과 1로만 구성된 문자열 → 이진 타겟으로 해석
            target = arg
        elif arg.isdigit():
            # 순수 숫자 (2~9 포함) → N 길이로 해석
            length = int(arg)
            random.seed(42)
            target = ''.join(str(random.randint(0, 1)) for _ in range(length))
            print(f"[설정] 길이 {length}의 랜덤 목표 해 생성")
        else:
            target = arg

    if len(sys.argv) > 2:
        alpha = float(sys.argv[2])

    if len(sys.argv) > 3:
        seed = int(sys.argv[3])

    n = len(target)
    print("=" * 60)
    print("Wishart Planted Ensemble QUBO 생성기")
    print("=" * 60)
    print(f"Target: {target} (N={n})")
    print(f"Alpha: {alpha} (M={int(alpha * n)})")
    if seed is not None:
        print(f"Seed: {seed}")

    # QUBO 생성
    Q = create_qubo_wishart(target, alpha=alpha, seed=seed)

    # 행렬 출력 (작을 때만)
    if n <= 12:
        print_q_matrix(Q, n)
        print_qubo_formula(Q)

    # 검증
    if n <= 20:
        print(f"\n[Brute Force 검증] (N={n})")
        bf_result = verify_brute_force(Q, target, n)
        print(f"  Target 에너지: {bf_result['target_energy']:.6f}")
        print(f"  최소 에너지:   {bf_result['best_energy']:.6f}")
        print(f"  최소 상태:     {bf_result['best_state']}")
        if bf_result['is_ground_state']:
            print("  ✓ Ground State 검증 성공!")
        else:
            print(f"  ✗ Ground State 검증 실패 (gap: {bf_result['energy_gap']:.6f})")
    else:
        print(f"\n[통계적 검증] (N={n})")
        is_local, is_global, stats = verify_ground_state(Q, target)
        print(f"  Target 에너지: {stats['target_energy']:.6f}")
        print(f"  최소 flip delta: {stats['min_flip_delta']:.6f}")
        print(f"  Local minimum: {'✓' if is_local else '✗'}")
        print(f"  랜덤 {stats['num_random_samples']}개 중 더 낮은 에너지: {stats['lower_count']}개")
        print(f"  Global minimum (추정): {'✓' if is_global else '✗'}")

    # 결과 저장
    output_dir = "qubo_results"
    os.makedirs(output_dir, exist_ok=True)
    filename_target = f"{target[:5]}_{n}"
    output_file = os.path.join(output_dir, f"wishart_{filename_target}_a{alpha}.txt")
    save_qubo_edgelist(Q, output_file, target)
    print(f"\n[저장 완료] {output_file}")
