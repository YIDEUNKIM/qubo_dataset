#!/usr/bin/env python3
"""
Signed Posiform QUBO 생성기

기존 posiform은 wrong tuple indicator에 양수 weight만 사용.
본 방법론은 음수 weight도 허용하되, GS 보장 범위 내로 제한.

핵심 아이디어:
  wrong tuple (x_i, x_j) = (w_i, w_j)에 weight w를 곱해서 f에 더할 때,
  - w > 0: 항상 안전 (기존 posiform)
  - w < 0: |w| < X - Y 범위에서만 안전
    X = 해당 wrong tuple 부분공간에서 f의 최솟값
    Y = global 최솟값 (target 에너지)

2단계 구성:
  Phase 1: 절의 일부를 양수 weight로 처리 → 기본 에너지 지형 수립
  Phase 2: 나머지를 음수 weight로 처리 → gap 기반 범위 내

장점:
  - Q 행렬에 양/음 혼합 → random QUBO와 유사 (planted solution 감지 어려움)
  - Frustration 도입 → SA에 더 어려운 문제
  - GS 수학적 보장 (전수검사 불필요, 한계치 준수로 보장)
"""

import random
import itertools
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from qubo_utils import calculate_energy, save_qubo_edgelist
from posiform.qubo_posiform import create_planted_2sat


def _clause_to_wrong_tuple(clause):
    """clause → (i, j, wrong_i, wrong_j) 변환. i < j 보장."""
    (v1, p1), (v2, p2) = clause
    wrong_1 = 0 if p1 else 1
    wrong_2 = 0 if p2 else 1
    if v1 <= v2:
        return v1, v2, wrong_1, wrong_2
    else:
        return v2, v1, wrong_2, wrong_1


def _add_posiform_term(Q, i, j, wi, wj, b):
    """
    wrong tuple (wi, wj) indicator에 weight b를 곱하여 Q에 추가.

    indicator 전개:
      (0,0): b*(1-x_i)*(1-x_j) = b - b*x_i - b*x_j + b*x_i*x_j
      (0,1): b*(1-x_i)*x_j     = b*x_j - b*x_i*x_j
      (1,0): b*x_i*(1-x_j)     = b*x_i - b*x_i*x_j
      (1,1): b*x_i*x_j

    Returns:
        constant: 상수항 변화량
    """
    constant = 0.0
    key_ii = (i, i)
    key_jj = (j, j)
    key_ij = (i, j)

    if wi == 0 and wj == 0:
        constant += b
        Q[key_ii] = Q.get(key_ii, 0) - b
        Q[key_jj] = Q.get(key_jj, 0) - b
        Q[key_ij] = Q.get(key_ij, 0) + b
    elif wi == 0 and wj == 1:
        Q[key_jj] = Q.get(key_jj, 0) + b
        Q[key_ij] = Q.get(key_ij, 0) - b
    elif wi == 1 and wj == 0:
        Q[key_ii] = Q.get(key_ii, 0) + b
        Q[key_ij] = Q.get(key_ij, 0) - b
    else:  # (1, 1)
        Q[key_ij] = Q.get(key_ij, 0) + b

    return constant


def _compute_energies_chunked(Q, n, constant, chunk_bits=22):
    """
    Q로부터 모든 2^n 상태의 에너지를 chunk 기반으로 계산.
    메모리 효율: chunk 단위로 bit 배열 생성 → 전체 n개 bit 배열을 동시에 올리지 않음.
    """
    size = 1 << n
    energies = np.full(size, constant, dtype=np.float64)
    chunk_size = min(size, 1 << chunk_bits)
    q_items = list(Q.items())

    for start in range(0, size, chunk_size):
        end = min(start + chunk_size, size)
        masks_chunk = np.arange(start, end, dtype=np.int32)
        # chunk 내 모든 bit 배열 계산 (chunk_size × n bytes)
        bits_chunk = [((masks_chunk >> bit) & 1).astype(np.int8) for bit in range(n)]

        for (i, j), w in q_items:
            if i == j:
                energies[start:end] += w * bits_chunk[i]
            else:
                energies[start:end] += w * bits_chunk[i] * bits_chunk[j]

        del bits_chunk, masks_chunk

    return energies


def _posiform_energy_delta_chunked(n, i, j, wi, wj, b, energies, chunk_bits=22):
    """posiform term 추가 시 에너지 증분 업데이트. chunk 기반."""
    size = 1 << n
    chunk_size = min(size, 1 << chunk_bits)

    for start in range(0, size, chunk_size):
        end = min(start + chunk_size, size)
        masks_chunk = np.arange(start, end, dtype=np.int32)
        xi = ((masks_chunk >> i) & 1).astype(np.int8)
        xj = ((masks_chunk >> j) & 1).astype(np.int8)

        if wi == 0 and wj == 0:
            delta = b * (1 - xi) * (1 - xj)
        elif wi == 0 and wj == 1:
            delta = b * (1 - xi) * xj
        elif wi == 1 and wj == 0:
            delta = b * xi * (1 - xj)
        else:
            delta = b * xi * xj

        energies[start:end] += delta


def _compute_gap_chunked(energies, n, target_mask, var_i, var_j, wrong_i, wrong_j,
                          chunk_bits=22):
    """
    wrong tuple 부분공간 최솟값과 target 에너지 차이(gap) 계산. chunk 기반.
    """
    target_energy = energies[target_mask]
    min_in_subspace = np.float64('inf')
    size = 1 << n
    chunk_size = min(size, 1 << chunk_bits)

    for start in range(0, size, chunk_size):
        end = min(start + chunk_size, size)
        masks_chunk = np.arange(start, end, dtype=np.int32)
        xi = (masks_chunk >> var_i) & 1
        xj = (masks_chunk >> var_j) & 1
        mask = (xi == wrong_i) & (xj == wrong_j)
        sub = energies[start:end][mask]
        if len(sub) > 0:
            chunk_min = np.min(sub)
            if chunk_min < min_in_subspace:
                min_in_subspace = chunk_min

    return min_in_subspace - target_energy


# ── n이 작을 때 빠른 전체 배열 버전 (n ≤ 23) ──

def _precompute_bits(n):
    """각 변수의 비트값 배열 사전 계산."""
    size = 1 << n
    masks = np.arange(size, dtype=np.int32)
    return [((masks >> bit) & 1).astype(np.int8) for bit in range(n)]


def _compute_energies_full(Q, n, constant, bits):
    """Q로부터 모든 2^n 상태의 에너지를 전체 배열로 계산."""
    size = 1 << n
    energies = np.full(size, constant, dtype=np.float64)
    for (i, j), w in Q.items():
        if i == j:
            energies += w * bits[i]
        else:
            energies += w * bits[i] * bits[j]
    return energies


def _posiform_energy_delta_full(n, i, j, wi, wj, b, energies, bits):
    """posiform term 추가 시 에너지 증분 업데이트. 전체 배열 버전."""
    xi = bits[i]
    xj = bits[j]
    if wi == 0 and wj == 0:
        energies += b * (1 - xi) * (1 - xj)
    elif wi == 0 and wj == 1:
        energies += b * (1 - xi) * xj
    elif wi == 1 and wj == 0:
        energies += b * xi * (1 - xj)
    else:
        energies += b * xi * xj


def _compute_gap_full(energies, target_mask, var_i, var_j, wrong_i, wrong_j, bits):
    """부분공간 gap 계산. 전체 배열 버전."""
    target_energy = energies[target_mask]
    subspace_mask = (bits[var_i] == wrong_i) & (bits[var_j] == wrong_j)
    min_in_subspace = np.min(energies[subspace_mask])
    return min_in_subspace - target_energy


def create_qubo_signed_posiform(target, coeff_range=(1.0, 3.0),
                                 negative_ratio=0.3, margin=0.01,
                                 max_clauses_factor=10,
                                 seed=None, verbose=True):
    """
    Signed Posiform QUBO 생성.

    Phase 1: 절의 (1 - negative_ratio) 비율을 양수 weight로 처리
    Phase 2: 나머지를 음수 weight로 처리 (gap 기반 범위)

    Args:
        target: 목표 비트스트링 (예: "10110")
        coeff_range: 양수 weight 범위 (lo, hi)
        negative_ratio: Phase 2 (음수) 비율 (0.0 ~ 1.0)
        margin: GS 보장 마진 (gap의 margin 비율만큼 여유)
        max_clauses_factor: 최대 clause 수 = factor * n
        seed: 난수 시드
        verbose: 출력 여부

    Returns:
        Q: QUBO 딕셔너리 {(i,j): weight}
        info: 메타정보
    """
    n = len(target)
    if n < 2:
        raise ValueError(f"변수 수가 2 미만입니다: n={n}")

    rng = random.Random(seed)
    target_bits = [int(b) for b in target]
    target_mask = sum(b << i for i, b in enumerate(target_bits))

    def log(msg):
        if verbose:
            print(msg)

    log(f"\n{'='*60}")
    log(f" Signed Posiform QUBO Generator")
    log(f"{'='*60}")
    log(f" Target: {target} (n={n})")
    log(f" 양수 범위: {coeff_range}, 음수 비율: {negative_ratio:.0%}")

    # ── 2-SAT 생성 ──
    sat_seed = seed if seed is not None else None
    clauses, is_unique = create_planted_2sat(
        target, max_clauses_factor=max_clauses_factor, seed=sat_seed)

    log(f"\n[2-SAT] 절 수: {len(clauses)}, 유일 해: {'OK' if is_unique else 'FAIL'}")

    # ── 절 분할: Phase 1 (양수) / Phase 2 (음수) ──
    indices = list(range(len(clauses)))
    rng.shuffle(indices)

    n_negative = max(0, int(len(clauses) * negative_ratio))
    n_positive = len(clauses) - n_negative

    positive_indices = indices[:n_positive]
    negative_indices = indices[n_positive:]

    log(f"[분할] Phase 1 (양수): {n_positive}절, Phase 2 (음수): {n_negative}절")

    # ── Phase 1: 양수 weight ──
    Q = {}
    constant = 0.0
    coeff_seed = (seed + 1000) if seed is not None else None
    rng_coeff = random.Random(coeff_seed)

    positive_weights = []
    for idx in positive_indices:
        clause = clauses[idx]
        i, j, wi, wj = _clause_to_wrong_tuple(clause)
        lo, hi = coeff_range
        b = rng_coeff.uniform(lo, hi)
        constant += _add_posiform_term(Q, i, j, wi, wj, b)
        positive_weights.append(b)

    log(f"\n[Phase 1] 양수 weight: avg={np.mean(positive_weights):.3f}")

    # ── Phase 2: 음수 weight (gap 기반) ──
    negative_weights = []
    skipped = 0
    use_chunked = n > 23  # n>23: chunk 기반 (메모리 절약), n≤23: 전체 배열 (빠름)

    if n_negative > 0:
        if use_chunked:
            # chunk 기반: 메모리 O(2^n * 8) + O(2^chunk_bits * n)
            energies = _compute_energies_chunked(Q, n, constant)

            for idx in negative_indices:
                clause = clauses[idx]
                i, j, wi, wj = _clause_to_wrong_tuple(clause)

                gap = _compute_gap_chunked(energies, n, target_mask,
                                           i, j, wi, wj)

                if gap <= margin:
                    lo, hi = coeff_range
                    b = rng_coeff.uniform(lo, hi)
                    constant += _add_posiform_term(Q, i, j, wi, wj, b)
                    positive_weights.append(b)
                    skipped += 1
                else:
                    max_neg = gap * (1.0 - margin)
                    b = -rng_coeff.uniform(margin, max_neg)
                    constant += _add_posiform_term(Q, i, j, wi, wj, b)
                    negative_weights.append(b)

                _posiform_energy_delta_chunked(n, i, j, wi, wj, b, energies)
        else:
            # 전체 배열: 빠르지만 메모리 O(2^n * n)
            bits = _precompute_bits(n)
            energies = _compute_energies_full(Q, n, constant, bits)

            for idx in negative_indices:
                clause = clauses[idx]
                i, j, wi, wj = _clause_to_wrong_tuple(clause)

                gap = _compute_gap_full(energies, target_mask,
                                        i, j, wi, wj, bits)

                if gap <= margin:
                    lo, hi = coeff_range
                    b = rng_coeff.uniform(lo, hi)
                    constant += _add_posiform_term(Q, i, j, wi, wj, b)
                    positive_weights.append(b)
                    skipped += 1
                else:
                    max_neg = gap * (1.0 - margin)
                    b = -rng_coeff.uniform(margin, max_neg)
                    constant += _add_posiform_term(Q, i, j, wi, wj, b)
                    negative_weights.append(b)

                _posiform_energy_delta_full(n, i, j, wi, wj, b, energies, bits)

    # 0에 가까운 항 제거
    Q = {k: v for k, v in Q.items() if abs(v) > 1e-15}

    if negative_weights:
        log(f"[Phase 2] 음수 weight: {len(negative_weights)}개, "
            f"avg={np.mean(negative_weights):.3f}, "
            f"min={np.min(negative_weights):.3f}")
    if skipped > 0:
        log(f"  gap 부족으로 양수 전환: {skipped}개")

    # ── 검증 ──
    if use_chunked:
        final_energies = _compute_energies_chunked(Q, n, constant)
    else:
        if 'bits' not in dir():
            bits = _precompute_bits(n)
        final_energies = _compute_energies_full(Q, n, constant, bits)
    target_energy = final_energies[target_mask]
    final_energies[target_mask] = float('inf')
    min_other = np.min(final_energies)
    energy_gap = min_other - target_energy
    is_gs = energy_gap > -1e-10

    log(f"\n[검증]")
    log(f"  Target 에너지: {target_energy:.6f}")
    log(f"  에너지 갭: {energy_gap:.6f}")
    log(f"  GS 보장: {'OK' if is_gs else 'FAIL'}")

    # Q 계수 부호 분석
    pos_count = sum(1 for v in Q.values() if v > 0)
    neg_count = sum(1 for v in Q.values() if v < 0)
    log(f"\n[Q 행렬]")
    log(f"  비영 항: {len(Q)}개 (양수: {pos_count}, 음수: {neg_count})")
    log(f"  양/음 비율: {neg_count}/{pos_count + neg_count} = "
        f"{neg_count/(pos_count+neg_count)*100:.1f}%")

    info = {
        'n': n,
        'num_clauses': len(clauses),
        'is_unique': is_unique,
        'n_positive_clauses': n_positive + skipped,
        'n_negative_clauses': len(negative_weights),
        'positive_weights': positive_weights,
        'negative_weights': negative_weights,
        'constant_offset': constant,
        'target_energy': target_energy,
        'energy_gap': energy_gap,
        'is_gs': is_gs,
        'q_pos_count': pos_count,
        'q_neg_count': neg_count,
    }

    return Q, info


if __name__ == "__main__":
    target = "10110"
    seed = None

    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if set(arg) <= {'0', '1'} and len(arg) >= 2:
            target = arg
        elif arg.isdigit():
            length = int(arg)
            if length < 2:
                length = 2
            random.seed(42)
            target = ''.join(str(random.randint(0, 1)) for _ in range(length))
            print(f"[설정] 길이 {length}의 랜덤 목표 해 생성: {target}")

    if len(sys.argv) > 2:
        seed = int(sys.argv[2])

    Q, info = create_qubo_signed_posiform(
        target, coeff_range=(1.0, 3.0), negative_ratio=0.3,
        seed=seed, verbose=True)

    # 결과 저장
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(output_dir, exist_ok=True)
    n = len(target)
    filename = f"signed_posiform_{target[:5]}_{n}.txt"
    output_file = os.path.join(output_dir, filename)
    save_qubo_edgelist(Q, output_file, target)
    print(f"\n[저장 완료] {output_file}")
