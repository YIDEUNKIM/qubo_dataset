#!/usr/bin/env python3
"""
Truth Table Concat QUBO Generator
====================================

k-bit Truth Table QUBO를 h개 생성하여 block-diagonal로 접합.
총 변수 수 = k*h (approx/random) 또는 k*h + aux (exact).

Landscape 모드:
  planted (기본): 블록별 target을 심은 진리표 -> QUBO. SA-easy bowl landscape.
  random: 블록별 무작위 진리표 -> unconstrained QP fit -> QUBO.
    Genuine local minima. Target은 brute force로 발견. Posiform hardening 필수.

Posiform Hardening (선택 / random일 때 필수):
  Q_final = Q_concat + alpha * Q_posiform
  posiform이 cross-block coupling을 추가하여 block-diagonal 분해 불가능.
  Ground state 유일성 보장.

파이프라인 (planted):
  1. 블록별로 preset_random_landscape -> 진리표 생성
  2. 각 블록을 approx 또는 exact 모드로 QUBO 변환
  3. 변수 인덱스를 shift하여 block-diagonal Q 구성
  4. (선택) posiform hardening: Q_final = Q_concat + alpha * Q_posiform
  5. 전체 target 검증

파이프라인 (random):
  1. 블록별로 random truth table -> unconstrained QP fit -> QUBO
  2. 각 블록의 ground state를 brute force로 탐색
  3. Block-diagonal 접합: Q_base
  4. Posiform hardening: Q_final = Q_base + alpha * Q_posiform
  5. 전체 target 검증

사용법:
    python truthtable_concat/qubo_truthtable_concat.py 1001111 10
    python truthtable_concat/qubo_truthtable_concat.py 1001111 10 --exact
    python truthtable_concat/qubo_truthtable_concat.py 1001111 10 --harden 0.1
    python truthtable_concat/qubo_truthtable_concat.py 1001111 10 --random --harden 0.01
    python truthtable_concat/qubo_truthtable_concat.py 1001111 10 --random --harden 0.01 --seed 42
"""

import sys
import os
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from etc.truthtable.qubo_truthtable import (
    create_qubo_truthtable, create_qubo_approx,
    preset_random_landscape, compute_aux_values,
)
from etc.posiform.qubo_posiform import create_qubo_posiform
from qubo_utils import calculate_energy, save_qubo_edgelist


def create_random_block_qubo(k, energy_range=(-5.0, 5.0), seed=None):
    """
    Random truth table -> unconstrained QP fit (with constant) -> block QUBO

    진리표 에너지를 전부 랜덤으로 생성 (target 심지 않음).
    상수항 포함 QP fit -> 상수 제거 -> QUBO. Brute force로 실제 GS 탐색.

    Args:
        k: 블록 변수 수
        energy_range: (lo, hi) uniform 분포 범위 (0 중심 권장)
        seed: 랜덤 시드

    Returns:
        Q: dict {(i,j): weight}
        gs_bitstring: QUBO의 ground state
        gs_energy: ground state 에너지
        block_info: landscape 메트릭
    """
    rng = np.random.RandomState(seed)
    n_states = 1 << k

    lo, hi = energy_range
    energies = rng.uniform(lo, hi, n_states)

    # Design matrix: 상수항 + QUBO 파라미터
    n_qubo_params = k + k * (k - 1) // 2
    A = np.zeros((n_states, 1 + n_qubo_params))

    for x_val in range(n_states):
        x_bits = [(x_val >> (k - 1 - j)) & 1 for j in range(k)]
        A[x_val, 0] = 1.0  # constant term
        col = 1
        for i in range(k):
            A[x_val, col] = x_bits[i]
            col += 1
        for i in range(k):
            for j in range(i + 1, k):
                A[x_val, col] = x_bits[i] * x_bits[j]
                col += 1

    # Unconstrained least squares (상수항 포함)
    q_full, _, _, _ = np.linalg.lstsq(A, energies, rcond=None)
    q_params = q_full[1:]

    # Build Q dict
    Q = {}
    col = 0
    for i in range(k):
        v = float(q_params[col])
        if abs(v) > 1e-10:
            Q[(i, i)] = v
        col += 1
    for i in range(k):
        for j in range(i + 1, k):
            v = float(q_params[col])
            if abs(v) > 1e-10:
                Q[(i, j)] = v
            col += 1

    # 모든 2^k 상태 평가
    qubo_energies = np.zeros(n_states)
    for x_val in range(n_states):
        x_str = format(x_val, f'0{k}b')
        qubo_energies[x_val] = calculate_energy(x_str, Q)

    # Ground state
    gs_idx = int(np.argmin(qubo_energies))
    gs_energy = float(qubo_energies[gs_idx])
    gs_bitstring = format(gs_idx, f'0{k}b')

    # Energy gap
    sorted_e = np.sort(qubo_energies)
    energy_gap = float(sorted_e[1] - sorted_e[0]) if n_states > 1 else 0.0

    # Local minima 수 (1-flip neighborhood)
    n_local_minima = 0
    for x_val in range(n_states):
        e = qubo_energies[x_val]
        is_local_min = True
        for j in range(k):
            neighbor = x_val ^ (1 << (k - 1 - j))
            if qubo_energies[neighbor] < e:
                is_local_min = False
                break
        if is_local_min:
            n_local_minima += 1

    # GS 유일성
    gs_count = int(np.sum(np.abs(qubo_energies - gs_energy) < 1e-10))

    # Fit quality
    fitted = A @ q_full
    rmse = float(np.sqrt(np.mean((fitted - energies) ** 2)))

    block_info = {
        'k': k,
        'gs_bitstring': gs_bitstring,
        'gs_energy': gs_energy,
        'energy_gap': energy_gap,
        'n_local_minima': n_local_minima,
        'gs_unique': gs_count == 1,
        'rmse': rmse,
        'n_params': n_qubo_params,
        'energy_range_actual': (float(np.min(qubo_energies)),
                                float(np.max(qubo_energies))),
    }

    return Q, gs_bitstring, gs_energy, block_info


def create_qubo_concat(target, h,
                       mode='approx', epsilon=0.01,
                       reduce_strategy='greedy',
                       posiform_scale=None,
                       posiform_coeff_range=(1.0, 1.0),
                       landscape='planted',
                       energy_range=(-5.0, 5.0),
                       seed=None, verbose=True):
    """
    k-bit Truth Table QUBO를 h개 block-diagonal로 접합 (+ 선택적 posiform hardening)

    Args:
        target: k-bit 문자열 (블록 크기 = len(target))
        h: 블록 수
        mode: 'approx' (보조변수 0, 기본값) / 'exact' (Rosenberg)
            landscape='random'일 때는 무시됨
        epsilon: 근사 모드 에너지 갭 하한
        reduce_strategy: exact 모드 차수축소 전략 ('original'/'cache'/'greedy')
        posiform_scale: alpha -- posiform 스케일링 계수 (None이면 hardening 안함)
        posiform_coeff_range: posiform clause 계수 범위
        landscape: 'planted' (기존) / 'random' (무작위 진리표, posiform 필수)
        energy_range: random landscape 에너지 범위 (planted일 때 무시)
        seed: 랜덤 시드 (블록 i는 seed+i 사용)
        verbose: 출력 여부

    Returns:
        Q: dict {(i,j): weight}
        info: dict (메타데이터)
    """
    k = len(target)

    if landscape == 'random':
        assert posiform_scale is not None, \
            "landscape='random'은 posiform hardening 필수 (posiform_scale 지정)"

    def log(msg):
        if verbose:
            print(msg)

    # Block target 결정
    if landscape == 'random':
        block_targets = []  # discovered in the loop
    elif posiform_scale is not None:
        target_rng = np.random.RandomState(
            (seed + 7777) if seed is not None else None)
        block_targets = []
        for _ in range(h):
            bt = ''.join(str(b) for b in target_rng.randint(0, 2, k))
            block_targets.append(bt)
    else:
        block_targets = [target] * h

    log(f"\n{'='*60}")
    log(f" Truth Table Concat QUBO Generator")
    log(f"{'='*60}")
    if landscape == 'random':
        log(f" Landscape: random (k={k})")
    elif posiform_scale is not None:
        log(f" 블록 target: 랜덤 (k={k})")
    else:
        log(f" 블록 target: {target} (k={k})")
    log(f" 블록 수: h={h}")
    if landscape != 'random':
        log(f" 모드: {mode}")
    if posiform_scale is not None:
        log(f" Posiform hardening: α={posiform_scale}")
    if landscape == 'random':
        log(f" Energy range: {energy_range}")
    else:
        log(f" 프리셋: random landscape (uniform 0.1~5.0)")
    if seed is not None:
        log(f" seed={seed}")

    Q_total = {}
    blocks_info = []
    all_aux_info = []
    total_local_minima = 0
    offset = 0

    for i in range(h):
        seed_i = (seed + i) if seed is not None else None

        if landscape == 'random':
            # Random landscape: unconstrained QP fit, GS by brute force
            Q_block, gs, gs_e, binfo = create_random_block_qubo(
                k, energy_range=energy_range, seed=seed_i)
            block_targets.append(gs)

            for (a, b), w in Q_block.items():
                Q_total[(a + offset, b + offset)] = w

            total_local_minima += binfo['n_local_minima']

            blocks_info.append({
                'block_idx': i, 'offset': offset,
                'n_total': k, 'n_aux': 0,
                'qubo_energy': gs_e, 'seed': seed_i,
                'target': gs,
                'energy_gap': binfo['energy_gap'],
                'n_local_minima': binfo['n_local_minima'],
                'rmse': binfo['rmse'],
            })

            log(f"  블록 {i}: target={gs}, E={gs_e:.4f}, "
                f"gap={binfo['energy_gap']:.4f}, LM={binfo['n_local_minima']}")

            offset += k

        else:
            # Planted landscape
            tt = preset_random_landscape(k, block_targets[i], seed=seed_i)

            if mode == 'approx':
                Q_block, info_block = create_qubo_approx(
                    tt, verbose=False, epsilon=epsilon)
            else:
                Q_block, info_block = create_qubo_truthtable(
                    tt, verbose=False, reduce_strategy=reduce_strategy)

            block_total_vars = info_block['n_total']

            for (a, b), w in Q_block.items():
                Q_total[(a + offset, b + offset)] = w

            for (y, a, b) in info_block['aux_info']:
                all_aux_info.append((y + offset, a + offset, b + offset))

            block_aux = info_block['aux_info']
            x_orig = [int(c) for c in block_targets[i]]
            x_block = compute_aux_values(x_orig, block_aux)
            block_target_str = ''.join(
                str(x_block[j]) for j in range(block_total_vars))
            block_qubo_energy = calculate_energy(block_target_str, Q_block)

            blocks_info.append({
                'block_idx': i, 'offset': offset,
                'n_total': block_total_vars, 'n_aux': info_block['n_aux'],
                'qubo_energy': block_qubo_energy, 'seed': seed_i,
                'target': block_targets[i],
            })

            log(f"  블록 {i}: target={block_targets[i]}, offset={offset}, "
                f"vars={block_total_vars}, aux={info_block['n_aux']}, "
                f"E_qubo={block_qubo_energy:.4f}")

            offset += block_total_vars

    n_total = offset
    n_original = k * h
    full_target = ''.join(block_targets)

    # exact 모드: 전체 target에 aux 비트 추가
    if landscape != 'random' and mode == 'exact' and all_aux_info:
        x_full = {}
        for bi_idx, bi in enumerate(blocks_info):
            block_n_total = bi['n_total']
            block_n_aux = bi['n_aux']
            block_k = block_n_total - block_n_aux
            bt = block_targets[bi_idx]
            for j in range(block_k):
                x_full[bi['offset'] + j] = int(bt[j])

        for (y, a, b) in all_aux_info:
            x_full[y] = x_full.get(a, 0) * x_full.get(b, 0)

        full_target_with_aux = ''.join(
            str(x_full.get(i, 0)) for i in range(n_total))
    else:
        full_target_with_aux = full_target

    # ─── Posiform Hardening ───
    posiform_info_out = None
    if posiform_scale is not None:
        log(f"\n[Posiform Hardening]")

        # 원래 변수 인덱스 매핑: posiform var j -> Q_total 실제 인덱스
        orig_var_map = []
        for bi in blocks_info:
            block_offset = bi['offset']
            for j in range(k):
                orig_var_map.append(block_offset + j)

        full_target_orig = full_target

        posiform_seed = (seed + 9999) if seed is not None else None
        # random landscape: targeted=True (블록이 target 안 심으므로 강한 유일성 필요)
        # planted landscape: targeted=False (블록이 이미 target 보장)
        posiform_targeted = (landscape == 'random')
        Q_posiform, posiform_info_out = create_qubo_posiform(
            full_target_orig,
            coeff_range=posiform_coeff_range,
            max_clauses_factor=20,
            seed=posiform_seed,
            targeted=posiform_targeted
        )

        # 유일성 실패 시 재시도
        if not posiform_info_out['is_unique']:
            for retry in range(5):
                posiform_seed = (seed + 10000 + retry) \
                    if seed is not None else None
                Q_posiform, posiform_info_out = create_qubo_posiform(
                    full_target_orig,
                    coeff_range=posiform_coeff_range,
                    max_clauses_factor=30,
                    seed=posiform_seed,
                    targeted=posiform_targeted
                )
                if posiform_info_out['is_unique']:
                    break

        if not posiform_info_out['is_unique']:
            warnings.warn(
                f"Posiform uniqueness 보장 실패 (n={len(full_target_orig)}): "
                f"ground state uniqueness가 보장되지 않습니다."
            )

        # 인덱스 리매핑 후 Q_total에 alpha * Q_posiform 추가
        for (a, b), w in Q_posiform.items():
            actual_a = orig_var_map[a]
            actual_b = orig_var_map[b]
            key = (min(actual_a, actual_b), max(actual_a, actual_b))
            Q_total[key] = Q_total.get(key, 0) + posiform_scale * w

        # 0에 가까운 항 제거
        Q_total = {kk: vv for kk, vv in Q_total.items() if abs(vv) > 1e-15}

        log(f"  Posiform unique: "
            f"{'OK' if posiform_info_out['is_unique'] else 'FAIL'}")
        log(f"  Posiform clauses: {posiform_info_out['num_clauses']}")
        posiform_target_e = posiform_info_out['target_energy']
        log(f"  P(x*) = {posiform_target_e:.6f}")
        log(f"  α·P(x*) = {posiform_scale * posiform_target_e:.6f}")

    # 검증: 전체 target으로 에너지 계산
    verify_energy = calculate_energy(full_target_with_aux, Q_total)
    block_energy_sum = sum(bi['qubo_energy'] for bi in blocks_info)

    if posiform_scale is not None and posiform_info_out is not None:
        posiform_contrib = posiform_scale * posiform_info_out['target_energy']
        expected_energy = block_energy_sum + posiform_contrib
    else:
        posiform_contrib = 0.0
        expected_energy = block_energy_sum

    log(f"\n 전체 target: {full_target[:50]}"
        f"{'...' if len(full_target) > 50 else ''} "
        f"(원래 변수 {n_original}개)")
    log(f" 전체 QUBO 크기: {n_total} (aux 포함)")
    log(f" Ground energy: {verify_energy:.4f}")
    log(f" 블록별 합계: {block_energy_sum:.4f}")
    if posiform_scale is not None:
        log(f" α·P(x*): {posiform_contrib:.4f}")
        log(f" 기대 합계: {expected_energy:.4f}")
    if landscape == 'random':
        log(f" 블록 local minima 합계: {total_local_minima}")

    energy_match = abs(verify_energy - expected_energy) < 1e-4
    log(f" 에너지 일치: {'YES' if energy_match else 'NO'}")
    if not energy_match:
        log(f"  *** 경고: 에너지 불일치! "
            f"diff={abs(verify_energy - expected_energy):.6f}")

    info = {
        'n_block': k,
        'h': h,
        'n_original': n_original,
        'n_total': n_total,
        'n_aux': n_total - n_original,
        'block_target': target,
        'block_targets': block_targets,
        'target': full_target,
        'target_with_aux': full_target_with_aux,
        'mode': mode,
        'landscape': landscape,
        'blocks': blocks_info,
        'aux_info': all_aux_info,
        'ground_energy': verify_energy,
        'posiform_scale': posiform_scale,
        'posiform_info': posiform_info_out,
        'total_local_minima': total_local_minima
        if landscape == 'random' else None,
        'energy_range': energy_range if landscape == 'random' else None,
    }

    return Q_total, info


def verify_ground_state(Q, target, n):
    """Brute force GS 검증 (n <= 25)"""
    assert n <= 25, f"n={n} too large for brute force"

    target_energy = calculate_energy(target, Q)
    best_energy = float('inf')
    best_state = None
    second_best = float('inf')

    for x_val in range(1 << n):
        x_str = format(x_val, f'0{n}b')
        e = calculate_energy(x_str, Q)
        if e < best_energy:
            second_best = best_energy
            best_energy = e
            best_state = x_str
        elif e < second_best:
            second_best = e

    return {
        'is_ground_state': abs(target_energy - best_energy) < 1e-8,
        'is_unique': (second_best - best_energy) > 1e-8,
        'target_energy': target_energy,
        'best_energy': best_energy,
        'best_state': best_state,
        'energy_gap': second_best - best_energy,
    }


# ─────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────

def main():
    results_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(results_dir, exist_ok=True)

    if len(sys.argv) < 3 or sys.argv[1] == '--help':
        print("사용법:")
        print("  python truthtable_concat/qubo_truthtable_concat.py "
              "TARGET H [옵션]")
        print()
        print("인자:")
        print("  TARGET  블록 target 비트스트링 (예: 1001111)")
        print("  H       블록 수")
        print()
        print("옵션:")
        print("  --exact          정확 모드 (Rosenberg, 기본: approx)")
        print("  --random         랜덤 landscape (posiform hardening 필수)")
        print("  --harden α       posiform hardening (α=0.1 등, 작을수록 어려움)")
        print("  --seed S         랜덤 시드")
        print("  --strategy S     차수축소 전략 "
              "(original/cache/greedy, 기본: greedy)")
        print()
        print("예시:")
        print("  python truthtable_concat/qubo_truthtable_concat.py "
              "1001111 10")
        print("  python truthtable_concat/qubo_truthtable_concat.py "
              "1001111 10 --exact")
        print("  python truthtable_concat/qubo_truthtable_concat.py "
              "1001111 10 --harden 0.1")
        print("  python truthtable_concat/qubo_truthtable_concat.py "
              "1001111 10 --random --harden 0.01")
        return

    target = sys.argv[1]
    h = int(sys.argv[2])

    use_exact = '--exact' in sys.argv
    use_random = '--random' in sys.argv

    def parse_opt(name, default):
        if name in sys.argv:
            return float(sys.argv[sys.argv.index(name) + 1])
        return default

    seed = int(parse_opt('--seed', -1))
    seed = seed if seed >= 0 else None

    reduce_strategy = 'greedy'
    if '--strategy' in sys.argv:
        idx = sys.argv.index('--strategy')
        reduce_strategy = sys.argv[idx + 1]

    posiform_scale = None
    if '--harden' in sys.argv:
        idx = sys.argv.index('--harden')
        posiform_scale = float(sys.argv[idx + 1])

    landscape = 'random' if use_random else 'planted'
    mode = 'exact' if use_exact else 'approx'

    Q, info = create_qubo_concat(
        target, h,
        mode=mode, reduce_strategy=reduce_strategy,
        posiform_scale=posiform_scale,
        landscape=landscape,
        seed=seed)

    k = info['n_block']
    n_total = info['n_total']

    print(f"\n{'='*60}")
    print(f" 결과 요약")
    print(f"{'='*60}")
    print(f"  블록 target: {target} (k={k})")
    print(f"  블록 수: h={h}")
    print(f"  Landscape: {landscape}")
    if landscape != 'random':
        print(f"  모드: {mode}")
    if posiform_scale is not None:
        print(f"  Posiform hardening: α={posiform_scale}")
        pi = info['posiform_info']
        if pi:
            print(f"  Posiform unique: "
                  f"{'OK' if pi['is_unique'] else 'FAIL'}")
            print(f"  Posiform clauses: {pi['num_clauses']}")
    print(f"  원래 변수: {info['n_original']}")
    print(f"  보조 변수: {info['n_aux']}")
    print(f"  QUBO 크기: {n_total} x {n_total}")
    print(f"  Q 비영 항: {len(Q)}개")
    print(f"  Ground energy: {info['ground_energy']:.4f}")
    print(f"  전체 target: {info['target'][:50]}"
          f"{'...' if len(info['target']) > 50 else ''}")
    if info.get('total_local_minima') is not None:
        print(f"  Local minima (blocks): {info['total_local_minima']}")

    # 저장
    if landscape == 'random':
        harden_tag = f"_a{posiform_scale}"
        filepath = os.path.join(
            results_dir,
            f"truthtable_concat_random_k{k}_h{h}"
            f"_total{n_total}{harden_tag}.csv")
    else:
        harden_tag = (f"_harden{posiform_scale}"
                      if posiform_scale is not None else "")
        filepath = os.path.join(
            results_dir,
            f"truthtable_concat_{mode}_k{k}_h{h}"
            f"_total{n_total}{harden_tag}.csv")
    save_qubo_edgelist(Q, filepath, info['target_with_aux'])
    print(f"\n  저장: {filepath}")


if __name__ == '__main__':
    main()
