#!/usr/bin/env python3
"""
Truth Table Concat QUBO Generator
====================================

k-bit Truth Table QUBO를 h개 생성하여 block-diagonal로 접합.
총 변수 수 = k*h (approx) 또는 k*h + aux (exact).

각 블록은 독립적이므로 전체 ground state = target 반복 h회 (유일성 보장).
블록마다 다른 seed로 진리표를 생성하여 landscape가 다르지만 target은 동일.

파이프라인:
  1. 블록별로 preset_random_landscape → 진리표 생성
  2. 각 블록을 approx 또는 exact 모드로 QUBO 변환
  3. 변수 인덱스를 shift하여 block-diagonal Q 구성
  4. 전체 target 검증

사용법:
    python truthtable_concat/qubo_truthtable_concat.py 1001111 10
    python truthtable_concat/qubo_truthtable_concat.py 1001111 10 --exact
    python truthtable_concat/qubo_truthtable_concat.py 1001111 10 --gap 2.0 --noise 1.0 --seed 42
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from truthtable.qubo_truthtable import (
    create_qubo_truthtable, create_qubo_approx,
    preset_random_landscape, compute_aux_values,
)
from qubo_utils import calculate_energy, save_qubo_edgelist


def create_qubo_concat(target, h,
                       mode='approx', epsilon=0.01,
                       reduce_strategy='greedy', seed=None, verbose=True):
    """
    k-bit Truth Table QUBO를 h개 block-diagonal로 접합
    block-diagonal: 큰 행렬이 대각선 방향으로 독립적인 작은 블록들로 이루어져 있고, 블록 바깥은 모두 0인 형태
    Args:
        target: k-bit 문자열 (블록 크기 = len(target))
        h: 블록 수
        mode: 'approx' (보조변수 0, 기본값) / 'exact' (Rosenberg)
        epsilon: 근사 모드 에너지 갭 하한
        reduce_strategy: exact 모드 차수축소 전략 ('original'/'cache'/'greedy')
        seed: 랜덤 시드 (블록 i는 seed+i 사용)
        verbose: 출력 여부

    Returns:
        Q: dict {(i,j): weight} — block-diagonal QUBO
        info: dict (메타데이터)
    """
    k = len(target)

    def log(msg):
        if verbose:
            print(msg)

    log(f"\n{'='*60}")
    log(f" Truth Table Concat QUBO Generator")
    log(f"{'='*60}")
    log(f" 블록 target: {target} (k={k})")
    log(f" 블록 수: h={h}")
    log(f" 모드: {mode}")
    log(f" 프리셋: random landscape (uniform 0.1~5.0)")
    if seed is not None:
        log(f" seed={seed}")

    Q_total = {}
    blocks_info = []
    all_aux_info = []
    offset = 0  # 현재 블록의 변수 인덱스 시작점

    for i in range(h):
        seed_i = (seed + i) if seed is not None else None
        tt = preset_random_landscape(k, target, seed=seed_i)

        if mode == 'approx':
            Q_block, info_block = create_qubo_approx(tt, verbose=False, epsilon=epsilon)
        else:
            Q_block, info_block = create_qubo_truthtable(
                tt, verbose=False, reduce_strategy=reduce_strategy)

        block_total_vars = info_block['n_total']

        # 변수 인덱스 shift
        for (a, b), w in Q_block.items():
            Q_total[(a + offset, b + offset)] = w

        # aux_info reindex
        for (y, a, b) in info_block['aux_info']:
            all_aux_info.append((y + offset, a + offset, b + offset))

        # 블록별 QUBO 에너지 계산 (target에서)
        block_aux = info_block['aux_info']
        x_orig = [int(c) for c in target]
        x_block = compute_aux_values(x_orig, block_aux)
        block_target_str = ''.join(str(x_block[j]) for j in range(block_total_vars))
        block_qubo_energy = calculate_energy(block_target_str, Q_block)

        blocks_info.append({
            'block_idx': i,
            'offset': offset,
            'n_total': block_total_vars,
            'n_aux': info_block['n_aux'],
            'qubo_energy': block_qubo_energy,
            'seed': seed_i,
        })

        log(f"  블록 {i}: offset={offset}, vars={block_total_vars}, "
            f"aux={info_block['n_aux']}, E_qubo={block_qubo_energy:.4f}")

        offset += block_total_vars

    n_total = offset
    n_original = k * h
    full_target = target * h

    # exact 모드: 전체 target에 aux 비트 추가
    if mode == 'exact' and all_aux_info:
        # 전체 target bitstring 구성 (원래 변수 + 보조변수)
        x_full = {}
        block_offset = 0
        for bi in blocks_info:
            block_n_total = bi['n_total']
            block_n_aux = bi['n_aux']
            block_k = block_n_total - block_n_aux
            # 원래 변수 값
            for j in range(block_k):
                x_full[block_offset + j] = int(target[j])
            block_offset += block_n_total

        # 보조변수 값 계산
        for (y, a, b) in all_aux_info:
            x_full[y] = x_full.get(a, 0) * x_full.get(b, 0)

        full_target_with_aux = ''.join(str(x_full.get(i, 0)) for i in range(n_total))
    else:
        full_target_with_aux = full_target

    # 검증: 전체 target으로 에너지 계산
    verify_energy = calculate_energy(full_target_with_aux, Q_total)
    block_energy_sum = sum(bi['qubo_energy'] for bi in blocks_info)
    log(f"\n 전체 target: {full_target} (원래 변수 {n_original}개)")
    log(f" 전체 QUBO 크기: {n_total} (aux 포함)")
    log(f" Ground energy: {verify_energy:.4f}")
    log(f" 블록별 합계: {block_energy_sum:.4f}")

    energy_match = abs(verify_energy - block_energy_sum) < 1e-4
    log(f" 에너지 일치: {'YES' if energy_match else 'NO'}")
    if not energy_match:
        log(f"  *** 경고: 에너지 불일치! diff={abs(verify_energy - block_energy_sum):.6f}")

    info = {
        'n_block': k,
        'h': h,
        'n_original': n_original,
        'n_total': n_total,
        'n_aux': n_total - n_original,
        'block_target': target,
        'target': full_target,
        'target_with_aux': full_target_with_aux,
        'mode': mode,
        'blocks': blocks_info,
        'aux_info': all_aux_info,
        'ground_energy': verify_energy,
    }

    return Q_total, info


# ─────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────

def main():
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(results_dir, exist_ok=True)

    if len(sys.argv) < 3 or sys.argv[1] == '--help':
        print("사용법:")
        print("  python truthtable_concat/qubo_truthtable_concat.py TARGET H [옵션]")
        print()
        print("인자:")
        print("  TARGET  블록 target 비트스트링 (예: 1001111)")
        print("  H       블록 수")
        print()
        print("옵션:")
        print("  --exact          정확 모드 (Rosenberg, 기본: approx)")
        print("  --seed S         랜덤 시드")
        print("  --strategy S     차수축소 전략 (original/cache/greedy, 기본: greedy)")
        print()
        print("예시:")
        print("  python truthtable_concat/qubo_truthtable_concat.py 1001111 10")
        print("  python truthtable_concat/qubo_truthtable_concat.py 1001111 10 --exact")
        print("  python truthtable_concat/qubo_truthtable_concat.py 1001111 10 --seed 42")
        return

    target = sys.argv[1]
    h = int(sys.argv[2])

    use_exact = '--exact' in sys.argv

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

    mode = 'exact' if use_exact else 'approx'

    Q, info = create_qubo_concat(
        target, h,
        mode=mode, reduce_strategy=reduce_strategy, seed=seed)

    k = info['n_block']
    n_total = info['n_total']

    print(f"\n{'='*60}")
    print(f" 결과 요약")
    print(f"{'='*60}")
    print(f"  블록 target: {target} (k={k})")
    print(f"  블록 수: h={h}")
    print(f"  모드: {mode}")
    print(f"  원래 변수: {info['n_original']}")
    print(f"  보조 변수: {info['n_aux']}")
    print(f"  QUBO 크기: {n_total} x {n_total}")
    print(f"  Q 비영 항: {len(Q)}개")
    print(f"  Ground energy: {info['ground_energy']:.4f}")
    print(f"  전체 target: {info['target'][:50]}{'...' if len(info['target']) > 50 else ''}")

    # 저장
    filepath = os.path.join(
        results_dir,
        f"truthtable_concat_{mode}_k{k}_h{h}_total{n_total}.csv")
    save_qubo_edgelist(Q, filepath, info['target_with_aux'])
    print(f"\n  저장: {filepath}")


if __name__ == '__main__':
    main()
