#!/usr/bin/env python3
"""
Truth Table QUBO Generator
===========================

진리표(비트스트링 → 에너지 매핑)로부터 QUBO를 구성한다.

파이프라인:
  1. Möbius 변환: 진리표 → 다중선형 다항식 (O(n·2^n))
  2. 항 분류: 상수 / 선형 / 이차 / 고차(≥3)
  3. Rosenberg 차수축소: 고차항 → 보조변수 + 패널티 → 2차화
  4. QUBO 조립 + 전수검사 검증

사용법:
    python truthtable/qubo_truthtable.py '{"000":3,"001":4,"010":4,"011":5,"100":3,"101":5,"110":2,"111":1}'
    python truthtable/qubo_truthtable.py --preset gap 8 10110011 --gap 2.0 --seed 42
    python truthtable/qubo_truthtable.py --preset valley 6 101010,010101 --barrier 5.0 --seed 42
"""

import numpy as np
import sys
import os
import json

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qubo_utils import calculate_energy, save_qubo_edgelist, print_q_matrix, print_qubo_formula


# ─────────────────────────────────────────────────
#  Step 1: Möbius 변환
# ─────────────────────────────────────────────────

def mobius_transform(truth_table, n):
    """
    진리표 → 다중선형 다항식 계수 (fast Möbius 역변환, O(n·2^n))

    비트스트링 "abc..." → x_0=a, x_1=b, x_2=c, ...
    mask의 bit i = variable x_i

    Args:
        truth_table: dict {bitstring: energy} / list [energy, ...] / callable
        n: 변수 개수

    Returns:
        coefficients: dict {frozenset(변수 인덱스): 계수}
    """
    size = 1 << n
    values = np.zeros(size, dtype=float)

    if isinstance(truth_table, dict):
        for bitstring, energy in truth_table.items():
            mask = 0
            for i, ch in enumerate(bitstring):
                if ch == '1':
                    mask |= (1 << i)
            values[mask] = energy
    elif callable(truth_table):
        for mask in range(size):
            bitstring = ''.join(str((mask >> i) & 1) for i in range(n))
            values[mask] = truth_table(bitstring)
    else:
        values[:] = truth_table

    # Fast Möbius 역변환 (포함-배제 원리)
    c = values.copy()
    for i in range(n):
        for s in range(size):
            if s & (1 << i):
                c[s] -= c[s ^ (1 << i)]

    # 0이 아닌 계수만 추출
    coefficients = {}
    for s in range(size):
        if abs(c[s]) > 1e-12:
            variables = frozenset(i for i in range(n) if s & (1 << i))
            coefficients[variables] = float(c[s])

    return coefficients


# ─────────────────────────────────────────────────
#  Step 2: 항 분류
# ─────────────────────────────────────────────────

def classify_terms(coefficients):
    """
    다중선형 계수를 차수별로 분류

    Returns:
        constant, linear {var: c}, quadratic {(i,j): c}, higher_order {frozenset: c}
    """
    constant = 0.0
    linear = {}
    quadratic = {}
    higher_order = {}

    for variables, coeff in coefficients.items():
        degree = len(variables)
        if degree == 0:
            constant = coeff
        elif degree == 1:
            (var,) = variables
            linear[var] = coeff
        elif degree == 2:
            i, j = sorted(variables)
            quadratic[(i, j)] = coeff
        else:
            higher_order[variables] = coeff

    return constant, linear, quadratic, higher_order


# ─────────────────────────────────────────────────
#  Step 3: Rosenberg 차수축소
# ─────────────────────────────────────────────────

def rosenberg_reduce(higher_order, n):
    """
    고차항을 Rosenberg 체이닝으로 2차화

    k차 단항식 c·x_{i1}...x_{ik}:
      y = x_{i1}·x_{i2} 도입, 패널티 M(x_{i1}x_{i2} - 2x_{i1}y - 2x_{i2}y + 3y) 추가
      → c·y·x_{i3}...x_{ik} 로 치환, 여전히 차수 > 2이면 반복

    Returns:
        reduced_linear:    dict {var: coeff}
        reduced_quadratic: dict {(i,j): coeff}
        aux_count:         보조변수 수
        aux_info:          list of (aux_idx, var_a, var_b)
    """
    reduced_linear = {}
    reduced_quadratic = {}
    aux_info = []
    next_aux = n

    pending = [(frozenset(vs), c) for vs, c in higher_order.items()]

    while pending:
        new_pending = []
        for variables, coeff in pending:
            var_list = sorted(variables)
            degree = len(var_list)

            if degree == 1:
                reduced_linear[var_list[0]] = reduced_linear.get(var_list[0], 0) + coeff
            elif degree == 2:
                key = (var_list[0], var_list[1])
                reduced_quadratic[key] = reduced_quadratic.get(key, 0) + coeff
            else:
                # 보조변수 y = var_list[0] * var_list[1]
                y = next_aux
                next_aux += 1
                a, b = var_list[0], var_list[1]
                aux_info.append((y, a, b))

                remaining = frozenset(var_list[2:]) | frozenset([y])
                new_pending.append((remaining, coeff))

        pending = new_pending

    return reduced_linear, reduced_quadratic, next_aux - n, aux_info


# ─────────────────────────────────────────────────
#  Step 4: 패널티 강도 결정
# ─────────────────────────────────────────────────

def compute_penalty_strength(truth_table, n):
    """
    최소 패널티 강도 M 계산

    M > max(E) - min(E) 이면 보조변수 위반 시 에너지가
    항상 정상 상태보다 높아짐 → ground state 보장
    """
    if isinstance(truth_table, dict):
        energies = list(truth_table.values())
    elif callable(truth_table):
        energies = [
            truth_table(''.join(str((m >> i) & 1) for i in range(n)))
            for m in range(1 << n)
        ]
    else:
        energies = list(truth_table)

    return float(max(energies) - min(energies) + 1.0)


# ─────────────────────────────────────────────────
#  Step 5: QUBO 조립
# ─────────────────────────────────────────────────

def assemble_qubo(n, constant, linear, quadratic,
                  reduced_linear, reduced_quadratic,
                  aux_count, aux_info, M):
    """
    Q 행렬 조립

    Returns:
        Q: dict {(i,j): weight}  (upper triangular, i ≤ j)
        offset: 상수항 (에너지 오프셋)
        total_vars: 총 변수 수
    """
    total_vars = n + aux_count
    Q = {}

    def add(i, j, val):
        key = (min(i, j), max(i, j))
        Q[key] = Q.get(key, 0) + val

    # 원래 선형항 → 대각
    for var, coeff in linear.items():
        add(var, var, coeff)

    # 원래 이차항
    for (i, j), coeff in quadratic.items():
        add(i, j, coeff)

    # 축소된 선형/이차항
    for var, coeff in reduced_linear.items():
        add(var, var, coeff)
    for (i, j), coeff in reduced_quadratic.items():
        add(i, j, coeff)

    # Rosenberg 패널티: M·(a·b - 2a·y - 2b·y + 3y)
    for (y, a, b) in aux_info:
        add(a, b, M)
        add(a, y, -2 * M)
        add(b, y, -2 * M)
        add(y, y, 3 * M)

    Q = {k: v for k, v in Q.items() if abs(v) > 1e-15}
    return Q, constant, total_vars


# ─────────────────────────────────────────────────
#  Step 6: 검증
# ─────────────────────────────────────────────────

def compute_aux_values(x_orig, aux_info):
    """원래 변수 값으로부터 올바른 보조변수 값 계산"""
    x = dict(enumerate(x_orig))
    for (y, a, b) in aux_info:
        x[y] = x[a] * x[b]
    return x


def verify_qubo(Q, truth_table, n, aux_info, offset):
    """
    모든 2^n 비트스트링에 대해 QUBO 에너지 == 진리표 에너지 검증

    Returns:
        errors: list of (bitstring, expected, actual)
        ground_state: (bitstring, energy)
    """
    errors = []
    min_energy = float('inf')
    min_bitstring = None

    for mask in range(1 << n):
        x_orig = [(mask >> i) & 1 for i in range(n)]
        bitstring = ''.join(str(b) for b in x_orig)

        x = compute_aux_values(x_orig, aux_info)

        energy = offset
        for (i, j), w in Q.items():
            if i == j:
                energy += w * x.get(i, 0)
            else:
                energy += w * x.get(i, 0) * x.get(j, 0)

        if isinstance(truth_table, dict):
            expected = truth_table[bitstring]
        elif callable(truth_table):
            expected = truth_table(bitstring)
        else:
            expected = truth_table[mask]

        if abs(energy - expected) > 1e-6:
            errors.append((bitstring, expected, energy))

        if energy < min_energy:
            min_energy = energy
            min_bitstring = bitstring

    return errors, (min_bitstring, min_energy)


# ─────────────────────────────────────────────────
#  에너지 함수 프리셋
# ─────────────────────────────────────────────────

def preset_energy_gap(n, target, gap=2.0, noise_scale=1.0, seed=None):
    """
    프리셋: 에너지 갭 제어

    E(target) = 0
    E(x ≠ target) = gap + |N(0, noise_scale)|

    gap이 작을수록 난이도 ↑ (양자 어닐링의 minimum spectral gap과 직결)
    """
    rng = np.random.default_rng(seed)
    truth_table = {}

    for mask in range(1 << n):
        bitstring = ''.join(str((mask >> i) & 1) for i in range(n))
        if bitstring == target:
            truth_table[bitstring] = 0.0
        else:
            truth_table[bitstring] = gap + abs(rng.normal(0, noise_scale))

    return truth_table


def preset_multi_valley(n, targets, gap=0.5, barrier_height=5.0, seed=None):
    """
    프리셋: 다중 계곡 (Multiple Local Minima)

    targets[0]: global minimum (energy = 0)
    targets[1:]: local minima (energy = gap)
    나머지: Hamming 거리 비례 에너지 + 장벽

    양자 tunneling 능력 벤치마크 (고전 SA는 local minima에 갇힘)
    """
    rng = np.random.default_rng(seed)
    truth_table = {}

    for mask in range(1 << n):
        bitstring = ''.join(str((mask >> i) & 1) for i in range(n))

        distances = [sum(a != b for a, b in zip(bitstring, t)) for t in targets]
        min_idx = int(np.argmin(distances))
        min_dist = distances[min_idx]

        if min_dist == 0:
            energy = 0.0 if min_idx == 0 else gap
        else:
            energy = gap + barrier_height * (min_dist / n) + rng.uniform(0, 0.5)

        truth_table[bitstring] = energy

    return truth_table


# ─────────────────────────────────────────────────
#  메인 진입점
# ─────────────────────────────────────────────────

def create_qubo_truthtable(truth_table, n=None, seed=None, verbose=True):
    """
    진리표로부터 QUBO 생성

    Args:
        truth_table: dict {bitstring: energy} / list / callable
        n: 변수 개수 (dict일 때 자동 추론)
        verbose: 진행 상황 출력 여부

    Returns:
        Q: dict {(i,j): weight}
        info: dict (메타데이터)
    """
    if n is None:
        if isinstance(truth_table, dict):
            n = len(next(iter(truth_table.keys())))
        elif isinstance(truth_table, (list, np.ndarray)):
            n = int(np.log2(len(truth_table)))
        else:
            raise ValueError("함수 기반 입력 시 n을 명시해야 합니다")

    def log(msg):
        if verbose:
            print(msg)

    log(f"\n{'='*60}")
    log(f" Truth Table QUBO Generator")
    log(f"{'='*60}")
    log(f" 변수 수: n = {n}")
    log(f" 진리표 크기: 2^{n} = {1 << n}")

    # Step 1
    log(f"\n[Step 1] Möbius 변환 (진리표 → 다중선형 다항식)...")
    coefficients = mobius_transform(truth_table, n)
    log(f"  비영 항 수: {len(coefficients)}")

    # Step 2
    log(f"\n[Step 2] 항 분류...")
    constant, linear, quadratic, higher_order = classify_terms(coefficients)
    log(f"  상수항: offset = {constant:.4f}")
    log(f"  선형항 (1차): {len(linear)}개")
    log(f"  이차항 (2차): {len(quadratic)}개")
    log(f"  고차항 (≥3차): {len(higher_order)}개")
    if higher_order:
        max_deg = max(len(v) for v in higher_order.keys())
        log(f"  최대 차수: {max_deg}")

    # Step 3
    if higher_order:
        log(f"\n[Step 3] Rosenberg 차수축소...")
        reduced_linear, reduced_quadratic, aux_count, aux_info = \
            rosenberg_reduce(higher_order, n)
        log(f"  보조변수: {aux_count}개")
    else:
        reduced_linear, reduced_quadratic, aux_count, aux_info = {}, {}, 0, []
        log(f"\n[Step 3] 고차항 없음 — 차수축소 불필요")

    # Step 4
    M = compute_penalty_strength(truth_table, n) if aux_count > 0 else 0
    if aux_count > 0:
        log(f"\n[Step 4] 패널티 강도: M = {M:.4f}")

    # Step 5
    log(f"\n[Step 5] QUBO 조립...")
    Q, offset, total_vars = assemble_qubo(
        n, constant, linear, quadratic,
        reduced_linear, reduced_quadratic,
        aux_count, aux_info, M
    )
    log(f"  총 변수: {total_vars} (원래 {n} + 보조 {aux_count})")
    log(f"  Q 비영 항: {len(Q)}개")

    # Step 6
    log(f"\n[Step 6] 전수검사 검증...")
    errors, ground_state = verify_qubo(Q, truth_table, n, aux_info, offset)

    if errors:
        log(f"  *** 검증 실패: {len(errors)}개 불일치 ***")
        for bs, exp, act in errors[:5]:
            log(f"    {bs}: 기대={exp:.4f}, 실제={act:.4f}")
    else:
        log(f"  모든 {1 << n}개 비트스트링 에너지 일치 확인")

    log(f"\n  Ground state: {ground_state[0]} (energy = {ground_state[1]:.4f})")

    info = {
        'n_original': n,
        'n_aux': aux_count,
        'n_total': total_vars,
        'offset': offset,
        'penalty_M': M,
        'ground_state': ground_state,
        'n_higher_order': len(higher_order),
        'aux_info': aux_info,
        'target': ground_state[0],
    }
    return Q, info


# ─────────────────────────────────────────────────
#  CLI
# ─────────────────────────────────────────────────

def main():
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(results_dir, exist_ok=True)

    if len(sys.argv) < 2 or sys.argv[1] == '--help':
        print("사용법:")
        print("  python qubo_truthtable.py '{\"000\":3,...}'")
        print("  python qubo_truthtable.py --preset gap N TARGET [--gap G] [--noise S] [--seed S]")
        print("  python qubo_truthtable.py --preset valley N T1,T2,... [--gap G] [--barrier B] [--seed S]")
        print()
        print("예시:")
        tt = '{"000":3,"001":4,"010":4,"011":5,"100":3,"101":5,"110":2,"111":1}'
        print(f"  python qubo_truthtable.py '{tt}'")
        print("  python qubo_truthtable.py --preset gap 8 10110011 --gap 2.0 --seed 42")
        print("  python qubo_truthtable.py --preset valley 6 101010,010101 --seed 42")
        return

    def parse_opt(name, default):
        if name in sys.argv:
            return float(sys.argv[sys.argv.index(name) + 1])
        return default

    if sys.argv[1] == '--preset':
        preset_name = sys.argv[2]
        seed = int(parse_opt('--seed', 0)) if '--seed' in sys.argv else None

        if preset_name == 'gap':
            n = int(sys.argv[3])
            target = sys.argv[4]
            gap = parse_opt('--gap', 2.0)
            noise = parse_opt('--noise', 1.0)
            print(f"프리셋: Energy Gap (n={n}, target={target}, gap={gap}, noise={noise})")
            truth_table = preset_energy_gap(n, target, gap=gap, noise_scale=noise, seed=seed)

        elif preset_name == 'valley':
            n = int(sys.argv[3])
            targets = sys.argv[4].split(',')
            gap = parse_opt('--gap', 0.5)
            barrier = parse_opt('--barrier', 5.0)
            print(f"프리셋: Multi-Valley (n={n}, targets={targets}, gap={gap}, barrier={barrier})")
            truth_table = preset_multi_valley(n, targets, gap=gap,
                                              barrier_height=barrier, seed=seed)
        else:
            print(f"알 수 없는 프리셋: {preset_name}")
            return
    else:
        truth_table = json.loads(sys.argv[1])
        truth_table = {k: float(v) for k, v in truth_table.items()}

    Q, info = create_qubo_truthtable(truth_table)

    n = info['n_original']
    total = info['n_total']

    print(f"\n{'='*60}")
    print(f" 결과 요약")
    print(f"{'='*60}")
    print(f"  원래 변수: {n}")
    print(f"  보조 변수: {info['n_aux']}")
    print(f"  QUBO 크기: {total} x {total}")
    print(f"  에너지 오프셋: {info['offset']:.4f}")
    if info['n_aux'] > 0:
        print(f"  패널티 M: {info['penalty_M']:.4f}")
    print(f"  Ground state: {info['ground_state'][0]} (energy = {info['ground_state'][1]:.4f})")

    if total <= 10:
        print_q_matrix(Q, total)
        print_qubo_formula(Q)

    # 저장: target에 보조변수 값도 포함
    x_orig = [int(info['target'][i]) for i in range(n)]
    x_full = compute_aux_values(x_orig, info['aux_info'])
    full_target = ''.join(str(x_full[i]) for i in range(total))

    filepath = os.path.join(results_dir, f"truthtable_n{n}_total{total}.csv")
    save_qubo_edgelist(Q, filepath, full_target)
    print(f"\n  저장: {filepath}")


if __name__ == '__main__':
    main()
