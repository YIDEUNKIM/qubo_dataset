#!/usr/bin/env python3
"""
Degree-Bounded Möbius QUBO Generator
=====================================

Target-centered 변수에서 degree ≤ d의 Möbius 계수를 직접 생성하여 QUBO 구성.
2^n 진리표 열거 없이 ground state를 수학적으로 보장.

핵심 아이디어:
  1. 변수 치환: z_i = x_i ⊕ t_i (target에서 z=0)
  2. z-공간에서 양수 1차항 + 제한된 음수 고차항 생성
  3. Ground State 정리로 유일 최소 보장 (Theorem 1, 아래 참조)
  4. z→x 변환 후 Rosenberg 차수축소 → QUBO

시간 복잡도: O(n^d) (d 고정 시 다항식, 2^n 불필요)
Ground state: 수학적 보장 (열거 불필요)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Theorem 1 (Ground State Guarantee):

  f(z) = Σ_{S⊆[n], |S|≥1} c_S · ∏_{i∈S} z_i  (z ∈ {0,1}^n) 에서

  조건:
    (C1) c_{i} > 0  ∀i ∈ [n]         (양수 1차항)
    (C2) Σ_{|S|≥2} |c_S| < min_i c_{i}  (고차항 바운드)

  결론: f(0) = 0 이고, f(z) > 0  ∀z ≠ 0  (target이 유일 ground state)

  증명:
    z ≠ 0이면 T = supp(z) = {i : z_i = 1} ≠ ∅.
    f(z) = Σ_{S⊆T, |S|≥1} c_S
         = Σ_{i∈T} c_{i} + Σ_{S⊆T, |S|≥2} c_S
         ≥ Σ_{i∈T} c_{i} - Σ_{S⊆T, |S|≥2} |c_S|
         ≥ min_i c_{i} - Σ_{|S|≥2} |c_S|
         > 0  (by C2)                                          □

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

참고 문헌:
  - Boros & Hammer (2002): Pseudo-Boolean optimization
  - Rosenberg (1975): Reduction of bivalent programming
  - Pelofske, Hahn, Djidjev (2024): posiform hardening 기법

사용법:
    python degree_bounded/qubo_degree_bounded.py TARGET [d] [--ratio γ] [--density ρ]
    python degree_bounded/qubo_degree_bounded.py 10110 3
    python degree_bounded/qubo_degree_bounded.py 50 4 --ratio 0.8 --density 0.5
    python degree_bounded/qubo_degree_bounded.py 100 3 --harden 0.1 --seed 42
"""

import random
import itertools
import sys
import os
import time
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from qubo_utils import calculate_energy, save_qubo_edgelist, print_q_matrix, print_qubo_formula


# ─────────────────────────────────────────────────
#  Step 1: Degree-Bounded Möbius 계수 생성
# ─────────────────────────────────────────────────

def generate_coefficients(n, d, linear_range=(1.0, 5.0),
                          higher_ratio=0.8, density=1.0,
                          seed=None):
    """
    z-공간 (target-centered) Möbius 계수 생성.

    Theorem 1의 조건 (C1), (C2)를 만족하도록 구성:
      - 1차항: uniform(linear_range) > 0
      - 고차항: |c_S| ≤ budget/num_terms, 총 합 < γ · min_i c_{i}

    Args:
        n: 변수 수
        d: 최대 차수 (2 ≤ d ≤ n)
        linear_range: 1차항 계수 범위 (양수)
        higher_ratio: γ — 고차항 바운드 비율 (0 < γ < 1)
                      Σ_{|S|≥2} |c_S| < γ · min_i c_{i}
        density: 고차항 밀도 (0~1, 가능한 항 중 생성할 비율)
        seed: 랜덤 시드

    Returns:
        coefficients: {frozenset: coeff} — z-공간 Möbius 계수
        info: 생성 메타정보
    """
    rng = random.Random(seed)

    coefficients = {}

    # 1차항: 양수 랜덤
    linear_coeffs = []
    for i in range(n):
        c = rng.uniform(*linear_range)
        coefficients[frozenset([i])] = c
        linear_coeffs.append(c)

    c_min = min(linear_coeffs)
    budget = c_min * higher_ratio

    # 고차항 후보 열거
    higher_terms = []
    for degree in range(2, d + 1):
        for combo in itertools.combinations(range(n), degree):
            higher_terms.append(frozenset(combo))

    # density에 따라 부분 선택
    if density < 1.0 and higher_terms:
        num_select = max(1, int(len(higher_terms) * density))
        rng.shuffle(higher_terms)
        higher_terms = higher_terms[:num_select]

    # 각 항에 균등 배분된 최대 크기 내에서 랜덤 생성
    num_higher = len(higher_terms)
    used_budget = 0.0

    if num_higher > 0:
        max_per_term = budget / num_higher
        for S in higher_terms:
            c = rng.uniform(-max_per_term, max_per_term)
            used_budget += abs(c)
            if abs(c) > 1e-12:
                coefficients[S] = c

    # 안전 검증 (수학적 보장)
    assert used_budget < c_min, (
        f"Ground state 보장 실패: used_budget={used_budget:.6f} >= c_min={c_min:.6f}")

    info = {
        'n': n,
        'd': d,
        'c_min': c_min,
        'budget': budget,
        'used_budget': used_budget,
        'num_linear': n,
        'num_higher_candidates': num_higher,
        'num_nonzero_higher': sum(1 for S, c in coefficients.items() if len(S) >= 2),
        'ground_state_margin': c_min - used_budget,
    }

    return coefficients, info


# ─────────────────────────────────────────────────
#  Step 2: z-공간 → x-공간 변환
# ─────────────────────────────────────────────────

def transform_z_to_x(z_coefficients, target):
    """
    z-공간 Möbius 계수를 x-공간으로 변환.

    치환: z_i = x_i (t_i=0), z_i = 1-x_i (t_i=1)

    각 항 c_S · Π_{i∈S} z_i 에서 t_i=1인 변수를 (1-x_i)로 대입:
      Π_{j∈S∩T1} (1-x_j) = Σ_{U⊆S∩T1} (-1)^|U| · Π_{j∈U} x_j

    → c_S · Σ_U (-1)^|U| · Π_{k∈(S∩T0)∪U} x_k

    시간 복잡도: O(|terms| · 2^d)  (d 고정 시 선형)

    Args:
        z_coefficients: {frozenset: coeff} — z-공간 계수
        target: target 비트스트링

    Returns:
        x_coefficients: {frozenset: coeff} — x-공간 계수
    """
    t1_set = frozenset(i for i in range(len(target)) if target[i] == '1')

    x_coefficients = {}

    for S, c_S in z_coefficients.items():
        flip_vars = list(S & t1_set)
        keep_vars = S - t1_set

        # Expand Π_{j∈flip_vars} (1-x_j)
        for mask in range(1 << len(flip_vars)):
            U = frozenset(
                flip_vars[b] for b in range(len(flip_vars)) if (mask >> b) & 1)
            sign = (-1) ** len(U)
            var_set = keep_vars | U

            coeff = c_S * sign
            if var_set in x_coefficients:
                x_coefficients[var_set] += coeff
            else:
                x_coefficients[var_set] = coeff

    x_coefficients = {k: v for k, v in x_coefficients.items() if abs(v) > 1e-12}

    return x_coefficients


# ─────────────────────────────────────────────────
#  Step 3: 항 분류
# ─────────────────────────────────────────────────

def classify_terms(coefficients):
    """
    Möbius 계수를 차수별로 분류.

    Returns:
        constant, linear, quadratic, higher_order
    """
    constant = 0.0
    linear = {}
    quadratic = {}
    higher_order = {}

    for S, c in coefficients.items():
        if len(S) == 0:
            constant += c
        elif len(S) == 1:
            var = next(iter(S))
            linear[var] = linear.get(var, 0) + c
        elif len(S) == 2:
            i, j = sorted(S)
            quadratic[(i, j)] = quadratic.get((i, j), 0) + c
        else:
            higher_order[S] = higher_order.get(S, 0) + c

    return constant, linear, quadratic, higher_order


# ─────────────────────────────────────────────────
#  Step 4: Rosenberg 차수축소 (greedy)
# ─────────────────────────────────────────────────

def rosenberg_reduce_greedy(higher_order, n):
    """
    고차항을 Rosenberg 차수축소 (빈도 기반 탐욕적 쌍 선택 + 보조변수 재활용).

    각 라운드:
      1. degree ≤ 2 항 분리
      2. 전체 pending에서 가장 빈번한 변수 쌍 (a,b) 선택
      3. 보조변수 y = x_a · x_b 도입 (캐시 확인)
      4. {a,b} → {y} 치환

    Returns:
        reduced_linear, reduced_quadratic, aux_count, aux_info
    """
    reduced_linear = {}
    reduced_quadratic = {}
    next_aux = n
    aux_info = []
    product_cache = {}

    pending = []
    for S, c in higher_order.items():
        if len(S) == 1:
            var = next(iter(S))
            reduced_linear[var] = reduced_linear.get(var, 0) + c
        elif len(S) == 2:
            pair = tuple(sorted(S))
            reduced_quadratic[pair] = reduced_quadratic.get(pair, 0) + c
        else:
            pending.append((S, c))

    while pending:
        # degree 1, 2 분리
        new_pending = []
        for S, c in pending:
            if len(S) == 1:
                var = next(iter(S))
                reduced_linear[var] = reduced_linear.get(var, 0) + c
            elif len(S) == 2:
                pair = tuple(sorted(S))
                reduced_quadratic[pair] = reduced_quadratic.get(pair, 0) + c
            else:
                new_pending.append((S, c))
        pending = new_pending

        if not pending:
            break

        # 쌍 빈도 계산
        pair_freq = {}
        for S, c in pending:
            vars_sorted = sorted(S)
            for i in range(len(vars_sorted)):
                for j in range(i + 1, len(vars_sorted)):
                    pair = (vars_sorted[i], vars_sorted[j])
                    pair_freq[pair] = pair_freq.get(pair, 0) + 1

        best_pair = max(pair_freq, key=pair_freq.get)
        a, b = best_pair

        if (a, b) in product_cache:
            y = product_cache[(a, b)]
        else:
            y = next_aux
            next_aux += 1
            product_cache[(a, b)] = y
            aux_info.append((y, a, b))

        new_pending = []
        for S, c in pending:
            if a in S and b in S:
                new_S = (S - {a, b}) | {y}
                new_pending.append((new_S, c))
            else:
                new_pending.append((S, c))
        pending = new_pending

    return reduced_linear, reduced_quadratic, next_aux - n, aux_info


# ─────────────────────────────────────────────────
#  Step 5: 패널티 강도 계산 (분석적)
# ─────────────────────────────────────────────────

def compute_penalty_strength(z_coefficients, reduced_linear,
                             reduced_quadratic, aux_info):
    """
    진리표 열거 없이 패널티 강도 M 계산.

    에너지 범위 상한 (z-공간):
      E_max ≤ Σ_{|S|≥1} |c_S|  (모든 계수 절댓값 합)
      E_min = 0  (target)
      → energy_range ≤ Σ |c_S|

    M = max(energy_range_bound, max_aux_coeff_sum) + 1.0
    """
    energy_range_bound = sum(abs(c) for c in z_coefficients.values())

    max_aux_coeff = 0.0
    if aux_info and (reduced_linear or reduced_quadratic):
        aux_vars = set(y for (y, a, b) in aux_info)
        aux_coeff_sum = {y: 0.0 for y in aux_vars}

        for var, coeff in reduced_linear.items():
            if var in aux_vars:
                aux_coeff_sum[var] += abs(coeff)

        for (i, j), coeff in reduced_quadratic.items():
            if i in aux_vars:
                aux_coeff_sum[i] += abs(coeff)
            if j in aux_vars:
                aux_coeff_sum[j] += abs(coeff)

        if aux_coeff_sum:
            max_aux_coeff = max(aux_coeff_sum.values())

    return max(energy_range_bound, max_aux_coeff) + 1.0


# ─────────────────────────────────────────────────
#  Step 6: QUBO 조립
# ─────────────────────────────────────────────────

def assemble_qubo(n, constant, linear, quadratic,
                  reduced_linear, reduced_quadratic,
                  aux_count, aux_info, M):
    """QUBO 딕셔너리 조립. offset = constant."""
    Q = {}
    total_vars = n + aux_count

    def add(i, j, val):
        key = (min(i, j), max(i, j))
        Q[key] = Q.get(key, 0) + val

    for var, coeff in linear.items():
        add(var, var, coeff)

    for (i, j), coeff in quadratic.items():
        add(i, j, coeff)

    for var, coeff in reduced_linear.items():
        add(var, var, coeff)

    for (i, j), coeff in reduced_quadratic.items():
        add(i, j, coeff)

    # Rosenberg 패널티: M(x_a·x_b - 2x_a·y - 2x_b·y + 3y)
    for (y, a, b) in aux_info:
        add(a, b, M)
        add(a, y, -2 * M)
        add(b, y, -2 * M)
        add(y, y, 3 * M)

    Q = {k: v for k, v in Q.items() if abs(v) > 1e-15}

    return Q, constant, total_vars


# ─────────────────────────────────────────────────
#  보조 함수
# ─────────────────────────────────────────────────

def compute_aux_values(x_orig, aux_info):
    """target의 보조변수 값 계산"""
    x = {}
    for i, v in enumerate(x_orig):
        x[i] = v
    for (y, a, b) in aux_info:
        x[y] = x[a] * x[b]
    return x


def verify_ground_state(Q, target, n, aux_info, offset):
    """
    Ground state 검증.
    n ≤ 20: brute force 완전 검증
    n > 20: single-flip + random sampling
    """
    x_orig = [int(c) for c in target]
    x_full = compute_aux_values(x_orig, aux_info)
    total_vars = max(x_full.keys()) + 1 if x_full else n
    target_str = ''.join(str(x_full.get(i, 0)) for i in range(total_vars))
    target_energy = calculate_energy(target_str, Q)

    if n <= 20:
        best_energy = float('inf')
        best_state = None
        num_degenerate = 0
        for bits in range(1 << n):
            orig = [0] * n
            for i in range(n):
                orig[i] = (bits >> i) & 1
            x_test = compute_aux_values(orig, aux_info)
            state_str = ''.join(str(x_test.get(i, 0)) for i in range(total_vars))
            energy = calculate_energy(state_str, Q)
            if energy < best_energy - 1e-10:
                best_energy = energy
                best_state = ''.join(str(orig[i]) for i in range(n))
                num_degenerate = 1
            elif abs(energy - best_energy) < 1e-10:
                num_degenerate += 1

        return {
            'method': 'brute_force',
            'target_energy': target_energy + offset,
            'best_energy': best_energy + offset,
            'best_state': best_state,
            'is_ground_state': abs(best_energy - target_energy) < 1e-10,
            'is_unique': num_degenerate == 1,
            'num_degenerate': num_degenerate,
        }

    # n > 20: single-flip verification
    min_flip_delta = float('inf')
    for i in range(n):
        flipped = list(x_orig)
        flipped[i] = 1 - flipped[i]
        x_test = compute_aux_values(flipped, aux_info)
        state_str = ''.join(str(x_test.get(j, 0)) for j in range(total_vars))
        delta = calculate_energy(state_str, Q) - target_energy
        min_flip_delta = min(min_flip_delta, delta)

    lower_count = 0
    num_samples = min(10000, 1 << n)
    rng = random.Random(12345)
    for _ in range(num_samples):
        rand_orig = [rng.randint(0, 1) for _ in range(n)]
        x_test = compute_aux_values(rand_orig, aux_info)
        state_str = ''.join(str(x_test.get(j, 0)) for j in range(total_vars))
        if calculate_energy(state_str, Q) < target_energy - 1e-10:
            lower_count += 1

    return {
        'method': 'statistical',
        'target_energy': target_energy + offset,
        'min_flip_delta': min_flip_delta,
        'is_local_min': min_flip_delta > -1e-10,
        'lower_count': lower_count,
        'is_likely_global': lower_count == 0,
    }


# ─────────────────────────────────────────────────
#  메인 진입점
# ─────────────────────────────────────────────────

def create_qubo_degree_bounded(target, d=3,
                                linear_range=(1.0, 5.0),
                                higher_ratio=0.8,
                                density=1.0,
                                posiform_scale=None,
                                posiform_coeff_range=(1.0, 1.0),
                                seed=None, verbose=True):
    """
    Degree-Bounded Möbius QUBO 생성.

    Pipeline:
      1. z-공간 Möbius 계수 생성 (degree ≤ d, Theorem 1 조건 충족)
      2. z→x 변환 (다항식 치환, O(|terms|·2^d))
      3. 항 분류 (constant/linear/quadratic/higher_order)
      4. Rosenberg 차수축소 (degree ≥ 3, greedy)
      5. 패널티 강도 계산 (분석적, 진리표 불필요)
      6. QUBO 조립
      7. (선택) Posiform hardening

    Args:
        target: target 비트스트링 또는 정수 (정수면 랜덤 target 생성)
        d: 최대 차수 (2~n)
        linear_range: 1차항 계수 범위
        higher_ratio: γ — 고차항 바운드 비율 (0 < γ < 1)
        density: 고차항 밀도 (0~1)
        posiform_scale: posiform hardening α (None이면 안함)
        posiform_coeff_range: posiform 계수 범위
        seed: 랜덤 시드
        verbose: 출력 여부

    Returns:
        Q: QUBO 딕셔너리 {(i,j): weight}
        info: 메타정보
    """
    # target이 정수이면 랜덤 생성
    if isinstance(target, int):
        rng = random.Random(seed)
        target = ''.join(str(rng.randint(0, 1)) for _ in range(target))

    n = len(target)
    d = min(d, n)

    def log(msg):
        if verbose:
            print(msg)

    log(f"\n{'='*60}")
    log(f" Degree-Bounded Möbius QUBO Generator")
    log(f"{'='*60}")
    log(f" target: {target[:60]}{'...' if n > 60 else ''} (n={n})")
    log(f" degree: d={d}")
    log(f" linear_range: {linear_range}")
    log(f" higher_ratio: γ={higher_ratio}")
    log(f" density: {density}")

    t0 = time.time()

    # Step 1
    z_coeffs, gen_info = generate_coefficients(
        n, d, linear_range=linear_range,
        higher_ratio=higher_ratio, density=density, seed=seed)

    log(f"\n[Step 1] z-공간 Möbius 계수 생성")
    log(f"  1차항: {gen_info['num_linear']}개, c_min={gen_info['c_min']:.4f}")
    log(f"  고차항: {gen_info['num_nonzero_higher']}/{gen_info['num_higher_candidates']}개")
    log(f"  Budget 사용: {gen_info['used_budget']:.4f} / {gen_info['budget']:.4f}")
    log(f"  GS margin: {gen_info['ground_state_margin']:.4f} > 0 ✓")

    # Step 2
    x_coeffs = transform_z_to_x(z_coeffs, target)
    log(f"\n[Step 2] z→x 변환: {len(z_coeffs)}항 → {len(x_coeffs)}항")

    # Step 3
    constant, linear, quadratic, higher_order = classify_terms(x_coeffs)
    log(f"\n[Step 3] 항 분류: 상수={constant:.4f}, "
        f"1차={len(linear)}, 2차={len(quadratic)}, 3차+={len(higher_order)}")

    # Step 4
    if higher_order:
        reduced_linear, reduced_quadratic, aux_count, aux_info = \
            rosenberg_reduce_greedy(higher_order, n)
        log(f"\n[Step 4] Rosenberg: 고차항 {len(higher_order)}개 → 보조변수 {aux_count}개")
    else:
        reduced_linear, reduced_quadratic, aux_count, aux_info = {}, {}, 0, []
        log(f"\n[Step 4] 고차항 없음 — 차수축소 불필요")

    # Step 5
    if aux_count > 0:
        M = compute_penalty_strength(
            z_coeffs, reduced_linear, reduced_quadratic, aux_info)
        log(f"\n[Step 5] 패널티 강도: M = {M:.4f}")
    else:
        M = 0

    # Step 6
    Q, offset, total_vars = assemble_qubo(
        n, constant, linear, quadratic,
        reduced_linear, reduced_quadratic,
        aux_count, aux_info, M)

    # target 에너지 계산
    x_orig = [int(c) for c in target]
    x_full = compute_aux_values(x_orig, aux_info)
    target_with_aux = ''.join(str(x_full.get(i, 0)) for i in range(total_vars))
    target_energy_qubo = calculate_energy(target_with_aux, Q)
    target_energy_total = target_energy_qubo + offset

    log(f"\n[Step 6] QUBO 조립")
    log(f"  변수: {total_vars} (원래 {n} + 보조 {aux_count})")
    log(f"  비영 항: {len(Q)}개")
    log(f"  에너지 오프셋: {offset:.4f}")
    log(f"  E_total(target) = {target_energy_total:.6f}")

    # Step 7: Posiform hardening
    posiform_info_out = None
    if posiform_scale is not None:
        try:
            from posiform.qubo_posiform import create_qubo_posiform
        except ImportError:
            log("  [경고] posiform 모듈 import 실패")
            posiform_scale = None

    if posiform_scale is not None:
        from posiform.qubo_posiform import create_qubo_posiform

        log(f"\n[Step 7] Posiform Hardening (α={posiform_scale})")

        posiform_seed = (seed + 9999) if seed is not None else None
        Q_posiform, posiform_info_out = create_qubo_posiform(
            target,
            coeff_range=posiform_coeff_range,
            max_clauses_factor=20,
            seed=posiform_seed,
            targeted=False
        )

        if not posiform_info_out['is_unique']:
            for retry in range(5):
                posiform_seed = (seed + 10000 + retry) if seed is not None else None
                Q_posiform, posiform_info_out = create_qubo_posiform(
                    target,
                    coeff_range=posiform_coeff_range,
                    max_clauses_factor=30,
                    seed=posiform_seed,
                    targeted=False
                )
                if posiform_info_out['is_unique']:
                    break

        if not posiform_info_out['is_unique']:
            warnings.warn(f"Posiform uniqueness 보장 실패 (n={n})")

        for (a, b), w in Q_posiform.items():
            key = (min(a, b), max(a, b))
            Q[key] = Q.get(key, 0) + posiform_scale * w

        Q = {k: v for k, v in Q.items() if abs(v) > 1e-15}

        target_energy_qubo = calculate_energy(target_with_aux, Q)
        target_energy_total = target_energy_qubo + offset

        log(f"  Posiform unique: {'OK' if posiform_info_out['is_unique'] else 'FAIL'}")
        log(f"  E_total(target) with posiform: {target_energy_total:.6f}")

    elapsed = time.time() - t0
    log(f"\n 생성 완료 ({elapsed:.3f}s)")

    info = {
        'n': n,
        'd': d,
        'n_original': n,
        'n_aux': aux_count,
        'n_total': total_vars,
        'target': target,
        'target_with_aux': target_with_aux,
        'target_energy': target_energy_total,
        'offset': offset,
        'penalty_M': M,
        'linear_range': linear_range,
        'higher_ratio': higher_ratio,
        'density': density,
        'c_min': gen_info['c_min'],
        'ground_state_margin': gen_info['ground_state_margin'],
        'num_higher_terms': gen_info['num_nonzero_higher'],
        'num_qubo_terms': len(Q),
        'aux_info': aux_info,
        'posiform_scale': posiform_scale,
        'posiform_info': posiform_info_out,
        'generation_time': elapsed,
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
        print("  python degree_bounded/qubo_degree_bounded.py TARGET [d] [옵션]")
        print()
        print("인자:")
        print("  TARGET  비트스트링 (예: 10110) 또는 정수 (랜덤 생성)")
        print("  d       최대 차수 (기본: 3)")
        print()
        print("옵션:")
        print("  --ratio γ      고차항 바운드 비율 (기본: 0.8)")
        print("  --density ρ    고차항 밀도 (기본: 1.0)")
        print("  --harden α     posiform hardening")
        print("  --seed S       랜덤 시드")
        print()
        print("예시:")
        print("  python degree_bounded/qubo_degree_bounded.py 10110 3")
        print("  python degree_bounded/qubo_degree_bounded.py 50 4 --ratio 0.8 --density 0.5")
        print("  python degree_bounded/qubo_degree_bounded.py 100 3 --harden 0.1 --seed 42")
        return

    # Target
    arg1 = sys.argv[1]
    if arg1.isdigit() and len(arg1) > 1 and all(c in '01' for c in arg1):
        target = arg1
    elif arg1.isdigit():
        target = int(arg1)
    else:
        target = arg1

    # Degree
    d = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 3

    # Options
    def parse_opt(name, default):
        if name in sys.argv:
            return float(sys.argv[sys.argv.index(name) + 1])
        return default

    higher_ratio = parse_opt('--ratio', 0.8)
    density = parse_opt('--density', 1.0)
    posiform_scale = None
    if '--harden' in sys.argv:
        posiform_scale = parse_opt('--harden', 0.1)

    seed = int(parse_opt('--seed', -1))
    seed = seed if seed >= 0 else None

    Q, info = create_qubo_degree_bounded(
        target, d=d,
        linear_range=(1.0, 5.0),
        higher_ratio=higher_ratio,
        density=density,
        posiform_scale=posiform_scale,
        seed=seed)

    n = info['n']
    total_vars = info['n_total']

    # 검증
    print(f"\n[Ground State 검증]")
    result = verify_ground_state(Q, info['target'], n, info['aux_info'], info['offset'])
    if result['method'] == 'brute_force':
        print(f"  방법: brute force (n={n})")
        print(f"  Ground state: {'OK' if result['is_ground_state'] else 'FAIL'}")
        print(f"  유일성: {'OK' if result['is_unique'] else 'FAIL'} "
              f"(축퇴도={result['num_degenerate']})")
        print(f"  Best state: {result['best_state']}")
        print(f"  Best energy: {result['best_energy']:.6f}")
    else:
        print(f"  방법: statistical (n={n})")
        print(f"  Local minimum: {'OK' if result['is_local_min'] else 'FAIL'}")
        print(f"  Min flip delta: {result['min_flip_delta']:.6f}")
        print(f"  Global (추정): {'OK' if result['is_likely_global'] else 'FAIL'}")

    if total_vars <= 10:
        print_q_matrix(Q, total_vars)
        print_qubo_formula(Q)

    # 저장
    target_str = info['target']
    harden_tag = f"_harden{posiform_scale}" if posiform_scale is not None else ""
    filepath = os.path.join(
        results_dir,
        f"degree_bounded_n{n}_d{d}{harden_tag}.csv")
    save_qubo_edgelist(Q, filepath, info['target_with_aux'])
    print(f"\n[저장 완료] {filepath}")


if __name__ == '__main__':
    main()
