"""
기댓값 0 QUBO 생성기

기본 모델(ZeroOffDiagonalModel): 각 off-diagonal 원소 q_ij의 기댓값이 0.
  double-flip ratio = single-flip ratio 합 → |계수| multiset {1, 1, 2}가 target 무관.
  diagonal E[q_ii]는 b_i 부호를 노출하나, SNR ~ √n으로 off-diagonal 누출(SNR ~ n)보다 양호.
"""

import itertools
import random
import warnings
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from abc import ABC, abstractmethod
from qubo_utils import calculate_energy, save_qubo_edgelist, print_q_matrix, print_qubo_formula

# -----------------------------------------------------------------------------
# Penalty Strategies (Strategy Pattern)
# -----------------------------------------------------------------------------

class PenaltyModel(ABC):
    """
    QUBO 페널티 비율을 정의하기 위한 추상 기본 클래스(Abstract Base Class)입니다.
    페널티 분배를 위한 다양한 전략 패턴을 허용합니다.
    """
    @abstractmethod
    def get_ratios(self, target_pair: tuple) -> dict:
        """
        penalty_state -> ratio 매핑 딕셔너리를 반환합니다.
        예시: 목표가 (0,0)인 경우, {(0,1): 1.0, (1,0): 1.0, (1,1): 1.65} 반환
        """
        pass

class DefaultZeroExpectationModel(PenaltyModel):
    """
    [DEPRECATED] LP 기반 비율 모델.

    문제점: LP 비율(1.64, 1.68 등)은 "3개 중 1개 선택" 카운팅에서 도출되었으나,
    현재 코드는 3개 penalty를 모두 추가하므로 전제가 불일치함.
    off-diagonal 기댓값이 0이 아님:
      target (0,0): E[q_ij] = μ(-1 - 1 + 1.64) = -0.36μ
      target (0,1): E[q_ij] = μ(+2.00 - 1.00 + 1.68) = +2.68μ
      target (1,0): E[q_ij] = μ(+2.00 - 1.00 + 1.68) = +2.68μ
      target (1,1): E[q_ij] = μ(+1.00 - 3.00 - 3.00) = -5.00μ

    대안: ZeroOffDiagonalModel 사용 (off-diagonal E[q_ij] = 0 수학적 보장).
    """
    def __init__(self):
        warnings.warn(
            "DefaultZeroExpectationModel은 deprecated입니다. "
            "off-diagonal E[q_ij] ≠ 0. ZeroOffDiagonalModel을 사용하세요.",
            DeprecationWarning,
            stacklevel=2,
        )
        # 정규화 기준: 최솟값 = 1.0 (LP 결과값 캡슐화)
        self._ratios_table = {
            # 정답 = (0, 0): 페널티
            (0, 0): {(0, 1): 1.00, (1, 0): 1.00, (1, 1): 1.64},
            # 정답 = (0, 1)
            (0, 1): {(0, 0): 2.00, (1, 0): 1.00, (1, 1): 1.68},
            # 정답 = (1, 0)
            (1, 0): {(0, 0): 2.00, (0, 1): 1.00, (1, 1): 1.68},
            # 정답 = (1, 1)
            (1, 1): {(0, 0): 1.00, (0, 1): 3.00, (1, 0): 3.00},
        }

    def get_ratios(self, target_pair: tuple) -> dict:
        return self._ratios_table.get(target_pair, {})

class ZeroOffDiagonalModel(PenaltyModel):
    """
    [기본 모델] 개별 q_{ij} 기댓값 0 모델.

    === 어떤 기댓값이 0인가 ===
    - off-diagonal E[q_ij] = 0: 모든 target pair에 대해 보장 ✓
    - diagonal E[q_ii] = 0: 보장 안 됨 ✗ (b_i=0이면 양수, b_i=1이면 음수)

    === 수학적 도출 ===
    q_ij 부호 패턴 (penalty state → q_ij 기여):
      (0,0)→+r, (0,1)→-r, (1,0)→-r, (1,1)→+r

    target (b_i,b_j)에 대해 3개 오답의 부호합이 0이려면:
      target (0,0): -R(01) - R(10) + R(11) = 0  →  R(11) = R(01) + R(10)
      target (0,1): +R(00) - R(10) + R(11) = 0  →  R(10) = R(00) + R(11)
      target (1,0): +R(00) - R(01) + R(11) = 0  →  R(01) = R(00) + R(11)
      target (1,1): +R(00) - R(01) - R(10) = 0  →  R(00) = R(01) + R(10)
    공통 규칙: double-flip ratio = single-flip ratio 합.

    대칭 single-flip (R_a = R_b = 1) → double-flip = 2. 비율 {1, 1, 2}.

    === Indistinguishability ===
    장점: |계수|의 multiset이 항상 {μ·1, μ·1, μ·2}이므로 q_ij 분포가
          target pair와 완전히 독립 (모든 모멘트 동일). Off-diagonal로는 구별 불가.
    단점: diagonal E[q_kk] ≠ 0 (부호가 b_k를 노출). 탐지 SNR ~ √n.

    === Minimax 최적성 ===
    off-diagonal zero 제약 하에서 비대칭 single-flip 비율(a ≠ b)을 허용하면
    diagonal bias가 변수 위치(i)에 의존하게 되어 더 나빠짐.
    따라서 대칭 {1, 1, 2}가 유일한 minimax 최적해.
    """
    def __init__(self):
        self._ratios_table = {
            # double-flip 상태에 2배, single-flip 상태에 1배
            (0, 0): {(0, 1): 1.0, (1, 0): 1.0, (1, 1): 2.0},
            (0, 1): {(0, 0): 1.0, (1, 0): 2.0, (1, 1): 1.0},
            (1, 0): {(0, 0): 1.0, (0, 1): 2.0, (1, 1): 1.0},
            (1, 1): {(0, 0): 2.0, (0, 1): 1.0, (1, 0): 1.0},
        }

    def get_ratios(self, target_pair: tuple) -> dict:
        return self._ratios_table.get(target_pair, {})


class BalancedModel(PenaltyModel):
    """
    [DEPRECATED] Minimax 균형 모델.

    원래 목적: max(|off_bias|, |diag_bias|)를 최소화.
    문제점: target (0,0)만 off-diagonal 0이고 나머지는 깨짐:
      target (0,0): E[q_ij] = μ(-1 - 1 + 2) = 0
      target (0,1): E[q_ij] = μ(+7/6 - 1 + 1) = +7μ/6
      target (1,0): E[q_ij] = μ(+7/6 - 1 + 1) = +7μ/6
      target (1,1): E[q_ij] = μ(+1 - 4/3 - 4/3) = -5μ/3

    off-diagonal zero 제약 하에서 minimax 최적해는 ZeroOffDiagonalModel과 동일함
    (대칭 single-flip 비율이 유일한 해). 별도 모델이 불필요.
    대안: ZeroOffDiagonalModel 사용.
    """
    def __init__(self):
        warnings.warn(
            "BalancedModel은 deprecated입니다. "
            "off-diagonal E[q_ij] ≠ 0 (target (0,1),(1,0),(1,1)). "
            "ZeroOffDiagonalModel을 사용하세요.",
            DeprecationWarning,
            stacklevel=2,
        )
        self._ratios_table = {
            (0, 0): {(0, 1): 1.0,  (1, 0): 1.0,  (1, 1): 2.0},
            (0, 1): {(0, 0): 7/6,  (1, 0): 1.0,  (1, 1): 1.0},
            (1, 0): {(0, 0): 7/6,  (0, 1): 1.0,  (1, 1): 1.0},
            (1, 1): {(0, 0): 1.0,  (0, 1): 4/3,  (1, 0): 4/3},
        }

    def get_ratios(self, target_pair: tuple) -> dict:
        return self._ratios_table.get(target_pair, {})


class SimpleUniformModel(PenaltyModel):
    """
    대체 모델: 비교를 위한 단순 균등 페널티 모델입니다.
    모든 오답에 대해 동일한 페널티를 부여합니다.
    """
    def get_ratios(self, target_pair: tuple) -> dict:
        # 모든 오답 상태에 1.0 부여
        all_states = [(0, 0), (0, 1), (1, 0), (1, 1)]
        return {s: 1.0 for s in all_states if s != target_pair}


def create_qubo_precise(target_binary_string, density=1.0, base_range=(1, 3), model: PenaltyModel = None,
                        penalty_mode='all', **kwargs):
    """
    정밀한 기댓값 0 QUBO 생성

    Args:
        target_binary_string: 목표 해 (예: "10110")
        density: 상호작용 추가 확률 (1.0 권장)
        base_range: 기본 r 샘플링 범위
        model: PenaltyModel 인스턴스 (None일 경우 ZeroOffDiagonalModel 사용)
        penalty_mode: 'all' = 3개 오답 전부에 페널티 (기본값)
                      'single' = 3개 중 1개를 랜덤 선택하여 페널티

    Returns:
        Q: QUBO 행렬 (딕셔너리)
    """
    if model is None:
        model = ZeroOffDiagonalModel()

    n = len(target_binary_string)
    Q = {}

    # 모든 큐빗 쌍(i, j)에 대해 반복 (i < j)
    for i in range(n):
        for j in range(i + 1, n):
            # density 확률로 상호작용 추가 여부 결정
            if random.random() > density:
                continue

            # 1. 목표 상태 확인 (이 쌍의 정답이 무엇이어야 하는지)
            bit_i = int(target_binary_string[i])
            bit_j = int(target_binary_string[j])
            target_state = (bit_i, bit_j)

            # 2. 모델을 통해 페널티 비율 가져오기 (Strategy Pattern 적용)
            ratios = model.get_ratios(target_state)

            # penalty_mode에 따라 적용할 오답 상태 결정
            if penalty_mode == 'single':
                # 3개 중 1개를 균등 확률로 선택
                chosen = random.choice(list(ratios.items()))
                penalty_items = [chosen]
            else:
                penalty_items = ratios.items()

            for penalty_state, ratio in penalty_items:
                # 3. 페널티 강도(r) 결정
                r = random.uniform(*base_range) * ratio
                
                s_i, s_j = penalty_state
                
                # 4. 페널티 항을 QUBO 행렬에 분배
                # 목표: 상태가 (s_i, s_j)일 때 에너지가 r만큼 증가해야 함
                
                if s_i == 0 and s_j == 0:
                    # Case 00: 페널티 식 = r * (1-x_i)(1-x_j)
                    # 전개: r(1 - x_i - x_j + x_i*x_j)
                    Q[(i, i)] = Q.get((i, i), 0) - r
                    Q[(j, j)] = Q.get((j, j), 0) - r
                    Q[(i, j)] = Q.get((i, j), 0) + r
                    
                elif s_i == 0 and s_j == 1:
                    # Case 01: 페널티 식 = r * (1-x_i)x_j
                    # 전개: r(x_j - x_i*x_j)
                    Q[(j, j)] = Q.get((j, j), 0) + r
                    Q[(i, j)] = Q.get((i, j), 0) - r
                    
                elif s_i == 1 and s_j == 0:
                    # Case 10: 페널티 식 = r * x_i(1-x_j)
                    # 전개: r(x_i - x_i*x_j)
                    Q[(i, i)] = Q.get((i, i), 0) + r
                    Q[(i, j)] = Q.get((i, j), 0) - r
                    
                else: # s_i == 1 and s_j == 1
                    # Case 11: 페널티 식 = r * x_i*x_j
                    # 전개: r*x_i*x_j
                    Q[(i, j)] = Q.get((i, j), 0) + r
    
    # 행/열 합 0으로 수렴하게 만들기 (옵션)
    if kwargs.get('balance_rows', False):
        # Ising 모델 기반으로 재생성 (기존 Q 덮어쓰기)
        # 이 방식은 E[RowSum] = 0을 만족하면서 Ground State를 완벽히 보존함
        Q = create_qubo_ising_derived(target_binary_string, density=density, base_range=base_range)
        
    return Q


def create_qubo_ising_derived(target, density=1.0, base_range=(1, 3)):
    """
    Ising 모델 (H = -sum J_ij s_i s_j, h_i=0)에서 유도된 QUBO.
    특징:
    1. Ground State가 Target과 일치함이 수학적으로 보장됨.
    2. QUBO 변환 시 Expected Row Sum이 0임.
    """
    n = len(target)
    Q = {}
    
    # Ising Spin 변환 (-1, 1)
    spins = [1 if b == '1' else -1 for b in target]
    
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() > density:
                continue
                
            # Ising Interaction J_ij
            # Target 상태에서 에너지가 최소화되려면 J_ij * s_i * s_j > 0 이어야 함 (H = - sum ...)
            # 즉 J_ij = alpha * s_i * s_j
            alpha = random.uniform(*base_range)
            J_ij = alpha * spins[i] * spins[j]
            
            # Convert Ising J to QUBO
            # H_term = - J_ij * s_i * s_j
            # s = 2x - 1
            # s_i s_j = (2x_i - 1)(2x_j - 1) = 4x_ix_j - 2x_i - 2x_j + 1
            # H_term = - J_ij (4x_ix_j - 2x_i - 2x_j + 1)
            #        = -4 J_ij x_ix_j + 2 J_ij x_i + 2 J_ij x_j - J_ij
            
            q_ij = -4 * J_ij
            q_i = 2 * J_ij
            q_j = 2 * J_ij
            
            Q[(i, j)] = Q.get((i, j), 0) + q_ij
            Q[(i, i)] = Q.get((i, i), 0) + q_i
            Q[(j, j)] = Q.get((j, j), 0) + q_j
            
    return Q


def create_qubo_2sat(target_binary_string, base_range=(1, 3), model=None, seed=None):
    """
    2-SAT 기반 incremental QUBO 생성 (single penalty mode).

    각 변수 쌍에 대해 3개 오답 중 1개만 선택하여 2-SAT clause로 추가.
    유일해가 보장될 때까지 반복 후, model ratio를 적용하여 QUBO 변환.

    E[q_ij] = 0 보존: ZeroOffDiagonalModel ratio {1,1,2}의 부호합이 0이므로,
    균등 선택 시 E[q_ij] = (1/3)·Σ(sign_k·ratio_k)·μ = 0.

    Args:
        target_binary_string: 목표 해 (예: "10110")
        base_range: 기본 r 샘플링 범위
        model: PenaltyModel 인스턴스 (None → ZeroOffDiagonalModel)
        seed: 난수 시드

    Returns:
        Q: QUBO 행렬 (딕셔너리)
        info: 메타정보 딕셔너리
    """
    from posiform.qubo_posiform import solve_2sat, check_2sat_uniqueness

    if model is None:
        model = ZeroOffDiagonalModel()

    n = len(target_binary_string)
    if n < 2:
        raise ValueError(f"변수 수가 2 미만입니다: n={n}")

    rng = random.Random(seed)
    target_bits = [int(b) for b in target_binary_string]
    target_sol = [b == 1 for b in target_bits]

    clauses = []       # 2-SAT clauses: [((var, is_positive), (var, is_positive)), ...]
    clause_info = []   # [(i, j, wrong_state), ...] — QUBO 변환용
    added = set()      # 중복 clause 방지

    max_iterations = 10 * n
    check_interval = max(1, n // 4)

    def wrong_to_clause(i, j, wrong):
        """Wrong state (s_i, s_j) → 2-SAT clause.
        wrong (s_i, s_j)를 배제: (x_i ≠ s_i) ∨ (x_j ≠ s_j)
        """
        lit1 = (i, wrong[0] == 0)  # s_i=0 → x_i(positive), s_i=1 → ¬x_i(negative)
        lit2 = (j, wrong[1] == 0)
        return (lit1, lit2)

    # --- Phase 1: 랜덤 clause 추가 ---
    is_unique = False
    for iteration in range(max_iterations):
        i, j = sorted(rng.sample(range(n), 2))
        ti, tj = target_bits[i], target_bits[j]
        target_pair = (ti, tj)

        # 3개 오답 중 1개 균등 선택
        wrong_states = [(vi, vj) for vi in (0, 1) for vj in (0, 1)
                        if (vi, vj) != target_pair]
        wrong = rng.choice(wrong_states)

        # 2-SAT clause 변환 + 중복 스킵
        clause = wrong_to_clause(i, j, wrong)
        if clause in added:
            continue
        added.add(clause)
        clauses.append(clause)
        clause_info.append((i, j, wrong))

        # 주기적 유일성 확인
        if len(clauses) % check_interval == 0:
            is_unique_check, flippable = check_2sat_uniqueness(n, clauses, target_sol)
            if is_unique_check:
                is_unique = True
                break

    # --- Phase 2: Targeted phase (flippable 변수 집중) ---
    if not is_unique:
        is_unique_check, flippable = check_2sat_uniqueness(n, clauses, target_sol)
        if is_unique_check:
            is_unique = True
        else:
            max_targeted = 3 * n
            for _ in range(max_targeted):
                if not flippable:
                    break
                fi = rng.choice(flippable)
                others = [k for k in range(n) if k != fi]
                partner = rng.choice(others)
                i, j = min(fi, partner), max(fi, partner)

                ti, tj = target_bits[i], target_bits[j]
                target_pair = (ti, tj)
                wrong_states = [(vi, vj) for vi in (0, 1) for vj in (0, 1)
                                if (vi, vj) != target_pair]
                rng.shuffle(wrong_states)

                for wrong in wrong_states:
                    clause = wrong_to_clause(i, j, wrong)
                    if clause not in added:
                        added.add(clause)
                        clauses.append(clause)
                        clause_info.append((i, j, wrong))
                        break

                if len(clauses) % max(1, len(flippable)) == 0:
                    is_unique_check, flippable = check_2sat_uniqueness(
                        n, clauses, target_sol)
                    if is_unique_check:
                        is_unique = True
                        break

    # --- QUBO 변환: clause_info → Q matrix (model ratio 적용) ---
    Q = {}
    for (i, j, wrong) in clause_info:
        ti, tj = target_bits[i], target_bits[j]
        target_pair = (ti, tj)
        ratios = model.get_ratios(target_pair)
        ratio = ratios[wrong]
        r = rng.uniform(*base_range) * ratio

        s_i, s_j = wrong
        if s_i == 0 and s_j == 0:
            # r * (1-x_i)(1-x_j) = r - r*x_i - r*x_j + r*x_i*x_j
            Q[(i, i)] = Q.get((i, i), 0) - r
            Q[(j, j)] = Q.get((j, j), 0) - r
            Q[(i, j)] = Q.get((i, j), 0) + r
        elif s_i == 0 and s_j == 1:
            # r * (1-x_i)*x_j = r*x_j - r*x_i*x_j
            Q[(j, j)] = Q.get((j, j), 0) + r
            Q[(i, j)] = Q.get((i, j), 0) - r
        elif s_i == 1 and s_j == 0:
            # r * x_i*(1-x_j) = r*x_i - r*x_i*x_j
            Q[(i, i)] = Q.get((i, i), 0) + r
            Q[(i, j)] = Q.get((i, j), 0) - r
        else:  # s_i == 1 and s_j == 1
            # r * x_i*x_j
            Q[(i, j)] = Q.get((i, j), 0) + r

    info = {
        'n': n,
        'num_clauses': len(clauses),
        'is_unique': is_unique,
        'model': model.__class__.__name__,
    }

    return Q, info


def balance_qubo_rows(Q, n):
    """Deprecated: 단순 보정은 Ground State를 깸. Ising Derived 사용 권장."""
    # 각 변수별 총 합 계산 (대각 + 비대각)
    # Q는 상삼각 형태지만 논리적으로 대칭이므로 i와 연결된 모든 j를 고려해야 함
    row_sums = {i: 0.0 for i in range(n)}
    
    for (i, j), weight in Q.items():
        if i == j:
            row_sums[i] += weight
        else:
            # i < j 형태
            row_sums[i] += weight
            row_sums[j] += weight
            
    # 대각 성분 보정
    for i in range(n):
        # 현재 합이 S_i라면, Q_ii에서 S_i를 빼면 새로운 합은 0이 됨
        # (Q_ii - S_i) + others = Q_ii + others - S_i = S_i - S_i = 0
        correction = row_sums[i]
        Q[(i, i)] = Q.get((i, i), 0) - correction
        
    return Q


def solve_brute_force(Q, n):
    """브루트포스 풀이"""
    if n > 20:
        return None, None, None

    best_energy = float('inf')
    best_solution = None
    all_results = []
    
    for bits in itertools.product([0, 1], repeat=n):
        current_state = "".join(map(str, bits))
        energy = calculate_energy(current_state, Q)
        all_results.append((current_state, energy))
        
        if energy < best_energy:
            best_energy = energy
            best_solution = current_state
    
    return best_solution, best_energy, all_results


def batch_test(num_tests=50, n_bits=10):
    """배치 테스트"""
    print(f"\n배치 테스트: {num_tests}회, {n_bits}비트")
    print("-" * 40)
    
    successes = 0
    
    for i in range(num_tests):
        target = ''.join(random.choice('01') for _ in range(n_bits))
        Q = create_qubo_precise(target, density=1.0)
        found, min_e, _ = solve_brute_force(Q, n_bits)
        
        if found == target:
            successes += 1
    
    print(f"성공률: {successes}/{num_tests} ({100*successes/num_tests:.1f}%)")
    return successes / num_tests


def large_scale_analysis(n_tests=500, n_bits=30, model=None):
    """
    대규모 계수 분석 (브루트포스 없이)

    diagonal과 off-diagonal 기댓값을 target 비트/pair 조건별로 분리 검증.
    off-diagonal이 target pair별로 동일한지 확인하는 것이 핵심.
    """
    if model is None:
        model = ZeroOffDiagonalModel()

    print(f"\n대규모 계수 분석: {n_tests}회, {n_bits}비트, 모델: {model.__class__.__name__}")
    print("-" * 40)

    # 대각 항: b_i 조건별
    xi_bi0 = []
    xi_bi1 = []
    # 비대각 항: 전체 + target pair 조건별
    xixj_all = []
    offdiag_by_target = {(0, 0): [], (0, 1): [], (1, 0): [], (1, 1): []}

    for _ in range(n_tests):
        target = ''.join(random.choice('01') for _ in range(n_bits))
        Q = create_qubo_precise(target, density=1.0, model=model)

        for (i, j), weight in Q.items():
            if i == j:
                if target[i] == '0':
                    xi_bi0.append(weight)
                else:
                    xi_bi1.append(weight)
            else:
                xixj_all.append(weight)
                target_pair = (int(target[i]), int(target[j]))
                offdiag_by_target[target_pair].append(weight)

    print(f"\n[대각 항] E[q_ii | b_i 조건별]:")
    print(f"  E[q_ii | b_i=0] = {np.mean(xi_bi0):.6f} (σ={np.std(xi_bi0):.4f})")
    print(f"  E[q_ii | b_i=1] = {np.mean(xi_bi1):.6f} (σ={np.std(xi_bi1):.4f})")

    print(f"\n[비대각 항] E[q_ij] 전체:")
    print(f"  E[q_ij] = {np.mean(xixj_all):.6f} (σ={np.std(xixj_all):.4f})")

    print(f"\n[비대각 항] E[q_ij | target pair별] (모두 ≈0이면 indistinguishable):")
    for pair in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        values = offdiag_by_target[pair]
        if values:
            print(f"  E[q_ij | target={pair}] = {np.mean(values):.6f} (σ={np.std(values):.4f}, n={len(values)})")
        else:
            print(f"  E[q_ij | target={pair}] = (데이터 없음)")


def compare_with_random_qubo(n_bits=20, n_tests=100):
    """
    생성된 QUBO와 순수 랜덤 QUBO의 계수 분포 비교
    """
    print(f"\n생성 QUBO vs 랜덤 QUBO 비교 ({n_bits}비트, {n_tests}회)")
    print("-" * 50)
    
    # 생성된 QUBO 계수 수집
    generated_diag = []
    generated_offdiag = []
    
    for _ in range(n_tests):
        target = ''.join(random.choice('01') for _ in range(n_bits))
        Q = create_qubo_precise(target, density=1.0)
        
        for (i, j), weight in Q.items():
            if i == j:
                generated_diag.append(weight)
            else:
                generated_offdiag.append(weight)
    
    # 순수 랜덤 QUBO 생성
    random_diag = []
    random_offdiag = []
    
    for _ in range(n_tests):
        for i in range(n_bits):
            # 대각: uniform(-10, 10)
            random_diag.append(random.uniform(-10, 10))
        
        for i in range(n_bits):
            for j in range(i+1, n_bits):
                # 비대각: uniform(-10, 10)
                random_offdiag.append(random.uniform(-10, 10))
    
    print("\n대각 항 (x_i 계수):")
    print(f"  생성 QUBO: 평균={np.mean(generated_diag):.4f}, 표준편차={np.std(generated_diag):.4f}")
    print(f"  랜덤 QUBO: 평균={np.mean(random_diag):.4f}, 표준편차={np.std(random_diag):.4f}")
    
    print("\n비대각 항 (x_ix_j 계수):")
    print(f"  생성 QUBO: 평균={np.mean(generated_offdiag):.4f}, 표준편차={np.std(generated_offdiag):.4f}")
    print(f"  랜덤 QUBO: 평균={np.mean(random_offdiag):.4f}, 표준편차={np.std(random_offdiag):.4f}")
    
    print("\n* 두 분포가 비슷할수록 구별하기 어려움")


def generate_dataset(num_problems=100, n_bits=50, output_dir='/home/claude'):
    """
    벤치마크 데이터셋 생성
    """
    import os
    
    print(f"\n데이터셋 생성: {num_problems}개, {n_bits}비트")
    print("-" * 40)
    
    dataset_path = os.path.join(output_dir, 'qubo_dataset')
    os.makedirs(dataset_path, exist_ok=True)
    
    metadata = []
    
    for i in range(num_problems):
        target = ''.join(random.choice('01') for _ in range(n_bits))
        Q = create_qubo_precise(target, density=1.0)
        
        filename = f'qubo_{i:04d}.csv'
        filepath = os.path.join(dataset_path, filename)
        save_qubo_edgelist(Q, filepath, target)
        
        metadata.append({
            'id': i,
            'filename': filename,
            'target': target,
            'n_bits': n_bits,
            'n_terms': len(Q)
        })
    
    # 메타데이터 저장
    meta_path = os.path.join(dataset_path, 'metadata.csv')
    with open(meta_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['id', 'filename', 'target', 'n_bits', 'n_terms'])
        writer.writeheader()
        writer.writerows(metadata)
    
    print(f"저장 완료: {dataset_path}")
    print(f"  - {num_problems}개의 QUBO 파일")
    print(f"  - metadata.csv")
    
    return dataset_path


if __name__ == "__main__":
    random.seed(42)
    
    print("=" * 60)
    print("기댓값 0 QUBO 생성기 (정밀 버전)")
    print("=" * 60)
    
    # 0. 단일 예제 상세 출력
    print("\n" + "=" * 60)
    print("단일 예제 상세 출력")
    print("=" * 60)
    
    target = "11000101010001101"
    balance_mode = False
    
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg.isdigit():
            # 숫자로 입력되면 해당 길이의 랜덤 타겟 생성
            length = int(arg)
            target = "".join(str(random.randint(0, 1)) for _ in range(length))
            print(f"[설정] 길이 {length}의 랜덤 목표 해를 생성했습니다.")
        elif arg == "balance":
            print("[설정] 행/열 평균 0 수렴 모드 활성화")
            balance_mode = True
        else:
            # 이진 문자열로 직접 입력된 경우
            target = arg
            
        # 추가 인자 확인 (순서 무관하게 balance 체크)
        if "balance" in sys.argv:
            balance_mode = True
            print("[설정] 행/열 평균 0 수렴 모드 활성화")
        
    print(f"\n목표 해: {target}")
    
    Q = create_qubo_precise(target, density=1.0, balance_rows=balance_mode)
    
    # Q 행렬 출력
    print_q_matrix(Q, len(target))
    
    # 수식 출력
    print_qubo_formula(Q)
    
    # 검증
    found, min_energy, all_results = solve_brute_force(Q, len(target))
    target_energy = calculate_energy(target, Q)
    
    print(f"\n검증 결과:")
    print(f"  목표 해 '{target}'의 에너지: {target_energy}")
    
    if found is None:
        print(f"  브루트포스 최소 해: (N > 20 이므로 생략됨)")
        print("  - 검증 생략됨")
    else:
        print(f"  브루트포스 최소 해 '{found}'의 에너지: {min_energy}")
        
        if found == target:
            print("  ✓ 성공: 목표 해가 최소값!")
        else:
            print("  ✗ 실패")

    # 결과 저장
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(output_dir, exist_ok=True)
    
    # 파일명 길이가 너무 길어지지 않도록 조정
    filename_target = f"{target[:5]}_{len(target)}"
    
    output_file = os.path.join(output_dir, f"qubo_{filename_target}.txt")
    
    save_qubo_edgelist(Q, output_file, target)
    print(f"\n[저장 완료] 결과 파일: {output_file}")
    
    # 에너지 순위
    if all_results:
        print("\n  에너지 순위 (상위 5개):")
        sorted_results = sorted(all_results, key=lambda x: x[1])
        for rank, (state, energy) in enumerate(sorted_results[:5], 1):
            marker = " <- 목표" if state == target else ""
            print(f"    {rank}. {state}: {energy:.4f}{marker}")
    else:
         print("\n  에너지 순위: (생략됨)")
    
    # 1. 배치 테스트
    batch_test(num_tests=30, n_bits=10)
    
    # 2. 대규모 계수 분석
    large_scale_analysis(n_tests=500, n_bits=30)
    
    # 3. 랜덤 QUBO와 비교
    compare_with_random_qubo(n_bits=20, n_tests=100)
    
    # 4. 데이터셋 생성 (예시)
    # generate_dataset(num_problems=10, n_bits=50)
