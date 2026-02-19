"""
TDD Step 1 (RED): 대각 항 기댓값 0 테스트

목표: E[Q_ii | b_i=0] ≈ 0 AND E[Q_ii | b_i=1] ≈ 0
또는 최소한: E[Q_ii | b_i=0] ≈ E[Q_ii | b_i=1] (target 비트 노출 방지)
"""

import random
import numpy as np
import sys

from qubo_zero_expectation import (
    create_qubo_precise,
    DefaultZeroExpectationModel,
    ZeroOffDiagonalModel,
    BalancedModel,
    PenaltyModel,
)


def measure_diagonal_bias(model, n=100, n_trials=200):
    """모델의 대각 편향을 측정."""
    diag_b0 = []
    diag_b1 = []
    offdiag = []

    for _ in range(n_trials):
        target = ''.join(str(random.randint(0, 1)) for _ in range(n))
        Q = create_qubo_precise(target, density=1.0, model=model)

        for i in range(n):
            val = Q.get((i, i), 0)
            if target[i] == '0':
                diag_b0.append(val)
            else:
                diag_b1.append(val)

        for (i, j), w in Q.items():
            if i != j:
                offdiag.append(w)

    return {
        'E_diag_b0': np.mean(diag_b0),
        'E_diag_b1': np.mean(diag_b1),
        'diag_bias': np.mean(diag_b1) - np.mean(diag_b0),
        'E_offdiag': np.mean(offdiag),
        'std_diag_b0': np.std(diag_b0),
        'std_diag_b1': np.std(diag_b1),
    }


def test_model(name, model):
    """모델 테스트: 대각 + 비대각 편향 측정."""
    random.seed(42)
    stats = measure_diagonal_bias(model, n=100, n_trials=200)

    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")
    print(f"  E[Q_ii | b=0] = {stats['E_diag_b0']:.4f} (±{stats['std_diag_b0']:.2f})")
    print(f"  E[Q_ii | b=1] = {stats['E_diag_b1']:.4f} (±{stats['std_diag_b1']:.2f})")
    print(f"  대각 편향 (b=1 - b=0) = {stats['diag_bias']:.4f}")
    print(f"  E[Q_ij] = {stats['E_offdiag']:.4f}")

    # 판정
    diag_ok = abs(stats['diag_bias']) < 0.5
    offdiag_ok = abs(stats['E_offdiag']) < 0.5

    print(f"\n  [판정]")
    print(f"    대각 편향 ≈ 0: {'PASS ✓' if diag_ok else 'FAIL ✗'}")
    print(f"    비대각 E[Q_ij] ≈ 0: {'PASS ✓' if offdiag_ok else 'FAIL ✗'}")

    return diag_ok, offdiag_ok


# ==============================================================
# 수학적 분석: 왜 불가능한가?
# ==============================================================
def prove_impossibility():
    """
    페널티 기반 구성에서 E[Q_ij]=0 AND E[Q_ii] 무편향이
    동시에 불가능함을 증명.
    """
    print("\n" + "#" * 60)
    print("# 수학적 분석: 동시 달성 가능성")
    print("#" * 60)

    print("""
각 페널티 상태가 Q_ii, Q_jj에 기여하는 패턴:

  페널티 상태 | Q_ii 기여 | Q_jj 기여 | Q_ij 기여
  (0,0)       | -r        | -r        | +r
  (0,1)       |  0        | +r        | -r
  (1,0)       | +r        |  0        | -r
  (1,1)       |  0        |  0        | +r

변수 정의 (12개 ratio):
  Target (0,0): r₁=(0,1), r₂=(1,0), r₃=(1,1)
  Target (0,1): r₄=(0,0), r₅=(1,0), r₆=(1,1)
  Target (1,0): r₇=(0,0), r₈=(0,1), r₉=(1,1)
  Target (1,1): r₁₀=(0,0), r₁₁=(0,1), r₁₂=(1,0)
""")

    print("=== 조건 1: E[Q_ij] = 0 (비대각 무편향) ===")
    print("  (E1) -r₁ - r₂ + r₃ = 0   → r₃ = r₁ + r₂")
    print("  (E2)  r₄ - r₅ + r₆ = 0   → r₅ = r₄ + r₆")
    print("  (E3)  r₇ - r₈ + r₉ = 0   → r₈ = r₇ + r₉")
    print("  (E4)  r₁₀ - r₁₁ - r₁₂ = 0 → r₁₀ = r₁₁ + r₁₂")

    print("\n=== 조건 2: E[Q_ii|b=0] = E[Q_ii|b=1] (대각 무편향) ===")

    print("""
변수 k의 Q_kk 기여 (이웃 m과의 페어에서):

  "i" 역할 (k < m), target = (b_k, b_m):
    (0,0): +r₂    (0,1): -r₄+r₅    (1,0): -r₇    (1,1): -r₁₀+r₁₂

  "j" 역할 (m < k), target = (b_m, b_k):
    (0,0): +r₁    (0,1): -r₄    (1,0): -r₇+r₈    (1,1): -r₁₀+r₁₁

  b_k별 합 (임의의 역할에서 동일):
    b_k=0, b_m=0: +r₂ (또는 +r₁)
    b_k=0, b_m=1: -r₄+r₅ (또는 -r₄)
    b_k=1, b_m=0: -r₇ (또는 -r₇+r₈... 역할에 따라 다름)

  "i" 역할 대각 무편향 조건:
    r₂ + (-r₄+r₅) = -r₇ + (-r₁₀+r₁₂)   ...(D1)
""")

    print("=== E[Q_ij]=0 조건을 D1에 대입 ===")
    print("""
  E2에서: r₅ = r₄ + r₆
  E4에서: r₁₀ = r₁₁ + r₁₂

  D1: r₂ - r₄ + (r₄ + r₆) = -r₇ - (r₁₁ + r₁₂) + r₁₂
      r₂ + r₆ = -r₇ - r₁₁

  즉: r₂ + r₆ + r₇ + r₁₁ = 0

  그런데 r₂, r₆, r₇, r₁₁ > 0 (양의 페널티 필수)

  ∴ 양수의 합 = 0은 불가능!
""")

    print("=" * 60)
    print("결론: E[Q_ij]=0과 대각 무편향을 동시에 달성하는 것은")
    print("      양의 페널티 제약 하에서 수학적으로 불가능하다.")
    print("=" * 60)


# ==============================================================
# 대안 탐색
# ==============================================================
def explore_alternatives():
    """가능한 대안들을 탐색."""
    print("\n\n" + "#" * 60)
    print("# 대안 탐색: 무엇이 가능한가?")
    print("#" * 60)

    print("""
불가능한 것:
  ✗ E[Q_ij] = 0 AND E[Q_ii|b=0] = E[Q_ii|b=1]
  ✗ E[Q_ij] = 0 AND E[Q_ii] = 0

가능한 것:
  ✓ E[Q_ij] = 0만 (현재 DefaultZeroExpectationModel) — 대각 편향 존재
  ✓ E[Q_ii|b=0] = E[Q_ii|b=1]만 — 비대각 편향 존재
  ✓ 대각 편향을 최소화하면서 비대각도 줄이는 타협안 (minimax)
""")


# ==============================================================
# 대안 1: 대각 무편향 모델 (비대각 편향 허용)
# ==============================================================
class DiagonalUnbiasedModel(PenaltyModel):
    """
    대각 무편향 모델: E[Q_ii|b=0] = E[Q_ii|b=1]을 보장.
    대신 E[Q_ij] ≠ 0 (비대각 편향 허용).

    도출: 대칭 가정 + D1 조건
      b = 1 + e + d 에서 e=1, d=1 → b=3
    """
    def __init__(self):
        self._ratios_table = {
            (0, 0): {(0, 1): 1.0, (1, 0): 1.0, (1, 1): 2.0},
            (0, 1): {(0, 0): 1.0, (1, 0): 1.0, (1, 1): 1.0},
            (1, 0): {(0, 0): 1.0, (0, 1): 1.0, (1, 1): 1.0},
            (1, 1): {(0, 0): 1.0, (0, 1): 3.0, (1, 0): 3.0},
        }

    def get_ratios(self, target_pair: tuple) -> dict:
        return self._ratios_table.get(target_pair, {})


# ==============================================================
# 대안 2: 양쪽 편향 최소화 (minimax 타협)
# ==============================================================
class MinimaxBiasModel(PenaltyModel):
    """
    대각 편향과 비대각 편향의 최대값을 최소화하는 타협 모델.
    완벽한 0은 달성 불가능하므로 양쪽 모두 작게 만드는 것이 목표.
    """
    def __init__(self):
        # BalancedModel과 유사하지만 대각 편향도 고려
        self._ratios_table = {
            (0, 0): {(0, 1): 1.0, (1, 0): 1.0, (1, 1): 2.0},
            (0, 1): {(0, 0): 7/6, (1, 0): 1.0, (1, 1): 1.0},
            (1, 0): {(0, 0): 7/6, (0, 1): 1.0, (1, 1): 1.0},
            (1, 1): {(0, 0): 1.0, (0, 1): 4/3, (1, 0): 4/3},
        }

    def get_ratios(self, target_pair: tuple) -> dict:
        return self._ratios_table.get(target_pair, {})


if __name__ == "__main__":
    # 불가능성 증명
    prove_impossibility()

    # 기존 모델 테스트
    print("\n\n" + "#" * 60)
    print("# 기존 모델 테스트")
    print("#" * 60)

    models = [
        ("DefaultZeroExpectation", DefaultZeroExpectationModel()),
        ("ZeroOffDiagonal", ZeroOffDiagonalModel()),
        ("Balanced", BalancedModel()),
        ("DiagonalUnbiased (새 모델)", DiagonalUnbiasedModel()),
    ]

    results = {}
    for name, model in models:
        diag_ok, offdiag_ok = test_model(name, model)
        results[name] = (diag_ok, offdiag_ok)

    # 결과 비교 테이블
    print("\n\n" + "=" * 70)
    print("모델별 비교 요약")
    print("=" * 70)
    print(f"{'모델':<30} | {'대각 무편향':>12} | {'비대각 E=0':>12}")
    print("-" * 60)
    for name, (d, o) in results.items():
        print(f"{name:<30} | {'PASS ✓' if d else 'FAIL ✗':>12} | {'PASS ✓' if o else 'FAIL ✗':>12}")

    print("\n결론: 두 조건을 동시에 PASS하는 모델은 존재할 수 없음 (수학적 증명)")
