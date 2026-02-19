"""
Zero-Expectation vs Wishart Q 행렬 구조 비교 분석
"""

import numpy as np
import random

from qubo_zero_expectation import create_qubo_precise, DefaultZeroExpectationModel
from qubo_wishart import create_qubo_wishart


def to_matrix(Q, n):
    """Q dict를 대칭 numpy 행렬로 변환."""
    M = np.zeros((n, n))
    for (i, j), v in Q.items():
        M[i][j] = v
        if i != j:
            M[j][i] = v
    return M


def analyze_q(name, Q, n, target):
    """Q 행렬의 구조적 특성 분석."""
    M = to_matrix(Q, n)
    spins = [1 if b == '1' else -1 for b in target]

    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    if n <= 10:
        print(f"\nQ 행렬:")
        print(np.round(M, 2))

    # 1. 대각 항 vs target 비트 상관관계
    diag = np.diag(M)
    diag_b0 = [diag[i] for i in range(n) if target[i] == '0']
    diag_b1 = [diag[i] for i in range(n) if target[i] == '1']

    print(f"\n[1] 대각 항 (Q_ii) vs target 비트")
    print(f"  target: {target}")
    print(f"  대각:   {np.round(diag, 2).tolist()}")
    print(f"  b_i=0인 Q_ii 평균: {np.mean(diag_b0):.4f}")
    print(f"  b_i=1인 Q_ii 평균: {np.mean(diag_b1):.4f}")
    print(f"  차이 (편향):        {np.mean(diag_b1) - np.mean(diag_b0):.4f}")

    # 2. 비대각 항 상관 구조
    offdiag = []
    for i in range(n):
        for j in range(i+1, n):
            offdiag.append(M[i][j])
    offdiag = np.array(offdiag)

    print(f"\n[2] 비대각 항 (Q_ij, i<j)")
    print(f"  평균: {np.mean(offdiag):.4f}")
    print(f"  표준편차: {np.std(offdiag):.4f}")

    # 3. 고유값 분석 → rank 구조
    eigvals = np.linalg.eigvalsh(M)
    eigvals_sorted = np.sort(np.abs(eigvals))[::-1]
    rank = np.sum(np.abs(eigvals) > 0.01)

    print(f"\n[3] 고유값 분석")
    print(f"  고유값: {np.round(sorted(eigvals), 3).tolist()}")
    print(f"  유효 rank (|λ|>0.01): {rank} / {n}")
    print(f"  최대 고유값: {max(eigvals):.4f}")
    print(f"  최소 고유값: {min(eigvals):.4f}")

    # 상위 고유값이 전체 에너지의 몇 %를 차지하는지
    total_var = np.sum(eigvals**2)
    top1_var = eigvals_sorted[0]**2 / total_var * 100
    top3_var = np.sum(eigvals_sorted[:min(3,n)]**2) / total_var * 100

    print(f"  상위 1개 고유값 에너지 비중: {top1_var:.1f}%")
    print(f"  상위 3개 고유값 에너지 비중: {top3_var:.1f}%")

    # 4. 행 합 분석
    row_sums = np.sum(M, axis=1)
    row_b0 = [row_sums[i] for i in range(n) if target[i] == '0']
    row_b1 = [row_sums[i] for i in range(n) if target[i] == '1']

    print(f"\n[4] 행 합 (row sum)")
    print(f"  전체 평균: {np.mean(row_sums):.4f}")
    print(f"  b_i=0 행 합 평균: {np.mean(row_b0):.4f}")
    print(f"  b_i=1 행 합 평균: {np.mean(row_b1):.4f}")

    return M


# ==============================================================
# 작은 예제 (N=5) - 구조를 직접 볼 수 있음
# ==============================================================
print("\n" + "#"*60)
print("# Part 1: 작은 예제 (N=5)")
print("#"*60)

random.seed(42)
np.random.seed(42)
target = "10110"

Q_ze = create_qubo_precise(target, density=1.0, model=DefaultZeroExpectationModel())
Q_w = create_qubo_wishart(target, alpha=0.7, seed=42)

analyze_q("Zero-Expectation (N=5)", Q_ze, 5, target)
analyze_q("Wishart alpha=0.7 (N=5)", Q_w, 5, target)


# ==============================================================
# 대규모 통계 (N=100, 50회 반복)
# ==============================================================
print("\n\n" + "#"*60)
print("# Part 2: 대규모 통계 (N=100, 50회 반복)")
print("#"*60)

n = 100
n_trials = 50

ze_diag_b0, ze_diag_b1 = [], []
ze_offdiag = []
ze_ranks = []

w_diag_b0, w_diag_b1 = [], []
w_offdiag = []
w_ranks = []

for trial in range(n_trials):
    target = ''.join(str(random.randint(0, 1)) for _ in range(n))

    # Zero-Expectation
    Q = create_qubo_precise(target, density=1.0, model=DefaultZeroExpectationModel())
    M = to_matrix(Q, n)
    diag = np.diag(M)
    eigvals = np.linalg.eigvalsh(M)

    for i in range(n):
        if target[i] == '0':
            ze_diag_b0.append(diag[i])
        else:
            ze_diag_b1.append(diag[i])

    for i in range(n):
        for j in range(i+1, n):
            ze_offdiag.append(M[i][j])

    ze_ranks.append(np.sum(np.abs(eigvals) > 0.01))

    # Wishart
    Q = create_qubo_wishart(target, alpha=0.7)
    M = to_matrix(Q, n)
    diag = np.diag(M)
    eigvals = np.linalg.eigvalsh(M)

    for i in range(n):
        if target[i] == '0':
            w_diag_b0.append(diag[i])
        else:
            w_diag_b1.append(diag[i])

    for i in range(n):
        for j in range(i+1, n):
            w_offdiag.append(M[i][j])

    w_ranks.append(np.sum(np.abs(eigvals) > 0.01))


print(f"\n{'='*70}")
print(f"  대규모 통계 비교 (N={n}, {n_trials}회)")
print(f"{'='*70}")

print(f"\n[1] 대각 항 편향 (target 비트 노출 여부)")
print(f"  {'':30s} {'E[Q_ii|b=0]':>12s} {'E[Q_ii|b=1]':>12s} {'차이(편향)':>12s}")
print(f"  {'Zero-Expectation':30s} {np.mean(ze_diag_b0):>12.4f} {np.mean(ze_diag_b1):>12.4f} {np.mean(ze_diag_b1)-np.mean(ze_diag_b0):>12.4f}")
print(f"  {'Wishart (alpha=0.7)':30s} {np.mean(w_diag_b0):>12.4f} {np.mean(w_diag_b1):>12.4f} {np.mean(w_diag_b1)-np.mean(w_diag_b0):>12.4f}")

print(f"\n[2] 비대각 항 통계")
print(f"  {'':30s} {'E[Q_ij]':>12s} {'Std[Q_ij]':>12s}")
print(f"  {'Zero-Expectation':30s} {np.mean(ze_offdiag):>12.4f} {np.std(ze_offdiag):>12.4f}")
print(f"  {'Wishart (alpha=0.7)':30s} {np.mean(w_offdiag):>12.4f} {np.std(w_offdiag):>12.4f}")

print(f"\n[3] 유효 Rank (|λ|>0.01)")
print(f"  {'Zero-Expectation':30s} 평균: {np.mean(ze_ranks):.1f} / {n}")
print(f"  {'Wishart (alpha=0.7)':30s} 평균: {np.mean(w_ranks):.1f} / {n}")

print(f"\n[4] 비대각 항 상관 (Q_ij와 Q_ik의 상관계수)")
# 같은 행에 속하는 비대각 원소들 간 상관 측정
ze_corrs = []
w_corrs = []

for trial in range(20):
    target = ''.join(str(random.randint(0, 1)) for _ in range(n))

    Q = create_qubo_precise(target, density=1.0, model=DefaultZeroExpectationModel())
    M = to_matrix(Q, n)
    # 행 0의 비대각 원소 vs 행 1의 비대각 원소 상관
    row0 = [M[0][j] for j in range(2, n)]
    row1 = [M[1][j] for j in range(2, n)]
    ze_corrs.append(np.corrcoef(row0, row1)[0, 1])

    Q = create_qubo_wishart(target, alpha=0.7)
    M = to_matrix(Q, n)
    row0 = [M[0][j] for j in range(2, n)]
    row1 = [M[1][j] for j in range(2, n)]
    w_corrs.append(np.corrcoef(row0, row1)[0, 1])

print(f"  {'Zero-Expectation':30s} 행간 상관: {np.mean(ze_corrs):.4f} (±{np.std(ze_corrs):.4f})")
print(f"  {'Wishart (alpha=0.7)':30s} 행간 상관: {np.mean(w_corrs):.4f} (±{np.std(w_corrs):.4f})")


# ==============================================================
# 결론
# ==============================================================
print(f"\n\n{'='*70}")
print("결론: 구별 가능 vs 구별 불가능의 차이")
print(f"{'='*70}")
print("""
1. 대각 편향: Wishart는 b_i=0과 b_i=1에서 Q_ii 평균이 크게 다름
   → 대각 항만 봐도 target 비트를 추측 가능
   → Zero-Expectation은 차이가 거의 0

2. Rank 구조: Wishart는 low-rank (M = alpha*N ≈ 70)
   → 고유값 분해로 즉시 탐지 가능
   → Zero-Expectation은 full-rank (≈ N)

3. 행간 상관: Wishart는 Q의 서로 다른 행이 상관됨
   → J = -(1/N) W W^T 에서 W가 공유되므로
   → Zero-Expectation은 각 (i,j) 페널티가 독립 → 행간 상관 ≈ 0

4. 핵심: 이 상관 구조가 metastable 상태를 만들어 SA를 어렵게 함
   → 상관 제거 = 구별 불가능 = SA에게 쉬움
   → 상관 유지 = 구별 가능 = SA에게 어려움
""")
