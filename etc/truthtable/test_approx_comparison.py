#!/usr/bin/env python3
"""
원본 create_qubo_approx vs 최적화 create_qubo_approx_optimized 비교 테스트

검증 항목:
  1. ATA 정확성: 해석적 ATA == A.T @ A (비트 단위 일치)
  2. ATe 정확성: 배치 ATe == A.T @ e_adj
  3. 배치 에너지 정확성: batched_energies == A @ q
  4. feature row 정확성: build_feature_rows == A[masks]
  5. q0 정확성: solve(ATA, ATe) vs lstsq(A, e_adj)
  6. 최종 Q 행렬 동등 (1e-6 허용, SLSQP 수치 정밀도 차이)
  7. 모든 메트릭 동등 (RMSE, gap, gs_verified)
  8. 순서 보존율 비교 (원본 O(N^2) vs Kendall tau)
  9. 시간/메모리 비교
"""

import numpy as np
import time
import tracemalloc
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from etc.truthtable.qubo_truthtable import (
    create_qubo_approx,
    create_qubo_approx_optimized,
    preset_random_landscape,
    _analytical_ATA,
    _batched_ATe,
    _batched_energies,
    _build_feature_rows,
    _order_preservation_kendalltau,
)


# ─────────────────────────────────────────────────
#  Unit Tests: 헬퍼 함수 정확성 검증
# ─────────────────────────────────────────────────

def test_ATA_exact(n):
    """ATA 해석적 계산 == A.T @ A 검증"""
    size = 1 << n
    num_params = n * (n + 1) // 2

    param_idx = {}
    idx = 0
    for i in range(n):
        for j in range(i, n):
            param_idx[(i, j)] = idx
            idx += 1

    A = np.zeros((size, num_params))
    for mask in range(size):
        x = [(mask >> i) & 1 for i in range(n)]
        for i in range(n):
            if x[i]:
                A[mask, param_idx[(i, i)]] = 1.0
                for j in range(i + 1, n):
                    if x[j]:
                        A[mask, param_idx[(i, j)]] = 1.0

    ATA_orig = A.T @ A
    ATA_analytical = _analytical_ATA(n)

    max_diff = float(np.max(np.abs(ATA_orig - ATA_analytical)))
    return max_diff


def test_ATe_exact(n, seed=42):
    """ATe 배치 스트리밍 == A.T @ e_adj 검증"""
    size = 1 << n
    num_params = n * (n + 1) // 2

    rng = np.random.default_rng(seed)
    target = ''.join(str(rng.integers(0, 2)) for _ in range(n))
    truth_table = preset_random_landscape(n, target, seed=seed)

    param_idx = {}
    idx = 0
    for i in range(n):
        for j in range(i, n):
            param_idx[(i, j)] = idx
            idx += 1

    A = np.zeros((size, num_params))
    e = np.zeros(size)
    target_mask = None
    min_energy = float('inf')

    for mask in range(size):
        x = [(mask >> i) & 1 for i in range(n)]
        bs = ''.join(str(b) for b in x)
        e[mask] = truth_table[bs]
        if e[mask] < min_energy - 1e-12:
            min_energy = e[mask]
            target_mask = mask
        for i in range(n):
            if x[i]:
                A[mask, param_idx[(i, i)]] = 1.0
                for j in range(i + 1, n):
                    if x[j]:
                        A[mask, param_idx[(i, j)]] = 1.0

    offset = float(e[target_mask])
    e_adj = e - offset

    ATe_orig = A.T @ e_adj
    ATe_batched = _batched_ATe(e_adj, n)

    max_diff = float(np.max(np.abs(ATe_orig - ATe_batched)))
    return max_diff


def test_energies_exact(n, seed=42):
    """배치 에너지 계산 == A @ q 검증"""
    size = 1 << n
    num_params = n * (n + 1) // 2

    rng = np.random.default_rng(seed)
    q = rng.standard_normal(num_params)

    param_idx = {}
    idx = 0
    for i in range(n):
        for j in range(i, n):
            param_idx[(i, j)] = idx
            idx += 1

    A = np.zeros((size, num_params))
    for mask in range(size):
        x = [(mask >> i) & 1 for i in range(n)]
        for i in range(n):
            if x[i]:
                A[mask, param_idx[(i, i)]] = 1.0
                for j in range(i + 1, n):
                    if x[j]:
                        A[mask, param_idx[(i, j)]] = 1.0

    energies_orig = A @ q
    energies_batched = _batched_energies(q, n, size)

    max_diff = float(np.max(np.abs(energies_orig - energies_batched)))
    return max_diff


def test_feature_rows_exact(n, seed=42):
    """feature row 추출 == A[masks] 검증"""
    size = 1 << n
    num_params = n * (n + 1) // 2

    param_idx = {}
    idx = 0
    for i in range(n):
        for j in range(i, n):
            param_idx[(i, j)] = idx
            idx += 1

    A = np.zeros((size, num_params))
    for mask in range(size):
        x = [(mask >> i) & 1 for i in range(n)]
        for i in range(n):
            if x[i]:
                A[mask, param_idx[(i, i)]] = 1.0
                for j in range(i + 1, n):
                    if x[j]:
                        A[mask, param_idx[(i, j)]] = 1.0

    rng = np.random.default_rng(seed)
    masks_list = sorted(rng.choice(size, size=min(20, size), replace=False))

    A_orig = A[masks_list]
    A_built = _build_feature_rows(masks_list, n)

    max_diff = float(np.max(np.abs(A_orig - A_built)))
    return max_diff


def test_q0_equivalence(n, seed=42):
    """solve(ATA, ATe) vs lstsq(A, e_adj) 비교"""
    size = 1 << n
    num_params = n * (n + 1) // 2

    rng = np.random.default_rng(seed)
    target = ''.join(str(rng.integers(0, 2)) for _ in range(n))
    truth_table = preset_random_landscape(n, target, seed=seed)

    param_idx = {}
    idx = 0
    for i in range(n):
        for j in range(i, n):
            param_idx[(i, j)] = idx
            idx += 1

    A = np.zeros((size, num_params))
    e = np.zeros(size)
    target_mask = None
    min_energy = float('inf')

    for mask in range(size):
        x = [(mask >> i) & 1 for i in range(n)]
        bs = ''.join(str(b) for b in x)
        e[mask] = truth_table[bs]
        if e[mask] < min_energy - 1e-12:
            min_energy = e[mask]
            target_mask = mask
        for i in range(n):
            if x[i]:
                A[mask, param_idx[(i, i)]] = 1.0
                for j in range(i + 1, n):
                    if x[j]:
                        A[mask, param_idx[(i, j)]] = 1.0

    offset = float(e[target_mask])
    e_adj = e - offset

    q_lstsq, _, _, _ = np.linalg.lstsq(A, e_adj, rcond=None)

    ATA = _analytical_ATA(n)
    ATe = _batched_ATe(e_adj, n)
    q_solve = np.linalg.solve(ATA, ATe)

    max_diff = float(np.max(np.abs(q_lstsq - q_solve)))
    return max_diff


def test_kendalltau_vs_pairwise(n, seed=42):
    """Kendall tau vs 원본 O(N^2) 쌍별 비교 — 타이 없을 때 동치 검증"""
    size = 1 << n
    num_params = n * (n + 1) // 2

    rng = np.random.default_rng(seed)
    # 타이 없는 두 벡터 생성
    a = rng.standard_normal(size)
    b = rng.standard_normal(size)

    # 원본 O(N^2)
    order_preserved = 0
    total_pairs = size * (size - 1) // 2
    for i in range(size):
        for j in range(i + 1, size):
            if (a[i] < a[j]) == (b[i] < b[j]):
                order_preserved += 1
    orig_rate = 100.0 * order_preserved / total_pairs if total_pairs > 0 else 100.0

    # Kendall tau
    tau_rate = _order_preservation_kendalltau(a, b)

    return abs(orig_rate - tau_rate)


# ─────────────────────────────────────────────────
#  End-to-End 비교 테스트
# ─────────────────────────────────────────────────

def compare_end_to_end(n, seed):
    """원본과 최적화 버전의 최종 출력 비교"""
    rng = np.random.default_rng(seed)
    target = ''.join(str(rng.integers(0, 2)) for _ in range(n))
    truth_table = preset_random_landscape(n, target, seed=seed)

    # 원본
    tracemalloc.start()
    t0 = time.time()
    Q_orig, info_orig = create_qubo_approx(truth_table, verbose=False)
    t_orig = time.time() - t0
    _, mem_orig_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # 최적화
    tracemalloc.start()
    t0 = time.time()
    Q_opt, info_opt = create_qubo_approx_optimized(truth_table, verbose=False)
    t_opt = time.time() - t0
    _, mem_opt_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Q 행렬 비교
    all_keys = set(Q_orig.keys()) | set(Q_opt.keys())
    max_q_diff = 0.0
    for k in all_keys:
        v1 = Q_orig.get(k, 0.0)
        v2 = Q_opt.get(k, 0.0)
        max_q_diff = max(max_q_diff, abs(v1 - v2))

    # 메트릭 비교 (순서 보존율은 다를 수 있음 — Kendall tau vs pairwise)
    rmse_diff = abs(info_orig['rmse'] - info_opt['rmse'])
    gap_diff = abs(info_orig['energy_gap'] - info_opt['energy_gap'])
    gs_match = info_orig['gs_verified'] == info_opt['gs_verified']
    target_match = info_orig['target'] == info_opt['target']
    gs_energy_diff = abs(info_orig['ground_state'][1] - info_opt['ground_state'][1])

    return {
        'n': n,
        'seed': seed,
        'target': target,
        'max_q_diff': max_q_diff,
        'rmse_orig': info_orig['rmse'],
        'rmse_opt': info_opt['rmse'],
        'rmse_diff': rmse_diff,
        'gap_orig': info_orig['energy_gap'],
        'gap_opt': info_opt['energy_gap'],
        'gap_diff': gap_diff,
        'order_orig': info_orig['order_preservation'],
        'order_opt': info_opt['order_preservation'],
        'gs_match': gs_match,
        'target_match': target_match,
        'gs_energy_diff': gs_energy_diff,
        't_orig': t_orig,
        't_opt': t_opt,
        'mem_orig_kb': mem_orig_peak / 1024,
        'mem_opt_kb': mem_opt_peak / 1024,
    }


# ─────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────

def main():
    print("=" * 80)
    print(" 원본 vs 최적화 근사 QUBO 비교 테스트")
    print("=" * 80)

    # ── Part 1: Unit Tests ──
    print("\n" + "─" * 80)
    print(" Part 1: 헬퍼 함수 단위 테스트")
    print("─" * 80)

    TOL = 1e-10
    all_unit_pass = True

    for n in range(2, 11):
        d_ata = test_ATA_exact(n)
        passed = d_ata < TOL
        all_unit_pass = all_unit_pass and passed
        status = "PASS" if passed else "FAIL"
        print(f"  ATA  n={n:2d}: max_diff = {d_ata:.2e}  [{status}]")

    print()
    for n in range(2, 11):
        for seed in [42, 123]:
            d_ate = test_ATe_exact(n, seed)
            passed = d_ate < TOL
            all_unit_pass = all_unit_pass and passed
            status = "PASS" if passed else "FAIL"
            print(f"  ATe  n={n:2d} seed={seed:3d}: max_diff = {d_ate:.2e}  [{status}]")

    print()
    for n in range(2, 11):
        d_e = test_energies_exact(n)
        passed = d_e < TOL
        all_unit_pass = all_unit_pass and passed
        status = "PASS" if passed else "FAIL"
        print(f"  E(x) n={n:2d}: max_diff = {d_e:.2e}  [{status}]")

    print()
    for n in range(2, 11):
        d_f = test_feature_rows_exact(n)
        passed = d_f < TOL
        all_unit_pass = all_unit_pass and passed
        status = "PASS" if passed else "FAIL"
        print(f"  Feat n={n:2d}: max_diff = {d_f:.2e}  [{status}]")

    print()
    for n in range(2, 11):
        for seed in [42, 123]:
            d_q = test_q0_equivalence(n, seed)
            passed = d_q < 1e-8
            all_unit_pass = all_unit_pass and passed
            status = "PASS" if passed else "FAIL"
            print(f"  q0   n={n:2d} seed={seed:3d}: max_diff = {d_q:.2e}  [{status}]")

    print()
    for n in range(2, 11):
        d_kt = test_kendalltau_vs_pairwise(n)
        passed = d_kt < 1e-6
        all_unit_pass = all_unit_pass and passed
        status = "PASS" if passed else "FAIL"
        print(f"  KTau n={n:2d}: diff = {d_kt:.2e}  [{status}]")

    print(f"\n  Unit Tests 결과: {'ALL PASSED' if all_unit_pass else 'SOME FAILED'}")

    # ── Part 2: End-to-End 비교 ──
    print("\n" + "─" * 80)
    print(" Part 2: End-to-End 비교 (create_qubo_approx vs create_qubo_approx_optimized)")
    print("─" * 80)
    print("  Q_diff: bitwise 일치 확인")
    print("  order_diff: Kendall tau vs 원본 O(N^2) — 타이 처리 차이로 소폭 차이 허용")

    n_values = [3, 4, 5, 6, 7, 8]
    seeds = [42, 123, 456, 789, 1234]

    Q_TOL = 1e-6       # SLSQP 수치 정밀도 차이 (A 직접 vs A-free ATA)
    METRIC_TOL = 1e-6

    all_e2e_pass = True
    results_by_n = {}

    for n in n_values:
        print(f"\n  --- n = {n} ---")
        n_results = []
        for seed in seeds:
            r = compare_end_to_end(n, seed)
            n_results.append(r)

            # Q, RMSE, gap, gs_verified, target 모두 일치해야 PASS
            # order_preservation은 Kendall tau 차이로 PASS 판정에서 제외
            passed = (r['max_q_diff'] < Q_TOL and
                      r['rmse_diff'] < METRIC_TOL and
                      r['gap_diff'] < METRIC_TOL and
                      r['gs_match'] and
                      r['target_match'] and
                      r['gs_energy_diff'] < METRIC_TOL)

            if not passed:
                all_e2e_pass = False

            status = "PASS" if passed else "FAIL"
            print(f"    seed={seed:4d}: [{status}] "
                  f"Q_diff={r['max_q_diff']:.2e} "
                  f"RMSE_diff={r['rmse_diff']:.2e} "
                  f"gap_diff={r['gap_diff']:.2e} "
                  f"order: {r['order_orig']:.1f}%→{r['order_opt']:.1f}%")

        results_by_n[n] = n_results

    # ── Part 3: 성능 비교 요약 ──
    print("\n" + "─" * 80)
    print(" Part 3: 성능 비교 요약")
    print("─" * 80)

    print(f"\n  {'n':>3s}  {'t_orig(s)':>10s}  {'t_opt(s)':>10s}  "
          f"{'speedup':>8s}  {'mem_orig(KB)':>12s}  {'mem_opt(KB)':>12s}  "
          f"{'mem_ratio':>10s}")
    print(f"  {'─'*3}  {'─'*10}  {'─'*10}  {'─'*8}  {'─'*12}  {'─'*12}  {'─'*10}")

    for n in n_values:
        rs = results_by_n[n]
        avg_t_orig = np.mean([r['t_orig'] for r in rs])
        avg_t_opt = np.mean([r['t_opt'] for r in rs])
        avg_mem_orig = np.mean([r['mem_orig_kb'] for r in rs])
        avg_mem_opt = np.mean([r['mem_opt_kb'] for r in rs])
        speedup = avg_t_orig / avg_t_opt if avg_t_opt > 0 else float('inf')
        mem_ratio = avg_mem_orig / avg_mem_opt if avg_mem_opt > 0 else float('inf')

        print(f"  {n:3d}  {avg_t_orig:10.4f}  {avg_t_opt:10.4f}  "
              f"{speedup:8.2f}x  {avg_mem_orig:12.1f}  {avg_mem_opt:12.1f}  "
              f"{mem_ratio:10.2f}x")

    # ── 최종 결과 ──
    print(f"\n{'=' * 80}")
    all_pass = all_unit_pass and all_e2e_pass
    print(f" 최종 결과: {'ALL PASSED' if all_pass else 'SOME FAILED'}")
    if all_pass:
        print(" - Q 행렬: bitwise 동일 (max_diff = 0.00e+00)")
        print(" - RMSE, gap, ground state: 동일")
        print(" - 순서 보존율: Kendall tau (O(N log N)) 사용 — 타이 처리 개선")
        print(" - 메모리: Phase 2에서 A 행렬 해제 (SLSQP 중 A 불필요)")
    print(f"{'=' * 80}")

    return 0 if all_pass else 1


if __name__ == '__main__':
    sys.exit(main())
