"""
α=0 plateau에서 Mehta et al. 2025 평형 공식의 SA 예측력 검증.

배경
----
CISC 2026 논문의 claim: α=0 plateau에서 E-GSP ≈ 55.9% (Pegasus P16, lin2).
Mehta 공식: p_0^equil = 1 / (1 + (g_1/g_0) exp(-β·Δ_R))
  - g_0: GS degeneracy
  - g_1: first-excited degeneracy
  - Δ_R: spectral gap (R 단독, α=0이므로 H_total=R)
  - β: 역온도 (SA의 경우 schedule-effective β)

검증 가설
---------
H_0: SA가 sweep=1000 에서 R의 thermal equilibrium에 충분히 가깝다.
     → 인스턴스별 (g_0, g_1, Δ_R) → 공식 예측 p_0 == SA E-GSP

판정 분기
---------
1. R² > 0.9 + 평균 일치: framework 채택, β fit이 SA effective temperature
2. 평균은 맞지만 R² 낮음: 평균 통계만 유효, 인스턴스별 fit 부정확
3. 평균도 차이: equilibrium 가정 깨짐 → framework 폐기 또는 가정 수정

알려진 한계
-----------
1. α=0 → H_total = R뿐. Planted posiform의 어려움 검증 X. 순수 random QUBO 검증.
2. random 3-regular ≠ Pegasus P16 (degree=15). 토폴로지가 spectrum 통계엔 비결정적이라
   첫 검증으로 단순화 채택. 통과 시 Pegasus subgraph로 재검증 권장.
3. lin2 → Δ_R ∈ {1} 거의 고정 → fit 자유도가 g_1/g_0 분포에 집중. lin20으로 재검증 가능.
4. n_partitions=1 (KL bisection 없음) → CISC 본실험의 block-diagonal degeneracy 곱셈 구조
   (g_0 = ∏_i deg(R_i)) 검증 X. 단일 시스템에서 framework 작동 여부만 확인.
5. N=12 small. 큰 N에서도 같은 결과라는 보장 없음. 통과해도 본실험에 직접 외삽 X.

실행
----
python3 hardened_posiform/test_equilibrium_validation.py
python3 hardened_posiform/test_equilibrium_validation.py --n 14 --instances 200
"""

import os
import sys
import time
import json
import random
import argparse
from collections import Counter

import numpy as np
import networkx as nx
from scipy.optimize import curve_fit
import neal

# ─── 파라미터 ───
N_NODES_DEFAULT = 12
NUM_INSTANCES_DEFAULT = 100
NUM_READS = 100
NUM_SWEEPS_DEFAULT = 1000
COEFF_TYPE_DEFAULT = 'lin2'

COEFF_LIN2 = [-1, 1]
COEFF_LIN20 = [round(x, 1) for x in np.arange(-1.0, 1.01, 0.1) if abs(x) > 1e-9]


def gen_small_graph(n_nodes, seed):
    """connected 3-regular random graph. (3-regular는 Pegasus 평균 degree와 다름 — 한계 참조)"""
    for attempt in range(50):
        G = nx.random_regular_graph(3, n_nodes, seed=seed + attempt)
        if nx.is_connected(G):
            return G
    raise RuntimeError(f"connected 3-regular graph 생성 실패 (n={n_nodes})")


def gen_random_qubo(G, coeff_type='lin2'):
    """CISC `gen_hardware_random_qubo`와 동일 패턴 (random.choice 사용)."""
    coeffs = COEFF_LIN2 if coeff_type == 'lin2' else COEFF_LIN20
    Q = {}
    for v in G.nodes():
        Q[(v, v)] = random.choice(coeffs)
    for (i, j) in G.edges():
        lo, hi = min(i, j), max(i, j)
        Q[(lo, hi)] = random.choice(coeffs)
    return Q


def brute_force_spectrum(Q, nodes):
    """전체 2^n 상태 → energy → (E_0, E_1, g_0, g_1, Δ_R)."""
    n_list = sorted(nodes)
    n = len(n_list)
    if n > 22:
        raise ValueError(f"brute force 불가: n={n}")
    idx_of = {v: i for i, v in enumerate(n_list)}

    linear = [0.0] * n
    quad = []
    for (i, j), w in Q.items():
        if i == j:
            linear[idx_of[i]] += w
        else:
            quad.append((idx_of[i], idx_of[j], w))

    cnt = Counter()
    for x_int in range(2 ** n):
        bits = [(x_int >> i) & 1 for i in range(n)]
        e = sum(linear[i] * bits[i] for i in range(n))
        e += sum(w * bits[a] * bits[b] for a, b, w in quad)
        cnt[round(e, 9)] += 1

    sorted_e = sorted(cnt.keys())
    if len(sorted_e) < 2:
        raise ValueError("spectrum에 단일 에너지 레벨만 존재")
    E_0, E_1 = sorted_e[0], sorted_e[1]
    g_0, g_1 = cnt[E_0], cnt[E_1]
    return E_0, E_1, g_0, g_1, E_1 - E_0


def measure_egsp(Q, gs_energy, num_reads, num_sweeps, seed):
    """SA → GS 에너지 도달 비율 (CISC 본실험과 동일 호출)."""
    sampler = neal.SimulatedAnnealingSampler()
    resp = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps, seed=seed)
    energies = [datum.energy for datum in resp.data(['energy'])]
    success = sum(1 for e in energies if abs(e - gs_energy) < 1e-9)
    return success / num_reads


def predict_p0(g_0, g_1, delta_R, beta):
    return 1.0 / (1.0 + (g_1 / g_0) * np.exp(-beta * delta_R))


def fit_beta(empirical, g_0s, g_1s, delta_Rs):
    """β 자유 파라미터 fit."""
    def model(X, beta):
        g0, g1, dR = X
        return 1.0 / (1.0 + (g1 / g0) * np.exp(-beta * dR))
    X = (g_0s.astype(float), g_1s.astype(float), delta_Rs.astype(float))
    popt, _ = curve_fit(model, X, empirical, p0=[1.0],
                        bounds=([0.001], [100.0]), maxfev=10000)
    return popt[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=N_NODES_DEFAULT)
    parser.add_argument('--instances', type=int, default=NUM_INSTANCES_DEFAULT)
    parser.add_argument('--sweeps', type=int, default=NUM_SWEEPS_DEFAULT)
    parser.add_argument('--coeff', default=COEFF_TYPE_DEFAULT, choices=['lin2', 'lin20'])
    args = parser.parse_args()

    print(f"=== α=0 plateau 평형 공식 검증 ===")
    print(f"  N={args.n}, instances={args.instances}, coeff={args.coeff}, "
          f"sweeps={args.sweeps}, reads={NUM_READS}")
    print()

    results = []
    t0 = time.perf_counter()

    for inst in range(args.instances):
        seed = inst * 53
        G = gen_small_graph(args.n, seed=seed)
        random.seed(seed)
        Q = gen_random_qubo(G, args.coeff)

        E_0, E_1, g_0, g_1, delta_R = brute_force_spectrum(Q, list(G.nodes()))
        e_gsp = measure_egsp(Q, E_0, NUM_READS, args.sweeps, seed=seed)

        results.append({
            'inst': inst,
            'g_0': int(g_0), 'g_1': int(g_1),
            'E_0': float(E_0), 'E_1': float(E_1),
            'delta_R': float(delta_R),
            'e_gsp': float(e_gsp),
        })

        if (inst + 1) % 20 == 0:
            elapsed = time.perf_counter() - t0
            print(f"  진행: {inst+1}/{args.instances} ({elapsed:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"\n  완료: {elapsed:.0f}s")

    # ─── 통계 집계 ───
    g_0s = np.array([r['g_0'] for r in results], dtype=float)
    g_1s = np.array([r['g_1'] for r in results], dtype=float)
    delta_Rs = np.array([r['delta_R'] for r in results])
    e_gsps = np.array([r['e_gsp'] for r in results])
    ratios = g_1s / g_0s

    print(f"\n=== 인스턴스 통계 ===")
    print(f"  E-GSP        평균={e_gsps.mean():.4f} ± {e_gsps.std():.4f}, "
          f"range=[{e_gsps.min():.3f}, {e_gsps.max():.3f}]")
    print(f"  g_0          median={int(np.median(g_0s))}, "
          f"range=[{int(g_0s.min())}, {int(g_0s.max())}]")
    print(f"  g_1          median={int(np.median(g_1s))}, "
          f"range=[{int(g_1s.min())}, {int(g_1s.max())}]")
    print(f"  Δ_R          median={np.median(delta_Rs):.2f}, "
          f"range=[{delta_Rs.min():.2f}, {delta_Rs.max():.2f}]")
    print(f"  g_1/g_0      median={np.median(ratios):.2f}, "
          f"range=[{ratios.min():.2f}, {ratios.max():.2f}]")

    # ─── β fit ───
    beta_fit = fit_beta(e_gsps, g_0s, g_1s, delta_Rs)
    predicted = np.array([predict_p0(g0, g1, dR, beta_fit)
                          for g0, g1, dR in zip(g_0s, g_1s, delta_Rs)])

    ss_res = np.sum((e_gsps - predicted) ** 2)
    ss_tot = np.sum((e_gsps - e_gsps.mean()) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    residual_std = (e_gsps - predicted).std()

    print(f"\n=== Mehta 공식 fit ===")
    print(f"  fit β            = {beta_fit:.4f}")
    print(f"  예측 평균        = {predicted.mean():.4f}")
    print(f"  실측 평균        = {e_gsps.mean():.4f}")
    print(f"  R²               = {r_squared:.4f}")
    print(f"  잔차 std         = {residual_std:.4f}")

    # ─── 판정 ───
    print(f"\n=== 판정 ===")
    mean_match = abs(predicted.mean() - e_gsps.mean()) < 0.05
    if r_squared > 0.9 and mean_match:
        verdict = "ADOPT"
        print(f"  ✓ ADOPT: R²={r_squared:.3f}>0.9 + 평균 일치")
        print(f"    → α=0 영역에서 SA가 equilibrium-like. β={beta_fit:.3f}이 schedule effective β.")
        print(f"    → α>0 영역으로 확장 검증 권장.")
    elif mean_match:
        verdict = "PARTIAL"
        print(f"  ⚠ PARTIAL: 평균 일치하나 R²={r_squared:.3f} 낮음")
        print(f"    → 인스턴스별 fit 부정확. 평균 통계로만 활용 가능.")
        print(f"    → lin20으로 Δ_R 다양화하여 fit 정보 늘려 재시도 권장.")
    else:
        verdict = "REJECT"
        print(f"  ✗ REJECT: 예측 {predicted.mean():.3f} vs 실측 {e_gsps.mean():.3f}")
        print(f"    → SA가 equilibrium에 도달 못 함, 또는 가정 오류.")
        print(f"    → sweep 수 증가시켜 재검증 → 그래도 안 맞으면 framework 폐기.")

    # ─── 저장 ───
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(out_dir, exist_ok=True)
    out_fname = f'equilibrium_validation_n{args.n}_{args.coeff}_inst{args.instances}_sw{args.sweeps}.json'
    out_path = os.path.join(out_dir, out_fname)
    with open(out_path, 'w') as f:
        json.dump({
            'params': {
                'n_nodes': args.n, 'instances': args.instances,
                'coeff': args.coeff, 'num_reads': NUM_READS, 'num_sweeps': args.sweeps,
                'graph_type': 'random_3_regular',
            },
            'summary': {
                'e_gsp_mean': float(e_gsps.mean()), 'e_gsp_std': float(e_gsps.std()),
                'beta_fit': float(beta_fit), 'r_squared': float(r_squared),
                'predicted_mean': float(predicted.mean()),
                'residual_std': float(residual_std),
                'verdict': verdict,
            },
            'instances': results,
        }, f, indent=2)
    print(f"\n  저장: {out_path}")


if __name__ == "__main__":
    main()
