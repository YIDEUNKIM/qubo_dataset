"""
SA fail mode 실증 실험 (Goal-driven).

가설:
  H0: SA의 인스턴스별 E-GSP가 Mehta thermal equilibrium 공식
      p_0 = 1 / (1 + (g_1/g_0) · exp(-β · Δ_R)) 로 충분히 예측됨.
  H1: SA에 (g_0, g_1, Δ_R)로 환원 불가 dynamic fail mode가 존재함.

Verifiable goals (G1~G5로 H1을 falsify 또는 confirm):
  G1: |mean(residual)| < 0.03            — 앙상블 평균 일치 (β fit 정상)
  G2: std(residual) > 0.05                — 인스턴스 변동 존재 (정보 부족)
  G3: max |Spearman ρ(r, X)| > 0.3        — 단일 dynamic 변수와 강한 신호
  G4: Mehta_pred 같은 bin 내 std(SA) > 0.1 — 같은 prediction, 다른 결과 (직접 증명)
  G5: Multiple regression R²(r ~ dyn) > 0.15 — 다변수 dynamic 신호 (G3 보완)

판정:
  G1 + G2 + (G3 or G5) + G4 통과 → H1 채택 (SA fail mode 증명)
  G2 통과 실패                    → H0 채택 (SA ≈ equilibrium)
  중간                            → 보류

Dynamic 변수 (brute force / Q dict 직계산):
  - n_local_minima   : gradient-descent local optima 수
  - gs_basin_fraction: random start → greedy → GS 도달 비율
  - max_abs_coupling : max|J_ij|
  - neg_coupling_ratio: 음수 quadratic edge 비율

실행:
  python3 hardened_posiform/test_sa_failmode.py
  python3 hardened_posiform/test_sa_failmode.py --instances 500 --n 12
"""

import os
import time
import json
import random
import argparse
from collections import Counter

import numpy as np
import networkx as nx
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
import neal

N_DEFAULT = 12
INST_DEFAULT = 1000
NUM_READS = 100
SWEEPS_DEFAULT = 1000
COEFF_LIN2 = [-1, 1]
BASIN_SAMPLES = 200


def gen_small_graph(n, seed):
    for k in range(50):
        G = nx.random_regular_graph(3, n, seed=seed + k)
        if nx.is_connected(G):
            return G
    raise RuntimeError(f"connected 3-regular 실패 n={n}")


def gen_random_qubo(G):
    Q = {}
    for v in G.nodes():
        Q[(v, v)] = random.choice(COEFF_LIN2)
    for (i, j) in G.edges():
        lo, hi = min(i, j), max(i, j)
        Q[(lo, hi)] = random.choice(COEFF_LIN2)
    return Q


def build_energy_table(Q, nodes):
    """2^n 상태 모두에 대한 에너지 list (index = bit-encoded x)."""
    n_list = sorted(nodes)
    n = len(n_list)
    idx_of = {v: i for i, v in enumerate(n_list)}
    linear = [0.0] * n
    quad = []
    for (i, j), w in Q.items():
        if i == j:
            linear[idx_of[i]] += w
        else:
            quad.append((idx_of[i], idx_of[j], w))

    energies = [0.0] * (1 << n)
    for x in range(1 << n):
        e = 0.0
        for i in range(n):
            if (x >> i) & 1:
                e += linear[i]
        for a, b, w in quad:
            if ((x >> a) & 1) and ((x >> b) & 1):
                e += w
        energies[x] = e
    return energies, n


def spectrum_stats(energies):
    cnt = Counter(round(e, 9) for e in energies)
    sorted_e = sorted(cnt.keys())
    if len(sorted_e) < 2:
        raise ValueError("단일 에너지 레벨")
    E_0, E_1 = sorted_e[0], sorted_e[1]
    return E_0, E_1, cnt[E_0], cnt[E_1], E_1 - E_0


def count_local_minima(energies, n):
    cnt = 0
    for x in range(1 << n):
        e = energies[x]
        is_lm = True
        for k in range(n):
            if energies[x ^ (1 << k)] < e:
                is_lm = False
                break
        if is_lm:
            cnt += 1
    return cnt


def gs_basin_fraction(energies, n, gs_energy, num_samples, rng):
    success = 0
    for _ in range(num_samples):
        x = rng.randrange(1 << n)
        while True:
            best_x, best_e = x, energies[x]
            for k in range(n):
                ne = energies[x ^ (1 << k)]
                if ne < best_e:
                    best_x, best_e = x ^ (1 << k), ne
            if best_x == x:
                break
            x = best_x
        if abs(energies[x] - gs_energy) < 1e-9:
            success += 1
    return success / num_samples


def measure_sa_gsp(Q, gs_energy, num_reads, num_sweeps, seed):
    sampler = neal.SimulatedAnnealingSampler()
    resp = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps, seed=seed)
    es = [d.energy for d in resp.data(['energy'])]
    return sum(1 for e in es if abs(e - gs_energy) < 1e-9) / num_reads


def predict_p0(g0, g1, dR, beta):
    return 1.0 / (1.0 + (g1 / g0) * np.exp(-beta * dR))


def fit_beta(p_sa, g0s, g1s, dRs):
    def model(X, beta):
        g0, g1, dR = X
        return 1.0 / (1.0 + (g1 / g0) * np.exp(-beta * dR))
    X = (g0s.astype(float), g1s.astype(float), dRs.astype(float))
    popt, _ = curve_fit(model, X, p_sa, p0=[1.0], bounds=([0.001], [100.0]), maxfev=10000)
    return popt[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=N_DEFAULT)
    parser.add_argument('--instances', type=int, default=INST_DEFAULT)
    parser.add_argument('--sweeps', type=int, default=SWEEPS_DEFAULT)
    args = parser.parse_args()

    print(f"=== SA fail mode 실증 (Goal-driven) ===")
    print(f"  N={args.n}, instances={args.instances}, sweeps={args.sweeps}, reads={NUM_READS}")
    print(f"  H0: SA = Mehta 평형 예측")
    print(f"  H1: SA에 dynamic fail mode 존재")
    print()

    results = []
    t0 = time.perf_counter()

    for inst in range(args.instances):
        seed = inst * 53
        G = gen_small_graph(args.n, seed=seed)
        random.seed(seed)
        Q = gen_random_qubo(G)

        energies, n = build_energy_table(Q, list(G.nodes()))
        E_0, E_1, g_0, g_1, dR = spectrum_stats(energies)
        n_lm = count_local_minima(energies, n)
        rng = random.Random(seed)
        basin = gs_basin_fraction(energies, n, E_0, BASIN_SAMPLES, rng)

        quad_items = [(k, v) for k, v in Q.items() if k[0] != k[1]]
        max_abs_J = max(abs(v) for _, v in quad_items)
        neg_ratio = sum(1 for _, v in quad_items if v < 0) / max(1, len(quad_items))

        p_sa = measure_sa_gsp(Q, E_0, NUM_READS, args.sweeps, seed=seed)

        results.append({
            'inst': inst,
            'g_0': int(g_0), 'g_1': int(g_1), 'delta_R': float(dR),
            'e_gsp': float(p_sa),
            'n_local_minima': int(n_lm),
            'gs_basin_fraction': float(basin),
            'max_abs_coupling': float(max_abs_J),
            'neg_coupling_ratio': float(neg_ratio),
        })

        if (inst + 1) % 100 == 0:
            elapsed = time.perf_counter() - t0
            print(f"  진행: {inst+1}/{args.instances} ({elapsed:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"\n  측정 완료: {elapsed:.0f}s")

    # ─── 분석 ───
    g0s = np.array([r['g_0'] for r in results], dtype=float)
    g1s = np.array([r['g_1'] for r in results], dtype=float)
    dRs = np.array([r['delta_R'] for r in results])
    p_sa = np.array([r['e_gsp'] for r in results])

    beta = fit_beta(p_sa, g0s, g1s, dRs)
    p_pred = np.array([predict_p0(g0, g1, dR, beta) for g0, g1, dR in zip(g0s, g1s, dRs)])
    r = p_sa - p_pred

    print(f"\n=== Mehta 앙상블 fit ===")
    print(f"  β = {beta:.4f}")
    print(f"  평균: 예측 {p_pred.mean():.4f}  실측 {p_sa.mean():.4f}")
    print(f"  잔차: mean = {r.mean():+.4f}, std = {r.std():.4f}, "
          f"range = [{r.min():+.3f}, {r.max():+.3f}]")

    G1 = abs(r.mean()) < 0.03
    G2 = r.std() > 0.05
    print(f"\n  G1 |mean(r)| < 0.03: {'✓' if G1 else '✗'}  ({abs(r.mean()):.4f})")
    print(f"  G2 std(r) > 0.05:    {'✓' if G2 else '✗'}  ({r.std():.4f})")

    # G3: Spearman corr with dynamic vars
    dyn = {
        'delta_R': dRs,
        'g_1/g_0': g1s / g0s,
        'n_local_minima': np.array([rr['n_local_minima'] for rr in results]),
        'gs_basin_fraction': np.array([rr['gs_basin_fraction'] for rr in results]),
        'max_abs_coupling': np.array([rr['max_abs_coupling'] for rr in results]),
        'neg_coupling_ratio': np.array([rr['neg_coupling_ratio'] for rr in results]),
    }
    print(f"\n=== Residual vs dynamic 변수 (Spearman ρ) ===")
    corrs = {}
    max_abs_rho = 0.0
    max_var = None
    for name, X in dyn.items():
        if np.std(X) < 1e-12:
            rho, pval = float('nan'), 1.0
        else:
            rho, pval = spearmanr(r, X)
        corrs[name] = float(rho) if not np.isnan(rho) else 0.0
        print(f"  {name:22s} ρ = {rho:+.4f}  (p = {pval:.2e})")
        if not np.isnan(rho) and abs(rho) > max_abs_rho:
            max_abs_rho = abs(rho)
            max_var = name
    G3 = max_abs_rho > 0.3
    print(f"\n  G3 max|ρ| > 0.3:     {'✓' if G3 else '✗'}  (max={max_abs_rho:.4f} on '{max_var}')")

    # G4: bin analysis
    n_bins = 10
    bin_edges = np.linspace(p_pred.min(), p_pred.max() + 1e-9, n_bins + 1)
    G4 = False
    max_bin_std = 0.0
    max_bin_idx = -1
    for b in range(n_bins):
        mask = (p_pred >= bin_edges[b]) & (p_pred < bin_edges[b + 1])
        if mask.sum() < 5:
            continue
        bs = p_sa[mask].std()
        if bs > max_bin_std:
            max_bin_std = bs
            max_bin_idx = b
        if bs > 0.1:
            G4 = True
    print(f"\n=== Same-prediction bin 분석 ===")
    if max_bin_idx >= 0:
        b = max_bin_idx
        mask = (p_pred >= bin_edges[b]) & (p_pred < bin_edges[b + 1])
        print(f"  최대 within-bin std: {max_bin_std:.4f}")
        print(f"    bin [{bin_edges[b]:.3f}, {bin_edges[b+1]:.3f}), n={mask.sum()}")
        print(f"    SA_GSP range: [{p_sa[mask].min():.3f}, {p_sa[mask].max():.3f}]")
    print(f"\n  G4 max within-bin std > 0.1: {'✓' if G4 else '✗'}  ({max_bin_std:.4f})")

    # G5: Multiple regression (residual ~ dynamic 변수들 standardized)
    valid_dyn = {name: X for name, X in dyn.items() if np.std(X) > 1e-12}
    X_mat = np.column_stack([(X - X.mean()) / X.std() for X in valid_dyn.values()])
    X_aug = np.column_stack([np.ones(len(r)), X_mat])
    coefs, *_ = np.linalg.lstsq(X_aug, r, rcond=None)
    r_hat = X_aug @ coefs
    ss_res_mr = float(np.sum((r - r_hat) ** 2))
    ss_tot_mr = float(np.sum((r - r.mean()) ** 2))
    R2_mr = 1.0 - ss_res_mr / ss_tot_mr if ss_tot_mr > 0 else float('nan')

    print(f"\n=== Multiple regression (잔차 ~ 모든 유효 dynamic 변수) ===")
    print(f"  R² = {R2_mr:.4f}")
    print(f"  Standardized coefficients:")
    for name, c in zip(valid_dyn.keys(), coefs[1:]):
        print(f"    {name:22s} = {c:+.4f}")
    G5 = R2_mr > 0.15
    print(f"\n  G5 multivariate R² > 0.15: {'✓' if G5 else '✗'} ({R2_mr:.4f})")

    # 종합 판정
    print(f"\n=== 종합 판정 ===")
    print(f"  G1 평균 일치          : {'✓' if G1 else '✗'}")
    print(f"  G2 인스턴스 변동      : {'✓' if G2 else '✗'}")
    print(f"  G3 단변수 signal      : {'✓' if G3 else '✗'}")
    print(f"  G4 bin within-std     : {'✓' if G4 else '✗'}")
    print(f"  G5 multivariate signal: {'✓' if G5 else '✗'}")

    h1_signal = G3 or G5
    if G1 and G2 and h1_signal and G4:
        verdict = "H1_ADOPTED"
        print(f"\n  ▶ H1 채택: SA에 dynamic fail mode 존재 — 증명 완료")
        print(f"    Mehta framework은 평균만 잡고 인스턴스별 hardness는 미흡.")
        if G3:
            print(f"    단변수 신호: '{max_var}' (ρ={max_abs_rho:+.3f})")
        if G5:
            top_dyn = sorted(zip(valid_dyn.keys(), coefs[1:]), key=lambda kv: -abs(kv[1]))[:3]
            print(f"    다변수 핵심 predictor (|coef| 상위 3):")
            for name, c in top_dyn:
                print(f"      {name}: {c:+.4f}")
        print(f"    G4: 같은 Mehta 예측에 다른 SA 결과 직접 증명 (within-bin std={max_bin_std:.3f})")
    elif not G2:
        verdict = "H0_ADOPTED"
        print(f"\n  ▶ H0 채택: SA는 평형에 매우 가까움 (fail mode 없음)")
    else:
        verdict = "INCONCLUSIVE"
        print(f"\n  ▶ 보류: 추가 검증 필요")

    # 저장
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f'sa_failmode_n{args.n}_inst{args.instances}_sw{args.sweeps}.json')
    with open(out_path, 'w') as f:
        json.dump({
            'params': {'n': args.n, 'instances': args.instances,
                       'num_reads': NUM_READS, 'num_sweeps': args.sweeps},
            'beta_fit': float(beta),
            'goals': {
                'G1_mean_match': bool(G1), 'mean_residual': float(r.mean()),
                'G2_variance': bool(G2), 'std_residual': float(r.std()),
                'G3_corr': bool(G3), 'max_corr_var': max_var, 'max_abs_rho': float(max_abs_rho),
                'G4_bin': bool(G4), 'max_within_bin_std': float(max_bin_std),
                'G5_multivariate': bool(G5), 'R2_multivariate': float(R2_mr),
                'multivariate_coefs': {n: float(c) for n, c in zip(valid_dyn.keys(), coefs[1:])},
                'verdict': verdict,
            },
            'correlations': corrs,
            'instances': results,
        }, f, indent=2)
    print(f"\n  저장: {out_path}")


if __name__ == "__main__":
    main()
