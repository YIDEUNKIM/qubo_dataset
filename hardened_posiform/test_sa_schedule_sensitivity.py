"""
SA schedule 민감도: ranking이 schedule-invariant인지 검증.

본 실험은 hardware_native_qubo.ipynb의 gen_hardware_native_qubo 파이프라인과 동일하게
구성된다 — target = R의 brute force GS, P는 alpha=1.0으로 한 번 생성,
Q(α) = R + α·P. test_sa_alpha_vulnerability.py와 동일한 인스턴스 셋을 재현.

배경
----
test_sa_alpha_vulnerability.py에서 α=0.1 transition zone에서 sweep=10 vs sweep=10000
ranking ρ=-0.36 (역전) 관측. 이게:
  (a) instance landscape의 본질적 다중모드 (다른 알고리즘도 같은 ranking) — 인스턴스 난이도
  (b) SA-specific 동역학 artifact (다른 SA schedule이 다른 ranking) — 알고리즘 취약점

검증 전략
---------
동일 인스턴스에 SA를 cooling schedule만 바꿔 여러 번 실행.
ranking이 schedule-invariant → 인스턴스 본질
ranking이 schedule-dependent → SA 알고리즘 artifact

schedule 변형:
  - geometric (neal 기본)
  - linear
  - β_range 좁게 / 넓게 변화
  - 동일 sweep 예산

만약 모든 schedule이 같은 ranking을 주면 → instance hardness가 ranking 결정 (H_inst)
schedule별로 ranking이 다르면 → SA의 schedule-specific bias (H_sched)

실행
----
python3 hardened_posiform/test_sa_schedule_sensitivity.py
python3 hardened_posiform/test_sa_schedule_sensitivity.py --alpha 0.1 --sweeps 100
"""

import os
import time
import json
import random
import argparse
import numpy as np
import networkx as nx
from scipy.stats import spearmanr
import neal
from pysat.solvers import Minisat22

N_DEFAULT = 12
INST_DEFAULT = 30
NUM_READS = 100
ALPHA_DEFAULT = 0.1   # transition zone (이전 실험에서 ρ 역전 관측)
SWEEPS_FOCUS = (10, 100, 1000)
COEFF_LIN2 = [-1, 1]


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


def posiform_planting_small(G, target_assignment, alpha):
    nodes = sorted(target_assignment.keys())
    n = len(nodes)
    node_to_idx = {v: i for i, v in enumerate(nodes)}
    edges = [(min(u, v), max(u, v)) for u, v in G.subgraph(nodes).edges()]
    if not edges:
        return {}
    Q_p = {}
    clauses = []
    all_tuples = [(0, 0), (0, 1), (1, 0), (1, 1)]
    max_steps = 200 * n

    def add_term(i, j):
        tt = (target_assignment[i], target_assignment[j])
        wt = [t for t in all_tuples if t != tt]
        wi, wj = random.choice(wt)
        lo, hi = min(i, j), max(i, j)
        if wi == 0 and wj == 0:
            Q_p[(i, i)] = Q_p.get((i, i), 0) - alpha
            Q_p[(j, j)] = Q_p.get((j, j), 0) - alpha
            Q_p[(lo, hi)] = Q_p.get((lo, hi), 0) + alpha
        elif wi == 0 and wj == 1:
            Q_p[(j, j)] = Q_p.get((j, j), 0) + alpha
            Q_p[(lo, hi)] = Q_p.get((lo, hi), 0) - alpha
        elif wi == 1 and wj == 0:
            Q_p[(i, i)] = Q_p.get((i, i), 0) + alpha
            Q_p[(lo, hi)] = Q_p.get((lo, hi), 0) - alpha
        else:
            Q_p[(lo, hi)] = Q_p.get((lo, hi), 0) + alpha
        lit_i = (node_to_idx[i] + 1) if wi == 0 else -(node_to_idx[i] + 1)
        lit_j = (node_to_idx[j] + 1) if wj == 0 else -(node_to_idx[j] + 1)
        clauses.append([lit_i, lit_j])

    def unique():
        with Minisat22() as s:
            for c in clauses:
                s.add_clause(c)
            blocking = []
            for v in nodes:
                idx = node_to_idx[v] + 1
                blocking.append(-idx if target_assignment[v] == 1 else idx)
            s.add_clause(blocking)
            return not s.solve()

    chk = max(1, n // 3)
    for step in range(max_steps):
        i, j = random.choice(edges)
        add_term(i, j)
        if (step + 1) % chk == 0 and unique():
            return {k: v for k, v in Q_p.items() if abs(v) > 1e-15}
    return {k: v for k, v in Q_p.items() if abs(v) > 1e-15}


def combine_qubo(R, P):
    Q = dict(R)
    for k, v in P.items():
        Q[k] = Q.get(k, 0) + v
    return Q


def full_energies(Q, nodes):
    """모든 2^n 상태의 에너지 (index = bit-encoded x)."""
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
    energies = np.zeros(1 << n)
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


def gs_energy_brute(Q, nodes):
    energies, _ = full_energies(Q, nodes)
    return float(energies.min())


def measure_egsp_schedule(Q, gs_energy, num_reads, num_sweeps, seed,
                          schedule_type='geometric', beta_range=None):
    """
    schedule_type: 'geometric' or 'linear'
    beta_range: (β_init, β_final) or None to use neal default
    """
    sampler = neal.SimulatedAnnealingSampler()
    kwargs = {
        'num_reads': num_reads,
        'num_sweeps': num_sweeps,
        'seed': seed,
        'beta_schedule_type': schedule_type,
    }
    if beta_range is not None:
        kwargs['beta_range'] = list(beta_range)
    resp = sampler.sample_qubo(Q, **kwargs)
    es = [d.energy for d in resp.data(['energy'])]
    return sum(1 for e in es if abs(e - gs_energy) < 1e-9) / num_reads


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=N_DEFAULT)
    parser.add_argument('--instances', type=int, default=INST_DEFAULT)
    parser.add_argument('--alpha', type=float, default=ALPHA_DEFAULT)
    parser.add_argument('--sweeps', default=','.join(str(s) for s in SWEEPS_FOCUS))
    args = parser.parse_args()

    sweeps_list = [int(s) for s in args.sweeps.split(',')]

    # Schedule 변형: geometric/linear × β_range 3종
    schedules = [
        ('geometric', None),                  # neal default
        ('geometric', (0.1, 5.0)),            # narrow range (cold start)
        ('geometric', (0.01, 20.0)),          # wide range (hot to cold)
        ('linear',    None),                  # linear cooling
        ('linear',    (0.1, 5.0)),
        ('linear',    (0.01, 20.0)),
    ]

    print(f"=== SA schedule 민감도 검증 ===")
    print(f"  N={args.n}, instances={args.instances}, α={args.alpha}, sweeps={sweeps_list}")
    print(f"  schedules: {len(schedules)}개")
    print()

    # 인스턴스 생성 — test_sa_alpha_vulnerability.py와 동일 파이프라인
    #   target = R의 brute force GS, P는 alpha=1.0으로 한 번 생성, Q = R + α·P
    instances = []
    t0 = time.perf_counter()
    for inst in range(args.instances):
        seed = inst * 53 + 7
        random.seed(seed)
        G = gen_small_graph(args.n, seed=seed)
        nodes = sorted(G.nodes())
        R = gen_random_qubo(G)
        energies_R, _ = full_energies(R, nodes)
        gs_idx = int(np.argmin(energies_R))
        target_asg = {v: ((gs_idx >> i) & 1) for i, v in enumerate(nodes)}
        random.seed(seed + 999983)
        P_unit = posiform_planting_small(G, target_asg, 1.0) if args.alpha > 0 else {}
        P_scaled = {k: args.alpha * v for k, v in P_unit.items()} if args.alpha > 0 else {}
        Q = combine_qubo(R, P_scaled)
        gs_e = gs_energy_brute(Q, nodes)
        instances.append({'inst': inst, 'seed': seed, 'Q': Q, 'gs_e': gs_e,
                          'target': target_asg})
    print(f"  인스턴스 생성: {time.perf_counter()-t0:.0f}s\n")

    # 측정: (instance, schedule, sweeps) → E-GSP
    # data[(sched_idx, sweeps)] = array of E-GSP per instance
    data = {}
    for si, (stype, br) in enumerate(schedules):
        sched_tag = f"{stype}|{'default' if br is None else f'β={br[0]},{br[1]}'}"
        print(f"  schedule {si}: {sched_tag}")
        ts = time.perf_counter()
        for s in sweeps_list:
            egsp_arr = np.zeros(len(instances))
            for ii, inst in enumerate(instances):
                egsp_arr[ii] = measure_egsp_schedule(
                    inst['Q'], inst['gs_e'], NUM_READS, s,
                    seed=inst['seed'], schedule_type=stype, beta_range=br
                )
            data[(si, s)] = egsp_arr
            print(f"    sweeps={s}: mean={egsp_arr.mean():.4f} std={egsp_arr.std():.4f}")
        print(f"    ({time.perf_counter()-ts:.0f}s)")

    # ─── 분석: schedule간 ranking 일치도 ───
    print(f"\n=== Schedule 간 ranking ρ (saturated 제외) ===")
    pairwise = {}
    for s in sweeps_list:
        print(f"\n  sweeps={s}:")
        rhos = []
        for i in range(len(schedules)):
            row = []
            for j in range(len(schedules)):
                ai, aj = data[(i, s)], data[(j, s)]
                m = (ai > 0.05) & (ai < 0.95) & (aj > 0.05) & (aj < 0.95)
                if m.sum() >= 5 and i != j:
                    rr, _ = spearmanr(ai[m], aj[m])
                    if np.isnan(rr):
                        rr = None
                else:
                    rr = 1.0 if i == j else None
                row.append(rr)
                if i < j and rr is not None:
                    rhos.append(rr)
            print(f"    sched[{i}]: " + ' '.join(
                f'{v:+.2f}' if v is not None else '  -  ' for v in row))
        if rhos:
            print(f"    median pairwise ρ: {float(np.median(rhos)):+.4f}")
            print(f"    min pairwise ρ:    {float(min(rhos)):+.4f}")
            pairwise[s] = {'median': float(np.median(rhos)), 'min': float(min(rhos)),
                          'count': len(rhos)}
        else:
            pairwise[s] = {'median': None, 'min': None, 'count': 0}

    # ─── 판정 ───
    print(f"\n=== 판정 ===")
    # 모든 sweep budget에서 schedule간 median ρ > 0.7 → schedule-invariant → instance 본질
    # 어느 budget에서든 median ρ < 0.5 → schedule-specific → SA artifact
    verdicts = {}
    for s in sweeps_list:
        m = pairwise[s].get('median')
        if m is None:
            verdicts[s] = 'INCONCLUSIVE'
            print(f"  sweeps={s:>5}: INCONCLUSIVE (saturated)")
        elif m > 0.7:
            verdicts[s] = 'INSTANCE'
            print(f"  sweeps={s:>5}: ρ_med={m:+.3f} > 0.7 → schedule-invariant (instance 본질)")
        elif m < 0.5:
            verdicts[s] = 'SCHEDULE'
            print(f"  sweeps={s:>5}: ρ_med={m:+.3f} < 0.5 → schedule-dependent (SA artifact)")
        else:
            verdicts[s] = 'PARTIAL'
            print(f"  sweeps={s:>5}: ρ_med={m:+.3f} (0.5~0.7) → partial schedule effect")

    if all(v in ('INSTANCE', 'INCONCLUSIVE') for v in verdicts.values()):
        final = "INSTANCE_DETERMINED"
        print(f"\n  ▶ 최종: schedule 무관 → ranking은 instance hardness 반영 (SA artifact 약함)")
    elif any(v == 'SCHEDULE' for v in verdicts.values()):
        final = "SCHEDULE_ARTIFACT"
        sw_bad = [s for s, v in verdicts.items() if v == 'SCHEDULE']
        print(f"\n  ▶ 최종: sweep {sw_bad}에서 ranking이 schedule 의존 → SA 알고리즘 취약점")
    else:
        final = "PARTIAL"
        print(f"\n  ▶ 최종: partial — 일부 sweep budget에서 schedule 효과 존재")

    # 저장
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f'sa_schedule_n{args.n}_alpha{args.alpha}_inst{args.instances}.json')
    save = {
        'params': {'n': args.n, 'instances': args.instances, 'alpha': args.alpha,
                   'sweeps': sweeps_list, 'num_reads': NUM_READS,
                   'schedules': [(t, list(b) if b else None) for t, b in schedules]},
        'pairwise_summary': {str(s): pairwise[s] for s in sweeps_list},
        'verdict_per_sweeps': {str(s): verdicts[s] for s in sweeps_list},
        'verdict_final': final,
        'data': {f'sched{i}_sw{s}': data[(i, s)].tolist()
                 for i in range(len(schedules)) for s in sweeps_list},
    }
    with open(out, 'w') as f:
        json.dump(save, f, indent=2)
    print(f"\n  저장: {out}")


if __name__ == "__main__":
    main()
