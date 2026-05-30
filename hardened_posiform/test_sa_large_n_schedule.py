"""
Large-n (Pegasus P16, n=5627) hardened posiform에서 SA schedule 민감도 검증.

목표: 같은 인스턴스에 6개 cooling schedule 적용 → ranking이 schedule 의존인지 검증.

기준:
  S1: ρ_med(schedule간 ranking) > 0.7 → schedule-invariant (instance 본질)
  S2: ρ_med < 0.5 → schedule-dependent (SA artifact)
  사이값 → partial

실행:
  python3 hardened_posiform/test_sa_large_n_schedule.py --inst 20 --alpha 0.01 --sweeps 10000
"""

import os
import time
import pickle
import json
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

import numpy as np
from scipy.stats import spearmanr


SCHEDULES = [
    ('geometric', None),
    ('geometric', (0.1, 5.0)),
    ('geometric', (0.01, 20.0)),
    ('linear',    None),
    ('linear',    (0.1, 5.0)),
    ('linear',    (0.01, 20.0)),
]


def build_Q(inst, alpha):
    Q = dict(inst['R_sum'])
    for k, v in inst['P'].items():
        Q[k] = Q.get(k, 0) + alpha * v
    return {k: v for k, v in Q.items() if v != 0.0}


def sa_worker(args):
    """neal SA → exact bit-match egsp + per-read mean HD/ΔE.
    energy 1e-9 비교는 부동소수점 false negative 양산 → bit-string 매치가 robust.
    HD/ΔE는 보조 측정 (egsp=0 영역의 valley 구조용).
    """
    import neal
    try:
        (Q, num_reads, num_sweeps, gs_energy, seed,
         schedule_type, beta_range, target_bits) = args
        sampler = neal.SimulatedAnnealingSampler()
        kwargs = dict(num_reads=num_reads, num_sweeps=num_sweeps,
                      seed=seed, beta_schedule_type=schedule_type)
        if beta_range is not None:
            kwargs['beta_range'] = list(beta_range)
        resp = sampler.sample_qubo(Q, **kwargs)
        sorted_nodes = sorted(target_bits.keys())
        hit = 0
        hd_sum = 0
        dE_sum = 0.0
        for sample, energy in zip(resp.record.sample, resp.record.energy):
            sample_dict = dict(zip(resp.variables, sample))
            hd = sum(1 for v in sorted_nodes if sample_dict[v] != target_bits[v])
            hd_sum += hd
            dE_sum += float(energy) - float(gs_energy)
            if hd == 0:
                hit += 1
        egsp = hit / num_reads
        hd_mean = hd_sum / num_reads
        dE_mean = dE_sum / num_reads
        return ('OK', egsp, hd_mean, dE_mean)
    except Exception as e:
        return ('ERR', str(e))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inst', type=int, default=20)
    parser.add_argument('--alpha', type=float, default=0.01)
    parser.add_argument('--sweeps', type=int, default=10000)
    parser.add_argument('--num_reads', type=int, default=100)
    parser.add_argument('--workers', type=int, default=max(1, mp.cpu_count() - 1))
    parser.add_argument('--pickle',
                        default='hardened_posiform/instances/instances_pegasus16_lin2_100.pkl')
    parser.add_argument('--out_tag', default='schedule')
    args = parser.parse_args()

    print(f"=== Large-n SA schedule 민감도 ===")
    print(f"  α={args.alpha}, sweeps={args.sweeps}, num_reads={args.num_reads}")
    print(f"  instances={args.inst}, schedules={len(SCHEDULES)}개, workers={args.workers}")
    print()

    with open(args.pickle, 'rb') as f:
        data = pickle.load(f)
    meta = data['meta']
    insts = data['instances'][:args.inst]
    print(f"  loaded: {meta['topology']} {meta['topo_size']}, n={meta['n']}, {len(insts)} insts\n")

    inst_Q_list = []
    for i, inst in enumerate(insts):
        Q = build_Q(inst, args.alpha)
        gs_e = inst['t_energy_r'] + args.alpha * inst['t_energy_p']
        inst_Q_list.append({'idx': i, 'seed': inst['seed'], 'Q': Q, 'gs_e': gs_e,
                             'target': inst['target']})

    # tasks: (inst_idx, sched_idx, args_tuple)
    tasks = []
    for d in inst_Q_list:
        for si, (stype, br) in enumerate(SCHEDULES):
            seed = (d['seed'] + 1) * 1009 + si * 7919
            tasks.append((d['idx'], si,
                          (d['Q'], args.num_reads, args.sweeps, d['gs_e'], seed, stype, br, d['target'])))
    print(f"  총 작업: {len(tasks)}")

    results = {}  # (inst_idx, si) -> (egsp, hd_mean, dE_mean)
    errors = []
    t0 = time.perf_counter()

    ctx = mp.get_context('fork') if 'fork' in mp.get_all_start_methods() else None
    pool_kwargs = {'max_workers': args.workers}
    if ctx is not None:
        pool_kwargs['mp_context'] = ctx

    with ProcessPoolExecutor(**pool_kwargs) as executor:
        fmap = {executor.submit(sa_worker, t[2]): (t[0], t[1]) for t in tasks}
        done = 0
        for fut in as_completed(fmap):
            inst_idx, si = fmap[fut]
            res = fut.result()
            if res[0] == 'OK':
                results[(inst_idx, si)] = (res[1], res[2], res[3])
            else:
                errors.append((inst_idx, si, res[1]))
            done += 1
            if done % 20 == 0 or done == len(tasks):
                el = time.perf_counter() - t0
                eta = el / done * (len(tasks) - done)
                print(f"    {done}/{len(tasks)} ({el:.0f}s, ETA {eta:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"\n  완료: {elapsed:.1f}s, errors={len(errors)}")

    # 결과 array (egsp / HD / ΔE)
    data_mat = np.full((len(SCHEDULES), len(insts)), np.nan)
    hd_mat = np.full((len(SCHEDULES), len(insts)), np.nan)
    dE_mat = np.full((len(SCHEDULES), len(insts)), np.nan)
    for inst_idx in range(len(insts)):
        for si in range(len(SCHEDULES)):
            r = results.get((inst_idx, si))
            if r is not None:
                data_mat[si, inst_idx] = r[0]
                hd_mat[si, inst_idx] = r[1]
                dE_mat[si, inst_idx] = r[2]

    print(f"\n  schedule별 평균 E-GSP / HD / ΔE:")
    for si, (stype, br) in enumerate(SCHEDULES):
        m = float(np.nanmean(data_mat[si]))
        s = float(np.nanstd(data_mat[si]))
        hd_m = float(np.nanmean(hd_mat[si]))
        dE_m = float(np.nanmean(dE_mat[si]))
        tag = f"{stype}|{'default' if br is None else f'β=({br[0]},{br[1]})'}"
        print(f"    [{si}] {tag:35s}: egsp={m:.4f}±{s:.4f}, HD={hd_m:.2f}, ΔE={dE_m:.4g}")

    # pairwise ρ
    print(f"\n  schedule간 pairwise ρ:")
    rho_matrix = np.zeros((len(SCHEDULES), len(SCHEDULES)))
    rho_nonsat = np.full((len(SCHEDULES), len(SCHEDULES)), np.nan)
    all_rhos = []
    nonsat_rhos = []
    for i in range(len(SCHEDULES)):
        for j in range(len(SCHEDULES)):
            a, b = data_mat[i], data_mat[j]
            r, _ = spearmanr(a, b)
            rho_matrix[i, j] = r if not np.isnan(r) else 0.0
            mask = (a > 0.05) & (a < 0.95) & (b > 0.05) & (b < 0.95)
            if mask.sum() >= 5:
                r_ns, _ = spearmanr(a[mask], b[mask])
                if not np.isnan(r_ns):
                    rho_nonsat[i, j] = r_ns
                    if i < j:
                        nonsat_rhos.append(r_ns)
            if i < j:
                all_rhos.append(rho_matrix[i, j])

    for i in range(len(SCHEDULES)):
        row = ' '.join(f"{rho_matrix[i,j]:+.2f}" for j in range(len(SCHEDULES)))
        print(f"    sched[{i}]: {row}")

    if all_rhos:
        rho_med_all = float(np.median(all_rhos))
        rho_min_all = float(np.min(all_rhos))
        print(f"\n  ρ_all (모든 인스턴스):    median={rho_med_all:+.4f}, min={rho_min_all:+.4f}")
    if nonsat_rhos:
        rho_med_ns = float(np.median(nonsat_rhos))
        rho_min_ns = float(np.min(nonsat_rhos))
        print(f"  ρ_nonsat (비포화):        median={rho_med_ns:+.4f}, min={rho_min_ns:+.4f}")

    # 판정
    print(f"\n=== 판정 ===")
    rho_med_for_verdict = rho_med_ns if nonsat_rhos else (rho_med_all if all_rhos else None)
    if rho_med_for_verdict is None:
        verdict = "INCONCLUSIVE — 모든 인스턴스 상수"
    elif rho_med_for_verdict > 0.7:
        verdict = "INSTANCE_DETERMINED — schedule-invariant ranking (인스턴스 본질이 ranking 결정)"
    elif rho_med_for_verdict < 0.5:
        verdict = "SCHEDULE_ARTIFACT — ranking이 cooling schedule에 의존"
    else:
        verdict = "PARTIAL — 부분적 schedule 효과 (0.5 < ρ_med ≤ 0.7)"
    print(f"  ▶ {verdict}")

    # save
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"sa_large_n_{args.out_tag}_alpha{args.alpha}_sw{args.sweeps}_inst{args.inst}.json")
    save = {
        'params': {
            'topology': meta['topology'], 'topo_size': meta['topo_size'],
            'n': meta['n'], 'coeff': meta['coeff'],
            'inst': args.inst, 'alpha': args.alpha, 'sweeps': args.sweeps,
            'num_reads': args.num_reads,
            'schedules': [(t, list(b) if b else None) for t, b in SCHEDULES],
        },
        'mean_egsp_per_sched': [float(np.nanmean(data_mat[i])) for i in range(len(SCHEDULES))],
        'mean_hd_per_sched':   [float(np.nanmean(hd_mat[i]))   for i in range(len(SCHEDULES))],
        'mean_dE_per_sched':   [float(np.nanmean(dE_mat[i]))   for i in range(len(SCHEDULES))],
        'rho_matrix': rho_matrix.tolist(),
        'rho_nonsat_matrix': [[None if np.isnan(v) else float(v)
                               for v in row] for row in rho_nonsat],
        'rho_med_all': float(np.median(all_rhos)) if all_rhos else None,
        'rho_med_nonsat': float(np.median(nonsat_rhos)) if nonsat_rhos else None,
        'verdict': verdict,
        'data_per_sched':    {str(i): data_mat[i].tolist() for i in range(len(SCHEDULES))},
        'data_hd_per_sched': {str(i): hd_mat[i].tolist()   for i in range(len(SCHEDULES))},
        'data_dE_per_sched': {str(i): dE_mat[i].tolist()   for i in range(len(SCHEDULES))},
        'elapsed_s': elapsed,
    }
    with open(out, 'w') as f:
        json.dump(save, f, indent=2)
    print(f"\n  저장: {out}")


if __name__ == "__main__":
    main()
