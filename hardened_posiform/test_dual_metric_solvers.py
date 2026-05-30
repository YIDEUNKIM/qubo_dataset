"""T3 — dual-metric solver-discrimination experiment for ieee_paper_2page.

Demonstrates that in the valley/pre-cliff regime, the *binary success rate* is
0% for BOTH a cheap and an expensive SA budget (so it cannot rank them), yet the
continuous metrics ΔE (and HD) clearly separate the two budgets. This is the
empirical justification for dual-metric reporting (defends novelty + the absence
of QPU runs: the metric is shown useful for classical solver comparison).

Two SA budgets on identical instances:
  - SA-1k : neal, 1000 sweeps  (cheap)
  - SA-5k : neal, 5000 sweeps  (5x compute)
α grid spans deep-valley → pre-cliff → cliff-onset → escape.

Pipeline (per-sample HD/energy/ΔE/is_target) is identical to extend_alphas.py /
the notebook cell-9 GSP pipeline, so results are directly comparable.

Output: results/t3_dual_metric_<stamp>/{summary.json, deltas_valley.npz, fig_dual_metric.{png,pdf}}
"""
import os
import csv
import json
import time
import pickle
import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

import numpy as np


def build_Q(inst, alpha):
    Q = dict(inst['R_sum'])
    for k, v in inst['P'].items():
        Q[k] = Q.get(k, 0) + alpha * v
    return {k: v for k, v in Q.items() if abs(v) > 1e-15}


# Forked workers inherit this read-only; Q is built per-worker (no 3 GB pre-build).
_INSTS = None


def sa_worker(args):
    """neal SA → per-sample (hd, delta, is_target). Q built inside the worker."""
    try:
        import neal
        (inst_id, alpha, num_reads, num_sweeps) = args
        inst = _INSTS[inst_id]
        Q_dict = build_Q(inst, alpha)
        target_energy = inst['t_energy_r'] + alpha * inst['t_energy_p']
        target_str = inst['target_str']
        sorted_nodes = inst['sorted_nodes']
        sampler = neal.SimulatedAnnealingSampler()
        resp = sampler.sample_qubo(Q_dict, num_reads=num_reads, num_sweeps=num_sweeps)
        target_bit_by_node = {nd: int(target_str[i]) for i, nd in enumerate(sorted_nodes)}
        hds, deltas, hits = [], [], []
        read_idx = 0
        for sample, energy, num_occ in resp.data(['sample', 'energy', 'num_occurrences']):
            sample = dict(sample)
            for nd in sorted_nodes:
                if nd not in sample:
                    sample[nd] = target_bit_by_node[nd]
            for _ in range(int(num_occ)):
                if read_idx >= num_reads:
                    break
                found_str = ''.join(str(sample[nd]) for nd in sorted_nodes)
                hd = sum(1 for a, b in zip(target_str, found_str) if a != b)
                delta = float(energy) - float(target_energy)
                hds.append(hd)
                deltas.append(delta)
                hits.append(int(found_str == target_str))
                read_idx += 1
            if read_idx >= num_reads:
                break
        return ('OK', inst_id, np.array(hds), np.array(deltas), np.array(hits))
    except Exception as e:
        return ('ERR', args[0], str(e), None, None)


def run_cell(insts, alpha, sweeps, num_reads, pool_kwargs):
    """Run one (alpha, sweeps) cell over all insts → flat per-sample arrays."""
    all_hd, all_dE, all_hit = [], [], []
    errors = 0
    with ProcessPoolExecutor(**pool_kwargs) as ex:
        futures = [ex.submit(sa_worker, (i, alpha, num_reads, sweeps))
                   for i in range(len(insts))]
        for fut in as_completed(futures):
            res = fut.result()
            if res[0] == 'OK':
                _, _, hd, dE, hit = res
                all_hd.append(hd); all_dE.append(dE); all_hit.append(hit)
            else:
                errors += 1
    hd = np.concatenate(all_hd); dE = np.concatenate(all_dE); hit = np.concatenate(all_hit)
    return hd, dE, hit, errors


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--pkl', default='instances/instances_pegasus16_lin2_500.pkl')
    p.add_argument('--alphas', default='1e-5,1e-3,1e-2,2e-2')
    p.add_argument('--sweeps', default='1000,5000')
    p.add_argument('--num_reads', type=int, default=100)
    p.add_argument('--max_inst', type=int, default=200)
    p.add_argument('--workers', type=int, default=12)  # capped: 15 GB box, swap full
    args = p.parse_args()

    alphas = [float(x) for x in args.alphas.split(',')]
    budgets = [int(x) for x in args.sweeps.split(',')]
    here = os.path.dirname(os.path.abspath(__file__))
    pkl_path = args.pkl if os.path.isabs(args.pkl) else os.path.join(here, args.pkl)

    print(f"=== T3 dual-metric solver discrimination ===")
    print(f"  pkl={pkl_path}")
    print(f"  αs={alphas}  budgets(sweeps)={budgets}  reads={args.num_reads}  "
          f"max_inst={args.max_inst}  workers={args.workers}")

    import gc
    with open(pkl_path, 'rb') as f:
        data = pickle.load(f)
    insts = list(data['instances'][:args.max_inst])
    meta = data['meta']
    del data  # free the unused instances (we only keep the first max_inst)
    # build target_str/sorted_nodes, then strip worker-unused fields (fork COW)
    _KEEP = {'R_sum', 'P', 't_energy_r', 't_energy_p', 'target_str', 'sorted_nodes'}
    for inst in insts:
        if 'sorted_nodes' not in inst or 'target_str' not in inst:
            sn = sorted(inst['target'].keys())
            inst['sorted_nodes'] = sn
            inst['target_str'] = ''.join(str(inst['target'][nd]) for nd in sn)
        for k in [k for k in inst if k not in _KEEP]:
            del inst[k]
    gc.collect()
    print(f"  loaded {meta['topology']}{meta['topo_size']} n={meta['n']} "
          f"coeff={meta['coeff']} using {len(insts)} insts")

    global _INSTS
    _INSTS = insts  # forked workers build their own Q from this

    ctx = mp.get_context('fork') if 'fork' in mp.get_all_start_methods() else None
    pool_kwargs = {'max_workers': args.workers}
    if ctx is not None:
        pool_kwargs['mp_context'] = ctx

    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    out_dir = os.path.join(here, 'results', f't3_dual_metric_{stamp}')
    os.makedirs(out_dir, exist_ok=True)

    rows = []          # summary per (alpha, sweeps)
    valley_deltas = {} # (alpha, sweeps) -> raw dE array (for distribution comparison)
    t0 = time.perf_counter()
    for alpha in alphas:
        for sw in budgets:
            tc = time.perf_counter()
            hd, dE, hit, errs = run_cell(insts, alpha, sw, args.num_reads, pool_kwargs)
            n = len(dE)
            row = {
                'alpha': alpha, 'sweeps': sw, 'n_samples': int(n), 'errors': errs,
                'success_rate': float(hit.mean()),
                'hd_mean': float(hd.mean()), 'hd_sem': float(hd.std(ddof=1) / np.sqrt(n)),
                'dE_mean': float(dE.mean()), 'dE_sem': float(dE.std(ddof=1) / np.sqrt(n)),
                'dE_median': float(np.median(dE)),
                'dE_q25': float(np.percentile(dE, 25)), 'dE_q75': float(np.percentile(dE, 75)),
            }
            rows.append(row)
            valley_deltas[f'{alpha:g}_{sw}'] = dE.astype(np.float32)
            el = time.perf_counter() - tc
            print(f"  α={alpha:<6g} sw={sw:<5d}: succ={row['success_rate']*100:6.2f}%  "
                  f"meanΔE={row['dE_mean']:.4f}±{row['dE_sem']:.4f}  medΔE={row['dE_median']:.4f}  "
                  f"meanHD={row['hd_mean']:.2f}  [{el:.0f}s]")

    elapsed = time.perf_counter() - t0
    summary = {'params': {'alphas': alphas, 'budgets': budgets, 'num_reads': args.num_reads,
                          'n_inst': len(insts), 'n': meta['n'], 'coeff': meta['coeff'],
                          'elapsed_s': elapsed}, 'rows': rows}
    with open(os.path.join(out_dir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    np.savez_compressed(os.path.join(out_dir, 'deltas.npz'), **valley_deltas)
    print(f"\n  done in {elapsed:.0f}s → {out_dir}")
    make_figure(summary, out_dir)


def make_figure(summary, out_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    rows = summary['rows']
    alphas = summary['params']['alphas']
    budgets = summary['params']['budgets']
    xpos = np.arange(len(alphas))
    palette = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red']
    colors = {sw: palette[i % len(palette)] for i, sw in enumerate(budgets)}

    def series(sw, key):
        d = {r['alpha']: r[key] for r in rows if r['sweeps'] == sw}
        return np.array([d[a] for a in alphas])

    # compact single-column (IEEE) stacked figure: success (top) + ΔE (bottom)
    fig, (axT, axB) = plt.subplots(2, 1, figsize=(3.5, 3.7), sharex=True)

    # (a) success rate — both pinned at 0 across the valley
    for sw in budgets:
        axT.plot(xpos, series(sw, 'success_rate') * 100, 'o-', color=colors[sw],
                 linewidth=1.6, markersize=5, label=f'{sw} sw')
    axT.set_yscale('symlog', linthresh=1)
    axT.set_ylim(bottom=0)
    axT.set_ylabel('success [%]', fontsize=9)
    axT.set_title('(a) target success rate', fontsize=9)
    axT.grid(True, alpha=0.3, which='both')
    axT.legend(fontsize=8, loc='upper left', title='SA budget', title_fontsize=8)
    axT.tick_params(labelsize=8)

    # (b) mean ΔE with SEM — clearly separates the two budgets in the valley
    for sw in budgets:
        m = series(sw, 'dE_mean'); s = series(sw, 'dE_sem')
        axB.errorbar(xpos, m, yerr=s, fmt='o-', color=colors[sw], linewidth=1.6,
                     markersize=5, capsize=2)
    axB.set_xticks(xpos); axB.set_xticklabels([f'{a:g}' for a in alphas],
                                              rotation=45, fontsize=8)
    axB.set_xlabel(r'$\alpha$', fontsize=9)
    axB.set_ylabel(r'mean $\Delta E\ (\pm$SEM)', fontsize=9)
    axB.set_title(r'(b) energy gap $\Delta E$', fontsize=9)
    axB.grid(True, alpha=0.3)
    axB.tick_params(labelsize=8)

    plt.tight_layout(h_pad=0.6)
    for ext in ('png', 'pdf'):
        out = os.path.join(out_dir, f'fig_dual_metric.{ext}')
        fig.savefig(out, dpi=200, bbox_inches='tight')
        print(f"  saved {out}")
    plt.close(fig)


if __name__ == '__main__':
    main()
