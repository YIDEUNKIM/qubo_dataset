"""
α=0.01 mixing failure가 SA-specific인지 algorithm-independent intrinsic hardness인지 검증.

배경: 정정 후 V3 — α=0.01에서 SA는 sweep=10⁵로도 mean E-GSP 0.30, plateau 미도달.
질문: 다른 알고리즘 (Tabu, Greedy)이 같은 30 인스턴스를 더 잘 푸는가?
  - 더 잘 풀면 → V3은 SA-specific algorithmic flaw (paper-worthy)
  - 동일 또는 더 못 풀면 → 인스턴스 본질적 어려움 (instance hardness lower bound)

알고리즘:
  - neal SA (geometric default, sweep=10^4, 10^5)
  - dwave-tabu (timeout 다양화)
  - dwave-greedy (deterministic local descent, 빠름)
  - SA + greedy 후속 (post-process)

측정: bit-string exact match metric.

실행:
  python3 hardened_posiform/test_sa_vs_alternative.py --inst 30 --alpha 0.01
"""

import os, time, pickle, json, argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import numpy as np


def build_Q(inst, alpha):
    Q = dict(inst['R_sum'])
    for k, v in inst['P'].items():
        Q[k] = Q.get(k, 0) + alpha * v
    return {k: v for k, v in Q.items() if v != 0.0}


def measure_target_match(resp, target_bits):
    """resp의 모든 read를 target bit-string과 exact match로 비교 → hit ratio."""
    sorted_nodes = sorted(target_bits.keys())
    hit = 0
    total = 0
    for sample in resp.record.sample:
        sample_dict = dict(zip(resp.variables, sample))
        if all(sample_dict[v] == target_bits[v] for v in sorted_nodes):
            hit += 1
        total += 1
    return hit / total if total > 0 else 0.0


def sa_worker(args):
    import neal
    Q, num_reads, num_sweeps, target_bits = args
    sampler = neal.SimulatedAnnealingSampler()
    resp = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps,
                               beta_schedule_type='geometric')
    return measure_target_match(resp, target_bits)


def tabu_worker(args):
    import tabu
    Q, num_reads, timeout_ms, target_bits = args
    sampler = tabu.TabuSampler()
    resp = sampler.sample_qubo(Q, num_reads=num_reads, timeout=timeout_ms)
    return measure_target_match(resp, target_bits)


def greedy_worker(args):
    import greedy
    Q, num_reads, target_bits = args
    sampler = greedy.SteepestDescentSampler()
    resp = sampler.sample_qubo(Q, num_reads=num_reads)
    return measure_target_match(resp, target_bits)


def run_experiment(args):
    print(f"=== α={args.alpha} 알고리즘 비교 (n={args.inst}, num_reads={args.num_reads}) ===\n")

    with open(args.pickle, 'rb') as f:
        data = pickle.load(f)
    meta = data['meta']
    insts = data['instances'][:args.inst]
    print(f"loaded: {meta['topology']} {meta['topo_size']}, n={meta['n']}, {len(insts)} insts\n")

    # Pre-build Q matrices
    inst_data = []
    for i, inst in enumerate(insts):
        Q = build_Q(inst, args.alpha)
        inst_data.append({'idx': i, 'Q': Q, 'target': inst['target']})

    workers = max(1, mp.cpu_count() - 1)
    ctx = mp.get_context('fork') if 'fork' in mp.get_all_start_methods() else None
    pool_kwargs = {'max_workers': workers}
    if ctx is not None:
        pool_kwargs['mp_context'] = ctx

    # Tasks
    configs = []
    # SA at multiple sweeps
    for sw in [10000, 100000]:
        configs.append(('SA_sw'+str(sw), 'sa', {'num_sweeps': sw}))
    # Tabu with various timeouts (in seconds)
    for to_ms in [1000, 5000, 20000]:
        configs.append(('Tabu_'+str(to_ms//1000)+'s', 'tabu', {'timeout_ms': to_ms}))
    # Greedy
    configs.append(('Greedy', 'greedy', {}))

    results = {}  # (config_name, inst_idx) -> egsp
    elapsed_per_config = {}

    for cfg_name, kind, kwargs in configs:
        print(f"\n--- {cfg_name} ---")
        t0 = time.perf_counter()
        tasks = []
        for d in inst_data:
            if kind == 'sa':
                tasks.append((d['idx'], (d['Q'], args.num_reads, kwargs['num_sweeps'], d['target'])))
            elif kind == 'tabu':
                tasks.append((d['idx'], (d['Q'], args.num_reads, kwargs['timeout_ms'], d['target'])))
            else:  # greedy
                tasks.append((d['idx'], (d['Q'], args.num_reads, d['target'])))

        with ProcessPoolExecutor(**pool_kwargs) as executor:
            if kind == 'sa':
                futures = {executor.submit(sa_worker, t[1]): t[0] for t in tasks}
            elif kind == 'tabu':
                futures = {executor.submit(tabu_worker, t[1]): t[0] for t in tasks}
            else:
                futures = {executor.submit(greedy_worker, t[1]): t[0] for t in tasks}
            for fut in as_completed(futures):
                idx = futures[fut]
                results[(cfg_name, idx)] = fut.result()

        elapsed = time.perf_counter() - t0
        elapsed_per_config[cfg_name] = elapsed
        egsp_arr = np.array([results[(cfg_name, i)] for i in range(args.inst)])
        print(f"  mean E-GSP = {egsp_arr.mean():.4f}, std = {egsp_arr.std():.4f}, "
              f"range = [{egsp_arr.min():.3f}, {egsp_arr.max():.3f}]  ({elapsed:.0f}s)")

    # Summary table
    print(f"\n=== 종합 비교 (α={args.alpha}, {args.inst} 인스턴스) ===")
    print(f"{'config':>15} | {'mean E-GSP':>10} | {'std':>8} | {'풀린 인스턴스 (≥0.5)':>20} | {'elapsed':>8}")
    print("-" * 80)
    by_cfg = {}
    for cfg_name, _, _ in configs:
        egsp_arr = np.array([results[(cfg_name, i)] for i in range(args.inst)])
        by_cfg[cfg_name] = egsp_arr
        n_solved = int((egsp_arr >= 0.5).sum())
        print(f"{cfg_name:>15} | {egsp_arr.mean():>10.4f} | {egsp_arr.std():>8.4f} | "
              f"{n_solved:>11}/{args.inst} ({n_solved*100/args.inst:.0f}%) | {elapsed_per_config[cfg_name]:>7.0f}s")

    # Per-instance comparison: which insts SA can't solve but alternatives can?
    print(f"\n=== per-instance comparison (mean E-GSP across reads) ===")
    print(f"{'inst':>4} | " + ' '.join(f"{n[:13]:>13}" for n, _, _ in configs))
    print("-" * (5 + 3 + 14 * len(configs)))
    sa_better_insts = []
    alt_better_insts = []
    for i in range(args.inst):
        row = [by_cfg[n][i] for n, _, _ in configs]
        line = f"{i:>4} | " + ' '.join(f"{v:>13.3f}" for v in row)
        print(line)
        sa_best = max(by_cfg[n][i] for n, _, _ in configs if 'SA' in n)
        alt_best = max(by_cfg[n][i] for n, _, _ in configs if 'SA' not in n)
        if alt_best > sa_best + 0.05:
            alt_better_insts.append((i, sa_best, alt_best))
        elif sa_best > alt_best + 0.05:
            sa_better_insts.append((i, sa_best, alt_best))

    print(f"\nAlternative > SA (>0.05 차이): {len(alt_better_insts)}/30")
    for i, sb, ab in alt_better_insts:
        print(f"  inst {i}: SA best={sb:.3f}, alt best={ab:.3f}")
    print(f"\nSA > Alternative (>0.05 차이): {len(sa_better_insts)}/30")
    for i, sb, ab in sa_better_insts:
        print(f"  inst {i}: SA best={sb:.3f}, alt best={ab:.3f}")

    # Verdict
    print(f"\n=== 판정 ===")
    sa_mean = max(by_cfg[n].mean() for n, _, _ in configs if 'SA' in n)
    alt_mean = max(by_cfg[n].mean() for n, _, _ in configs if 'SA' not in n)
    if alt_mean > sa_mean + 0.10:
        verdict = f"V3_SA_SPECIFIC — 대체 알고리즘이 SA보다 mean {alt_mean-sa_mean:+.3f} 더 잘 풀음. V3는 SA-specific."
    elif sa_mean > alt_mean + 0.10:
        verdict = f"V3_INSTANCE_INTRINSIC — SA가 대체 알고리즘보다 더 잘 풀음. V3는 instance가 hardware-specific 또는 SA-favorable."
    else:
        verdict = f"V3_ALGORITHMS_TIED — SA와 대체 알고리즘 거의 동일 수준 (|Δmean|={abs(sa_mean-alt_mean):.3f}). 인스턴스 본질적 어려움 강한 증거."
    print(f"  ▶ {verdict}")

    # Save
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"sa_vs_alt_alpha{args.alpha}_inst{args.inst}.json")
    save = {
        'params': {'alpha': args.alpha, 'inst': args.inst, 'num_reads': args.num_reads,
                   'configs': [(n, k, kw) for n, k, kw in configs]},
        'mean_egsp_per_config': {n: float(by_cfg[n].mean()) for n, _, _ in configs},
        'std_per_config': {n: float(by_cfg[n].std()) for n, _, _ in configs},
        'n_solved_per_config': {n: int((by_cfg[n] >= 0.5).sum()) for n, _, _ in configs},
        'elapsed_per_config': elapsed_per_config,
        'per_inst_per_config': {n: by_cfg[n].tolist() for n, _, _ in configs},
        'alt_better_insts': [{'inst': i, 'sa_best': sb, 'alt_best': ab}
                              for i, sb, ab in alt_better_insts],
        'verdict': verdict,
    }
    with open(out, 'w') as f:
        json.dump(save, f, indent=2)
    print(f"\n  저장: {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--inst', type=int, default=30)
    parser.add_argument('--alpha', type=float, default=0.01)
    parser.add_argument('--num_reads', type=int, default=100)
    parser.add_argument('--pickle',
                        default='hardened_posiform/instances/instances_pegasus16_lin2_100.pkl')
    args = parser.parse_args()
    run_experiment(args)
