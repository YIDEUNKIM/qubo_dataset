"""
α=0.1 universal trap 13개에 극단적 SA schedule 적용 — SA family escape 가능성 최종 점검.

배경: 6 standard schedules에서 모두 trap → universal trap.
       sweep=10^4 ~ 10^6 모두 같은 결과 → mixing time 부족 아님.

극단적 schedule 옵션:
  - 매우 wide β_range: β=(0.001, 100) — 거의 무한 온도에서 매우 차가운 끝까지
  - 매우 wide + 긴 sweep: 10^5 sweep
  - β=(0.0001, 200) — 더 극단

기준:
  - 풀림 → SA family 내 escape 가능 (단지 default schedule이 부적절)
  - 0% 그대로 → SA Markov chain의 absolute structural trap
"""

import os, time, pickle, json, argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import numpy as np


TRAP_INSTS = [0, 5, 10, 12, 13, 15, 16, 18, 20, 21, 23, 27, 29]

EXTREME_SCHEDULES = [
    ('geometric', (0.001, 100.0)),
    ('geometric', (0.0001, 200.0)),
    ('linear',    (0.001, 100.0)),
    ('linear',    (0.0001, 200.0)),
    # 기존 default도 reference로
    ('geometric', None),
]


def build_Q(inst, alpha):
    Q = dict(inst['R_sum'])
    for k, v in inst['P'].items():
        Q[k] = Q.get(k, 0) + alpha * v
    return {k: v for k, v in Q.items() if v != 0.0}


def sa_worker(args):
    """neal SA → exact bit-match egsp + per-read mean HD/ΔE (부동소수점 안전).
    HD/ΔE는 egsp=0 영역에서 valley 깊이 진단용 보조 측정.
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
    parser.add_argument('--alpha', type=float, default=0.1)
    parser.add_argument('--sweeps', type=int, default=100000)
    parser.add_argument('--num_reads', type=int, default=100)
    parser.add_argument('--workers', type=int, default=max(1, mp.cpu_count() - 1))
    parser.add_argument('--pickle',
                        default='hardened_posiform/instances/instances_pegasus16_lin2_100.pkl')
    args = parser.parse_args()

    print(f"=== Extreme schedules on universal trap 13 ===")
    print(f"  α={args.alpha}, sweeps={args.sweeps}, schedules={len(EXTREME_SCHEDULES)}, trap insts={len(TRAP_INSTS)}")

    with open(args.pickle, 'rb') as f:
        data = pickle.load(f)
    meta = data['meta']
    insts = [data['instances'][i] for i in TRAP_INSTS]
    print(f"  loaded {len(insts)} trap insts (indices: {TRAP_INSTS})\n")

    # Build Q once per inst
    inst_Q_list = []
    for orig_idx, inst in zip(TRAP_INSTS, insts):
        Q = build_Q(inst, args.alpha)
        gs_e = inst['t_energy_r'] + args.alpha * inst['t_energy_p']
        inst_Q_list.append({'orig_idx': orig_idx, 'seed': inst['seed'],
                            'Q': Q, 'gs_e': gs_e, 'target': inst['target']})

    tasks = []
    for d in inst_Q_list:
        for si, (stype, br) in enumerate(EXTREME_SCHEDULES):
            seed = (d['seed'] + 1) * 1009 + si * 7919 + args.sweeps
            tasks.append((d['orig_idx'], si,
                          (d['Q'], args.num_reads, args.sweeps, d['gs_e'],
                           seed, stype, br, d['target'])))
    print(f"  총 작업: {len(tasks)} (= 13 trap × {len(EXTREME_SCHEDULES)} schedules)")

    results = {}  # (orig_idx, si) -> (egsp, hd_mean, dE_mean)
    t0 = time.perf_counter()
    ctx = mp.get_context('fork') if 'fork' in mp.get_all_start_methods() else None
    pool_kwargs = {'max_workers': args.workers}
    if ctx is not None:
        pool_kwargs['mp_context'] = ctx

    with ProcessPoolExecutor(**pool_kwargs) as executor:
        fmap = {executor.submit(sa_worker, t[2]): (t[0], t[1]) for t in tasks}
        done = 0
        for fut in as_completed(fmap):
            orig_idx, si = fmap[fut]
            res = fut.result()
            if res[0] == 'OK':
                results[(orig_idx, si)] = (res[1], res[2], res[3])
            done += 1
            if done % 10 == 0 or done == len(tasks):
                el = time.perf_counter() - t0
                eta = el / done * (len(tasks) - done)
                print(f"    {done}/{len(tasks)} ({el:.0f}s, ETA {eta:.0f}s)")

    elapsed = time.perf_counter() - t0
    print(f"\n  완료: {elapsed:.1f}s")

    # 결과 표
    sched_labels = [f"{s}|{'def' if br is None else f'β={br[0]},{br[1]}'}" for s, br in EXTREME_SCHEDULES]
    print(f"\n  trap inst × schedule heatmap (E-GSP):")
    print(f"{'inst':>5} | " + ' '.join(f"{l[:22]:>22}" for l in sched_labels))
    print("-" * (5 + 3 + 23 * len(EXTREME_SCHEDULES)))
    any_solved = {i: False for i in TRAP_INSTS}
    for orig_idx in TRAP_INSTS:
        row = [results.get((orig_idx, si), (float('nan'),))[0] for si in range(len(EXTREME_SCHEDULES))]
        if any(v > 0.05 for v in row):
            any_solved[orig_idx] = True
        line = f"{orig_idx:>5} | " + ' '.join(f"{v:>22.4f}" for v in row)
        print(line)

    # HD / ΔE 표 (보조)
    print(f"\n  trap inst × schedule (HD mean):")
    print(f"{'inst':>5} | " + ' '.join(f"{l[:22]:>22}" for l in sched_labels))
    for orig_idx in TRAP_INSTS:
        row = [results.get((orig_idx, si), (float('nan'),)*3)[1] for si in range(len(EXTREME_SCHEDULES))]
        line = f"{orig_idx:>5} | " + ' '.join(f"{v:>22.2f}" for v in row)
        print(line)
    print(f"\n  trap inst × schedule (ΔE mean):")
    print(f"{'inst':>5} | " + ' '.join(f"{l[:22]:>22}" for l in sched_labels))
    for orig_idx in TRAP_INSTS:
        row = [results.get((orig_idx, si), (float('nan'),)*3)[2] for si in range(len(EXTREME_SCHEDULES))]
        line = f"{orig_idx:>5} | " + ' '.join(f"{v:>22.4g}" for v in row)
        print(line)

    n_solved = sum(any_solved.values())
    print(f"\n  Trap 13개 중 어떤 extreme schedule이라도 풀어낸 인스턴스: {n_solved}/13")
    if n_solved == 0:
        verdict = "ABSOLUTE_SA_TRAP — 13/13 모든 extreme schedule에서도 풀이 불가. SA family 내 escape 불가능"
    elif n_solved < 13:
        verdict = f"PARTIAL_ESCAPE — {n_solved}/13만 extreme schedule로 풀림. 나머지는 SA absolute trap"
    else:
        verdict = "ALL_ESCAPED — 13/13 모두 풀림 → 기존 default schedule의 한계, SA family 내 escape 가능"
    print(f"\n  ▶ {verdict}")

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"sa_extreme_trap_alpha{args.alpha}_sw{args.sweeps}.json")
    save = {
        'params': {'alpha': args.alpha, 'sweeps': args.sweeps, 'num_reads': args.num_reads,
                   'trap_insts': TRAP_INSTS,
                   'schedules': [(t, list(b) if b else None) for t, b in EXTREME_SCHEDULES]},
        'results': {f"{i}_{si}": (results[(i, si)][0] if (i, si) in results else None)
                    for i in TRAP_INSTS for si in range(len(EXTREME_SCHEDULES))},
        'results_hd': {f"{i}_{si}": (results[(i, si)][1] if (i, si) in results else None)
                       for i in TRAP_INSTS for si in range(len(EXTREME_SCHEDULES))},
        'results_dE': {f"{i}_{si}": (results[(i, si)][2] if (i, si) in results else None)
                       for i in TRAP_INSTS for si in range(len(EXTREME_SCHEDULES))},
        'any_solved_per_inst': {str(i): any_solved[i] for i in TRAP_INSTS},
        'n_solved': n_solved,
        'verdict': verdict,
        'elapsed_s': elapsed,
    }
    with open(out, 'w') as f:
        json.dump(save, f, indent=2)
    print(f"\n  저장: {out}")


if __name__ == "__main__":
    main()
