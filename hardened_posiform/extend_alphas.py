"""α ∈ {0.02, 0.04, 0.07} 추가 측정 (lin2, Pegasus P16, 200 inst × 100 reads × 1000 sweeps).

기존 노트북 cell 9 GSP-vs-α 파이프라인을 standalone으로 재현.
- inst pkl: instances_pegasus16_lin2_200.pkl
- per-sample (hd, energy, delta_energy) 수집
- per-instance summary
- per-α aggregate (avg_hd, avg_dE, gsp_*)
- 결과 dir: results/gsp_extend_alphas_lin2_inst200_<timestamp>/
"""
import os
import sys
import csv
import json
import time
import pickle
import gc
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


# Forked workers inherit this read-only (copy-on-write) so we never pickle the
# ~6 MB-per-instance QUBOs across the pool nor materialise all 500 Q dicts at
# once in the parent (that 3 GB spike was OOM-killing workers on this 15 GB box).
_INSTS = None


def sa_worker(args):
    """neal SA → per-sample records (read_idx, hd, energy, delta, is_target, is_R_gs, eng_ok).

    Q is built inside the worker from the forked-global instance list.
    """
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
        # lin2엔 0이 없어서 모든 노드가 살아있지만, 방어적으로 missing 처리
        target_bit_by_node = {nd: int(target_str[i]) for i, nd in enumerate(sorted_nodes)}
        records = []
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
                is_target = (found_str == target_str)
                delta = float(energy) - float(target_energy)
                eng_ok = abs(delta) < 1e-7
                # α > 0이면 target이 unique GS → is_R_gs = is_target
                is_R_gs = is_target
                records.append((read_idx, hd, float(energy), delta,
                                bool(is_target), bool(is_R_gs), bool(eng_ok)))
                read_idx += 1
            if read_idx >= num_reads:
                break
        return ('OK', inst_id, records)
    except Exception as e:
        try:
            return ('ERR', args[0], str(e))
        except Exception:
            return ('ERR', -1, str(e))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pkl', default='hardened_posiform/instances/instances_pegasus16_lin2_200.pkl')
    parser.add_argument('--alphas', default='0.02,0.04,0.07')
    parser.add_argument('--num_reads', type=int, default=100)
    parser.add_argument('--sweeps', type=int, default=1000)
    parser.add_argument('--workers', type=int, default=max(1, mp.cpu_count() - 1))
    parser.add_argument('--out_tag', default='extend')
    args = parser.parse_args()

    alphas = [float(x) for x in args.alphas.split(',')]
    print(f"=== α 확장 실험 (HD+ΔE) ===")
    print(f"  pkl: {args.pkl}")
    print(f"  αs: {alphas}, num_reads={args.num_reads}, sweeps={args.sweeps}, workers={args.workers}")

    t_load = time.perf_counter()
    with open(args.pkl, 'rb') as f:
        data = pickle.load(f)
    meta = data['meta']
    insts = data['instances']
    print(f"  loaded: {meta['topology']}{meta['topo_size']}, n={meta['n']}, coeff={meta['coeff']}, "
          f"{len(insts)} insts ({time.perf_counter()-t_load:.1f}s)")

    # 각 inst의 target_str / sorted_nodes 사전 추출 후, 워커가 안 쓰는 큰 필드 제거
    # (target/all_block_gs/partitions … → fork COW 메모리 절감)
    _KEEP = {'R_sum', 'P', 't_energy_r', 't_energy_p', 'target_str', 'sorted_nodes'}
    for inst in insts:
        if 'sorted_nodes' not in inst or 'target_str' not in inst:
            sorted_nodes = sorted(inst['target'].keys())
            target_str = ''.join(str(inst['target'][nd]) for nd in sorted_nodes)
            inst['sorted_nodes'] = sorted_nodes
            inst['target_str'] = target_str
        for k in [k for k in inst if k not in _KEEP]:
            del inst[k]
    gc.collect()

    # 출력 dir
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    out_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'results',
        f'gsp_{args.out_tag}_alphas_{meta["coeff"]}_inst{len(insts)}_{stamp}'
    )
    os.makedirs(out_dir, exist_ok=True)
    print(f"  결과 dir: {out_dir}\n")

    # 누적 row containers
    all_samples_rows = []  # alpha, instance_id, read_idx, hd, energy, delta, is_target, is_R_gs, eng_ok
    all_inst_rows = []     # alpha, instance_id, hd_mean, hd_min, hd_max, hd_std, eng_c, bit_c, gs_c, best_e, best_d, target_density
    gsp_results = {}

    ctx = mp.get_context('fork') if 'fork' in mp.get_all_start_methods() else None
    pool_kwargs = {'max_workers': args.workers}
    if ctx is not None:
        pool_kwargs['mp_context'] = ctx

    global _INSTS
    _INSTS = insts  # forked workers build their own Q from this (no 3 GB pre-build)

    t0_all = time.perf_counter()
    for ai, alpha in enumerate(alphas):
        # 병렬 실행 — Q는 워커 내부에서 build (inst_id만 전달)
        records_by_inst = {}
        errors = []
        t_alpha = time.perf_counter()
        with ProcessPoolExecutor(**pool_kwargs) as executor:
            futures = [
                executor.submit(sa_worker, (ii, alpha, args.num_reads, args.sweeps))
                for ii in range(len(insts))
            ]
            for fut in as_completed(futures):
                res = fut.result()
                if res[0] == 'OK':
                    _, inst_id, records = res
                    records_by_inst[inst_id] = records
                else:
                    _, inst_id, msg = res
                    errors.append((inst_id, msg))
                    if len(errors) <= 3:
                        print(f"  [worker error] α={alpha} inst={inst_id}: {msg}")

        # α별 집계
        eng_total = bit_total = gs_total = hd_total = reads_total = 0
        dE_total = 0.0
        for inst_id in sorted(records_by_inst.keys()):
            records = records_by_inst[inst_id]
            inst = insts[inst_id]
            target_density = sum(int(c) for c in inst['target_str']) / len(inst['target_str'])
            hd_arr = []
            best_energy = float('inf')
            best_delta = float('inf')
            eng_c = bit_c = gs_c = 0
            for (read_idx, hd, energy, delta, is_target, is_R_gs, eng_ok) in records:
                all_samples_rows.append((alpha, inst_id, read_idx, hd, energy, delta,
                                          int(is_target), int(is_R_gs), int(eng_ok)))
                hd_arr.append(hd)
                if energy < best_energy:
                    best_energy = energy
                    best_delta = delta
                if eng_ok:
                    eng_c += 1
                if is_target:
                    bit_c += 1
                if is_R_gs:
                    gs_c += 1
                hd_total += hd
                dE_total += float(delta)
                reads_total += 1
                eng_total += int(eng_ok)
                bit_total += int(is_target)
                gs_total += int(is_R_gs)
            hd_np = np.array(hd_arr, dtype=np.int32) if hd_arr else np.array([0], dtype=np.int32)
            all_inst_rows.append((alpha, inst_id,
                                   float(hd_np.mean()), int(hd_np.min()),
                                   int(hd_np.max()), float(hd_np.std()),
                                   eng_c, bit_c, gs_c,
                                   float(best_energy), float(best_delta),
                                   target_density))

        gsp_eng = eng_total / reads_total if reads_total else 0
        gsp_bit = bit_total / reads_total if reads_total else 0
        gsp_gs  = gs_total  / reads_total if reads_total else 0
        avg_hd  = hd_total  / reads_total if reads_total else 0
        avg_dE  = dE_total  / reads_total if reads_total else 0.0
        gsp_results[alpha] = {
            'gsp_eng': gsp_eng, 'gsp_bit': gsp_bit, 'gsp_gs': gsp_gs,
            'eng_ok': eng_total, 'bit_ok': bit_total, 'gs_ok': gs_total,
            'hd_sum': hd_total, 'avg_hd': avg_hd,
            'dE_sum': dE_total, 'avg_dE': avg_dE,
            'total': reads_total, 'eng_found': eng_total,
        }
        el = time.perf_counter() - t_alpha
        el_all = time.perf_counter() - t0_all
        eta = el_all / (ai + 1) * (len(alphas) - ai - 1)
        print(f"  α={alpha:<7}: GS={gsp_gs:.3f} Bit={gsp_bit:.3f} HD={avg_hd:.2f} ΔE={avg_dE:.4g} "
              f"[{ai+1}/{len(alphas)}, α={el:.0f}s, total={el_all:.0f}s, ETA {eta:.0f}s]")
        del records_by_inst
        gc.collect()

    elapsed = time.perf_counter() - t0_all
    print(f"\n  완료: {elapsed:.1f}s, 총 records: {len(all_samples_rows):,}")

    # CSV 저장
    samples_path = os.path.join(out_dir, 'samples.csv')
    with open(samples_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['alpha', 'instance_id', 'read_idx', 'hd', 'energy', 'delta_energy',
                    'is_target', 'is_R_gs', 'eng_ok'])
        w.writerows(all_samples_rows)
    print(f"  samples.csv: {len(all_samples_rows):,} rows → {samples_path}")

    instances_path = os.path.join(out_dir, 'instances.csv')
    with open(instances_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['alpha', 'instance_id', 'hd_mean', 'hd_min', 'hd_max', 'hd_std',
                    'eng_count', 'bit_count', 'gs_count', 'best_energy', 'best_delta',
                    'target_density'])
        w.writerows(all_inst_rows)
    print(f"  instances.csv: {len(all_inst_rows):,} rows → {instances_path}")

    data_json = {
        'params': {
            'topology': meta['topology'], 'topo_size': meta['topo_size'],
            'coeff': meta['coeff'], 'num_instances': len(insts),
            'num_reads': args.num_reads, 'fixed_sweep': args.sweeps,
            'alphas': alphas, 'elapsed_s': elapsed,
            'n': meta['n'],
        },
        'results': {f'{a:g}': gsp_results[a] for a in alphas},
    }
    with open(os.path.join(out_dir, 'data.json'), 'w') as f:
        json.dump(data_json, f, indent=2)
    print(f"  data.json → {os.path.join(out_dir, 'data.json')}")


if __name__ == '__main__':
    main()
