"""QPU 파일럿 — 본 실험 전 소규모 타당성/탐색 (Advantage_system6.4).

확인 목표:
  (a) active-set 정합: 인스턴스 Q_α의 노드/커플러가 현재 칩에 다 있는가
  (b) 분해 threshold α*: QPU가 planted αP를 α 얼마부터 느끼는가 (wrong-only ΔE가
      α=0 baseline에서 벗어나는 첫 α)
  (c) QPU 실패율: SA로 고른 집단이 QPU에서도 갇히는가 (안 갇히면 생존편향 경고)

⚠ 실제 QPU 제출 = Leap 시간 소모. 토큰/Leap 설정 필요 (`dwave config` 또는 DWAVE_API_TOKEN).
  파일럿 비용 ≈ N_INST × len(ALPHAS) × 37ms access (기본 8×6 ≈ 1.8s). 매우 작음.

산출: results/qpu_pilot_<stamp>/{summary.json, samples.csv} + 콘솔 리포트.
"""
import os, csv, json, pickle, datetime
from collections import defaultdict
import numpy as np

PKL      = 'instances/instances_pegasus16_lin2_hard300.pkl'
N_INST   = 30
ALPHAS   = [0.0, 1e-9, 1e-5, 1e-4, 1e-3, 3e-3, 5e-3, 1e-2, 1.5e-2, 2e-2]   # SA와 동일 그리드
NUM_READS = 100
ANNEAL_US = 20
GAUGES   = 2                                      # spin-reversal 게이지 평균


def build_Q(inst, alpha):
    Q = dict(inst['R_sum'])
    for k, v in inst['P'].items():
        Q[k] = Q.get(k, 0) + alpha * v
    return {k: v for k, v in Q.items() if abs(v) > 1e-15}


def qubo_energy(sample, Q):
    e = 0.0
    for (i, j), w in Q.items():
        e += w * sample[i] if i == j else w * sample[i] * sample[j]
    return e


def main():
    data = pickle.load(open(PKL, 'rb'))
    insts = data['instances'][:N_INST]
    meta = data['meta']
    solver_name = meta.get('qpu_solver')
    print(f'pkl: {PKL}  | insts: {len(insts)}  | target solver: {solver_name}')

    from dwave.system import DWaveSampler
    qpu = DWaveSampler(solver={'name': solver_name}) if solver_name \
        else DWaveSampler(solver={'topology__type': 'pegasus'})
    print(f'connected solver: {qpu.solver.name}')
    if solver_name and qpu.solver.name != solver_name:
        raise RuntimeError(f'solver mismatch: got {qpu.solver.name}, want {solver_name}')

    qpu_nodes = set(qpu.nodelist)
    qpu_edges = {frozenset(e) for e in qpu.edgelist}

    # ---- (a) active-set 검증 (전 인스턴스, 제출 전) ----
    valid = []
    for idx, inst in enumerate(insts):
        keys = set(inst['R_sum']) | set(inst['P'])           # 모든 α의 노드/커플러 합집합
        nodes = {i for k in keys for i in k}
        edges = {frozenset(k) for k in keys if k[0] != k[1]}
        miss_n = nodes - qpu_nodes
        miss_e = edges - qpu_edges
        if miss_n or miss_e:
            print(f'  inst {idx}: INVALID — missing {len(miss_n)} nodes, {len(miss_e)} edges')
        else:
            valid.append((idx, inst))
    print(f'(a) active-set OK: {len(valid)}/{len(insts)} insts valid')
    if not valid:
        raise RuntimeError('valid instance 0개 — calibration drift 또는 solver 불일치. 중단.')

    # 고정 스케일 (auto_scale=False 대비) — 본 실험 검증용으로 계산만, 파일럿은 auto_scale 사용
    cmax = max(abs(v) for _, inst in valid for a in ALPHAS for v in build_Q(inst, a).values())
    print(f'    (참고) 전 α 최대 |계수| = {cmax:.3f} → 고정 스케일 후보 SCALE={1/cmax:.3f}')

    # ---- 제출 ----
    gauge_ok = [True]
    def sample(Q):
        kw = dict(num_reads=NUM_READS, annealing_time=ANNEAL_US, auto_scale=True,
                  answer_mode='raw')
        if gauge_ok[0]:
            try:
                return qpu.sample_qubo(Q, num_spin_reversal_transforms=GAUGES, **kw)
            except (TypeError, ValueError, KeyError):        # 이 솔버/버전이 미지원
                gauge_ok[0] = False
                print('  [warn] num_spin_reversal_transforms 미지원 — 게이지 평균 없이 진행 '
                      '(본 실험은 SpinReversalTransformComposite로 추가 필요)')
        return qpu.sample_qubo(Q, **kw)

    rows = []          # samples.csv
    acc = defaultdict(lambda: {'dE': [], 'hd': [], 'n': 0, 'wrong': 0})
    for idx, inst in valid:
        tnode = {nd: int(inst['target_str'][k]) for k, nd in enumerate(inst['sorted_nodes'])}
        for a in ALPHAS:
            Q = build_Q(inst, a)
            te = inst['t_energy_r'] + a * inst['t_energy_p']
            resp = sample(Q)
            ri = 0
            for s, _e, occ in resp.data(['sample', 'energy', 'num_occurrences']):
                s = dict(s)
                for _ in range(int(occ)):
                    if ri >= NUM_READS:
                        break
                    e = qubo_energy(s, Q)
                    dE = e - te
                    hd = sum(s[nd] != tnode[nd] for nd in inst['sorted_nodes'])
                    is_t = int(hd == 0)
                    is_R_gs = int(abs(qubo_energy(s, inst['R_sum']) - inst['t_energy_r']) < 1e-7)
                    eng_ok = int(abs(dE) < 1e-7)
                    rows.append([a, idx, ri, hd, e, dE, is_t, is_R_gs, eng_ok])
                    acc[a]['n'] += 1
                    if not is_t:
                        acc[a]['wrong'] += 1; acc[a]['dE'].append(dE); acc[a]['hd'].append(hd)
                    ri += 1
        print(f'  inst {idx} done')

    # ---- 리포트 ----
    print(f'\n{"α":>8} {"실패율":>7} {"meanΔE":>8} {"medΔE":>7} {"meanHD":>8} {"medHD":>7}')
    summ = {}
    for a in ALPHAS:
        d = acc[a]; fr = d['wrong'] / d['n'] if d['n'] else 0
        dE = np.array(d['dE']); hd = np.array(d['hd'])
        mE = float(dE.mean()) if len(dE) else float('nan')
        meE = float(np.median(dE)) if len(dE) else float('nan')
        mH = float(hd.mean()) if len(hd) else float('nan')
        meH = float(np.median(hd)) if len(hd) else float('nan')
        summ[f'{a:g}'] = dict(fail_rate=fr, n_wrong=d['wrong'], wrong_dE=mE,
                              wrong_dE_med=meE, wrong_HD=mH, wrong_HD_med=meH)
        print(f'{a:>8g} {fr*100:>6.1f}% {mE:>8.2f} {meE:>7.2f} {mH:>8.1f} {meH:>7.0f}')

    # (b) α* 대략: α=0 baseline 대비 wrongΔE가 1.2배↑ 되는 첫 α (>0)
    base = summ['0']['wrong_dE']
    astar = next((a for a in ALPHAS if a > 0 and summ[f'{a:g}']['wrong_dE'] > 1.2 * base
                  and not np.isnan(summ[f'{a:g}']['wrong_dE'])), None)
    print(f'\n(b) 대략 lift-off α* ≈ {astar}  (baseline wrongΔE@α=0 = {base:.3f})')
    print(f'(c) QPU 실패율: 위 표 — 낮으면 집단이 QPU에선 안 갇힘(생존편향 경고)')

    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    out = f'results/qpu_pilot_{stamp}'; os.makedirs(out, exist_ok=True)
    with open(f'{out}/samples.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['alpha','instance_id','read_idx','hd','energy','delta_energy','is_target','is_R_gs','eng_ok'])
        w.writerows(rows)
    json.dump({'solver': qpu.solver.name, 'n_inst': len(valid), 'alphas': ALPHAS,
               'num_reads': NUM_READS, 'anneal_us': ANNEAL_US, 'gauges': GAUGES,
               'cmax': cmax, 'alpha_star_approx': astar, 'rows': summ},
              open(f'{out}/summary.json', 'w'), indent=2)
    print(f'\nsaved → {out}/')


if __name__ == '__main__':
    main()
