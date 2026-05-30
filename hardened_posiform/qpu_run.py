"""QPU 본런 (단계별/resumable) — SA 300 동일 집단 직접 비교용 (Advantage_system6.4).

100개는 300개의 부분집합이라 단계별로 돌릴 수 있다:
  단계1: INST_START=0,   INST_END=100   (먼저 빠르게)
  단계2: INST_START=100, INST_END=300   (이어서, 앞 100개 재실행 X)
각 단계는 results/qpu_run_i{START}-{END}_<stamp>/ 에 저장. instance_id는 전역 인덱스라
단계들을 그대로 합칠 수 있다(make_overlay.py가 qpu_run_* 전부 합산).

프로토콜: 고정 스케일(auto_scale=False, **전체 300 기준 공통 SCALE** — 단계 무관 일관)
+ spin-reversal 게이지 평균 + ΔE 원래 Q_α로 numpy 재계산 + SA 동일 9칼럼.
⚠ Leap 토큰 필요. wall-clock은 제출 수(=인스턴스×α×GAUGES)에 비례.
"""
import os, csv, json, pickle, datetime
from collections import defaultdict
import numpy as np

PKL        = 'instances/instances_pegasus16_lin2_hard300.pkl'
COHORT_N   = 300                 # 전체 집단 크기(고정 SCALE 계산 기준 — 단계 무관 일관)
INST_START = 100                 # ★ 이 단계 인스턴스 범위 [START, END)
INST_END   = 300
ALPHAS     = [0.0, 1e-9, 1e-5, 1e-4, 1e-3, 3e-3, 5e-3, 1e-2, 1.5e-2, 2e-2]  # SA와 동일
NUM_READS  = 100
ANNEAL_US  = 20    # ★ 어닐링 시간(µs). 20=기본·단일 budget. T-sweep은 이 값만 바꿔 재실행.
GAUGES     = 2     # 1=게이지 없음(제출 1회/문제, 빠름). ≥2=composite(제출 G배, wall-clock↑)
FIXED_SCALE = None # None=auto_scale(문제별 range 채움, 권장·본런). 숫자=고정 스케일(정밀도 바닥 대조용, 예 0.641)


def build_Q(inst, a):
    Q = dict(inst['R_sum'])
    for k, v in inst['P'].items():
        Q[k] = Q.get(k, 0) + a * v
    return {k: v for k, v in Q.items() if abs(v) > 1e-15}


def arrs(D, nidx):
    ii = np.fromiter((nidx[i] for (i, j) in D), int)
    jj = np.fromiter((nidx[j] for (i, j) in D), int)
    ww = np.fromiter(D.values(), float)
    return ii, jj, ww


def main():
    data = pickle.load(open(PKL, 'rb'))
    cohort = data['instances'][:COHORT_N]
    stage = data['instances'][INST_START:INST_END]
    # 솔버: 6.4(pkl meta) 우선해 일관성 유지, 없으면 Advantage_system6로 폴백
    # (2026-05 분기 점검 rename — 동일 칩·동일 그래프, active-set 검증이 안전망)
    pref = data['meta']['qpu_solver']                        # 'Advantage_system6.4'
    candidates = [pref] + [c for c in ['Advantage_system6'] if c != pref]
    from dwave.system import DWaveSampler
    qpu = None
    for nm in candidates:
        try:
            qpu = DWaveSampler(solver={'name': nm}); break
        except Exception as e:
            print(f'  solver {nm} 불가 ({type(e).__name__}) → 다음 후보')
    if qpu is None:
        raise RuntimeError(f'후보 솔버 모두 불가: {candidates}')
    sname = qpu.solver.name
    print(f'solver {sname} (후보 {candidates}) | 단계 인스턴스 [{INST_START},{INST_END}) = {len(stage)}개 '
          f'| T={ANNEAL_US}us | gauges={GAUGES}')

    # active-set 검증 (이 단계 인스턴스, 제출 전)
    qn = set(qpu.nodelist); qe = {frozenset(e) for e in qpu.edgelist}
    bad = 0
    for inst in stage:
        keys = set(inst['R_sum']) | set(inst['P'])
        if ({i for k in keys for i in k} - qn) or ({frozenset(k) for k in keys if k[0] != k[1]} - qe):
            bad += 1
    print(f'active-set: {len(stage)-bad}/{len(stage)} valid')
    if bad:
        raise RuntimeError(f'{bad} insts invalid — calibration drift. 중단.')

    # 스케일: 기본 auto_scale(문제별 hardware range 채움, 본런). FIXED_SCALE 지정 시 고정(정밀도 바닥 대조).
    if FIXED_SCALE is None:
        SCALE = None
        print('scale: auto_scale (문제별 range 채움)')
    else:
        SCALE = float(FIXED_SCALE)
        cmax = max(abs(v) for inst in cohort for a in ALPHAS for v in build_Q(inst, a).values())
        print(f'scale: 고정 SCALE={SCALE:.4f} (참고 cmax over {COHORT_N}={cmax:.3f})')

    if GAUGES >= 2:
        from dwave.preprocessing import SpinReversalTransformComposite
        sampler = SpinReversalTransformComposite(qpu); skw = dict(num_spin_reversal_transforms=GAUGES)
    else:
        sampler = qpu; skw = {}

    def submit(Q):
        if SCALE is None:
            return sampler.sample_qubo(Q, num_reads=NUM_READS, annealing_time=ANNEAL_US,
                                       auto_scale=True, answer_mode='raw', **skw)
        return sampler.sample_qubo({k: SCALE * v for k, v in Q.items()}, num_reads=NUM_READS,
                                   annealing_time=ANNEAL_US, auto_scale=False,
                                   answer_mode='raw', **skw)

    submit(build_Q(stage[0], ALPHAS[-1]))       # 긴 런 전 1문제 사전 테스트
    print('  사전 제출 테스트 OK (게이지·고정스케일·범위 통과)')

    rows = []
    acc = defaultdict(lambda: {'dE': [], 'hd': [], 'n': 0, 'w': 0})
    for li, inst in enumerate(stage):
        gid = INST_START + li                    # 전역 instance_id
        nodes = inst['sorted_nodes']; nidx = {nd: k for k, nd in enumerate(nodes)}
        Rii, Rjj, Rww = arrs(inst['R_sum'], nidx)
        Pii, Pjj, Pww = arrs(inst['P'], nidx)
        tgt = np.fromiter((int(inst['target_str'][k]) for k in range(len(nodes))), int)
        tr, tp = inst['t_energy_r'], inst['t_energy_p']
        for a in ALPHAS:
            resp = submit(build_Q(inst, a))
            ri = 0
            for s, _e, occ in resp.data(['sample', 'energy', 'num_occurrences']):
                x = np.fromiter((s[nd] for nd in nodes), int)
                Re = float((Rww * x[Rii] * x[Rjj]).sum())
                Pe = float((Pww * x[Pii] * x[Pjj]).sum())
                dE = (Re - tr) + a * (Pe - tp)
                hd = int((x != tgt).sum())
                e = Re + a * Pe
                is_t = int(hd == 0); is_R = int(abs(Re - tr) < 1e-7); ok = int(abs(dE) < 1e-7)
                for _ in range(int(occ)):
                    if ri >= NUM_READS:
                        break
                    rows.append([a, gid, ri, hd, e, dE, is_t, is_R, ok])
                    acc[a]['n'] += 1
                    if not is_t:
                        acc[a]['w'] += 1; acc[a]['dE'].append(dE); acc[a]['hd'].append(hd)
                    ri += 1
                if ri >= NUM_READS:
                    break
        if (li + 1) % 20 == 0:
            print(f'  {li+1}/{len(stage)} done')

    print(f'\n{"α":>8} {"실패율":>7} {"meanΔE":>8} {"medΔE":>7} {"meanHD":>8} {"medHD":>7}')
    summ = {}
    for a in ALPHAS:
        d = acc[a]; fr = d['w'] / d['n'] if d['n'] else 0
        dE = np.array(d['dE']); hd = np.array(d['hd'])
        g = lambda v, fn: float(fn(v)) if len(v) else float('nan')
        summ[f'{a:g}'] = dict(fail_rate=fr, n_wrong=d['w'],
                              wrong_dE=g(dE, np.mean), wrong_dE_med=g(dE, np.median),
                              wrong_HD=g(hd, np.mean), wrong_HD_med=g(hd, np.median))
        r = summ[f'{a:g}']
        print(f'{a:>8g} {fr*100:>6.1f}% {r["wrong_dE"]:>8.2f} {r["wrong_dE_med"]:>7.2f} '
              f'{r["wrong_HD"]:>8.1f} {r["wrong_HD_med"]:>7.0f}')

    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    out = f'results/qpu_run_i{INST_START}-{INST_END}_{stamp}'; os.makedirs(out, exist_ok=True)
    with open(f'{out}/samples.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['alpha', 'instance_id', 'read_idx', 'hd', 'energy', 'delta_energy',
                    'is_target', 'is_R_gs', 'eng_ok'])
        w.writerows(rows)
    json.dump({'solver': qpu.solver.name, 'inst_range': [INST_START, INST_END],
               'alphas': ALPHAS, 'num_reads': NUM_READS, 'anneal_us': ANNEAL_US,
               'gauges': GAUGES, 'scale_mode': 'auto' if FIXED_SCALE is None else 'fixed',
               'scale': SCALE, 'rows': summ},
              open(f'{out}/summary.json', 'w'), indent=2)
    print(f'\nsaved → {out}/   (다음 단계: INST_START={INST_END}, INST_END=...로 재실행)')


if __name__ == '__main__':
    main()
