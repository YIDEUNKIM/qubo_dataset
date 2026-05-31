"""SA vs QPU overlay. QPU는 results/qpu_run_i*/ 의 모든 단계를 합쳐(오답만) 집계.
(qpu_run_* 없으면 qpu_pilot_* 로 fallback.) → ieee_paper_2page/figures/fig_sa_vs_qpu.{png,pdf}
ΔE는 각 솔버의 작은 α floor로 정규화해 공통 단일축, HD는 공통 로그축.
"""
import os, glob, csv
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, ScalarFormatter, NullFormatter

# QPU: 단계 dir 전부 합산 (전역 instance_id라 disjoint)
dirs = sorted(glob.glob('results/qpu_run_i*')) or sorted(glob.glob('results/qpu_pilot_*'))
agg = defaultdict(lambda: {'dE': [], 'hd': []}); seen = set(); seenkey = set()
for d in dirs:
    p = os.path.join(d, 'samples.csv')
    if not os.path.exists(p):
        continue
    with open(p) as f:
        r = csv.reader(f); h = next(r); ci = {n: i for i, n in enumerate(h)}
        for row in r:
            key = (row[ci['instance_id']], row[ci['alpha']], row[ci['read_idx']])
            if key in seenkey:                         # (인스턴스,α,read) 중복 제거
                continue
            seenkey.add(key)
            if int(row[ci['is_target']]) != 0:         # 오답만
                continue
            a = float(row[ci['alpha']])
            agg[a]['dE'].append(float(row[ci['delta_energy']]))
            agg[a]['hd'].append(int(row[ci['hd']]))
            seen.add(int(row[ci['instance_id']]))
alphas = sorted(agg)
qE = [float(np.mean(agg[a]['dE'])) for a in alphas]
qH = [float(np.mean(agg[a]['hd'])) for a in alphas]
print(f'QPU dirs: {len(dirs)} | 인스턴스 {len(seen)}개 | α {len(alphas)}점')

# SA, QPU와 동일한 100개 인스턴스(고정 최난해 집단의 [0,100))에 대한 오답만 평균
SA = {0.0: (458.8, 0.576), 1e-9: (73.5, 1.102), 1e-5: (73.8, 1.104), 1e-4: (73.2, 1.161),
      1e-3: (68.3, 1.661), 3e-3: (58.4, 2.541), 5e-3: (49.2, 3.079), 1e-2: (30.7, 3.345),
      1.5e-2: (18.4, 2.792), 2e-2: (11.2, 2.160)}
def sa(a):
    return SA[min(SA, key=lambda k: abs(k - a))]
sH = [sa(a)[0] for a in alphas]; sE = [sa(a)[1] for a in alphas]

x = np.arange(len(alphas)); lab = [f'{a:g}' for a in alphas]
BLUE, RED = 'tab:blue', 'tab:red'
fig, (axH, axE) = plt.subplots(2, 1, figsize=(3.5, 2.6), sharex=True)

axH.plot(x, sH, 's--', color=BLUE, ms=4, lw=1.4, label='SA')
axH.plot(x, qH, 'o-', color=RED, ms=4, lw=1.4, label='QPU')
axH.set_yscale('log'); axH.set_ylim(8, 600)
axH.set_yticks([10, 100, 500])
axH.yaxis.set_major_formatter(ScalarFormatter()); axH.yaxis.set_minor_formatter(NullFormatter())
axH.set_ylabel('HD', fontsize=9)
axH.set_title('(a) HD  (mean)', fontsize=9)
axH.grid(alpha=0.3, which='both')
axH.legend(fontsize=6, loc='lower left', ncol=2, markerscale=0.8,
           handlelength=1.2, handletextpad=0.4, columnspacing=0.8, borderpad=0.3)
axH.tick_params(labelsize=8)

# ΔE를 각 솔버의 작은 α floor(평탄구간 α∈[1e-9,1e-4] 평균; α=0은 축퇴라 제외)로
# 정규화 → 공통 단일축. floor 대비 정점비가 데이터로 결정되어 모양 비교가 정직.
plat = [i for i, a in enumerate(alphas) if 1e-9 <= a <= 1e-4]
saf = float(np.mean([sE[i] for i in plat])); qpf = float(np.mean([qE[i] for i in plat]))
sEn = np.asarray(sE) / saf; qEn = np.asarray(qE) / qpf
print(f'floor: SA={saf:.3f}  QPU={qpf:.3f}  | peak: SA={sEn.max():.2f}x  QPU={qEn.max():.2f}x')
axE.axhline(1.0, color='gray', ls=':', lw=0.8)
axE.plot(x, sEn, 's--', color=BLUE, ms=4, lw=1.4)
axE.plot(x, qEn, 'o-', color=RED, ms=4, lw=1.4)
axE.set_ylim(0.4, 3.3); axE.set_yticks([1, 2, 3]); axE.set_ylabel(r'$\Delta E$ / floor', fontsize=9)
i1 = int(np.argmax(sEn)); i2 = int(np.argmax(qEn))
axE.annotate(f'{sEn[i1]:.1f}' + r'$\times$', (x[i1], sEn[i1]), textcoords='offset points',
             xytext=(0, 3), fontsize=7, color=BLUE, ha='center')
axE.annotate(f'{qEn[i2]:.1f}' + r'$\times$', (x[i2], qEn[i2]), textcoords='offset points',
             xytext=(2, 4), fontsize=7, color=RED, ha='left')
axE.set_title(r'(b) $\Delta E$ / floor (mean)', fontsize=9)
axE.set_xticks(x); axE.set_xticklabels(lab, rotation=45, fontsize=8)
axE.set_xlabel(r'$\alpha$', fontsize=9); axE.grid(alpha=0.3); axE.tick_params(labelsize=8)

plt.tight_layout(h_pad=0.5)
out = 'ieee_paper_2page/figures/fig_sa_vs_qpu'
for ext in ('png', 'pdf'):
    fig.savefig(out + '.' + ext, dpi=200, bbox_inches='tight')
plt.close(fig)
print('saved', out + '.png/pdf')
print(f'{"α":>8} {"SA_ΔE":>7} {"QPU_ΔE":>8} | {"SA_HD":>7} {"QPU_HD":>7}')
for i, a in enumerate(alphas):
    print(f'{a:>8g} {sE[i]:>7.2f} {qE[i]:>8.2f} | {sH[i]:>7.1f} {qH[i]:>7.1f}')
