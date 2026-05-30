"""Fixed roster of 300 hardest instances, wrong-only ΔE/HD across α (excl 0.1).

Roster = the 300 instances with the most total wrong samples (the ones that stay
trapped longest) — a CONSTANT instance set across α (satisfies the "일정" req).
At each α we take only these 300's wrong (is_target==0) samples. The wrong-sample
count within the fixed roster shrinks at high α (reported) — that residual tail is
the failures of the hardest instances.
"""
import csv, os
from collections import defaultdict, Counter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, ScalarFormatter, NullFormatter

RUN = 'results/gsp_fullrun500_alphas_lin2_inst500_20260521_173814'
SAMP = os.path.join(RUN, 'samples.csv')
N = 5612
NCOH = 300

# pass 1: per-instance total wrong count (excl 0.1) → rank hardest
tot_wrong = Counter()
rows_cache = []  # (a, iid, hd, dE, is_tgt) — cache to avoid 2nd disk pass
with open(SAMP) as f:
    r = csv.reader(f); h = next(r); ci = {n: i for i, n in enumerate(h)}
    for row in r:
        a = float(row[ci['alpha']])
        if a == 0.1:
            continue
        iid = int(row[ci['instance_id']]); tgt = int(row[ci['is_target']])
        hd = int(row[ci['hd']]); dE = float(row[ci['delta_energy']])
        rows_cache.append((a, iid, hd, dE, tgt))
        if tgt == 0:
            tot_wrong[iid] += 1

roster = set(i for i, _ in tot_wrong.most_common(NCOH))
print(f'roster: {len(roster)} hardest instances (by total wrong samples)')

# merge extra α from extend run (cohort ranking already fixed on original data)
import glob
for ed in sorted(glob.glob('results/gsp_extend500peak_*')):
    ep = os.path.join(ed, 'samples.csv')
    if not os.path.exists(ep):
        continue
    with open(ep) as f:
        r = csv.reader(f); h = next(r); ci2 = {n: i for i, n in enumerate(h)}
        for row in r:
            a = float(row[ci2['alpha']])
            rows_cache.append((a, int(row[ci2['instance_id']]), int(row[ci2['hd']]),
                               float(row[ci2['delta_energy']]), int(row[ci2['is_target']])))
    print(f'  merged extra αs from {ed}')

# pass 2: per-α wrong-only stats over the fixed roster
ens = defaultdict(lambda: {'hd': [], 'dE': []})
fully_wrong = defaultdict(Counter)  # alpha -> {iid: wrong_count}
for a, iid, hd, dE, tgt in rows_cache:
    if iid not in roster:
        continue
    if tgt == 0:
        ens[a]['hd'].append(hd); ens[a]['dE'].append(dE); fully_wrong[a][iid] += 1

alphas = sorted(ens)
def col(key, fn): return np.array([fn(np.array(ens[a][key])) for a in alphas])
hd_mean, hd_med = col('hd', np.mean), col('hd', np.median)
dE_mean, dE_med = col('dE', np.mean), col('dE', np.median)
dE_sem = col('dE', lambda x: x.std(ddof=1)/np.sqrt(len(x)))
nwrong = [len(ens[a]['hd']) for a in alphas]
n_inst_anyfail = [len(fully_wrong[a]) for a in alphas]
n_inst_allfail = [sum(1 for c in fully_wrong[a].values() if c == 100) for a in alphas]

print(f'\n{"α":>8} {"오답샘플":>8} {"실패인스턴스":>10} {"100%오답":>8} {"meanΔE":>8} {"medΔE":>6} {"meanHD":>8} {"medHD":>6}')
for i, a in enumerate(alphas):
    print(f'{a:>8g} {nwrong[i]:>8d} {n_inst_anyfail[i]:>10d}/{NCOH} {n_inst_allfail[i]:>8d} '
          f'{dE_mean[i]:>8.3f} {dE_med[i]:>6.2f} {hd_mean[i]:>8.1f} {hd_med[i]:>6.0f}')

# figure — only α where the fixed cohort is genuinely trapped (≥~90% reads fail)
PLOT_MAX = 0.02
pa = [a for a in alphas if a <= PLOT_MAX + 1e-15]
pi = [alphas.index(a) for a in pa]
hd_mean, hd_med = hd_mean[pi], hd_med[pi]
dE_mean, dE_med, dE_sem = dE_mean[pi], dE_med[pi], dE_sem[pi]
alphas = pa
x = np.arange(len(alphas)); lab = [f'{a:g}' for a in alphas]
fig, (axL, axR) = plt.subplots(2, 1, figsize=(3.6, 2.6), sharex=True)
axL.plot(x, hd_mean, 's--', color='tab:blue', lw=1.5, ms=4, label='mean')
axL.plot(x, hd_med, 'o-', color='red', lw=1.5, ms=4, label='median')
axL.set_yscale('log'); axL.set_ylim(8, 600)
axL.yaxis.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0)))   # 10,20,50,100,200,500
axL.yaxis.set_major_formatter(ScalarFormatter())                # plain integers, not 10^x
axL.yaxis.set_minor_formatter(NullFormatter())
axL.set_ylabel('HD', fontsize=9)
axL.set_title('(a) HD of failures', fontsize=8.5)
axL.grid(alpha=0.3, which='both'); axL.legend(fontsize=4.4, ncol=2, loc='upper right')
axL.tick_params(labelsize=7)
axR.errorbar(x, dE_mean, yerr=dE_sem, fmt='s--', color='tab:blue', lw=1.5, ms=4, capsize=2, label='mean')
axR.plot(x, dE_med, 'o-', color='red', lw=1.5, ms=4, label='median')
axR.axhline(1.0, color='gray', ls=':', lw=0.8)               # ΔR floor (1 R-gap; labelled in caption)
pk = int(np.nanargmax(dE_mean))
axR.annotate(f'peak {dE_mean[pk]:.2f}', xy=(x[pk], dE_mean[pk]), xytext=(x[pk]-2.4, 0.55),
             fontsize=7, ha='left', va='center', arrowprops=dict(arrowstyle='->', lw=0.7))
axR.set_ylim(0, dE_mean[pk] * 1.12)                          # LINEAR y — peak now prominent
axR.set_yticks([0, 1, 2, 3])
axR.set_xticks(x); axR.set_xticklabels(lab, rotation=45, fontsize=7)
axR.set_xlabel(r'$\alpha$', fontsize=9); axR.set_ylabel(r'$\Delta E$', fontsize=9)
axR.set_title(r'(b) $\Delta E$ of failures', fontsize=8.5)
axR.grid(alpha=0.3, which='both'); axR.tick_params(labelsize=7)
plt.tight_layout(h_pad=0.5)
fig_out = 'ieee_paper_2page/figures/fig_traj'
for ext in ('png', 'pdf'):
    fig.savefig(fig_out + '.' + ext, dpi=200, bbox_inches='tight')
plt.close(fig)
# dump cohort roster (original instance_ids) for reuse (pkl, future runs)
import json
with open('results/hard300_roster.json', 'w') as f:
    json.dump({'roster_ids': sorted(roster), 'ncohort': NCOH,
               'selection': 'top-300 by total wrong samples (excl α=0.1) in gsp_fullrun500',
               'src_run': RUN, 'plot_alphas': [float(a) for a in alphas]}, f)
print(f'\nsaved → {fig_out}.png/pdf ; roster → results/hard300_roster.json')
