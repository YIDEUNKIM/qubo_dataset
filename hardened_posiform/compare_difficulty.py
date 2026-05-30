"""Three ΔE definitions compared, to show that including correct samples
distorts (hides) the difficulty of the hard cases.

(A) all-samples mean ΔE      = current fig_traj — conflates success-rate & trap depth
(B) wrong-only, varying set  = conditional trap depth | failure (population shifts with α)
(C) wrong-only, fixed set = (B) on a fixed always-trapped set (α≤0.02 only)

α=0.1 excluded (per goal). Reads (A,B) from fig_traj_summary.json, (C) from the
latest fixed_cohort_*/summary.json.
"""
import glob, json, os
import numpy as np
import matplotlib.pyplot as plt

SUM = 'ieee_paper_2page/figures/fig_traj_summary.json'
js = json.load(open(SUM))
rows = [r for r in js['rows'] if r['alpha'] != 0.1]          # exclude 0.1
alphas = [r['alpha'] for r in rows]

allmean = [r['dE_mean'] for r in rows]
allsem  = [r['dE_sem'] for r in rows]
trmean  = [r['dE_mean_trapped'] for r in rows]
trfrac  = [r['trapped_frac'] for r in rows]

# fixed cohort (C) — latest run
fc_dir = sorted(glob.glob('results/fixed_cohort_*'))[-1]
fc = json.load(open(os.path.join(fc_dir, 'summary.json')))
fc381 = fc['C381']                                            # {alpha_str: {...}}
def fcval(a):
    k = f'{a:g}'
    return fc381[k]['dE_mean'] if k in fc381 else np.nan
fcmean = [fcval(a) for a in alphas]

# ---- table ----
print(f'{"α":>8} {"trap%":>6} | {"(A)all":>8} {"(B)wrong":>9} {"(C)fixed":>9}')
print('-'*48)
for i, a in enumerate(alphas):
    c = f'{fcmean[i]:.3f}' if not np.isnan(fcmean[i]) else '   —'
    print(f'{a:>8g} {trfrac[i]*100:>6.1f} | {allmean[i]:>8.3f} {trmean[i]:>9.3f} {c:>9}')

# ---- figure ----
x = np.arange(len(alphas))
fig, ax = plt.subplots(figsize=(4.2, 2.8))
ax.errorbar(x, allmean, yerr=allsem, fmt='s--', color='gray', capsize=2,
            markersize=4, lw=1.3, label='(A) all samples (current)')
ax.plot(x, trmean, 'o-', color='tab:red', markersize=4, lw=1.6,
        label='(B) wrong-only (varying set)')
ax.plot(x, fcmean, 'D-', color='tab:blue', markersize=4, lw=1.6,
        label='(C) wrong-only (fixed set)')
ax.set_yscale('symlog', linthresh=1e-2)
ax.set_ylim(bottom=0)
ax.set_xticks(x); ax.set_xticklabels([f'{a:g}' for a in alphas], rotation=45, fontsize=7)
ax.set_xlabel(r'$\alpha$', fontsize=9)
ax.set_ylabel(r'mean $\Delta E$ ($\Delta_R$)', fontsize=9)
ax.set_title('Difficulty of the hard cases is hidden by correct samples', fontsize=8.5)
ax.grid(alpha=0.3, which='both')
ax.legend(fontsize=6.8, loc='lower left')
# annotate trapped fraction on the (B) curve at the cliff
for i, a in enumerate(alphas):
    if a in (0.04, 0.07):
        ax.annotate(f'{trfrac[i]*100:.1f}% wrong', (x[i], trmean[i]),
                    textcoords='offset points', xytext=(0, 8), fontsize=6,
                    ha='center', color='tab:red')
plt.tight_layout()
out = os.path.join(fc_dir, 'fig_compare_difficulty')
for ext in ('png', 'pdf'):
    fig.savefig(out + '.' + ext, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'\nsaved → {out}.png/pdf')
