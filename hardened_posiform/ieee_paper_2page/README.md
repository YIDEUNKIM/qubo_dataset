# Posiform-Planted QUBO Difficulty — 2-page IEEE conference paper

## Build

```bash
cd hardened_posiform/cisc_paper_2page
pdflatex main.tex && pdflatex main.tex   # 2-pass for refs (Korean: main_ko.tex)
```

## Files

- `main.tex` / `main_ko.tex` — paper source (EN / KO, kept in sync), IEEEtran conference, 2 pages
- `figures/fig_traj.pdf` — main figure: wrong-only HD + ΔE over the fixed set of 300 hardest instances, 10 α (0–0.02)
- `figures/fig_dual_metric.pdf` — SA budget comparison (1000 vs 5000 sweeps), referenced in Discussion
- `figures/fig_compare_difficulty.pdf` — all-sample vs wrong-only vs fixed-set ΔE (rationale; used in EXPERIMENT_SUMMARY)
- `figures/fig_traj_summary.json` — frozen all-sample / varying-trapped stats (input to `compare_difficulty.py`)
- `EXPERIMENT_SUMMARY.md` / `.pdf` — plain-language experiment writeup
- `QPU_EXPERIMENT_PLAN.md` — plan to replicate on D-Wave QPU under the same fixed set and α

## Approach

Difficulty is measured **conditional on failure** (wrong-only ΔE/HD) over a **fixed set of the 300 hardest instances** (removes survivorship bias). All-sample ΔE collapses to ~0 at large α — a success-rate effect, not shallower traps — and hides the difficulty. Valid for α ≤ 0.02; beyond it these instances escape.

## Key findings

1. α=0 degeneracy fingerprint: HD maximal while ΔE median 0 (56.6% land on a degenerate R ground state)
2. ΔE = ΔR + αΔP: flat ΔR floor (≈1) through the valley → non-monotone peak (3.28, 3.3× the floor) at α=0.01
3. HD monotone descent vs ΔE peak — failures grow bit-closer yet energetically deeper
4. ΔE ranks two SA budgets (1000 vs 5000 sweeps) that the success rate scores identically (0%)

## Data & regeneration

- SA runs: `../results/gsp_fullrun500_alphas_lin2_inst500_20260521_173814/` (10 original α) +
  `../results/gsp_extend500peak_alphas_lin2_inst500_20260527_162048/` (3 peak-region α)
- Fixed set: `../instances/instances_pegasus16_lin2_hard300.pkl`; roster `../results/hard300_roster.json`

```bash
python3 ../analyze_hard300.py        # → fig_traj.{png,pdf} + hard300_roster.json
python3 ../compare_difficulty.py     # → fig_compare_difficulty.{png,pdf}
# add new α (SA): python3 ../extend_alphas.py --pkl ../instances/instances_pegasus16_lin2_500.pkl --alphas <list> --workers 4
```
