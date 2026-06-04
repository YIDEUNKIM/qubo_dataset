# Hardened Posiform Planting — instances, experiment code, QPU data

This directory holds the dataset, experiment code, and raw QPU output behind the paper.
For the full reproduction walkthrough (figure + instance usage), see the
[top-level README](../README.md).

## Construction

Following Pelofske, Hahn & Djidjev (npj Unconventional Computing **2**:17, 2025), a random
±1 spin-glass QUBO is combined with a posiform-planted QUBO and scaled by a coefficient `α`:

```
Q_α(x) = Σ Rᵢ(x) + α·P(x)
```

* **`Rᵢ`** — the Pegasus P16 graph (`n = 5612`) is split by recursive bisection into blocks of
  ≤ 16 qubits; each block gets a random `lin2` (±1) QUBO that is solved exactly by brute force.
  The planted target `x*` is the concatenation of the block ground states, so it is the global
  minimum of `Σ Rᵢ`.
* **`P`** — the posiform of a planted 2-SAT whose only solution is `x*`.

Hence `x*` is the **unique** global minimum of `Q_α` for every `α > 0`. Smaller `α` lets the
rugged random landscape dominate, making the instance harder for simulated annealing.

## Contents

| Path | What it is |
|------|------------|
| `instances/instances_pegasus16_lin2_100.pkl` | The fixed set of 100 hardest instances (Pegasus P16, `n = 5612`, ±1 coefficients). Stores `R_sum` and `P` separately so `Q_α = R_sum + α·P` can be assembled at any `α` |
| `make_overlay.py` | Builds the paper figure (SA vs. QPU) from `results/` plus the embedded SA lookup table, and prints the floor/peak numbers |
| `qpu_run.py`, `qpu_pilot.py` | D-Wave Advantage (`Advantage_system6.4`) sampling scripts |
| `hardware_native_qubo.ipynb` | Notebook driving the SA + QPU experiment on these instances |
| `results/qpu_run_i0-100_20260527_221806/` | Raw per-sample QPU output (`samples.csv`) + `summary.json` consumed by `make_overlay.py` |

## Reproduce the figure

```bash
mkdir -p ieee_paper_2page/figures   # output directory make_overlay.py writes to
python3 make_overlay.py             # prints the SA/QPU ΔE & HD table and saves the figure
```

See the [top-level README](../README.md) for the expected output and for loading the instance
set. Re-running `qpu_run.py` or the notebook requires D-Wave Leap access; the included QPU
samples are sufficient to regenerate the figure offline.
