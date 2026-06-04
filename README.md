# Posiform-Planted QUBO Benchmarks — Paper Data & Code

Reproduction package for the paper

> **Difficulty Analysis of Posiform-Planted QUBO Benchmarks via Energy and Hamming Metrics**
> Yideun Kim, Dongjae Lee — Kangwon National University

This repository is intentionally minimal: it contains **only the data and code behind the
paper's figure and reported numbers**, so the results can be independently checked. It also
publishes the benchmark instance set used in the experiments for reuse by other solvers.

The paper studies the *hardened posiform planting* construction of Pelofske, Hahn & Djidjev
(npj Unconventional Computing **2**:17, 2025), where a random ±1 spin-glass QUBO
`R = Σ Rᵢ` is combined with a posiform-planted QUBO `P` (unique minimum `x*`) scaled by a
coefficient `α`:

```
Q_α(x) = Σ Rᵢ(x) + α·P(x)
```

Difficulty is measured **conditional on failure** (over reads that miss `x*`) on a **fixed set
of the 100 hardest instances**, via the energy gap `ΔE = E(x) − E(x*)` and the Hamming
distance `HD(x, x*)`. These resolve the uninformative 0 % success rate into a continuous
`ΔE = ΔR + α·ΔP` curve with a non-monotone peak.

---

## Repository contents

| Path | What it is |
|------|------------|
| `hardened_posiform/instances/instances_pegasus16_lin2_100.pkl` | **Published instance set** — 100 hardened-posiform instances on the D-Wave Pegasus P16 graph (`n = 5612`, `lin2` = ±1 coefficients), ≈98 MB |
| `hardened_posiform/make_overlay.py` | Regenerates the paper figure (Fig. 1, SA vs. QPU) and prints the reported floor/peak numbers |
| `hardened_posiform/results/qpu_run_i0-100_20260527_221806/` | Raw per-sample D-Wave Advantage output (`samples.csv`) + `summary.json` consumed by the figure |
| `hardened_posiform/qpu_run.py`, `qpu_pilot.py` | D-Wave Advantage sampling scripts (the QPU experiment) |
| `hardened_posiform/hardware_native_qubo.ipynb` | Notebook driving the SA + QPU experiment on these instances |

---

## 1. Reproduce the paper figure

The figure overlays the D-Wave Advantage (`Advantage_system6.4`) run on the 100 fixed
instances against the simulated-annealing (SA) baseline. The QPU per-sample data is included;
the SA means (a sweep over the same 100-hardest cohort) are embedded as a lookup table inside
`make_overlay.py`.

```bash
cd hardened_posiform
mkdir -p ieee_paper_2page/figures        # output directory the script writes to
python3 make_overlay.py
```

Expected output (matches Fig. 1 and Sec. III of the paper):

```
QPU dirs: 1 | instances 100 | α 10 points
floor: SA=1.12  QPU=23.94  | peak: SA=2.98x  QPU=1.82x
       α   SA_ΔE   QPU_ΔE |   SA_HD  QPU_HD
       0    0.58    24.28 |   458.8   501.0
   1e-09    1.10    23.85 |    73.5   499.8
    ...
    0.01    3.35    42.88 |    30.7   285.9
    0.02    2.16    42.05 |    11.2   167.1
```

* `ΔE` floor ≈ 1·Δ_R (SA) — the quantitative form of the local-minimum conjecture.
* `ΔE` peak at `α = 10⁻²`: **3.0×** the floor for SA, **1.8×** for the QPU.
* `HD` falls monotonically (SA 73 → 11, QPU 501 → 167) while `ΔE` peaks.

The figure is written to `ieee_paper_2page/figures/fig_sa_vs_qpu.{png,pdf}`.

---

## 2. Use the published instance set

Each instance stores the random part `R_sum` and the posiform part `P` separately, so any
difficulty level `α` can be assembled on demand as `Q_α = R_sum + α·P`
(smaller `α` is harder for SA). **`Q` is keyed by D-Wave Pegasus node IDs** (not `0…n−1`),
so energies are computed against the `target` map `{node: bit}`.

```python
import pickle

with open("hardened_posiform/instances/instances_pegasus16_lin2_100.pkl", "rb") as f:
    data = pickle.load(f)

print(data["meta"])           # topology, n, qpu_solver, ...
inst = data["instances"][0]   # instance 0 of 100

alpha = 0.01                  # paper's reference value
Q = dict(inst["R_sum"])
for k, v in inst["P"].items():
    Q[k] = Q.get(k, 0) + alpha * v

# Planted optimum (node -> bit) and its energy
target    = inst["target"]                               # {node: 0/1}
gs_energy = inst["t_energy_r"] + alpha * inst["t_energy_p"]

# Integrity check: energy of the target under Q_α must equal gs_energy
def energy(bits, Q):
    return sum(w * bits[i] * bits[j] for (i, j), w in Q.items())

assert abs(energy(target, Q) - gs_energy) < 1e-6
```

**Instance keys:** `R_sum`, `P` (each `{(node_i, node_j): weight}`, Pegasus node IDs) ·
`target` (`{node: bit}`, use this for energies) · `target_str` (same optimum as a compact
bitstring in `sorted_nodes` order) · `sorted_nodes` (bit position ↔ Pegasus node ID) ·
`t_energy_r`, `t_energy_p` (target energy of the `R` and `P` parts) · `seed`.
The ground state is unique **by construction** (each block `Rᵢ` solved exactly by brute force;
`P` is the posiform of a planted 2-SAT with a unique solution).

### Download

```bash
git clone https://github.com/YIDEUNKIM/qubo_dataset.git   # pkl is a normal blob (~98 MB)

# or fetch just the instance file
curl -L -o instances_pegasus16_lin2_100.pkl \
  https://raw.githubusercontent.com/YIDEUNKIM/qubo_dataset/main/hardened_posiform/instances/instances_pegasus16_lin2_100.pkl
```

---

## Dependencies

```bash
pip install numpy neal dimod matplotlib
# QPU sampling (qpu_run.py / notebook) additionally needs: pip install dwave-system
```

`make_overlay.py` and the instance loader need only `numpy` / `matplotlib`. Running
`qpu_run.py` or the notebook requires D-Wave Leap access.

## Scope and notes

* **Minimal by design.** This public repository carries only the dataset and code that the
  paper's results rest on. The broader generator/experiment code is kept outside it.
* **SA data** behind Fig. 1 is the fixed-cohort SA sweep, stored as the lookup table in
  `make_overlay.py`; the QPU raw samples are included in full under `results/`.
* The full 300/500-instance pickles used to *select* the hardest 100 are not redistributed
  (GitHub size limits); the published 100-instance subset is the fixed evaluation set and is
  sufficient to re-run any solver against a known optimum.

## Citation

```bibtex
@inproceedings{kim_posiform_difficulty,
  title     = {Difficulty Analysis of Posiform-Planted QUBO Benchmarks
               via Energy and Hamming Metrics},
  author    = {Kim, Yideun and Lee, Dongjae},
  booktitle = {IEEE International Conference on Quantum Computing and Engineering (QCE)},
  year      = {2026}
}
```
