# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository generates **QUBO (Quadratic Unconstrained Binary Optimization) benchmark datasets** with known optimal solutions. The core research question (from the PDFs): given a target binary string x*, can we construct a QUBO problem where x* is the ground state, yet the problem is statistically indistinguishable from a randomly generated QUBO? This enables evaluation of QUBO solvers (quantum annealers, QAOA, simulated annealing) by providing problems with verified answers.

## Key Concepts

- **QUBO**: Minimize E(x) = x^T Q x over binary variables x_i in {0,1}. NP-hard in general.
- **Zero Expectation**: The generated Q matrix coefficients should have E[q_ii] = 0 and E[q_ij] = 0, making the problem indistinguishable from random QUBO. Achieved by LP-optimized penalty ratios (e.g., sampling E[r0], E[r0'] at 3x other values).
- **Penalty Construction (from PDF Step 1 & 2)**: For each non-optimal (x_i, x_j) pair, add r * indicator_function as penalty. Three penalty terms per qubit pair, each with independent random r > 0.
- **Ising Conversion**: QUBO variables x_i in {0,1} map to spins s_i in {-1,+1} via s = 2x - 1.

## Architecture

### QUBO Generators (core)
- **`qubo_generator.py`** - Basic dense QUBO generator. Adds pairwise interactions using (x_i - x_j)^2 or (x_i + x_j - 1)^2 penalty terms. Good for small problems (n <= 20).
- **`qubo_generator_claude.py`** - PDF-faithful implementation. Assigns independent penalty per non-optimal state per qubit pair. Supports `balance_expectation` flag for 3x ratio sampling.
- **`qubo_zero_expectation.py`** - **Primary generator.** Uses Strategy Pattern with `PenaltyModel` ABC:
  - `DefaultZeroExpectationModel` - LP-derived ratios (e.g., target=(0,0): penalties {(0,1):1.0, (1,0):1.0, (1,1):1.65})
  - `SimpleUniformModel` - Equal penalty for all wrong states (baseline comparison)
  - `create_qubo_precise(target, density, model)` - Main entry point
  - `create_qubo_ising_derived(target)` - Alternative: derives QUBO from Ising model (J_ij = alpha * s_i * s_j), guarantees E[row sum] = 0
- **`qubo_hard_mode.py`** - Backbone + frustration model. Star graph (W_STRONG=20.0) ensures ground state; random weak edges (W_WEAK=0.2) with noise_ratio probability of frustration.
- **`qubo_wishart.py`** - **Wishart Planted Ensemble generator.** Constructs SA-hard QUBO via orthogonal Gaussian projection (W^T t = 0). Difficulty controlled by alpha=M/N parameter. Phase transition at alpha_c ≈ 0.95 (N=100). See `WISHART_EXPERIMENT.md` for details.

### Solvers & Converters
- **`dwave_simple_test.py`** - Loads saved QUBO from `qubo_results/` edge-list files, solves via `neal.SimulatedAnnealingSampler`.
- **`qubo_to_qaoa.py`** - Converts QUBO to Ising Hamiltonian (h, J, offset), generates Qiskit QAOA script.
- **`run_qaoa_n5.py`** - Generated Qiskit QAOA script for a 5-qubit example.

### Analysis & Verification
- **`check_zero_expectation.py`** - Validates that row/column means of Q matrix converge to 0.
- **`compare_experiment.py`** - Batch comparison: generates QUBO, solves with SA, checks exact/energy/symmetry match. Uses `neal`.
- **`test_hard_mode.py`** - Same as compare_experiment but for hard mode QUBO.
- **`test_wishart.py`** - SA experiment framework for Wishart ensemble: alpha sweep, N scaling, Wishart vs Hard Mode comparison, hardness metrics (TTS, spectral gap).
- **`visualize_qubo.py`** - 3D surface plot of Q matrix using matplotlib.

### Data
- **`qubo_results/`** - Saved QUBO problems as CSV edge-lists (format: `# target,<bitstring>\ni,j,weight\n...`).

## Commands

```bash
# Generate a QUBO with default target "11000101010001101"
python3 qubo_zero_expectation.py

# Generate with specific binary target
python3 qubo_zero_expectation.py 10110

# Generate with random target of length N
python3 qubo_zero_expectation.py 50

# Generate with Ising-derived row-balanced mode
python3 qubo_zero_expectation.py 10110 balance

# Generate hard mode QUBO (args: target_or_length [noise_ratio])
python3 qubo_hard_mode.py 11001 0.1

# Solve a saved QUBO file with simulated annealing (requires neal)
python3 dwave_simple_test.py qubo_results/qubo_00100_50.txt

# Run SA comparison experiment (args: n_bits [num_runs] [balance])
python3 compare_experiment.py 500 4
python3 compare_experiment.py 500 4 balance

# Run hard mode experiment (args: noise_ratio [density])
python3 test_hard_mode.py 0.1 1.0

# Generate Wishart planted ensemble QUBO (args: target alpha [seed])
python3 qubo_wishart.py 10110 0.7
python3 qubo_wishart.py 10110 0.7 42

# Run Wishart alpha sweep experiment (args: n_bits [num_runs])
python3 test_wishart.py 100 10

# Run Wishart scaling experiment (args: --scaling [alpha])
python3 test_wishart.py --scaling 0.7

# Run Wishart vs Hard Mode comparison
python3 test_wishart.py --compare

# Check zero-expectation property of a saved QUBO
python3 check_zero_expectation.py

# Convert QUBO to QAOA/Ising and generate Qiskit script
python3 qubo_to_qaoa.py 10110

# Visualize QUBO matrix (args: n_bits)
python3 visualize_qubo.py 50
```

## Dependencies

- **Core**: `numpy` (coefficient analysis)
- **Solvers**: `neal` (D-Wave simulated annealing), `dimod` (D-Wave ecosystem)
- **QAOA**: `qiskit`, `qiskit-algorithms` (optional, for quantum circuit simulation)
- **Visualization**: `matplotlib`

## Q Matrix Format

Q is stored as a Python dict `{(i, j): weight}` where i <= j (upper triangular). Edge-list files use CSV: `i,j,weight` with a `# target,<bitstring>` header. The `calculate_energy(x_str, Q)` function is shared across modules.

## Language

Code comments and print statements are primarily in Korean. Variable names and docstrings mix Korean and English.
