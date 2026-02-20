# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository generates **QUBO (Quadratic Unconstrained Binary Optimization) benchmark datasets** with known optimal solutions. The core research question (from the PDFs): given a target binary string x*, can we construct a QUBO problem where x* is the ground state, yet the problem is statistically indistinguishable from a randomly generated QUBO? This enables evaluation of QUBO solvers (quantum annealers, QAOA, simulated annealing) by providing problems with verified answers.

## Key Concepts

- **QUBO**: Minimize E(x) = x^T Q x over binary variables x_i in {0,1}. NP-hard in general.
- **Zero Expectation**: The generated Q matrix coefficients should have E[q_ii] = 0 and E[q_ij] = 0, making the problem indistinguishable from random QUBO. Achieved by LP-optimized penalty ratios (e.g., sampling E[r0], E[r0'] at 3x other values).
- **Penalty Construction (from PDF Step 1 & 2)**: For each non-optimal (x_i, x_j) pair, add r * indicator_function as penalty. Three penalty terms per qubit pair, each with independent random r > 0.
- **Ising Conversion**: QUBO variables x_i in {0,1} map to spins s_i in {-1,+1} via s = 2x - 1.
- **Quiet Planting**: Planted random 3-SAT where target satisfies all clauses. Rosenberg reduction linearizes cubic penalty z1*z2*z3 into QUBO with auxiliary variable y=z1*z2. Clause density alpha=m/n; alpha < 3.86 ensures statistical indistinguishability from random 3-SAT (Krzakala & Zdeborova, 2009). **Degeneracy problem**: all SAT solutions have identical QUBO energy (penalty=0), so clause weights alone cannot distinguish the planted target. A **planted field** (small linear bias per variable) is required to break degeneracy and make the target uniquely recoverable.

## Directory Structure

```
qubo_dataset/
├── qubo_utils.py                           ← 공유 유틸리티 (calculate_energy, save_qubo_edgelist, print_q_matrix, print_qubo_formula)
├── CLAUDE.md
├── README.md
├── .gitignore
├── docs/                                   ← 실험 결과 문서
│
├── zero_expectation/                       ← Zero Expectation 방법론
│   ├── qubo_zero_expectation.py            ← 생성기 (Strategy Pattern + PenaltyModel)
│   ├── test_zero_expectation.py            ← SA 스케일링 실험
│   ├── test_diagonal_zero.py              ← 대각 편향 분석 + 불가능성 증명
│   ├── analyze_q_structure.py             ← Q 행렬 구조 비교 (ZeroExp vs Wishart)
│   └── results/                           ← 생성된 QUBO 파일
│
├── hard_mode/                              ← Hard Mode (Backbone + Frustration)
│   └── qubo_hard_mode.py                  ← 생성기
│
├── wishart/                                ← Wishart Planted Ensemble
│   ├── qubo_wishart.py                    ← 생성기
│   ├── test_wishart.py                    ← SA 실험 (alpha sweep, scaling, comparison)
│   └── results/                           ← 생성된 QUBO 파일
│
├── quiet_planting/                         ← Quiet Planting (3-SAT → Rosenberg)
│   ├── qubo_quiet_planted.py              ← 생성기
│   ├── test_quiet_planted.py              ← SA 실험 (alpha sweep, scaling, 4-way comparison)
│   └── results/                           ← 생성된 QUBO 파일
│
├── posiform/                               ← Posiform Planting (2-SAT → Posiform)
│   ├── qubo_posiform.py                   ← 생성기
│   ├── test_posiform.py                   ← SA 실험 (scaling, coeff sweep, 5-way comparison)
│   └── results/                           ← 생성된 QUBO 파일
│
└── posiform_hardened/                      ← Hardened Posiform (Random QUBO + Posiform overlay)
    ├── qubo_posiform_hardened.py           ← 생성기 (Pelofske 2024)
    ├── test_posiform_hardened.py           ← SA 실험 (sweep, scaling, comparison)
    └── results/                           ← 생성된 QUBO 파일
```

## Architecture

### Shared Utilities
- **`qubo_utils.py`** — 모든 방법론에서 공통으로 사용하는 함수:
  - `calculate_energy(x, Q)` — QUBO 에너지 계산
  - `save_qubo_edgelist(Q, filepath, target)` — Edge List CSV 저장
  - `print_q_matrix(Q, n)` — Q 행렬 그리드 출력
  - `print_qubo_formula(Q)` — QUBO 수식 출력

### QUBO Generators (core)
- **`zero_expectation/qubo_zero_expectation.py`** - **Primary generator.** Uses Strategy Pattern with `PenaltyModel` ABC:
  - `DefaultZeroExpectationModel` - LP-derived ratios (e.g., target=(0,0): penalties {(0,1):1.0, (1,0):1.0, (1,1):1.65})
  - `SimpleUniformModel` - Equal penalty for all wrong states (baseline comparison)
  - `create_qubo_precise(target, density, model)` - Main entry point
  - `create_qubo_ising_derived(target)` - Alternative: derives QUBO from Ising model (J_ij = alpha * s_i * s_j), guarantees E[row sum] = 0
- **`hard_mode/qubo_hard_mode.py`** - Backbone + frustration model. Star graph (W_STRONG=20.0) ensures ground state; random weak edges (W_WEAK=0.2) with noise_ratio probability of frustration.
- **`wishart/qubo_wishart.py`** - **Wishart Planted Ensemble generator.** Constructs SA-hard QUBO via orthogonal Gaussian projection (W^T t = 0). Difficulty controlled by alpha=M/N parameter. Phase transition at alpha_c ≈ 0.95 (N=100). See `docs/WISHART_EXPERIMENT.md` for details.
- **`quiet_planting/qubo_quiet_planted.py`** - **Quiet Planting generator.** Planted random 3-SAT → Rosenberg reduction → QUBO. Clause density alpha=m/n controls difficulty. QUBO size = n(1+alpha) (auxiliary variables). Key parameters:
  - `alpha`: clause density (3.86=condensation, 4.27=SAT threshold)
  - `clause_weight_range`: per-clause random weights (does NOT break SAT-solution degeneracy — all SAT solutions get penalty=0 regardless of weights)
  - `field_strength`: **planted field** — adds small linear bias toward target. This IS what breaks degeneracy. Recommended 0.1~1.0. SA scaling with field=0.5: N=100(100%) → N=300(70%) → N=500(20%) → N=750+(0%). See `docs/QUIET_PLANTING_EXPERIMENT.md`

- **`posiform/qubo_posiform.py`** - **Posiform Planting generator.** Hahn, Pelofske, Djidjev (2023) 기반: planted 2-SAT → posiform → QUBO. 보조변수 없음 (QUBO 크기 = n). Tarjan SCC 기반 2-SAT solver + uniqueness 검사로 target이 유일한 ground state임을 보장. Key functions:
  - `create_planted_2sat(target)`: 2단계 (random + targeted) clause 추가로 유일한 2-SAT 해 생성
  - `posiform_to_qubo(n, clauses, coeff_range)`: wrong tuple 배제 → posiform → QUBO 변환
  - `create_qubo_posiform(target, coeff_range, seed)`: 메인 진입점
  - `coeff_range`: posiform 계수 범위 (default (1.0, 3.0))
- **`posiform_hardened/qubo_posiform_hardened.py`** - **Hardened Posiform Planting generator.** Pelofske, Hahn, Djidjev (2024) 기반: 변수를 disjoint subgraph로 분할 → 각 subgraph에 random QUBO 생성 → subproblem ground state concatenate → posiform overlay. `Q_final = Σ R_i + α × P`. Key parameters:
  - `posiform_scale` (α): 작을수록 어려움 (0.01 = hardest)
  - `coeff_type`: 'lin2' ({-1,+1}) vs 'lin20' (21단계)
  - `subgraph_size`: 각 subgraph 변수 수 (default 15, brute force 가능)

### Analysis & Verification
- **`zero_expectation/test_zero_expectation.py`** - SA scaling experiment for Zero-Expectation QUBO.
- **`zero_expectation/test_diagonal_zero.py`** - Diagonal bias analysis + impossibility proof.
- **`zero_expectation/analyze_q_structure.py`** - Q matrix structure comparison (ZeroExp vs Wishart).
- **`wishart/test_wishart.py`** - SA experiment framework for Wishart ensemble: alpha sweep, N scaling, Wishart vs Hard Mode comparison, hardness metrics (TTS, spectral gap).
- **`quiet_planting/test_quiet_planted.py`** - SA experiment framework for Quiet Planting: alpha sweep (2.0–5.0), N scaling, 4-way comparison (Quiet vs Wishart vs ZeroExp vs HardMode).
- **`posiform/test_posiform.py`** - SA experiment framework for Posiform Planting: N scaling, coefficient range sweep, 5-way comparison (Posiform vs Quiet vs Wishart vs ZeroExp vs HardMode).
- **`posiform_hardened/test_posiform_hardened.py`** - SA experiment framework for Hardened Posiform: sweep transition, N scaling, hardened vs plain comparison.

### Data
- **`<method>/results/`** - 각 방법론별 생성된 QUBO 파일. CSV edge-list 형식 (`# target,<bitstring>\ni,j,weight\n...`). 각 생성기가 자신의 `results/` 디렉토리에 자동 저장.

## Commands

```bash
# Generate a QUBO with default target "11000101010001101"
python3 zero_expectation/qubo_zero_expectation.py

# Generate with specific binary target
python3 zero_expectation/qubo_zero_expectation.py 10110

# Generate with random target of length N
python3 zero_expectation/qubo_zero_expectation.py 50

# Generate with Ising-derived row-balanced mode
python3 zero_expectation/qubo_zero_expectation.py 10110 balance

# Generate hard mode QUBO (args: target_or_length [noise_ratio])
python3 hard_mode/qubo_hard_mode.py 11001 0.1

# Generate Wishart planted ensemble QUBO (args: target alpha [seed])
python3 wishart/qubo_wishart.py 10110 0.7
python3 wishart/qubo_wishart.py 10110 0.7 42

# Run Wishart alpha sweep experiment (args: n_bits [num_runs])
python3 wishart/test_wishart.py 100 10

# Run Wishart scaling experiment (args: --scaling [alpha])
python3 wishart/test_wishart.py --scaling 0.7

# Run Wishart vs Hard Mode comparison
python3 wishart/test_wishart.py --compare

# Generate Quiet Planting QUBO (args: target alpha [seed])
python3 quiet_planting/qubo_quiet_planted.py 10110 4.2
python3 quiet_planting/qubo_quiet_planted.py 10110 4.2 42

# Run Quiet Planting alpha sweep experiment (args: n_bits [num_runs])
python3 quiet_planting/test_quiet_planted.py 50 10

# Run Quiet Planting scaling experiment (args: --scaling [alpha])
python3 quiet_planting/test_quiet_planted.py --scaling 4.2

# Run 4-way comparison (Quiet vs Wishart vs ZeroExp vs HardMode)
python3 quiet_planting/test_quiet_planted.py --compare

# Generate Posiform Planting QUBO (args: target [seed])
python3 posiform/qubo_posiform.py 10110
python3 posiform/qubo_posiform.py 10110 42

# Run Posiform scaling experiment (args: --scaling [num_runs])
python3 posiform/test_posiform.py --scaling 10

# Run Posiform coefficient range sweep (args: --coeff [num_runs])
python3 posiform/test_posiform.py --coeff 10

# Run 5-way comparison (Posiform vs Quiet vs Wishart vs ZeroExp vs HardMode)
python3 posiform/test_posiform.py --compare

# Run Zero-Expectation SA scaling experiment
python3 zero_expectation/test_zero_expectation.py 10,20,50,100 10

# Run diagonal bias analysis
python3 zero_expectation/test_diagonal_zero.py

# Run Q matrix structure comparison
python3 zero_expectation/analyze_q_structure.py
```

## Dependencies

- **Core**: `numpy` (coefficient analysis)
- **Solvers**: `neal` (D-Wave simulated annealing), `dimod` (D-Wave ecosystem)
- **QAOA**: `qiskit`, `qiskit-algorithms` (optional, for quantum circuit simulation)
- **Visualization**: `matplotlib`

## Q Matrix Format

Q is stored as a Python dict `{(i, j): weight}` where i <= j (upper triangular). Edge-list files use CSV: `i,j,weight` with a `# target,<bitstring>` header. The `calculate_energy(x_str, Q)` function is in `qubo_utils.py` and shared across all modules.

## Language

Code comments and print statements are primarily in Korean. Variable names and docstrings mix Korean and English.
