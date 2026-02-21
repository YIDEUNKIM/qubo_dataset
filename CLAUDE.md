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
│   ├── qubo_posiform.py                   ← 생성기 (기본 posiform)
│   ├── test_posiform.py                   ← SA 실험 (scaling, coeff sweep, 5-way comparison)
│   └── results/                           ← 생성된 QUBO 파일
│
├── hardened_posiform/                      ← Hardened Posiform (Random QUBO + Posiform)
│   ├── qubo_posiform_hardened.py          ← 생성기 (Pelofske, Hahn, Djidjev 2024)
│   ├── test_posiform_hardened.py          ← SA 실험 (sweep 전이, N-scaling, 비교)
│   ├── papers/                            ← 참고 논문 PDF
│   └── results/                           ← 생성된 QUBO 파일
│
└── mceliece/                               ← McEliece Cryptographic QUBO
    ├── qubo_mceliece.py                   ← 생성기 (Mandrà et al. 2024)
    ├── test_mceliece.py                   ← SA 실험 (m-scaling, t-sweep, sweep 전이, 6-way 비교)
    ├── papers/                            ← 참고 논문 PDF
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
- **`hardened_posiform/qubo_posiform_hardened.py`** - **Hardened Posiform Planting generator.** Pelofske, Hahn, Djidjev (2024) 기반: discrete coefficient random QUBO + posiform QUBO 결합 (Q_final = Σ R_i + α × P). 기본 posiform의 SA-easy 한계를 극복. Ground state 유일성 수학적 보장. Key parameters:
  - `posiform_scale` (α): 작을수록 어려움 (0.01 hardest, 0.1 easier)
  - `coeff_type`: 'lin2' ({-1,+1}) vs 'lin20' ({-1,-0.9,...,0.9,1}) — lin2가 더 어려움
  - `max_subgraph_size`: brute force 가능 크기 (≤23, default 15)
  - SA sweep 전이 (N=500): lin2,α=0.01은 5000 sweeps에서도 29.8% 성공률
- **`mceliece/qubo_mceliece.py`** - **McEliece Cryptographic QUBO generator.** Mandrà et al. (arXiv:2308.09704) 기반: McEliece 공개키를 Ising 스핀 시스템으로 캐스팅 → p-body 상호작용 → Rosenberg 차수축소 → QUBO. 암호학적 보안에 기반한 hardness. Key parameters:
  - `m`: GF(2^m) 확장 차수 (N=2^m). m≥5는 차수축소 비용으로 생성이 매우 느림
  - `t`: 에러 정정 능력 (k ≈ N - m*t, d ≥ 2t+1)
  - QUBO 크기: total_vars = k + num_aux (보조변수). m=4,t=2일 때 k=8, aux≈25~40
  - `create_qubo_mceliece(target, m, t, seed)`: 메인 진입점. Returns (Q, info)
  - `extract_original_solution(sample, k)`: SA 결과에서 원래 k개 변수 추출

### Analysis & Verification
- **`zero_expectation/test_zero_expectation.py`** - SA scaling experiment for Zero-Expectation QUBO.
- **`zero_expectation/test_diagonal_zero.py`** - Diagonal bias analysis + impossibility proof.
- **`zero_expectation/analyze_q_structure.py`** - Q matrix structure comparison (ZeroExp vs Wishart).
- **`wishart/test_wishart.py`** - SA experiment framework for Wishart ensemble: alpha sweep, N scaling, hardness metrics (TTS, spectral gap).
- **`quiet_planting/test_quiet_planted.py`** - SA experiment framework for Quiet Planting: alpha sweep (2.0–5.0), N scaling, 3-way comparison (Quiet vs Wishart vs ZeroExp).
- **`posiform/test_posiform.py`** - SA experiment framework for Posiform Planting: N scaling, coefficient range sweep, 4-way comparison (Posiform vs Quiet vs Wishart vs ZeroExp).
- **`hardened_posiform/test_posiform_hardened.py`** - SA experiment framework for Hardened Posiform: sweep transition (S-curve), N-scaling, hardened vs plain comparison.
- **`mceliece/test_mceliece.py`** - SA experiment framework for McEliece Cryptographic QUBO: m-scaling (GF(2^m) 차수 vs 난이도), t-sweep (에러 정정 능력 vs 난이도), sweep transition (S-curve), 6-way comparison. 주의: m≥5는 Rosenberg 차수축소 비용으로 QUBO 생성이 매우 느림 → m=3,4로 제한.

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

# Generate Wishart planted ensemble QUBO (args: target alpha [seed])
python3 wishart/qubo_wishart.py 10110 0.7
python3 wishart/qubo_wishart.py 10110 0.7 42

# Run Wishart alpha sweep experiment (args: n_bits [num_runs])
python3 wishart/test_wishart.py 100 10

# Run Wishart scaling experiment (args: --scaling [alpha])
python3 wishart/test_wishart.py --scaling 0.7

# Generate Quiet Planting QUBO (args: target alpha [seed])
python3 quiet_planting/qubo_quiet_planted.py 10110 4.2
python3 quiet_planting/qubo_quiet_planted.py 10110 4.2 42

# Run Quiet Planting alpha sweep experiment (args: n_bits [num_runs])
python3 quiet_planting/test_quiet_planted.py 50 10

# Run Quiet Planting scaling experiment (args: --scaling [alpha])
python3 quiet_planting/test_quiet_planted.py --scaling 4.2

# Run 3-way comparison (Quiet vs Wishart vs ZeroExp)
python3 quiet_planting/test_quiet_planted.py --compare

# Generate Posiform Planting QUBO (args: target [seed])
python3 posiform/qubo_posiform.py 10110
python3 posiform/qubo_posiform.py 10110 42

# Run Posiform scaling experiment (args: --scaling [num_runs])
python3 posiform/test_posiform.py --scaling 10

# Run Posiform coefficient range sweep (args: --coeff [num_runs])
python3 posiform/test_posiform.py --coeff 10

# Run 4-way comparison (Posiform vs Quiet vs Wishart vs ZeroExp)
python3 posiform/test_posiform.py --compare

# Generate Hardened Posiform QUBO (args: n [coeff_type] [posiform_scale] [seed])
python3 hardened_posiform/qubo_posiform_hardened.py 500 lin2 0.01 42
python3 hardened_posiform/qubo_posiform_hardened.py 100 lin20 0.1

# Run Hardened Posiform sweep transition experiment
python3 hardened_posiform/test_posiform_hardened.py --sweep 10

# Run Hardened Posiform N-scaling experiment
python3 hardened_posiform/test_posiform_hardened.py --scaling 10

# Run Hardened vs Plain Posiform comparison
python3 hardened_posiform/test_posiform_hardened.py --compare 10

# Generate McEliece Cryptographic QUBO (args: target [m] [t] [seed])
python3 mceliece/qubo_mceliece.py 10110
python3 mceliece/qubo_mceliece.py 10110 4 2 42

# Run McEliece m-scaling experiment (args: --m-scaling [num_runs])
python3 mceliece/test_mceliece.py --m-scaling 10

# Run McEliece t-parameter sweep (args: --t-sweep [num_runs])
python3 mceliece/test_mceliece.py --t-sweep 10

# Run McEliece SA sweep transition (args: --sweep [num_instances])
python3 mceliece/test_mceliece.py --sweep 10

# Run 6-way comparison (McEliece vs Hardened vs Posiform vs Quiet vs Wishart vs ZeroExp)
python3 mceliece/test_mceliece.py --compare 10

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
