# QUBO Benchmark Dataset Generator

양자 어닐링 및 QUBO 솔버 벤치마크를 위한 **정답이 알려진 QUBO 문제 생성기**.

핵심 기능: 임의의 목표 비트스트링이 ground state임이 수학적으로 보장된 QUBO를 생성하여, 솔버의 정확도를 정량적으로 측정할 수 있다.

---

## 프로젝트 구조

```
qubo_dataset/
├── qubo_utils.py                 # 공유 유틸리티 (calculate_energy, save_qubo_edgelist 등)
├── docs/                         # 실험 보고서
│
├── zero_expectation/             # Zero Expectation (E[q_ij]=0 보장)
│   ├── qubo_zero_expectation.py
│   ├── test_zero_expectation.py
│   ├── test_diagonal_zero.py
│   ├── analyze_q_structure.py
│   └── results/
│
├── wishart/                      # Wishart Planted Ensemble (SA-hard)
│   ├── qubo_wishart.py
│   ├── test_wishart.py
│   └── results/
│
├── quiet_planting/               # Quiet Planting (3-SAT → Rosenberg)
│   ├── qubo_quiet_planted.py
│   ├── test_quiet_planted.py
│   └── results/
│
├── posiform/                     # Posiform Planting (2-SAT → Posiform)
│   ├── qubo_posiform.py
│   ├── test_posiform.py
│   └── results/
│
├── hardened_posiform/            # Hardened Posiform (Random QUBO + Posiform)
│   ├── qubo_posiform_hardened.py
│   ├── test_posiform_hardened.py
│   ├── papers/
│   └── results/
│
├── signed_posiform/              # Signed Posiform (양/음 혼합 posiform)
│   ├── qubo_signed_posiform.py
│   ├── test_signed_posiform.py
│   ├── SIGNED_POSIFORM_EXPERIMENT.md
│   └── results/
│
├── truthtable/                   # Truth Table QUBO (Möbius + Rosenberg / 근사 QP)
│   ├── qubo_truthtable.py
│   ├── test_truthtable.py
│   ├── test_approx_comparison.py
│   ├── papers/
│   └── results/
│
├── truthtable_concat/            # Truth Table Concat (block-diagonal 접합)
│   ├── qubo_truthtable_concat.py
│   ├── test_truthtable_concat.py
│   └── results/
│
├── mceliece/                     # McEliece Cryptographic QUBO
│   ├── qubo_mceliece.py
│   ├── test_mceliece.py
│   ├── papers/
│   └── results/
│
└── docs/                         # 실험 보고서 및 분석 문서
    ├── METHODOLOGY_COMPARISON.md
    ├── POSIFORM_EXPERIMENT.md
    ├── QUIET_PLANTING_EXPERIMENT.md
    ├── SIGNAL_HARDNESS_DILEMMA.md
    └── ...
```

각 방법론 디렉토리에 `results/`(생성된 QUBO 파일), `papers/`(참조 논문 PDF)가 포함됨.

## 생성 방식 비교

| 생성기 | QUBO 크기 | Ground State | SA 난이도 | 구별 불가능 | 핵심 논문 |
|--------|:---------:|:-----------:|:---------:|:----------:|----------|
| **Wishart** | n | 수학적 (유한정밀도 제외) | **SA-hard** | X (low-rank) | Hamze et al. 2020 |
| **Hardened Posiform** | n | **수학적 (유일)** | SA-moderate | X (block-diagonal) | Pelofske et al. 2024 |
| **McEliece** | k+aux | 조건부 (M 의존) | **SA-hard** (m=4) | 미분석 | Mandrà et al. 2024 |
| **Quiet Planting** | n(1+α) | 조건부 (field 필요) | SA-medium (f=0.5) | **O** (α<3.86) | Krzakala & Zdeborova 2009 |
| **Truth Table Concat** | k×h (+aux) | **수학적** | 설정 가변 | X (block-diagonal) | 내부 연구 |
| **Truth Table** | n+aux | **수학적** | 설정 가변 | X | 내부 연구 |
| **Signed Posiform** | n | **수학적 (유일)** | SA-easy~moderate | X (양/음 비율) | 내부 연구 |
| **Posiform** | n | **수학적 (유일)** | SA-easy | X (sparse, 대각편향) | Hahn et al. 2023 |
| **Zero Expectation** | n | 수학적 | **SA-trivial** | **O** (E[q_ij]=0) | 내부 연구 |

> **Ground State "수학적"**: 구성에 의해 target이 최소 에너지임이 증명됨. "조건부"는 파라미터에 따라 GS가 깨질 수 있음.

## Quick Start

```bash
# Posiform: 보조변수 없이 유일한 ground state 보장
python3 posiform/qubo_posiform.py 10110

# Hardened Posiform: posiform + random QUBO 결합으로 난이도 증가
python3 hardened_posiform/qubo_posiform_hardened.py 500 lin2 0.01 42

# Wishart: SA-hard QUBO (alpha=0.7이 가장 어려움)
python3 wishart/qubo_wishart.py 10110 0.7

# Quiet Planting: 통계적 은닉성 (alpha<3.86)
python3 quiet_planting/qubo_quiet_planted.py 10110 4.2

# McEliece: 암호학적 hardness 기반 QUBO (m=4, t=2)
python3 mceliece/qubo_mceliece.py 10110

# Signed Posiform: 양/음 혼합 weight로 Q 행렬 자연스러움 향상
python3 signed_posiform/qubo_signed_posiform.py 10110

# Truth Table: 진리표 기반 에너지 landscape 직접 설계
python3 truthtable/qubo_truthtable.py --preset random 8 10110011

# Truth Table Concat: block-diagonal 접합으로 큰 QUBO 생성
python3 truthtable_concat/qubo_truthtable_concat.py 1001111 10 --random --harden 0.01

# Zero Expectation: E[q_ij]=0 보장
python3 zero_expectation/qubo_zero_expectation.py 10110
```

## SA 실험

```bash
# Posiform N 스케일링 (10 runs)
python3 posiform/test_posiform.py --scaling 10

# Hardened Posiform sweep 전이 실험
python3 hardened_posiform/test_posiform_hardened.py --sweep 10

# Wishart alpha sweep (N=100, 10 runs)
python3 wishart/test_wishart.py 100 10

# Quiet Planting N 스케일링
python3 quiet_planting/test_quiet_planted.py --scaling 4.2

# McEliece m-scaling 실험
python3 mceliece/test_mceliece.py --m-scaling 10

# Signed Posiform negative ratio sweep
python3 signed_posiform/test_signed_posiform.py --neg-sweep 20

# Truth Table sweep 전이 (S-curve)
python3 truthtable/test_truthtable.py --sweep 100

# Truth Table Concat h-scaling (6 configs)
python3 truthtable_concat/test_truthtable_concat.py --scaling 10

# Truth Table Concat 14-way comparison
python3 truthtable_concat/test_truthtable_concat.py --compare 10

# 전체 비교 (최대 규모)
python3 mceliece/test_mceliece.py --compare 10
```

## 주요 실험 결과 요약

### SA 성공률 비교 (num_reads=200, num_sweeps=max(1000, 10×QUBO_vars))

| N | Posiform | Hardened (lin2,α=0.01) | Quiet (α=4.2, f=0.5) | Wishart (α=0.7) | ZeroExp |
|---:|:--------:|:---------------------:|:--------------------:|:-------------------:|:-------:|
| 100 | **100%** | ~100% | **100%** | ~10% | **100%** |
| 200 | **100%** | ~100% | **100%** | ~0% | **100%** |
| 300 | **100%** | ~100% | **40%** | ~0% | **100%** |
| 500 | **100%** | ~90% | **0%** | ~0% | **100%** |
| 1000 | **100%** | ~40% | 0% | ~0% | **100%** |

**SA 난이도 분류 기준** (표준 조건: num_sweeps=10×QUBO_vars, num_reads=200):
- **SA-trivial**: p(N=1000) = 100%, local minima 1개 — **Zero Expectation**
- **SA-easy**: p(N=1000) ≥ 90% — Posiform
- **SA-moderate**: 10% ≤ p(N=500) < 90% — Hardened Posiform (lin2, α=0.01)
- **SA-hard**: p(N=500) < 10% — Wishart (α=0.7), Quiet Planting (N>300), McEliece (m=4,t=2: ~0%@k=8)

**SA 실패 양상의 차이**:
- **Wishart**: Hamming ~ N/2 → SA가 metastable trap에 빠져 target과 완전히 다른 위치에 도달
- **Quiet Planting**: Hamming 1~3 → SA가 target 근처까지 접근하나 마지막 몇 비트를 못 뒤집음 (glassy)
- **Zero Expectation**: frustration = 0, 에너지 지형이 단일 funnel → SA가 항상 target 도달

## 각 생성기 상세

### Zero Expectation
Off-diagonal 기대값 E[q_ij]=0을 보장하여 무작위 QUBO와 통계적으로 구별 불가능. Strategy Pattern: `ZeroOffDiagonalModel` (default, 비율 {1,1,2}). SA-trivial (에너지 지형이 단일 funnel).

### Wishart Planted Ensemble
가우시안 직교 투영(W^T t = 0)을 통해 SA-hard QUBO 구성. α=M/N으로 난이도 조절. 상전이 α_c ≈ 0.95 (N=100).

### Quiet Planting (3-SAT → Rosenberg)
Planted random 3-SAT에서 Rosenberg 차수축소로 QUBO 변환. α < 3.86이면 random 3-SAT과 통계적으로 구별 불가능. QUBO 크기 = n(1+α) (보조변수). **Planted field** (f=0.1~1.0)이 SAT 해 축퇴를 깨뜨리는 핵심.

### Posiform Planting (2-SAT → Posiform)
Hahn, Pelofske, Djidjev (2023). Planted 2-SAT → posiform → QUBO. 보조변수 없음 (QUBO 크기 = n). Tarjan SCC 기반 유일성 보장. SA-easy이지만 GS 보장이 가장 강력.

### Hardened Posiform
Pelofske, Hahn, Djidjev (2024). Discrete coefficient random QUBO + posiform QUBO 결합 (Q = Σ R_i + α × P). α 작을수록 어려움. `coeff_type`: `lin2` ({-1,+1}, 더 어려움) vs `lin20`.

### Signed Posiform
기존 posiform에 음수 weight 허용. Phase 1 (양수) → Phase 2 (gap 기반 음수). Q 행렬 양/음 비율 ~50%로 자연스러운 QUBO 생성. 생성 한계: n ≤ 25 (gap 계산에 O(2^n)).

### Truth Table QUBO
진리표(비트스트링 → 에너지)에서 Möbius 변환 → Rosenberg 차수축소 → QUBO. 에너지 landscape를 직접 설계 가능. 차수축소 전략: `greedy`(기본, aux 95.8% 절감), `cache`, `original`. 근사 모드 (`create_qubo_approx_optimized`): 보조변수 0개, n≤23 실용적.

### Truth Table Concat
k-bit Truth Table QUBO를 h개 block-diagonal 접합. `planted`/`random` landscape + 선택적 posiform hardening. 큰 QUBO (N=k×h)를 블록별로 생성하여 GS 보장.

### McEliece Cryptographic QUBO
Mandrà et al. (arXiv:2308.09704). McEliece 공개키 → Ising 스핀 → Rosenberg 차수축소 → QUBO. 암호학적 보안에 기반한 hardness. m=GF(2^m) 차수, t=에러 정정 능력. m≥5는 생성이 매우 느림.

## 출력 형식

Q 행렬은 Python dict `{(i, j): weight}` (upper triangular, i ≤ j)로 저장. Edge-list CSV 파일: `# target,<bitstring>` 헤더 + `i,j,weight` 행.

## 의존성

```bash
pip install numpy neal dimod matplotlib
```

- `numpy`: 행렬 연산, 계수 분석
- `neal`: D-Wave Simulated Annealing Sampler
- `dimod`: D-Wave 에코시스템
- `matplotlib`: 실험 결과 시각화
- `qiskit`, `qiskit-algorithms`: (선택) QAOA 양자 회로 시뮬레이션

## 문서

| 문서 | 내용 |
|------|------|
| [방법론 비교](docs/METHODOLOGY_COMPARISON.md) | 전체 방법론 SA 벤치마크 종합 비교 |
| [Posiform 실험 보고서](docs/POSIFORM_EXPERIMENT.md) | N=1000까지 100% 성공, SA-easy 분석 |
| [Quiet Planting 실험 보고서](docs/QUIET_PLANTING_EXPERIMENT.md) | 축퇴 문제, planted field, SA 상전이 |
| [Signal-Hardness Dilemma](docs/SIGNAL_HARDNESS_DILEMMA.md) | 은닉성과 난이도의 근본적 딜레마 분석 |
| [Signed Posiform 실험](signed_posiform/SIGNED_POSIFORM_EXPERIMENT.md) | 음수 weight 비율에 따른 난이도 변화 |
| [참고문헌](docs/REFERENCES.md) | 이론적 배경 논문 목록 |
