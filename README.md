# QUBO Benchmark Dataset Generator

양자 어닐링 및 QUBO 솔버 벤치마크를 위한 **정답이 알려진 QUBO 문제 생성기**.

핵심 기능: 임의의 목표 비트스트링이 ground state임이 수학적으로 보장된 QUBO를 생성하여, 솔버의 정확도를 정량적으로 측정할 수 있다.

---

## 프로젝트 구조

```
qubo_dataset/
├── qubo_utils.py                 # 공유 유틸리티 (calculate_energy, save_qubo_edgelist 등)
├── docs/                         # 실험 보고서
│   ├── POSIFORM_EXPERIMENT.md
│   ├── QUIET_PLANTING_EXPERIMENT.md
│   └── ...
│
├── zero_expectation/             # Zero Expectation (E[q_ij]=0 보장)
│   ├── qubo_zero_expectation.py
│   ├── test_zero_expectation.py
│   ├── results/
│   └── README.md
│
├── wishart/                      # Wishart Planted Ensemble (SA-hard)
│   ├── qubo_wishart.py
│   ├── test_wishart.py
│   ├── results/
│   ├── papers/
│   └── README.md
│
├── quiet_planting/               # Quiet Planting (3-SAT → Rosenberg)
│   ├── qubo_quiet_planted.py
│   ├── test_quiet_planted.py
│   ├── results/
│   ├── papers/
│   └── README.md
│
├── posiform/                     # Posiform Planting (2-SAT → Posiform)
│   ├── qubo_posiform.py
│   ├── test_posiform.py
│   ├── results/
│   ├── papers/
│   └── README.md
│
└── hardened_posiform/            # Hardened Posiform (Random QUBO + Posiform)
    ├── qubo_posiform_hardened.py
    ├── test_posiform_hardened.py
    ├── results/
    ├── papers/
    └── README.md
```

각 방법론 디렉토리에 `results/`(생성된 QUBO 파일), `papers/`(참조 논문 PDF), `README.md`(상세 문서)가 포함됨.

## 생성 방식 비교

| 생성기 | QUBO 크기 | Ground State | SA 난이도 | 구별 불가능 | 벤치마크 적합 | 핵심 논문 |
|--------|:---------:|:-----------:|:---------:|:----------:|:----------:|----------|
| **Wishart** | n | 수학적 (유한정밀도 제외) | **SA-hard** | X (low-rank) | 조건부 | Hamze et al. 2020 |
| **Hardened Posiform** | n | **수학적 (유일)** | SA-moderate | 미분석 | **O** | Pelofske et al. 2024 |
| **Quiet Planting** | n(1+alpha) | 조건부 (field 필요) | field 의존적 | **O** (alpha<3.86) | 조건부 | Krzakala & Zdeborova 2009 |
| **Posiform** | n | **수학적 (유일)** | SA-easy | 미분석 | **O** | Hahn et al. 2023 |
| **Zero Expectation** | n(n-1)/2 | 수학적 | SA-easy | **O** (E[q_ij]=0) | **O** | 내부 연구 |

> **벤치마크 적합**: 주장하는 ground state가 실제로 맞는지의 신뢰도. "조건부"는 파라미터에 따라 GS가 깨질 수 있음을 의미.

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

# 4-way 비교 (Posiform vs Quiet vs Wishart vs ZeroExp)
python3 posiform/test_posiform.py --compare
```

## 주요 실험 결과 요약

### SA 성공률 비교 (num_reads=1~200)

| N | Posiform | Hardened (lin2,α=0.01) | Quiet (field=0.5) | Wishart (alpha=0.7) | ZeroExp |
|---:|:--------:|:---------------------:|:-----------------:|:-------------------:|:-------:|
| 100 | **100%** | ~100% | 100% | ~10% | ~100% |
| 500 | **100%** | **~13%** (500 sweeps) | 20% | ~0% | ~100% |
| 1000 | **100%** | - | 0% | ~0% | ~100% |

**SA 난이도 분류 기준** (표준 조건: num_sweeps=10N, num_reads=200):
- **SA-hard**: p(N=500) < 10% — Wishart (alpha=0.7), Quiet Planting (N>300)
- **SA-moderate**: 10% ≤ p(N=500) < 90% — Hardened Posiform (lin2, α=0.01)
- **SA-easy**: p(N=1000) ≥ 90% — Posiform, Zero Expectation

자세한 분석: [`docs/POSIFORM_EXPERIMENT.md`](docs/POSIFORM_EXPERIMENT.md), [`docs/QUIET_PLANTING_EXPERIMENT.md`](docs/QUIET_PLANTING_EXPERIMENT.md)

## 문서

| 문서 | 내용 |
|------|------|
| [Posiform 실험 보고서](docs/POSIFORM_EXPERIMENT.md) | N=1000까지 100% 성공, SA-easy 분석 |
| [Quiet Planting 실험 보고서](docs/QUIET_PLANTING_EXPERIMENT.md) | 축퇴 문제, planted field, SA 상전이 |
| [각 방법론 README](posiform/README.md) | 이론, 구현, 파라미터, 참고문헌 |

## 의존성

```bash
pip install numpy neal dimod
```

- `numpy`: 행렬 연산
- `neal`: D-Wave Simulated Annealing Sampler
- `dimod`: D-Wave 에코시스템
