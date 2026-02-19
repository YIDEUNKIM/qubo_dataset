# QUBO Benchmark Dataset Generator

양자 어닐링 및 QUBO 솔버 벤치마크를 위한 **정답이 알려진 QUBO 문제 생성기**.

핵심 기능: 임의의 목표 비트스트링이 ground state임이 수학적으로 보장된 QUBO를 생성하여, 솔버의 정확도를 정량적으로 측정할 수 있다.

---

## 프로젝트 구조

```
qubo_dataset/
├── qubo_wishart.py            # Wishart Planted Ensemble QUBO 생성기 (SA-hard)
├── qubo_zero_expectation.py   # 기댓값 0 QUBO 생성기 (구별 불가능) + 공용 유틸리티
├── qubo_hard_mode.py          # Backbone + Frustration 모델 (비교용)
├── test_wishart.py            # SA 실험 프레임워크
├── qubo_results/              # 생성된 QUBO 파일
├── docs/
│   ├── WISHART_EXPERIMENT.md  # 실험 결과 보고서
│   └── DESIGN_RATIONALE.md   # 설계 근거 및 개념 정리
├── 기댓값0만들기_2트.pdf        # 참고 논문
└── 읽어봐줘!.pdf               # 참고 논문
```

## 생성 방식 비교

| 생성기 | 설계 공간 | Ground State 보장 | SA-hard | 구별 불가능 |
|--------|----------|------------------|---------|-----------|
| `qubo_wishart.py` | Ising → QUBO | 수학적 보장 ($W^T t = 0$) | O (alpha < 1.0) | X |
| `qubo_zero_expectation.py` | QUBO 직접 | LP 최적화 페널티 | X | O ($E[Q_{ij}] = 0$) |
| `qubo_hard_mode.py` | Ising → QUBO | 조건부 (W_STRONG >> W_WEAK) | X (SA 100% 성공) | X |

## Quick Start

### Wishart Planted Ensemble (SA-hard QUBO 생성)

```bash
# 기본 사용: target=10110, alpha=0.7
python3 qubo_wishart.py 10110 0.7

# 랜덤 target (길이 100), alpha=0.7, seed=42
python3 qubo_wishart.py 100 0.7 42
```

### SA 실험

```bash
# Alpha sweep (N=100, 10 runs per alpha)
python3 test_wishart.py 100 10

# N 스케일링 (alpha=0.7 고정)
python3 test_wishart.py --scaling 0.7

# Wishart vs Hard Mode 비교
python3 test_wishart.py --compare

# Hardness metrics (TTS, 에너지비)
python3 test_wishart.py --metrics
```

## 주요 실험 결과

### Alpha Sweep (N=100, SA num_reads=200)

| Alpha | SA 성공률 | 에너지 정확도 | 구간 |
|-------|----------|-------------|------|
| 0.50  | 0%       | 94.6%       | Hard |
| 0.70  | 0%       | 92.7%       | **최고 난이도** |
| 0.85  | 50%      | 95.3%       | 전이 구간 |
| 0.95  | 100%     | 100%        | Easy |

### Scaling (alpha=0.7, SA num_reads=200)

| N   | SA 성공률 | 에너지 정확도 |
|-----|----------|-------------|
| 50  | 100%     | 100%        |
| 100 | 0%       | 92.0%       |
| 300 | 0%       | 91.6%       |
| 500 | 0%       | 91.3%       |

- **N=50→100 사이에서 상전이**: N=100부터 SA 완전 실패
- **에너지 정확도 ~91%로 포화**: SA가 metastable 상태에 갇혀 ground state 대비 ~9% gap 유지
- **해밍거리 ≈ N/2**: SA가 찾은 해는 정답과 무관 (랜덤 수준)

### Wishart vs Hard Mode (N=100)

| 방식 | SA 성공률 |
|------|----------|
| Wishart (alpha=0.7) | **10%** |
| Hard Mode (noise=0.1) | 100% |

## 문서

- **[실험 결과 보고서](docs/WISHART_EXPERIMENT.md)** — 알고리즘 상세, 전체 실험 데이터, 분석
- **[설계 근거](docs/DESIGN_RATIONALE.md)** — Ising vs QUBO 설계, 구별 가능성 vs 계산 난이도, QPU 전망

## 의존성

```bash
pip install numpy neal dimod
```

- `numpy`: 행렬 연산
- `neal`: D-Wave Simulated Annealing Sampler
- `dimod`: D-Wave 에코시스템
