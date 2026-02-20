# QUBO Benchmark Dataset Generator

양자 어닐링 및 QUBO 솔버 벤치마크를 위한 **정답이 알려진 QUBO 문제 생성기**.

핵심 기능: 임의의 목표 비트스트링이 ground state임이 보장된 QUBO를 생성하여, 솔버의 정확도를 정량적으로 측정할 수 있다.

---

## 방법론 비교

### 핵심 특성 비교

| | Wishart | McEliece | Quiet Planting | Posiform | Posiform Hardened | Zero Expectation | Hard Mode |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **Ground State 보장** | 수학적 | 조건부 | 조건부 | 수학적 (유일) | 수학적 (유일) | 수학적 | 조건부 |
| **SA 난이도** | Hard | 미측정 | 중간 | Easy | 조절 가능 | Easy | Easy |
| **보조변수** | 없음 | 있음 (대량) | 있음 | 없음 | 없음 | 없음 | 없음 |
| **QUBO 크기** | n | k + O(kN) | n(1+alpha) | n | n | n | n |
| **난이도 조절** | alpha | m, t | alpha, field | coeff_range | alpha, coeff_type | density | noise_ratio |
| **랜덤 QUBO와 구별** | 가능 (low-rank) | 미분석 | 불가 (alpha<3.86) | 미분석 | 가능 (block-diagonal) | 불가 (E[q]=0) | 가능 (backbone) |

### 상세 비교

| 특성 | Wishart | McEliece | Quiet Planting | Posiform | Posiform Hardened | Zero Expectation | Hard Mode |
|------|---------|----------|----------------|----------|-------------------|------------------|-----------|
| **원리** | 직교 Gaussian 투영 W^T t = 0 | Goppa 코드 + McEliece 암호 프로토콜 | Planted 3-SAT + Rosenberg 축소 | Planted 2-SAT + Posiform 변환 | Random QUBO + Posiform overlay | LP 최적화 penalty ratio | Star backbone + frustration |
| **GS 보장 근거** | W^T t = 0 → t가 Ising ground state | p-local: H(target)=t 증명됨. 2-local: M 의존 | SAT 해 = 에너지 0. field로 축퇴 해소 | Posiform 양수 계수 → 유일 최솟값 | ΣR_i(x*) 최적 + P(x*)=0 → 유일 | 비최적 상태에 양의 penalty | Backbone이 target 강제 |
| **GS 보장 한계** | 유한 정밀도 수치 오차 | Rosenberg 페널티 M 부족 시 실패 (w>=10) | field=0이면 모든 SAT해가 동일 에너지 | 없음 (수학적 유일) | 없음 (수학적 유일) | diagonal bias 존재 | noise_ratio 높으면 깨짐 |
| **난이도 근거** | 상전이 (alpha_c ~ 0.95) | 암호학적 보안 (O(2^k) PT) | SAT 상전이 (alpha ~ 4.27) | 없음 (SA-trivial) | alpha 작을수록 rugged landscape | 없음 | 없음 |
| **핵심 파라미터** | alpha = M/N (0.5~0.8 = hard) | m (GF 크기), t (에러 가중치) | alpha (절 밀도), field (편향 강도) | coeff_range (계수 범위) | posiform_scale, coeff_type | density (간선 밀도) | noise_ratio (좌절 비율) |
| **QUBO 변수 수** | n | k + 보조변수 (m=4: 71개) | n + ceil(alpha*n) | n | n | n | n |
| **기반 논문** | Hamze et al. 2020 | Mandra et al. 2025 | Krzakala & Zdeborova 2009 | Hahn, Pelofske & Djidjev 2023 | Pelofske, Hahn & Djidjev 2024 | 자체 연구 | - |

### SA 성공률 비교

| N | Wishart (a=0.7) | Quiet (field=0.5) | Posiform | Hardened (a=0.01) | Zero Expectation | Hard Mode |
|---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 100 | ~10% | 100% | **100%** | 미측정 | ~100% | ~100% |
| 500 | ~0% | 20% | **100%** | 미측정 | ~100% | ~100% |
| 1000 | ~0% | 0% | **100%** | 미측정 | ~100% | ~100% |

McEliece, Posiform Hardened은 SA 스케일링 실험 미수행.

### 벤치마크 적합성 평가

| 평가 기준 | Wishart | McEliece | Quiet Planting | Posiform | Posiform Hardened | Zero Expectation | Hard Mode |
|-----------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| GS 보장 신뢰도 | A | B+ | B | **A+** | **A+** | A | C |
| 솔버 난이도 | **A+** | ? | B+ | F | B (조절 가능) | F | F |
| 보조변수 없음 | O | **X** | X | O | O | O | O |
| 확장성 (N>1000) | O | X (크기 폭발) | O | O | O (subgraph brute force 제한) | O | O |
| **종합** | **1순위** | 연구 가치 | 3순위 | GS 검증용 | **2순위 후보** | 참고용 | 참고용 |

**권장 사용 시나리오**:
- **솔버 성능 비교**: Wishart (유일하게 SA-hard + 수학적 GS 보장)
- **난이도 조절 가능한 벤치마크**: Posiform Hardened (수학적 유일 GS + alpha로 난이도 제어)
- **정확도 검증**: Posiform (100% 유일 GS, SA가 반드시 찾음)
- **암호학적 난이도 연구**: McEliece (p-local Ising에서 O(2^k))
- **통계적 은닉성 연구**: Quiet Planting (alpha < 3.86)

---

## 프로젝트 구조

```
qubo_dataset/
├── qubo_utils.py                 # 공유 유틸리티 (calculate_energy, save_qubo_edgelist 등)
├── docs/                         # 실험 보고서
│
├── wishart/                      # Wishart Planted Ensemble (SA-hard)
│   ├── qubo_wishart.py
│   ├── test_wishart.py
│   └── papers/
│
├── mceliece/                     # McEliece Cryptographic QUBO (암호학적 난이도)
│   ├── qubo_mceliece.py
│   ├── README.md                 # 상세 문서 (논문 대조, 페널티 M 분석)
│   └── papers/
│
├── quiet_planting/               # Quiet Planting (3-SAT → Rosenberg)
│   ├── qubo_quiet_planted.py
│   ├── test_quiet_planted.py
│   └── papers/
│
├── posiform/                     # Posiform Planting (2-SAT → Posiform, SA-easy)
│   ├── qubo_posiform.py
│   ├── test_posiform.py
│   └── papers/
│
├── posiform_hardened/            # Hardened Posiform (Random QUBO + Posiform overlay)
│   ├── qubo_posiform_hardened.py
│   ├── test_posiform_hardened.py
│   ├── run_experiments.py
│   ├── README.md                 # 상세 문서 (방법론, GS 증명, 한계점)
│   ├── papers/
│   └── results/
│
├── zero_expectation/             # Zero Expectation (E[q_ij]=0 보장)
│   ├── qubo_zero_expectation.py
│   └── test_zero_expectation.py
│
└── hard_mode/                    # Backbone + Frustration (기준선)
    └── qubo_hard_mode.py
```

## Quick Start

```bash
# Wishart: SA-hard QUBO (alpha=0.7이 가장 어려움)
python3 wishart/qubo_wishart.py 10110 0.7

# McEliece: 암호학적 QUBO (m=3, t=1)
python3 mceliece/qubo_mceliece.py 1011 3 1 42

# Quiet Planting: 통계적 은닉성 (alpha<3.86)
python3 quiet_planting/qubo_quiet_planted.py 10110 4.2

# Posiform: 유일한 ground state 보장 (SA-easy)
python3 posiform/qubo_posiform.py 10110

# Posiform Hardened: 난이도 조절 가능 (alpha=0.01이 가장 어려움)
python3 posiform_hardened/qubo_posiform_hardened.py 50 lin2 0.01 42

# Zero Expectation: E[q_ij]=0 보장
python3 zero_expectation/qubo_zero_expectation.py 10110

# Hard Mode: 기준선
python3 hard_mode/qubo_hard_mode.py 10110
```

## SA 실험

```bash
# Wishart alpha sweep (N=100, 10 runs)
python3 wishart/test_wishart.py 100 10

# Wishart N 스케일링
python3 wishart/test_wishart.py --scaling 0.7

# Quiet Planting N 스케일링
python3 quiet_planting/test_quiet_planted.py --scaling 4.2

# Posiform N 스케일링
python3 posiform/test_posiform.py --scaling 10

# Posiform Hardened: sweep 전이 실험
python3 posiform_hardened/test_posiform_hardened.py --sweep

# Posiform Hardened: N 스케일링
python3 posiform_hardened/test_posiform_hardened.py --scaling

# Posiform Hardened: Hardened vs Plain 비교
python3 posiform_hardened/test_posiform_hardened.py --compare

# 다중 비교 (Posiform vs Quiet vs Wishart vs ZeroExp vs HardMode)
python3 posiform/test_posiform.py --compare
```

## 의존성

```bash
pip install numpy neal dimod
```

- `numpy`: 행렬 연산
- `neal`: D-Wave Simulated Annealing Sampler
- `dimod`: D-Wave 에코시스템

## 각 방법론 상세 문서

| 방법론 | 문서 |
|--------|------|
| McEliece | [mceliece/README.md](mceliece/README.md) — 논문 대조, 페널티 M 분석, 한계점 |
| Posiform 실험 | [docs/POSIFORM_EXPERIMENT.md](docs/POSIFORM_EXPERIMENT.md) |
| Posiform Hardened | [posiform_hardened/README.md](posiform_hardened/README.md) — 방법론, GS 유일성 증명, 난이도 조절 |
| Quiet Planting 실험 | [docs/QUIET_PLANTING_EXPERIMENT.md](docs/QUIET_PLANTING_EXPERIMENT.md) |
| Wishart 실험 | [docs/WISHART_EXPERIMENT.md](docs/WISHART_EXPERIMENT.md) |
