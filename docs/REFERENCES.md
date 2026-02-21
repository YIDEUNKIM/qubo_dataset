# 참고 문헌 (References)

## 1. 핵심 논문

### Wishart Planted Ensemble

- **논문**: "Wishart Planted Ensemble: A Tunably Rugged Pairwise Ising Model with a First-Order Phase Transition"
- **저자**: Firas Hamze, Jack Raymond, Christopher A. Pattison, Katja Biswas, Helmut G. Katzgraber
- **출판**: Physical Review E, 2020
- **arXiv**: [arXiv:2009.11217](https://arxiv.org/abs/2009.11217)
- **사용처**: `wishart/qubo_wishart.py`, `wishart/test_wishart.py`
- **핵심 내용**:
  - 목표 비트스트링 t가 ground state임이 수학적으로 보장되는 Ising 모델 구성법
  - alpha = M/N 단일 파라미터로 난이도 제어
  - alpha_c 근처에서 1차 상전이 발생, SA가 metastable 상태에 갇힘
  - 직교 Gaussian 투영 (W^T t = 0)으로 ground state 보장

### Posiform Planting (2-SAT → QUBO)

- **논문**: "Using 2-SAT to Generate QUBO Instances with Known Optimal Solutions"
- **저자**: Georg Hahn, Elijah Pelofske, Hristo N. Djidjev
- **출판**: IEEE International Conference on Quantum Computing and Engineering (QCE), 2023
- **사용처**: `posiform/qubo_posiform.py`, `posiform/test_posiform.py`
- **핵심 내용**:
  - Planted 2-SAT → posiform → QUBO 파이프라인
  - 보조변수 없이 QUBO 크기 = N 유지
  - Tarjan SCC 기반 2-SAT 유일성 검사로 GS 유일성 수학적 보장
  - SA에 대해 trivially easy (2-SAT이 P에 속하므로 smooth landscape)

### Hardened Posiform Planting

- **논문**: "Increasing the Hardness of Posiform Planting Using Random QUBOs for Programmable Quantum Annealer Benchmarking"
- **저자**: Elijah Pelofske, Georg Hahn, Hristo N. Djidjev
- **출판**: npj Unconventional Computing, 2025
- **사용처**: `hardened_posiform/qubo_posiform_hardened.py`, `hardened_posiform/test_posiform_hardened.py`
- **핵심 내용**:
  - Random discrete-coefficient QUBO + posiform QUBO 결합 (Q = ΣR_i + α×P)
  - 기본 posiform의 SA-trivial 한계를 극복
  - GS 유일성 수학적 보장 유지
  - alpha_scale(α)과 coeff_type으로 난이도 연속 조절

### McEliece Cryptographic QUBO

- **논문**: "Ising formulation of the McEliece problem"
- **저자**: Salvatore Mandrà, Gianluca Passarelli, Gian Giacomo Guerreschi
- **출판**: Future Generation Computer Systems (FGCS), 2025
- **arXiv**: [arXiv:2308.09704](https://arxiv.org/abs/2308.09704)
- **사용처**: `mceliece/qubo_mceliece.py`, `mceliece/test_mceliece.py`
- **핵심 내용**:
  - McEliece 암호 프로토콜의 공개키를 Ising 스핀 시스템으로 캐스팅
  - GF(2^m) 위의 Goppa 코드 → McEliece 키 생성 → p-local Ising 상호작용
  - Rosenberg 보조변수로 p-body → 2-body 차수 축소 → QUBO
  - 암호학적 보안에 기반한 computational hardness
  - 파라미터 (m, t)로 난이도 조절: m=GF(2^m) 확장 차수, t=에러 정정 능력
  - **제한**: m≥5에서 Rosenberg 차수축소의 지수적 비용으로 QUBO 생성 불가. Eq.14 exact decomposition으로 해결 가능 (미구현)

### Quiet Planting (Planted 3-SAT)

- **논문**: "Hiding Quiet Solutions in Random Constraint Satisfaction Problems"
- **저자**: Florent Krzakala, Lenka Zdeborova
- **출판**: Physical Review Letters, 102(23), 238701, 2009
- **사용처**: `quiet_planting/qubo_quiet_planted.py`, `quiet_planting/test_quiet_planted.py`
- **핵심 내용**:
  - Planted random 3-SAT에서 alpha < 3.86이면 랜덤과 통계적 구별 불가
  - Rosenberg reduction으로 3-SAT → QUBO (보조변수 m개)
  - 축퇴 해소를 위한 planted field 필요

---

## 2. 프로젝트 내부 참고 자료

### 기댓값 0 QUBO 설계

- **파일**: `기댓값0만들기_2트.pdf`
- **내용**: Zero-expectation QUBO 구성 방법론. Q 행렬의 계수가 E[q_ii] = 0, E[q_ij] = 0을 만족하도록 LP 최적화된 페널티 비율을 설계
- **사용처**: `qubo_zero_expectation.py`

### QUBO 벤치마크 개요

- **파일**: `읽어봐줘!.pdf`
- **내용**: 프로젝트 배경 및 QUBO 벤치마크 데이터셋 생성의 동기. 정답이 알려진 QUBO를 생성하여 솔버 성능을 정량적으로 평가하는 방법론
- **사용처**: 프로젝트 전반

---

## 3. 이론적 배경

### QUBO (Quadratic Unconstrained Binary Optimization)

- **정의**: $E(x) = x^T Q x = \sum Q_{ii}x_i + \sum_{i<j} Q_{ij}x_ix_j$, $x_i \in \{0, 1\}$
- **특성**: NP-hard (일반적인 경우)
- **참고 문헌**:
  - Kochenberger, G. A., et al. "The unconstrained binary quadratic programming problem: a survey." *Journal of Combinatorial Optimization*, 28(1), 58-81, 2014
  - Glover, F., Kochenberger, G., & Du, Y. "Quantum Bridge Analytics I: a tutorial on formulating and using QUBO models." *4OR*, 17(4), 335-371, 2019

### Ising Model

- **정의**: $H(s) = -\sum_{i<j} J_{ij} s_i s_j - \sum_i h_i s_i$, $s_i \in \{-1, +1\}$
- **QUBO-Ising 등가성**: $s_i = 2x_i - 1$ 치환으로 상호 변환 가능
- **사용처**: `qubo_wishart.py`, `qubo_hard_mode.py`, `qubo_zero_expectation.py` (Ising-derived 모드)

### Simulated Annealing (SA)

- **개요**: Metropolis-Hastings 기반 확률적 최적화. 온도를 점진적으로 낮추며 에너지 최소 상태 탐색
- **구현**: D-Wave `neal.SimulatedAnnealingSampler`
- **참고 문헌**:
  - Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. "Optimization by simulated annealing." *Science*, 220(4598), 671-680, 1983

### QAOA (Quantum Approximate Optimization Algorithm)

- **개요**: 변분 양자 알고리즘. 양자 중첩과 변분 최적화를 결합하여 조합 최적화 문제를 근사적으로 해결
- **구현**: Qiskit (`qubo_to_qaoa.py`, `run_qaoa_n5.py`)
- **참고 문헌**:
  - Farhi, E., Goldstone, J., & Gutmann, S. "A Quantum Approximate Optimization Algorithm." [arXiv:1411.4028](https://arxiv.org/abs/1411.4028), 2014

### Quantum Annealing

- **개요**: 양자 터널링 효과를 이용하여 에너지 장벽을 관통. SA가 열적 요동으로 넘어야 하는 장벽을 양자역학적으로 통과
- **플랫폼**: D-Wave QPU
- **프로젝트 관련**: Wishart 벤치마크로 "QPU가 SA 대비 얼마나 더 잘 푸는가"를 정량 측정하는 것이 궁극적 목표
- **참고 문헌**:
  - Kadowaki, T. & Nishimori, H. "Quantum annealing in the transverse Ising model." *Physical Review E*, 58(5), 5355, 1998
  - Johnson, M. W., et al. "Quantum annealing with manufactured spins." *Nature*, 473(7346), 194-198, 2011

---

## 4. 소프트웨어 & 라이브러리

| 라이브러리 | 용도 | 참고 |
|-----------|------|------|
| `numpy` | 행렬 연산 (Wishart 구성, 에너지 계산) | [numpy.org](https://numpy.org) |
| `neal` | D-Wave Simulated Annealing Sampler | [docs.ocean.dwavesys.com](https://docs.ocean.dwavesys.com/en/stable/docs_neal/) |
| `dimod` | D-Wave 에코시스템 (BQM, 샘플러 인터페이스) | [docs.ocean.dwavesys.com](https://docs.ocean.dwavesys.com/en/stable/docs_dimod/) |
| `qiskit` | QAOA 양자 회로 시뮬레이션 | [qiskit.org](https://qiskit.org) |
| `matplotlib` | Q 행렬 시각화 | [matplotlib.org](https://matplotlib.org) |

---

## 5. 프로젝트 내 개념-파일 매핑

| 개념 | 관련 파일 | 참고 논문/이론 |
|------|----------|--------------|
| Wishart Planted Ensemble | `wishart/qubo_wishart.py` | Hamze et al. (2020) |
| Zero-Expectation QUBO | `zero_expectation/qubo_zero_expectation.py` | `기댓값0만들기_2트.pdf` |
| Quiet Planting | `quiet_planting/qubo_quiet_planted.py` | Krzakala & Zdeborova (2009) |
| Posiform Planting | `posiform/qubo_posiform.py` | Hahn, Pelofske, Djidjev (2023) |
| Hardened Posiform | `hardened_posiform/qubo_posiform_hardened.py` | Pelofske, Hahn, Djidjev (2025) |
| McEliece QUBO | `mceliece/qubo_mceliece.py`, `mceliece/test_mceliece.py` | Mandrà et al. (2025) |
| Ising-QUBO 변환 | `wishart/qubo_wishart.py` | Ising model 표준 변환 |
| SA 벤치마크 | 각 방법론별 `test_*.py` | Kirkpatrick et al. (1983) |
| **전체 비교** | **`docs/METHODOLOGY_COMPARISON.md`** | **6-way SA 벤치마크 결과** |
