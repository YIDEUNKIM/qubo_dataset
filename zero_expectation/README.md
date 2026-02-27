# Zero-Expectation QUBO 생성기

## 개요

**랜덤 QUBO와 통계적으로 구별 불가능한** planted QUBO를 생성하는 방법론. Q 행렬의 비대각 계수가 **E[q_ij] = 0** 을 만족하도록, 각 큐빗 쌍의 3가지 오답 상태에 부여하는 페널티 비율을 수학적으로 도출한다.

기본 모델(**ZeroOffDiagonalModel**): double-flip ratio = single-flip ratio 합 → 비율 {1, 1, 2}. 모든 target pair에서 off-diagonal 기댓값 0이 수학적으로 보장됨.

## 이론적 배경

### QUBO와 Planted 문제

QUBO: E(x) = x^T Q x, x_i in {0,1}을 최소화하는 문제.

**Planted QUBO**: 목표 비트스트링 x\*가 주어졌을 때, x\*가 ground state인 Q를 구성하되, Q의 계수 분포가 랜덤 QUBO와 구별 불가능하도록 만드는 것이 핵심 질문이다.

### 왜 구별 불가능해야 하는가?

#### 1. 블라인드 벤치마크에서의 부정행위 방지

QUBO 벤치마크의 목적은 **솔버의 순수한 최적화 능력**을 측정하는 것이다. 그런데 Q 행렬에서 target 정보가 누출되면, 솔버 개발자가 최적화 없이 통계 분석만으로 정답을 추출할 수 있다:

```
1. Q 행렬의 대각 항 부호를 분석
2. q_ii가 양수면 target bit = 0, 음수면 target bit = 1로 추정
3. 최적화를 수행하지 않고도 높은 정확도로 정답 추측 가능
4. 벤치마크 무의미
```

이를 방지하려면, Q 행렬의 계수 분포가 target에 대한 정보를 누출하지 않아야 한다.

#### 2. 암호학적 비유: Semantic Security

이 요구사항은 암호학의 **semantic security**와 정확히 동일한 구조이다:

| 암호학 | QUBO 벤치마크 |
|--------|--------------|
| ciphertext가 plaintext 정보를 누출하면 안 됨 | Q 행렬이 target 정보를 누출하면 안 됨 |
| 랜덤 문자열과 구별 불가능한 ciphertext | 랜덤 QUBO와 구별 불가능한 planted QUBO |
| 구별 가능 → 복호화 가능성 | 구별 가능 → target 추출 가능성 |

#### 3. 구체적 누출 경로와 Zero-Expectation의 대응

Q 행렬에서 target 정보가 누출되는 두 가지 경로:

| 누출 경로 | 공격 방법 | Zero-Expectation의 대응 |
|----------|----------|----------------------|
| **대각 편향**: E[q_ii\|b=0] ≠ E[q_ii\|b=1] | 대각 항의 부호/크기로 target bit 방향 추정 | LP로 비대각 무편향 달성. 대각 편향은 잔존 (동시 달성 불가) |
| **비대각 편향**: E[q_ij] ≠ 0 | 비대각 항의 부호 패턴에서 target pair 추론 | **E[q_ij] = 0** 보장 → 비대각 항에서 정보 누출 차단 |

다른 방법론이 구별 가능한 이유:

| 방법론 | 누출 경로 | 공격 방법 |
|--------|----------|----------|
| Wishart | Q가 low-rank 구조 (rank M) | 고유값 분해로 rank-M 탐지 + 대각 편향으로 target 추측 |
| Hardened Posiform | Q가 block-diagonal 구조 | 구조 분석으로 planted 여부 탐지 가능 |

#### 4. 주의: "구별 불가능 ≠ SA-hard"

두 성질은 **독립적**이다. Zero-Expectation은 구별은 불가능하지만 SA-trivial이고, Wishart는 구별 가능하지만 SA-hard이다.

- **구별 불가능** = Q 행렬 분석으로 target 정보를 추출할 수 없음 (정보 보안)
- **SA-hard** = 에너지 지형에 metastable trap이 있어 SA가 ground state를 찾지 못함 (계산 난이도)

SA-hard하려면 Q 계수 간 상관 구조가 필요하고, 구별 불가능하려면 상관 구조가 없어야 한다. 근본적으로 충돌하는 요구사항이므로, 두 성질의 완전한 결합은 열린 연구 문제이다. 자세한 논의는 [`docs/DESIGN_RATIONALE.md`](../docs/DESIGN_RATIONALE.md) Section 4 참조.

### 페널티 기반 QUBO 구성법

각 큐빗 쌍 (i, j)에서 목표 상태 (b_i, b_j) 외의 3가지 오답 상태에 독립적인 양의 페널티 r > 0를 부여한다. 오답 상태 (s_i, s_j)에 대한 페널티 r은 Q 행렬에 다음과 같이 분배된다:

| 오답 상태 (s_i, s_j) | 페널티 함수 | Q_ii 기여 | Q_jj 기여 | Q_ij 기여 |
|:-------------------:|:----------:|:---------:|:---------:|:---------:|
| (0, 0) | r(1-x_i)(1-x_j) | -r | -r | +r |
| (0, 1) | r(1-x_i)x_j | 0 | +r | -r |
| (1, 0) | rx_i(1-x_j) | +r | 0 | -r |
| (1, 1) | rx_ix_j | 0 | 0 | +r |

목표 상태 (b_i, b_j)에서는 모든 페널티 함수가 0이므로, target에서의 에너지 기여 = 0. 오답 상태에서는 해당 페널티 > 0이므로 에너지 증가 → **target이 ground state**.

### Off-Diagonal 기댓값 0 조건의 도출

Q 행렬의 비대각 항 q_ij는 3개 페널티의 q_ij 기여 합이다. 위 표에서 부호 패턴:

```
(0,0) → +r,  (0,1) → -r,  (1,0) → -r,  (1,1) → +r
```

E[q_ij] = 0이 되려면, 3개 오답의 부호 가중합이 0이어야 한다. 각 target pair에 대해:

```
target (0,0): -R(01) - R(10) + R(11) = 0  →  R(11) = R(01) + R(10)
target (0,1): +R(00) - R(10) + R(11) = 0  →  R(10) = R(00) + R(11)
target (1,0): +R(00) - R(01) + R(11) = 0  →  R(01) = R(00) + R(11)
target (1,1): +R(00) - R(01) - R(10) = 0  →  R(00) = R(01) + R(10)
```

**공통 규칙: double-flip ratio = single-flip ratio 합.** 대칭 single-flip (R_a = R_b = 1)이면 double-flip = 2. 비율 **{1, 1, 2}**.

### 대각 편향과 비대각 편향의 동시 달성 불가능성

#### 대수적 불가능성 증명

`test_diagonal_zero.py`에서 수학적으로 증명:

- **E[Q_ij] = 0** (비대각 무편향)과 **E[Q_ii|b=0] = E[Q_ii|b=1]** (대각 무편향)을 **동시에** 달성하는 것은, 양의 페널티 제약 하에서 **불가능**하다.
- 증명: 동시 달성 조건을 전개하면 r_2 + r_6 + r_7 + r_11 = 0이 필요하나, 모두 양수이므로 모순.

따라서 **off-diagonal 무편향을 우선**하고 diagonal 편향은 허용하는 것이 합리적 선택이다 (diagonal 탐지 SNR ~ √n < off-diagonal 탐지 SNR ~ n).

#### 심층 분석: 대각 편향 = 에너지 갭 (analyze_diagonal_deep.py)

대수적 증명을 넘어, 대각 편향이 **왜 제거 불가능한지**에 대한 물리적/정보이론적 분석:

**1. 대각 편향 ≡ 에너지 갭 (동일 메커니즘)**

변수 k의 single-bit flip 에너지 변화:

```
ΔE_k = Q_kk + Σ_{m≠k} Q_km · b_m
```

E[Q_km] = 0 (ZeroOffDiagonal 보장)이므로:

```
E[ΔE_k] = E[Q_kk]   (교차항 기댓값 소멸)
```

실험 검증 (N=20, 1000회 반복): `|E[ΔE_k] - E[Q_kk]|` < 5×10⁻¹³ (기계 정밀도).

**핵심**: 대각 편향을 제거하면 E[ΔE_k] = 0이 되어 single-flip 에너지 장벽이 사라지고, ground state가 파괴된다. 실제로 대각 항을 0으로 설정하면 **100/100 (100%)** 인스턴스에서 GS가 깨짐.

정확한 스케일링:

| b_k 값 | E[Q_kk] | 페어당 기여 |
|:------:|:-------:|:---------:|
| 0 | +(n-1) × μ | +μ |
| 1 | -(n-1) × μ | -μ |

여기서 μ = E[r] = 2.0 (base_range=(1,3)). 대각 편향의 크기 ~ 2n은 에너지 갭의 크기 ~ 2n과 **정확히 일치**.

**2. 노이즈 주입의 비효과성**

대각에 N(0, σ²) 노이즈를 추가하여 편향을 은닉하려는 시도:

| σ | GS 보존율 (N=50) | 대각 SNR | SNR 감소율 |
|--:|:---------------:|:--------:|:---------:|
| 0 | 100% | 0.9975 | — |
| 10 | 100% | 0.9876 | 1.0% |
| 20 | 98% | 0.9604 | 3.7% |
| **31.5** (최적) | **95%** | **0.9500** | **4.8%** |
| 50 | 36% | 0.8491 | 14.9% |

**최적 σ ≈ 31.5에서도 SNR은 0.95 (5%만 감소)**. 편향 크기(~2n)가 에너지 갭 크기(~2n)와 같으므로, 갭을 유지하면서 편향만 숨기는 것이 원리적으로 불가능하다.

**3. Ising-Derived와의 비교: 반대 편향 패턴**

| 모드 | 대각 SNR | 비대각 SNR | 비대각 부호 탐지 |
|------|:--------:|:---------:|:--------------:|
| **ZeroOffDiagonal** | 0.999 | **0.003** | 불가 |
| **Ising-derived** | **0.105** | 0.962 | **100%** (742,500/742,500) |

Ising-derived 모드는 E[Q_kk|b_k] = 0 (대각 무편향)을 달성하지만, 대신 sign(Q_ij) = -sign(s_i·s_j)로 비대각 부호가 target을 100% 노출한다. **두 모드의 편향이 정확히 상보적** — 한쪽을 0으로 만들면 다른 쪽이 완전히 노출되는 구조.

**4. 하이브리드 접근: 한계적 개선**

Q_hybrid = (1-λ) × Q_zero + λ × Q_ising으로 두 모드를 혼합:

| λ | 대각 SNR | 비대각 SNR | max(대각, 비대각) |
|--:|:--------:|:---------:|:----------------:|
| 0.0 | 0.999 | 0.003 | 0.999 |
| 0.3 | 0.907 | 0.420 | 0.907 |
| **0.5** | **0.830** | **0.830** | **0.830** (minimax) |
| 0.7 | 0.420 | 0.907 | 0.907 |
| 1.0 | 0.105 | 0.962 | 0.962 |

Minimax 최적점 λ=0.5에서 max_SNR = 0.83. 순수 ZeroOffDiagonal(0.999) 대비 17% 개선이지만 여전히 **높은 탐지율**. GS는 모든 λ에서 100% 보존.

#### 정보이론적 결론

```
대각 편향 제거 불가능의 근본 원인:

1. Off-diagonal E[q_ij] = 0이 완벽하면
   → GS를 보장하는 모든 정보가 대각에 집중
   → 대각 편향 = 에너지 갭 (동일 물리량)

2. 대각 편향 크기 ~ O(n) = 에너지 갭 크기 ~ O(n)
   → 편향을 O(n)만큼 줄이면 갭도 O(n)만큼 줄어듦
   → "편향만 숨기고 갭은 유지"가 원리적 불가능

3. Ising-derived로 대각을 고치면 비대각이 100% 노출
   → 두 채널이 정보를 보완적으로 분담
   → 총 정보량은 보존 (채널 간 재분배만 가능)
```

> **한 줄 요약**: 양의 페널티 제약 하에서 ground state를 보장하려면 반드시 탐지 가능한 편향이 존재해야 하며, off-diagonal과 diagonal 중 하나를 선택할 수밖에 없다. ZeroOffDiagonalModel은 더 위험한 off-diagonal 누출(SNR~n)을 완벽 차단하고, 덜 위험한 diagonal 누출(SNR~√n)만 허용하는 **최적 전략**이다.

## 구현 방식

### Strategy Pattern

`PenaltyModel` 추상 기본 클래스(ABC)를 통해 다양한 페널티 전략을 플러그인 방식으로 교체 가능:

```python
class PenaltyModel(ABC):
    @abstractmethod
    def get_ratios(self, target_pair: tuple) -> dict:
        """penalty_state -> ratio 매핑 반환"""
        pass
```

### 제공 모델

#### 1. ZeroOffDiagonalModel (기본)

모든 target pair에서 **E[q_ij] = 0** 보장. double-flip ratio = single-flip 합.

```python
(0, 0): {(0, 1): 1.0, (1, 0): 1.0, (1, 1): 2.0},  # 2.0 = 1.0 + 1.0
(0, 1): {(0, 0): 1.0, (1, 0): 2.0, (1, 1): 1.0},
(1, 0): {(0, 0): 1.0, (0, 1): 2.0, (1, 1): 1.0},
(1, 1): {(0, 0): 2.0, (0, 1): 1.0, (1, 0): 1.0},
```

- **어떤 기댓값이 0인가**: off-diagonal E[q_ij] = 0 ✓ / diagonal E[q_ii] = 0 ✗
- **장점**: |q_ij|의 multiset이 항상 {1, 1, 2}로 target pair와 **완전히 독립** (모든 모멘트 동일). Off-diagonal로는 구별 불가.
- **단점**: diagonal E[q_kk] ≠ 0 (부호가 b_k를 노출). 탐지 SNR ~ √n.
- **Minimax 최적성**: off-diagonal zero 제약 하에서 비대칭 single-flip(a ≠ b)을 허용하면 diagonal bias가 변수 위치(i)에 의존하여 더 나빠짐. 대칭 {1, 1, 2}가 유일한 minimax 최적해.

#### 2. SimpleUniformModel (기준선)

모든 오답에 동일 페널티 (ratio = 1.0). 편향 최적화 없음. 비교 실험용.

#### 3. DefaultZeroExpectationModel [DEPRECATED]

LP로 최적화된 비율 (1.64, 1.68 등).

```python
(0, 0): {(0, 1): 1.00, (1, 0): 1.00, (1, 1): 1.64},
(0, 1): {(0, 0): 2.00, (1, 0): 1.00, (1, 1): 1.68},
...
```

**Deprecated 사유**: LP 비율은 "3개 중 1개 선택" 카운팅에서 도출되었으나, 현재 코드는 3개 penalty를 모두 추가하므로 전제가 불일치. off-diagonal 기댓값이 0이 아님:

| target pair | E[q_ij] (이론) | E[q_ij] (실험, μ=2) |
|:-----------:|:--------------:|:-------------------:|
| (0,0) | -0.36μ | -0.73 |
| (0,1) | +2.68μ | +5.36 |
| (1,0) | +2.68μ | +5.35 |
| (1,1) | -5.00μ | -9.99 |

#### 4. BalancedModel [DEPRECATED]

Minimax 균형 모델.

**Deprecated 사유**: target (0,0)만 off-diagonal 0이고 나머지 target pair는 깨짐:

| target pair | E[q_ij] (이론) | E[q_ij] (실험, μ=2) |
|:-----------:|:--------------:|:-------------------:|
| (0,0) | 0 | +0.01 |
| (0,1) | +7μ/6 | +2.33 |
| (1,0) | +7μ/6 | +2.34 |
| (1,1) | -5μ/3 | -3.34 |

off-diagonal zero 제약 하에서 minimax 최적해는 ZeroOffDiagonalModel과 동일 (대칭 single-flip이 유일해). 별도 모델이 불필요.

### Ising-Derived 모드

`create_qubo_ising_derived(target)`: Ising 모델에서 유도된 QUBO.

```
J_ij = alpha * s_i * s_j,  alpha ~ Uniform(1, 3)
s_i = 2*x_i - 1 치환으로 QUBO 변환:
  Q_ij = -4 * J_ij
  Q_ii += 2 * J_ij, Q_jj += 2 * J_ij
```

- **E[Row Sum] = 0** 보장 (Ising 구조에 의해 행 합 기댓값이 대칭)
- Ground state가 target과 일치함이 수학적으로 보장 (J_ij * s_i * s_j > 0)

### QUBO 생성 파이프라인

```
1. 모든 큐빗 쌍 (i, j) 순회  (i < j)
   ├── density 확률로 상호작용 추가 여부 결정
   ├── target pair (b_i, b_j) 확인
   ├── PenaltyModel.get_ratios(target_pair) → 3개 오답 비율
   └── 각 오답 상태에 대해:
       ├── r = Uniform(1, 3) * ratio
       └── 페널티 함수 전개 → Q_ii, Q_jj, Q_ij에 기여 누적

2. Q 행렬 반환 (상삼각 딕셔너리 형태)
```

### 핵심 파라미터

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `density` | 1.0 | 큐빗 쌍 상호작용 추가 확률. 1.0이면 모든 쌍 사용 |
| `base_range` | (1, 3) | 기본 페널티 r의 균등분포 범위 |
| `model` | ZeroOffDiagonalModel | 페널티 비율 전략 |
| `balance_rows` | False | True이면 Ising-derived 모드 사용 |

### 주요 함수

| 함수 | 설명 |
|------|------|
| `create_qubo_precise(target, density, model)` | **메인 진입점**. 밀도와 모델 지정 |
| `create_qubo_ising_derived(target, density)` | Ising 모델 기반 행-균형 QUBO |
| `solve_brute_force(Q, n)` | N <= 20일 때 완전 탐색 검증 |
| `batch_test(num_tests, n_bits)` | 배치 정확도 테스트 |
| `large_scale_analysis(n_tests, n_bits, model)` | 대규모 계수 분포 분석 (target pair별 off-diagonal 검증 포함) |
| `compare_with_random_qubo(n_bits, n_tests)` | 생성 QUBO vs 순수 랜덤 QUBO 비교 |

### 반환값

```python
Q = create_qubo_precise(target, density=1.0)
# Q: QUBO 딕셔너리 {(i,j): weight} (i <= j, 상삼각)
```

## SA 난이도 특성

### SA에 대해 Trivially 쉬운 문제

ZeroOffDiagonalModel 기준, **num_reads=1** (단 1회 SA 시도), num_sweeps=1000:

| N | 성공률 | 에너지비 | 해밍거리 | 시간 |
|---:|:------:|:-------:|:-------:|:----:|
| 10 | **100%** (10/10) | 1.0000 | 0 | 0.00s |
| 20 | **100%** (10/10) | 1.0000 | 0 | 0.00s |
| 50 | **100%** (10/10) | 1.0000 | 0 | 0.00s |
| 100 | **100%** (10/10) | 1.0000 | 0 | 0.01s |
| 200 | **100%** (10/10) | 1.0000 | 0 | 0.02s |
| 300 | **100%** (10/10) | 1.0000 | 0 | 0.04s |
| 500 | **100%** (10/10) | 1.0000 | 0 | 0.11s |
| 1000 | **100%** (10/10) | 1.0000 | 0 | 0.42s |

N=1000에서도 단 1회 시도로 100% 성공. 7개 방법론 중 가장 안정적으로 SA-trivial.

### 원인 분석: 왜 SA가 즉시 푸는가

#### 1. 로컬 미니마가 1개뿐 (= target 자체)

N=16에서 전수 탐색(2^16 = 65,536개 상태) 결과:

| 방식 | 로컬 미니마 수 | 비율 |
|------|:------------:|:----:|
| **Zero-Expectation** | **1개** | 0.0015% |
| Hardened Posiform (lin2, α=0.01) | **4개** | 0.0061% |
| Wishart (alpha=0.7) | **224개** | 0.3418% |

Zero-Expectation은 에너지 랜드스케이프 전체에 **함정(trap)이 하나도 없다**. SA가 아무 상태에서 출발해도 1-bit flip으로 에너지가 낮아지는 방향이 항상 존재하고, 그 방향을 따라가면 반드시 target에 도달한다.

#### 2. Frustration(조건 간 상충)이 0

> **Frustration**은 스핀 글래스 물리학의 정식 용어(Toulouse, 1977)로, **"모든 조건을 동시에 만족시킬 수 없는 상황"**을 뜻한다. QUBO에서는 **서로 다른 pair의 페널티가 한 변수에 대해 모순되는 방향을 요구하는 것**이다.

예시 — 변수 x1에 대해:

```
페어(1,2)의 Q_12가 "x1=1이 좋다"고 지시
페어(1,3)의 Q_13가 "x1=0이 좋다"고 지시
→ x1을 어떻게 놓아도 한쪽은 손해 = frustration
```

**Zero-Expectation에서 frustration이 불가능한 이유 — 페어별 독립 구성**:

```
페어(1,2): target[1], target[2]만 보고 독립적으로 페널티 생성
페어(1,3): target[1], target[3]만 보고 독립적으로 페널티 생성
페어(2,3): target[2], target[3]만 보고 독립적으로 페널티 생성

→ 모든 페어가 같은 target 방향을 가리킴
→ x1에 대해 모든 페어가 "target[1]로 가라"고 일관된 지시
→ 충돌 불가능
```

반면 **Wishart**는 공유 행렬 W에서 모든 Q_ij가 파생(J_ij = -(1/M) Σ_μ W_iμ W_jμ)되므로, 커플링이 상관되고 자연스럽게 상충이 발생한다.

N=20에서 실측:

| 방식 | Frustrated pairs |
|------|:---------------:|
| **Zero-Expectation** | **0 / 190** (0%) |
| Wishart (alpha=0.7) | 25 / 190 (~13%) |

**핵심**: 페어별 독립 구성 → 모든 페널티가 target 방향으로 정렬 → frustration = 0 → 로컬 미니마 없음.

#### 3. 에너지가 해밍거리에 선형 비례 (Funnel 구조)

N=16 전수 탐색, 해밍거리별 에너지 프로파일:

**Zero-Expectation:**

| 해밍거리 d | 평균 에너지 | 표준편차 | 신호/노이즈 |
|:---------:|:---------:|:------:|:---------:|
| 0 (target) | -109.41 | 0.00 | - |
| 1 | -72.38 | 8.90 | - |
| 2 | -37.34 | 10.11 | - |
| 3 | -4.29 | 9.60 | - |

에너지 증가 기울기: **~37/bit flip**, 셸 내 표준편차: **~9** → **신호(37) >> 노이즈(9)**. SA가 어디서든 "target 방향이 어딘지" 명확히 감지.

**Wishart (alpha=0.7):**

| 해밍거리 d | 평균 에너지 | 표준편차 | 신호/노이즈 |
|:---------:|:---------:|:------:|:---------:|
| 0 (target) | -2.43 | 0.00 | - |
| 1 | -1.16 | 0.49 | - |
| 2 | -0.06 | 1.04 | - |
| 3 | +0.87 | 1.50 | - |

에너지 증가 기울기: **~1.3/bit flip**, 셸 내 표준편차: **~1.5** → **신호(1.3) < 노이즈(1.5)**. SA가 방향을 잡을 수 없음. d=3에서 이미 d=0보다 낮은 에너지 상태가 존재.

#### 4. 에너지 갭이 거대함

N=20 전수 탐색:

| 방식 | Ground (E0) | 1st Excited (E1) | 절대 갭 | 상대 갭 |
|------|:----------:|:---------------:|:------:|:------:|
| **Zero-Expectation** | -329.22 | -290.30 | **38.92** | **11.8%** |
| Wishart (alpha=0.7) | -4.71 | -4.71 | **~0** | **~0%** |

#### 5. Single-bit flip 에너지 장벽

N=20, target으로부터 1비트 뒤집었을 때:

| 방식 | 최소 Delta E | 평균 Delta E | 최대 Delta E |
|------|:----------:|:----------:|:----------:|
| **Zero-Expectation** | **+39.67** | +54.21 | +93.25 |
| Wishart (alpha=0.7) | +0.61 | +1.21 | +2.37 |

ZeroExp의 최소 delta(39.67)가 Wishart의 최대 delta(2.37)보다 **17배** 큼.

### 구조적 원인 요약

**비유: 공을 굴린다고 생각하면**

```
볼록 이차함수 (매끄러운 그릇):     Zero-Expectation (자갈 깔린 깔때기):    Wishart (계란판):

에너지                              에너지                                 에너지
  |                                   |                                      |
  |\                                  |\.                                    |\   /\   /\
  | \                                 | '.:                                  | \_/  \_/  \
  |  \                                |   ':.                                |  ^    ^    \___
  |   \                               |    .:'.                              | trap  trap
  |    \___                            |      :._                             |
  +--------                           +--------                              +--------
  공이 미끄러져서                     자갈에 잠깐 걸려도                    움푹 파인 곳에 빠져서
  바닥에 도달                         옆으로 밀면 또 굴러감                 진짜 바닥에 못 감
```

- **볼록 이차함수**: 모든 방향이 내리막. 연속적이고 매끄러움.
- **Zero-Expectation**: 거시적으로는 깔때기지만, 같은 해밍거리 내에서 에너지가 ±σ만큼 출렁거림 (σ~9 vs 기울기~37). 개별 pair 수준에서는 wrong→correct 뒤집기가 에너지를 올릴 수도 있지만, n-1개 pair를 합산하면 **대수의 법칙**으로 항상 내리막 방향이 존재.
- **Wishart**: 커플링이 공유 벡터를 통해 상관되어 frustration(조건 간 상충) 발생. 곳곳에 local minima(trap)가 있어서 SA가 탈출 불가.

**왜 깔때기인가 — 한 줄 요약**:

```
"둘 다 틀린" penalty(ratio 2) > "하나만 틀린" penalty(ratio 1)
→ 틀린 bit 하나를 고치면, 관련된 모든 pair에서 평균적으로 에너지 감소
→ n이 크면 대수의 법칙으로 "평균적으로"가 "거의 확실히"로
→ local minimum이 target 외에 존재할 수 없음 (= frustration-free)
```

**핵심**: 페널티 구조에서 모든 pair가 독립적으로 target 방향을 가리키므로 frustration = 0, local minima = 1. Q 행렬의 1차 통계량(평균, 분산)을 위장하는 것은 에너지 **함수** 구조(상관 구조, frustration, local minima)를 위장하는 것과 전혀 다르다.

### 다른 방법론과의 비교

| 방법론 | SA 난이도 | GS 보장 | 구별 불가능 | 난이도 조절 | N=500 성공률 |
|--------|:---------:|:-------:|:----------:|:----------:|:----------:|
| **Zero Expectation** | **trivial** | 수학적 | **E[q_ij]=0** | X | **100%** |
| Posiform | trivial | 수학적 (유일) | 미분석 | X | 100% |
| Hard Posiform (α=0.01) | moderate | 수학적 (유일) | X | **O** | 90% |
| Quiet Planting (f=0.5) | medium | 조건부 | **O** (α<3.86) | **O** | **0%** |
| Wishart (α=0.7) | **hard** | 수학적 (유한정밀도) | X | **O** | 0% |
| McEliece (m=4,t=2) | **hard** | 조건부 | 미분석 | **O** | ~0% (k=8) |

> 전체 비교: [`docs/METHODOLOGY_COMPARISON.md`](../docs/METHODOLOGY_COMPARISON.md)

### 벤치마크로서의 위치

**장점**:
- E[q_ij] = 0으로 Q 행렬 비대각 항이 랜덤 QUBO와 통계적으로 동일 (ZeroOffDiagonalModel)
- 수학적 도출이 간결하고 검증 가능 (double-flip = single-flip 합)
- Strategy Pattern으로 다양한 페널티 전략 실험 가능
- 솔버 정확성 검증(sanity check)에 적합 — 이 문제를 못 풀면 솔버 구현에 문제가 있는 것

**한계**:
- SA가 N=500에서도 100% 성공 → 난이도 벤치마크로 부적합
- 비대각 무편향(E[q_ij]=0)과 대각 무편향을 동시에 달성 불가 → 완전한 통계적 위장 아님
- 솔버 간 성능 변별력이 전혀 없음

## 파일 구성

| 파일 | 역할 |
|------|------|
| `qubo_zero_expectation.py` | 생성기 (Strategy Pattern + PenaltyModel + Ising-derived) |
| `test_zero_expectation.py` | SA N-scaling 실험 |
| `test_diagonal_zero.py` | 대각 편향 분석 + 불가능성 증명 + 모델 간 비교 |
| `analyze_diagonal_deep.py` | 대각 편향 심층 분석 (편향=갭 동치, 노이즈 주입, Ising 비교, 하이브리드) |
| `analyze_q_structure.py` | Q 행렬 구조 비교 (ZeroExp vs Wishart: 고유값, rank, 행간 상관) |
| `results/` | 생성된 QUBO 파일 |

## 참고 문헌

### 직접 참고

1. **프로젝트 내부 PDF**: `기댓값0만들기_2트.pdf` — Zero-expectation QUBO 구성의 LP 최적화 방법론. 페널티 비율 도출 과정, 비대각 기댓값 0 조건의 수학적 유도.
2. **프로젝트 내부 PDF**: `읽어봐줘!.pdf` — QUBO 벤치마크 데이터셋 생성 동기 및 요구사항.

### QUBO 일반

3. **Kochenberger, G. A., et al.** "The unconstrained binary quadratic programming problem: a survey." *Journal of Combinatorial Optimization*, 28(1), 58-81, 2014.
4. **Glover, F., Kochenberger, G., & Du, Y.** "Quantum Bridge Analytics I: a tutorial on formulating and using QUBO models." *4OR*, 17(4), 335-371, 2019.

### 의사 불 최적화 (Pseudo-Boolean Optimization)

5. **Boros, E. & Hammer, P. L.** "Pseudo-Boolean optimization." *Discrete Applied Mathematics*, 123(1-3), 155-225, 2002. — 페널티 함수의 이론적 기반. 이진 변수 다항식의 최적화 프레임워크.

### Ising-QUBO 변환

6. 표준 치환: s_i = 2x_i - 1. Ising H = -sum J_ij s_i s_j → QUBO Q_ij = -4J_ij, Q_ii += 2J_ij.

## 사용법

```bash
# 기본 실행 (기본 목표 해)
python3 zero_expectation/qubo_zero_expectation.py

# 이진 목표 지정
python3 zero_expectation/qubo_zero_expectation.py 10110

# 랜덤 목표 (길이 50)
python3 zero_expectation/qubo_zero_expectation.py 50

# Ising-derived 행-균형 모드
python3 zero_expectation/qubo_zero_expectation.py 10110 balance

# SA 스케일링 실험
python3 zero_expectation/test_zero_expectation.py 10,20,50,100 10

# 대각 편향 분석 + 불가능성 증명
python3 zero_expectation/test_diagonal_zero.py

# Q 행렬 구조 분석 (ZeroExp vs Wishart)
python3 zero_expectation/analyze_q_structure.py
```
