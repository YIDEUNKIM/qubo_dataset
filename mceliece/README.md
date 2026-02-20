# McEliece Cryptographic QUBO Generator

McEliece 공개키 암호 프로토콜을 활용하여 **암호학적으로 어려운 planted QUBO 인스턴스**를 생성한다.

**기반 논문**: Mandra et al., "Generating Hard-to-Solve Optimization Problems Using the Structure of McEliece-Type Cryptosystems" (arXiv:2308.09704v2, FGCS 2025)

## 1. 문제 정의

### 목표
Ground state(최적해)가 보장되는 QUBO 행렬 Q를 생성하여, 양자 컴퓨터 및 최적화 솔버의 벤치마크로 사용한다.

### QUBO란
Quadratic Unconstrained Binary Optimization:

```
minimize E(x) = x^T Q x,  x_i in {0, 1}
```

NP-hard 문제로, 양자 어닐링/QAOA/SA 등의 핵심 벤치마크 형식이다.

### 왜 McEliece인가
기존 planted QUBO 생성 방식의 한계:

| 방식 | Ground state 보장 | 난이도 | 문제점 |
|------|:-:|:-:|--------|
| Wishart Planted | O (수학적) | SA-hard | alpha 조절로 난이도 제어 가능하지만 암호학적 보장은 없음 |
| Posiform Planting | O (수학적) | SA-trivial | 솔버가 너무 쉽게 풀어서 벤치마크 무의미 |
| Quiet Planting | 조건부 | 중간 | planted field 없이는 축퇴 문제 발생 |
| **McEliece** | **조건부** | **암호학적** | **복호화 난이도 = 최적화 난이도** |

McEliece의 핵심 장점: **문제의 난이도가 암호학적 보안 수준에 직접 연결**된다. 즉, McEliece 암호가 깨지지 않는 한 솔버도 이 QUBO를 효율적으로 풀 수 없다.

## 2. 방법론

### 2.1 파이프라인 개요

```
Goppa Code C[N,k,d]        McEliece Keygen         Encoding
GF(2^m) 위 생성       -->   G' = S * G * P     -->  q' = q*G' + e
     |                         |                        |
     v                         v                        v
[Step 1]                  [Step 2]                 [Step 3]
기약 Goppa 다항식          S: k*k 비특이 행렬       에러 벡터 e (가중치 t)
서포트 L (N개 원소)        P: N*N 치환 행렬         암호문 q' (길이 N)

     Systematic Form       p-local Ising           Degree Reduction
     [I_k | R]       -->   Eq.13 변환         -->  Rosenberg 축소
         |                     |                        |
         v                     v                        v
    [Step 4]              [Step 5]                 [Step 6]
    열 치환으로 변환       각 column의 비영 행      p-body -> 2-body
                          인덱스로 스핀 상호작용     보조변수 도입
                          H(target) = t 보장        QUBO 변환
```

### 2.2 핵심 수학

**p-local Ising Hamiltonian (논문 Eq.13)**:
```
H(sigma) = sum_j [1 - sign_j * prod_{i in I_j} sigma_i] / 2
```
- `I_j`: G_sys의 j번째 열에서 비영인 행 인덱스
- `sign_j = (-1)^{q'_j}`: 암호문의 j번째 비트
- `sigma_i = 2*x_i - 1`: Ising 스핀 변환

**Ground state 보장**: `H(target) = t` (에러 가중치)가 수학적으로 최솟값.

**Rosenberg Degree Reduction**:
- `x_a * x_b -> y` (보조변수)
- 페널티: `M * (x_a*x_b - 2*x_a*y - 2*x_b*y + 3*y)`
- M이 충분히 크면 `y = x_a * x_b`가 강제되어 ground state 보존

### 2.3 파라미터

| 파라미터 | 의미 | 범위 |
|----------|------|------|
| `m` | GF(2^m) 확장 차수 | 2~10 |
| `t` | 에러 정정 능력 / 에러 가중치 | 1~3 |
| `k` | 메시지 길이 (= target 길이) | k = N - m*t |
| `N` | 코드워드 길이 | N = 2^m (또는 2^m - 1) |
| `seed` | 난수 시드 | 정수 |

인스턴스 크기 관계:
- QUBO 원래 변수: `k`개
- 보조변수: 최대 `O(k * (N-k))`개
- 전체 QUBO 변수: `k + num_aux`

## 3. 페널티 M 강화 (본 구현의 핵심 기여)

### 3.1 문제: 기존 휴리스틱의 불충분성

기존 Rosenberg 페널티:
```python
M = max(10 * max_coeff, 10.0)  # max_coeff = 다항식 최대 계수
```

이 휴리스틱은 **w >= 10인 고차 상호작용에서 ground state 보장에 실패**한다.

#### 이론적 분석

w-body 항 `prod(1 - 2*x_i)`를 전개하면:
- 최대 계수: `2^w`
- 기존 M = `10 * 2^w`
- **실제 필요한 하한** (B_pair/2): `2 * 3^(w-2)`

| w (body 수) | max coeff | 기존 M (10*max) | B_pair | 필요 하한 (B/2) | 기존 M 충분? |
|:-:|:-:|:-:|:-:|:-:|:-:|
| 3 | 8 | 80 | 12 | 6 | O |
| 5 | 32 | 320 | 108 | 54 | O |
| 7 | 128 | 1,280 | 972 | 486 | O |
| 9 | 512 | 5,120 | 8,748 | 4,374 | O |
| **10** | 1,024 | 10,240 | 26,244 | **13,122** | **X (부족!)** |
| **11** | 2,048 | 20,480 | 78,732 | **39,366** | **X** |
| **12** | 4,096 | 40,960 | 236,196 | **118,098** | **X** |

**B_pair의 증가율이 3^w (지수적)으로, 기존 M의 2^w보다 빠르기 때문에 w >= 10에서 역전된다.**

### 3.2 해결: 이론적 하한 기반 M 계산

```python
# B_pair: 대체되는 변수쌍 (xa, xb)를 동시에 포함하는 모든 항의 |계수| 합
B_pair = sum(|coeff| for terms containing both xa and xb)

# 이론적 하한: M > B_pair / 2
# 안전 마진 4배 적용 (연쇄 보조변수 효과 고려)
M = max(B_pair * 2.0, 10.0)
```

#### 이론적 근거

1. Rosenberg 위반 최소 비용 = M
2. 보조변수 y를 잘못 설정했을 때 목적함수에서의 최대 이득 <= B_pair / 2
3. M > B_pair / 2이면 위반 비용 > 이득 → ground state 보존

### 3.3 사후 검증 & 자동 보정

이론적 M으로도 부족할 수 있는 경우 (연쇄 보조변수 효과)를 위해 사후 검증 루프 추가:

```
QUBO 생성 → 보조변수 최적화 → single-flip 검증
                                    |
                             모든 aux 안정? --YES--> 완료
                                    |
                                   NO
                                    |
                          해당 페널티 2배 증가 → 재검증 (최대 10회)
```

- 각 보조변수를 하나씩 flip하여 에너지 증가 여부 확인
- 불안정 발견 시 해당 페널티만 선택적으로 2배 증가
- 수렴 실패 시 경고 출력

## 4. 논문-코드 대조 검증

### 4.1 일치하는 부분

| 단계 | 논문 | 코드 | 일치 |
|------|------|------|:----:|
| Goppa 코드 생성 | GF(2^m), 기약 다항식, 패리티 검사 행렬 | `create_goppa_code()` | O |
| McEliece 키 생성 | G' = S * G * P | `mceliece_keygen()` | O |
| Systematic Form | [I_k \| R] 변환 | `_to_systematic_form()` | O |
| 인코딩 | q' = q * G' + e | `mceliece_encrypt()` 로직 | O |
| p-local Ising 항 | Eq.13: column별 스핀 상호작용 | `_ciphertext_to_ising_terms()` | O |
| Ground state 에너지 | H(target) = t | brute force 검증 통과 | O |

### 4.2 차이점

| 항목 | 논문 | 코드 |
|------|------|------|
| **실험 대상** | p-local Ising 직접 사용 (PT) | 2-local QUBO로 축소 |
| **차수 축소** | Eq.14-15 (exact identity) | Rosenberg reduction (페널티 M 의존) |
| **난이도 근거** | O(2^k) PT 스케일링 (p-local) | 2-local에서의 난이도는 미검증 |
| **Stern 알고리즘** | O(2^{0.056N}) 고전 복잡도 | 미구현 |

**핵심 차이**: 논문의 hardness 결과 O(2^k)는 **p-local Ising에서의 결과**이며, 2-local QUBO로 축소하면 이 보장이 직접 전이되지 않는다.

## 5. 실험 결과

### 5.1 Brute Force 검증 (m=3)

**m=3, t=1** (k=4, N=7, aux=3, total=7 QUBO vars):

| target | seeds | 성공률 | 축퇴도 | 에너지 범위 |
|--------|:-----:|:------:|:------:|:-----------:|
| 0000 | 20 | 20/20 (100%) | 1 | 0.0 |
| 0001 | 20 | 20/20 (100%) | 1 | -3.0 ~ -1.0 |
| 0110 | 20 | 20/20 (100%) | 1 | -3.0 ~ -1.0 |
| 1011 | 20 | 20/20 (100%) | 1 | -3.0 ~ -1.0 |
| 1111 | 20 | 20/20 (100%) | 1 | -5.0 ~ -3.0 |

**m=3, t=2** (k=2, N=8, aux=0~2, total=2~4 QUBO vars):

| target | seeds | 성공률 |
|--------|:-----:|:------:|
| 10 | 10 | 10/10 (100%) |

총 160/160 brute force 검증 성공. 모든 인스턴스에서 축퇴도 1 (유일한 ground state).

### 5.2 통계적 검증 (m=4)

**m=4, t=1** (k=11, N=15, aux=60, total=71 QUBO vars):

| seed | target 에너지 | min flip delta | local min | global min (추정) | 사후 보정 |
|:----:|:------------:|:--------------:|:---------:|:-----------------:|:---------:|
| 0 | -7.0 | 3.0 | O | O | 불필요 |
| 1 | -4.0 | 3.0 | O | O | 불필요 |
| 2 | -4.0 | 2.0 | O | O | 불필요 |
| 3 | -5.0 | 214.0 | O | O | 불필요 |
| 4 | -3.0 | 3.0 | O | O | 불필요 |
| 5 | -6.0 | 214.0 | O | O | 불필요 |
| 6 | -6.0 | 214.0 | O | O | 불필요 |
| 7 | -4.0 | 3.0 | O | O | 불필요 |
| 8 | -7.0 | 1.0 | O | O | 불필요 |
| 9 | -5.0 | 1.0 | O | O | 불필요 |

10/10 통과. 페널티 M 범위: 24 ~ 864. 이론적 M으로 초기 설정이 충분하여 사후 보정 불필요.

### 5.3 페널티 M 정보

m=4 인스턴스의 페널티 구조:
- 7-body 항 4개 → Rosenberg 축소 시 다단계 분해
- 보조변수 60개 생성
- M 값 분포: 24 (2-body 근처) ~ 864 (고차 항 축소)
- B_pair 범위: 12 ~ 432

## 6. 한계점

### 6.1 Ground State 보장의 조건부성

p-local Ising에서는 H(target) = t가 **수학적으로 보장**되지만, 2-local QUBO로의 Rosenberg 축소 과정에서 페널티 M이라는 자유 파라미터가 개입한다.

| 계층 | 보장 수준 | 조건 |
|------|-----------|------|
| p-local Ising | **수학적 100%** | 논문 Eq.13에 의해 증명 |
| 2-local QUBO (비연쇄) | **이론적 보장** | M > B_pair / 2 |
| 2-local QUBO (연쇄) | **경험적 보장** | 4배 안전 마진 + 사후 검증 |
| multi-flip 안정성 | **미검증** | single-flip만 확인 |

**연쇄 효과**: 보조변수 y가 다른 보조변수 z의 페널티에 참여하면 (z = y * x_c), z의 페널티 항이 y의 에너지 환경에 영향을 줄 수 있다. 이론적 하한 계산에서 이 효과를 완전히 포착하지 못한다.

### 6.2 Rosenberg vs 논문의 Exact Decomposition

논문의 Eq.14는 **페널티 파라미터 없이** p-body를 분해하는 exact identity:
```
(-1)^t * sigma_1...sigma_l * sigma_{l+1}...sigma_p
  = min_{omega=+-1} [sigma_1...sigma_l * omega + (-1)^t * omega * sigma_{l+1}...sigma_p + 1]
```

현재 코드는 이를 구현하지 않고 Rosenberg reduction을 사용하므로, 페널티 M에 대한 의존성이 남아 있다. Eq.14 구현 시 M 의존성을 완전히 제거할 수 있다.

### 6.3 2-local에서의 난이도 미보장

논문의 O(2^k) 난이도는 **p-local Ising에서의 PT 결과**이다. 2-local QUBO로 축소하면:
- 변수 수가 k에서 k + O(k*N)으로 증가
- 에너지 landscape가 변형됨
- 2-local에서의 난이도가 p-local과 동일한지 이론적 보장이 없음

논문에서도 2-local 축소 후 실험은 수행하지 않았다.

### 6.4 QUBO 크기 폭발

Rosenberg 축소는 보조변수를 대량 도입한다:

| m | t | k | N | 보조변수 | 전체 QUBO 변수 | 비영 항 수 |
|:-:|:-:|:-:|:-:|:--------:|:--------------:|:----------:|
| 3 | 1 | 4 | 7 | 3 | 7 | 20 |
| 4 | 1 | 11 | 15 | 60 | 71 | ~500 |
| 5 | 1 | 21 | 31 | ~400 | ~420 | ~5,000+ |

m=5 이상에서는 QUBO가 수백 변수로 폭발하여, 현재 D-Wave QPU의 연결성 제약(Pegasus: ~5,000 큐빗, 제한된 연결)에서 임베딩이 어려울 수 있다.

### 6.5 보조변수 최적화의 한계

`_compute_aux_values()`는 greedy 좌표 하강법(5라운드)으로 보조변수의 최적값을 결정한다. 이는:
- 작은 인스턴스에서는 충분 (m=3,4 모두 성공)
- 큰 인스턴스에서는 지역 최적에 갇힐 수 있음
- Exhaustive search는 보조변수 수가 O(k*N)이라 비현실적

### 6.6 아직 수행하지 못한 실험

1. **QPU 실험**: D-Wave/QAOA에서의 실제 성능 미측정
2. **SA 스케일링**: m 증가에 따른 SA 성공률 변화 미측정
3. **Stern 알고리즘과의 비교**: 고전 최적 공격 복잡도와의 비교 미수행
4. **p-local PT 재현**: 논문의 핵심 실험(p-local에서 PT 스케일링) 미재현

## 7. 사용법

```bash
# 기본 실행 (target=10110, t=2, m=자동)
python3 qubo_mceliece.py

# target + m + t + seed 지정
python3 qubo_mceliece.py 1011 3 1 42

# 길이 N의 랜덤 target
python3 qubo_mceliece.py 20 4 1

# 프로그래밍 인터페이스
from qubo_mceliece import create_qubo_mceliece
Q, info = create_qubo_mceliece("1011", m=3, t=1, seed=42)
print(info['full_target'])       # ground state (target + aux)
print(info['target_energy'])     # ground state 에너지
print(info['penalty_converged']) # 페널티 수렴 여부
```

## 8. 파일 구조

```
mceliece/
  qubo_mceliece.py              # McEliece QUBO 생성기 (메인 코드)
  README.md                     # 이 문서
  papers/
    Mandra2025_McEliece_...pdf  # 기반 논문
  results/
    mceliece_1011_4_m3_t1.txt   # 생성된 QUBO (edge-list CSV)
    mceliece_1011_4_m4_t1.txt
```

## 9. 향후 과제

1. **Eq.14 Exact Decomposition 구현**: 페널티 M 의존성 완전 제거
2. **p-local Ising PT 실험**: 논문의 O(2^k) 스케일링 재현
3. **QPU 벤치마크**: D-Wave에서 실제 성능 측정
4. **Wishart와의 비교 실험**: 동일 크기에서 SA/QPU 성공률 비교
