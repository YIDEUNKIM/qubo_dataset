# McEliece Cryptographic QUBO 생성기

## 개요

**Mandrà, Passarelli, Guerreschi (arXiv:2308.09704, FGCS 2025)**에 기반한 암호학적 QUBO 생성 방법론. McEliece 암호 프로토콜의 공개키를 Ising 스핀 시스템으로 캐스팅하여, **암호학적 보안에 기반한 computational hardness**를 가진 planted QUBO 인스턴스를 생성한다. Ground state 찾기가 McEliece 복호화와 동치이므로, 코드 이론의 난이도가 QUBO에 직접 반영된다.

## 이론적 배경

### McEliece 암호 시스템

1978년 McEliece가 제안한 코드 기반 공개키 암호 시스템. NIST 포스트양자 암호 표준화 후보 중 하나인 Classic McEliece의 기반이다.

**핵심 구성요소**:
- **비밀키**: Goppa 코드의 생성 행렬 G, 스크램블 행렬 S, 치환 행렬 P
- **공개키**: G' = S·G·P (구조가 숨겨진 생성 행렬)
- **인코딩**: q' = target·G' + ε (메시지에 에러 추가)
- **복호화**: 비밀키 없이 에러 벡터 ε를 제거하는 것이 NP-hard

### Goppa 코드

GF(2^m) 위의 이진 Goppa 코드 C[N, k, d]:

| 파라미터 | 의미 | 관계 |
|---------|------|------|
| N | 코드워드 길이 | 2^m (또는 2^m - 1, t에 따라) |
| k | 메시지 차원 | ≥ N - m·t |
| d | 최소 거리 | ≥ 2t + 1 |
| m | GF(2^m) 확장 차수 | 유한체 크기 결정 |
| t | 에러 정정 능력 | 정정 가능한 최대 에러 수 |

**Goppa 다항식** g(z): GF(2^m) 위의 degree-t 기약다항식. 서포트 L = {α ∈ GF(2^m) : g(α) ≠ 0}에 대해 패리티 검사 행렬 H를 구성.

### QUBO 변환 파이프라인

```
1. Goppa 코드 C[N, k, d] 생성 (GF(2^m) 위)
     ↓
2. McEliece 키 생성: G' = S·G·P
     ↓
3. 체계적 형태 변환: G_sys = [I_k | R]
     ↓
4. 인코딩: q' = target·G_sys + ε
     ↓
5. p-local Ising 항 생성
   각 컬럼 j: H_j = [1 - (-1)^{q'_j} · Π_{i∈I_j} σ_i] / 2
   I_j = 컬럼 j에서 비영인 행 인덱스
     ↓
6. Rosenberg 차수 축소: p-body → 2-body + 보조변수
   x_a·x_b → y, 패널티: M·(x_a·x_b - 2·x_a·y - 2·x_b·y + 3·y)
     ↓
7. Ising → QUBO 변환
```

### 체계적 형태와 Ising 항

공개키를 체계적 형태 [I_k | R]로 변환하면:

- **처음 k개 컬럼** (I_k 부분): 가중치 1 → **1-body 항** (차수 축소 불필요)
- **나머지 N-k개 컬럼** (R 부분): 가중치 w ≥ 2 → **w-body 항** → Rosenberg 축소 필요

체계적 형태에서 첫 k비트가 메시지(target)에 대응하므로, ground state의 첫 k비트가 정확히 target이 된다.

### 차수 축소 (Rosenberg Reduction)

3차 이상의 Ising 상호작용을 2차(QUBO)로 변환:

1. 가장 빈번한 변수 쌍 (x_a, x_b) 선택
2. 보조변수 y 도입: x_a·x_b → y
3. 패널티 항 추가: M·(x_a·x_b - 2·x_a·y - 2·x_b·y + 3·y)
4. M이 충분히 크면 y = x_a·x_b가 최적에서 보장됨
5. 차수가 2 이하가 될 때까지 반복

**비용**: 컬럼 가중치 w인 항의 축소에 보조변수 O(w-1)개 필요. w는 R 행렬의 컬럼별 1의 개수로, m이 커지면 지수적으로 증가할 수 있다.

### 다른 방법론과의 차이

| 특성 | McEliece | Wishart | Posiform | Quiet Planting |
|------|:--------:|:-------:|:--------:|:--------------:|
| 난이도 근거 | **암호학적 보안** | Metastable 상태 | 없음 (P 문제) | NP-completeness (3-SAT) |
| 보조변수 | **있음** (대량) | 없음 | 없음 | 있음 (Rosenberg) |
| QUBO 크기 | k + aux | n | n | n(1+alpha) |
| GS 보장 | 조건부 (M 의존) | 수학적 | **수학적 (유일)** | 조건부 (field 의존) |
| 난이도 조절 | m, t | alpha | coeff_range | alpha, field |
| 확장성 | **m≤4로 제한** | 임의 N | 임의 N | 임의 N |

## 구현 방식

### 전체 구조

```
qubo_mceliece.py 내부 구조:

1. GF(2^m) 유한체 클래스 (GF2m)
   - log/antilog 룩업 테이블로 곱셈/나눗셈 O(1)
   - 기약다항식 테이블: m=2~10 지원

2. Goppa 코드 생성 (create_goppa_code)
   - 랜덤 기약 Goppa 다항식 탐색
   - 패리티 검사 행렬 H 구성 (GF(2^m) → 이진 전개)
   - H의 영공간으로 생성 행렬 G 계산

3. McEliece 키 생성 (mceliece_keygen)
   - S: 랜덤 가역 k×k 이진 행렬
   - P: 랜덤 N×N 치환 행렬
   - G_pub = S·G·P

4. 체계적 형태 변환 (_to_systematic_form)
   - Gaussian elimination with column pivoting
   - G_sys = [I_k | R]

5. p-local Ising 항 생성 (_ciphertext_to_ising_terms)
   - 각 컬럼의 비영 행 인덱스로 상호작용 항 구성

6. Rosenberg 차수 축소 (_reduce_multibody_to_qubo)
   - 반복적 쌍 축소: 가장 빈번한 변수쌍 우선
   - 패널티 강도 M = 10 × max(|coeff|)
```

### 핵심 파라미터

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `m` | 자동 선택 | GF(2^m) 확장 차수. N=2^m. **m≥5는 생성 비용 O(2^w)로 비실용적** |
| `t` | 2 | 에러 정정 능력. k ≈ N - m·t. t↑ → k↓, 보조변수 수 변화 |
| `seed` | None | 재현성을 위한 난수 시드 (Goppa, 키, 인코딩에 분리 적용) |

### 주요 함수

| 함수 | 설명 |
|------|------|
| `GF2m(m)` | GF(2^m) 유한체 (log/antilog 테이블, 곱셈/나눗셈/다항식 연산) |
| `create_goppa_code(m, t, seed)` | Goppa 코드 생성 (G, goppa_poly, support 반환) |
| `mceliece_keygen(G, seed)` | McEliece 키 생성 (G_pub, S, P, perm 반환) |
| `_to_systematic_form(G_pub)` | 공개키 → 체계적 형태 [I_k \| R] 변환 |
| `_ciphertext_to_ising_terms(G_sys, ciphertext, k)` | p-local Ising 항 생성 |
| `_reduce_multibody_to_qubo(terms, k)` | Rosenberg 차수 축소 → QUBO 변환 |
| `create_qubo_mceliece(target, m, t, seed)` | **메인 진입점** |
| `extract_original_solution(sample, k)` | SA 결과에서 원래 k개 변수 추출 |
| `verify_ground_state(Q, target, info)` | 통계적 GS 검증 (single-flip + 랜덤 샘플) |
| `verify_brute_force(Q, target, info)` | Brute force GS 완전 검증 (total_vars ≤ 20) |
| `_compute_aux_values(Q, target, k, num_aux)` | 보조변수 최적값 탐욕적 결정 |

### 반환값

```python
Q, info = create_qubo_mceliece(target, m=4, t=2, seed=42)
# Q: QUBO 딕셔너리 {(i,j): weight} (i ≤ j, upper triangular)
# info: {
#   'n':              k (메시지 차원),
#   'N':              코드워드 길이 (support 크기),
#   'm':              GF(2^m) 확장 차수,
#   't':              에러 정정 능력,
#   'num_aux':        보조변수 수,
#   'total_vars':     k + num_aux (전체 QUBO 변수 수),
#   'target_energy':  target + 최적 aux에서의 QUBO 에너지,
#   'constant_offset': Ising → QUBO 변환 상수항,
#   'error_weight':   에러 벡터의 해밍 가중치,
#   'goppa_poly':     Goppa 다항식 계수 리스트,
#   'full_target':    target + aux 전체 비트스트링,
# }
```

### 보조변수 처리

McEliece QUBO는 원래 변수(메시지 k비트)와 보조변수(Rosenberg 축소에서 도입)로 구성된다. SA 결과에서 원래 변수만 추출할 때:

```python
# SA 실행 후
best_sample = ss.first.sample
found_orig = extract_original_solution(best_sample, info['n'])  # 처음 k개 변수
expected_orig = info['full_target'][:info['n']]  # 기대되는 원래 변수값

# 성공 판정
if found_orig == expected_orig:
    print("Ground state 발견!")
```

## SA 실험 결과

### 실험 1: m-Scaling (GF(2^m) 확장 차수 vs SA 성공률)

t=2 고정, num_reads=200, sweeps=max(1000, 10×total_vars):

| m | N(=2^m) | k | total_vars | aux | SA 성공률 | 비고 |
|---|:-------:|:-:|:----------:|:---:|:---------:|------|
| 3 | 8 | 2 | 2 | 0 | **100%** | 변수 2개, trivial |
| 4 | 16 | 8 | 33~46 | 25~38 | **~0%** | 보조변수 25+개, SA-hard |
| 5+ | 32+ | — | — | — | 생성 불가 | Rosenberg O(2^w) 비용 |

**핵심 관찰**: m=3→4에서 SA 난이도가 급격히 증가. m=3에서는 Goppa 다항식의 degree가 2이므로 R 행렬의 컬럼 가중치가 낮아 보조변수가 거의 없다. m=4에서는 k=8, R 행렬의 컬럼이 복잡해져 보조변수가 25~38개로 급증한다.

### 실험 2: t-Parameter Sweep (에러 정정 능력 vs SA 성공률)

m=4 고정, num_reads=200:

| t | k(≈N-m·t) | total_vars | aux | SA 성공률 | 해석 |
|---|:---------:|:----------:|:---:|:---------:|------|
| 1 | 11 | 71 | 60 | **0%** | 가장 많은 보조변수, 가장 어려움 |
| 2 | 8 | 33~46 | 25~38 | **~0-50%** | 여전히 어려움 |
| 3 | 4 | 8~11 | 4~7 | **~100%** | 보조변수 소수, SA-easy |

**역설적 결과**: t↑ → k↓ + 에러 정정 능력↑이므로 암호학적 보안은 증가하지만, SA 난이도는 **감소**한다. 이유: t가 커지면 k가 줄어 메시지 비트가 적어지고, R 행렬도 작아져 보조변수 수가 급감하기 때문이다. SA 난이도는 암호학적 난이도가 아닌 **보조변수 수(QUBO 크기)**에 의해 결정된다.

### 실험 3: SA Sweep Count 전이

m=4, t=2, reads=50/sweep, 인스턴스 10개:

| sweeps | 1 | 5 | 10 | 50 | 100 | 500 | 1000 | 5000 | 10000 |
|:------:|---|---|----|----|-----|-----|------|------|-------|
| GS rate | 0% | 0% | 0% | 0% | 0% | 0% | 0% | ~2% | ~5% |

10000 sweeps에서도 ~5%에 불과. Wishart(alpha=0.7)과 유사한 수준의 SA-hard 행동을 보인다.

다른 config들의 비교:

| Config | Sweep 전이 패턴 |
|--------|----------------|
| m=3, t=1 (k=4, aux=3) | sweeps≥50에서 ~60-80% 성공 (적당히 어려움) |
| m=3, t=2 (k=2, aux=0) | 모든 sweep에서 100% (trivial — 변수 2개) |
| m=4, t=2 (k=8, aux=25+) | sweeps=10000에서도 ~5% (SA-hard) |
| m=4, t=3 (k=4, aux=4~7) | sweeps≥100에서 ~100% (보조변수 소수) |

### 실험 4: 6-Way 방법론 비교 (n_bits=8)

n_bits=8, num_reads=200, num_sweeps=1000:

| Method | 성공률 | 설명 |
|--------|:------:|------|
| **Posiform** | 100% | SA-trivial, 유일한 GS 보장 |
| **Zero Expectation** | 100% | SA-trivial, 구별 불가능 |
| **Hardened (lin2, α=0.1)** | ~100% | α=0.1에서 N=8은 easy |
| **Quiet Planting (α=4.2)** | ~60-100% | N=8은 임계점 미만 |
| **Wishart (α=0.7)** | ~30-50% | N=8은 아직 어려운 구간 진입 전 |
| **McEliece (m=4, t=2)** | **~0%** | k=8이지만 total_vars=33+, SA-hard |

**주목할 점**: McEliece는 동일한 "원래 변수 수" n_bits=8에서 다른 방법론보다 압도적으로 어렵다. 이는 보조변수가 탐색 공간을 2^8에서 2^33+으로 확장하기 때문이다. 단, 이 비교는 QUBO 크기가 동일하지 않으므로 공정한 비교가 아님에 유의.

### SA-Hard인 이유

1. **대규모 보조변수**: m=4, t=2에서 원래 변수 8개에 보조변수 25~38개가 추가되어 총 33~46개 QUBO 변수. 탐색 공간이 2^33 이상.

2. **페널티 계층 구조**: Rosenberg 차수 축소에서 생성된 패널티 항들이 에너지 지형에 복잡한 계곡과 장벽을 형성. 패널티 강도 M이 크면 보조변수가 올바른 값에 고정되지만, SA가 원래 변수와 보조변수를 동시에 올바르게 맞춰야 하므로 어렵다.

3. **상관된 보조변수**: 보조변수 y = x_a·x_b의 값이 원래 변수 값에 의존. SA가 원래 변수를 뒤집으면 연관된 모든 보조변수도 함께 업데이트되어야 하지만, SA는 한 번에 1비트만 뒤집는다.

4. **암호학적 난이도와의 관계**: 논문에서는 p-local Ising 시스템에서 O(2^k) parallel tempering 스케일링을 보였다. 현재 구현의 Rosenberg 방식은 이 암호학적 난이도에 **추가로** 보조변수 비용까지 더하므로, 실제 암호학적 난이도보다 더 어려운 문제를 만든다.

## 제한사항 및 미래 작업

### 현재 제한사항

1. **m≥5 생성 불가**: Rosenberg 차수축소의 비용이 O(2^w) (w = R 행렬의 컬럼 가중치). m=5에서 w가 최대 22까지 증가하여 단일 항의 축소에 수백만 연산 필요. 60초+ 소요로 실험 불가.

2. **페널티 M 의존적 GS 보장**: Rosenberg 패널티 강도 M이 충분히 크지 않으면 보조변수가 올바른 값에 고정되지 않아 ground state가 깨질 수 있다. 현재 M = 10 × max(|coeff|)로 설정하며, brute force 검증(m=3)과 통계적 검증(m=4)으로 확인.

3. **보조변수에 의한 공정성 문제**: McEliece QUBO의 SA 난이도는 암호학적 난이도보다 보조변수 수에 의해 결정됨. 다른 방법론과 동일 조건 비교가 어렵다.

### Eq.14 Exact Decomposition (미구현)

논문 Section 2의 Eq.14는 p-body Ising 상호작용을 **보조변수 없이** 2-body로 분해하는 exact identity:

```
Π_{i∈S} σ_i = Σ_{T⊆S} (-2)^{|S|-|T|} · Σ_{i∈T} σ_i / C
```

이를 구현하면:
- 보조변수 완전 제거 → QUBO 크기 = k (원래 변수 수만)
- 페널티 M 의존성 제거 → GS 보장이 무조건적
- m≥5 확장 가능 → 대규모 실험 가능
- **순수한 암호학적 난이도**만 반영하는 QUBO 생성

이는 현재 가장 높은 우선순위의 미구현 기능이다.

## 파일 구성

| 파일 | 역할 |
|------|------|
| `qubo_mceliece.py` | 생성기 (GF(2^m) 유한체 → Goppa 코드 → McEliece 키 → Ising → Rosenberg → QUBO) |
| `test_mceliece.py` | SA 실험 (m-scaling, t-sweep, sweep 전이, 6-way 비교) |
| `results/` | 생성된 QUBO 파일 (edge-list CSV 형식) |

## 참고 문헌

### 핵심 논문

1. **Mandrà, S., Passarelli, G., & Guerreschi, G. G.** "Ising formulation of the McEliece problem." *Future Generation Computer Systems (FGCS)*, 2025. [arXiv:2308.09704](https://arxiv.org/abs/2308.09704)
   - McEliece 공개키를 Ising 스핀 시스템으로 캐스팅
   - p-local Ising 상호작용에서 O(2^k) parallel tempering 스케일링
   - Eq.14 exact decomposition (보조변수 없는 2-body 변환)

### 코드 이론

2. **McEliece, R. J.** "A public-key cryptosystem based on algebraic coding theory." *DSN Progress Report*, 42-44, 1978.
   - McEliece 암호 시스템 원논문

3. **Goppa, V. D.** "A new class of linear correcting codes." *Problemy Peredachi Informatsii*, 6(3), 24-30, 1970.
   - Goppa 코드의 정의

### 차수 축소

4. **Rosenberg, I. G.** "Reduction of bivalent maximization to the quadratic case." *Cahiers du Centre d'Etudes de Recherche Operationnelle*, 17, 71-74, 1975.
   - 고차 의사불 함수의 2차화 (x_a·x_b → y + 패널티)

### 관련 방법론 (이 프로젝트 내)

5. **Hahn, G., Pelofske, E., & Djidjev, H.** "Using 2-SAT to generate QUBO instances with known optimal solutions." *IEEE QCE*, 2023. — Posiform planting
6. **Pelofske, E., Hahn, G., & Djidjev, H.** "Increasing the Hardness of Posiform Planting Using Random QUBOs." *npj Unconventional Computing*, 2025. — Hardened Posiform
7. **Hamze, F. et al.** "Wishart Planted Ensemble." *Physical Review E*, 2020. — Wishart planted ensemble

## 사용법

```bash
# 기본 생성 (target="10110", m/t 자동 선택)
python3 mceliece/qubo_mceliece.py 10110

# m, t, seed 지정
python3 mceliece/qubo_mceliece.py 10110 4 2 42

# 길이로 랜덤 목표 생성
python3 mceliece/qubo_mceliece.py 8

# SA m-scaling 실험 (m=3,4에서 SA 성공률)
python3 mceliece/test_mceliece.py --m-scaling 10

# SA t-parameter sweep (m=4, t=1,2,3)
python3 mceliece/test_mceliece.py --t-sweep 10

# SA sweep 전이 실험 (S-curve)
python3 mceliece/test_mceliece.py --sweep 10

# 6-way 방법론 비교
python3 mceliece/test_mceliece.py --compare 10
```

### 출력 예시 (m=4, t=2, target="10110")

```
============================================================
McEliece Cryptographic QUBO 생성기
============================================================
Target: 10110 (k=5)
t=2 (m은 자동 선택)

[코드 정보]
  GF(2^4), N=16, k=8
  에러 가중치 t=2
  보조변수 수: 28
  전체 QUBO 변수: 36
  QUBO 비영 항 수: 203

[에너지 검증]
  Target 에너지: -6.000000
  상수 오프셋: 8.000000

[통계적 검증] (전체 변수 36개)
  Target 에너지: -6.000000
  최소 flip delta: 1.000000
  Local minimum: OK
  랜덤 5000개 중 더 낮은 에너지: 0개
  Global minimum (추정): OK
```
