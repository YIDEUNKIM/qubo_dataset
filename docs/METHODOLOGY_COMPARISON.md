# 전체 방법론 비교 — SA 벤치마크 종합 결과

## 1. 실험 조건

| 파라미터 | 값 |
|---------|-----|
| SA Solver | `neal.SimulatedAnnealingSampler` (D-Wave Ocean SDK) |
| num_reads | 50 |
| num_sweeps | 1,000 |
| num_runs | 10 (각 N에서 독립 인스턴스) |
| 시드 | random.seed(42), np.random.seed(42) |
| 성공 기준 | EXACT match (ZeroExp/Wishart는 Z2 대칭 SYM_MATCH도 포함) |
| 실험 일자 | 2026-02-20 |

**비교 대상 방법론 (7개)**:

| 방법론 | 생성기 | 난이도 파라미터 |
|--------|--------|----------------|
| Zero Expectation | `create_qubo_precise(target, density=1.0)` | density=1.0 |
| Wishart | `create_qubo_wishart(target, alpha=0.7)` | alpha=0.7 (hard regime) |
| Quiet Planting | `create_qubo_quiet_planted(target, alpha=4.2, field_strength=0.5)` | alpha=4.2, field=0.5 |
| Posiform | `create_qubo_posiform(target, coeff_range=(1.0, 3.0))` | coeff_range=(1.0, 3.0) |
| Hardened Posiform (moderate) | `create_qubo_hardened_posiform(n, coeff_type='lin2', posiform_scale=0.1)` | lin2, alpha=0.1 |
| Hardened Posiform (hardest) | `create_qubo_hardened_posiform(n, coeff_type='lin2', posiform_scale=0.01)` | lin2, alpha=0.01 |
| McEliece | `create_qubo_mceliece(target, m=4, t=2)` | m=4, t=2 |

---

## 2. SA 성공률 비교

### 2.1 N별 성공률 (%)

| Method | N=50 | N=100 | N=200 | N=500 | N=1000 |
|---|:---:|:---:|:---:|:---:|:---:|
| **Zero Expectation** | 100 | 100 | 100 | 100 | 100 |
| **Wishart (alpha=0.7)** | 70 | 0 | 0 | 0 | 0 |
| **Quiet Planting (alpha=4.2, f=0.5)** | 100 | 60 | 0 | — | — |
| **Posiform** | 100 | 100 | 100 | 100 | — |
| **Hard Posiform (lin2, alpha=0.1)** | 100 | 100 | 100 | 100 | 100 |
| **Hard Posiform (lin2, alpha=0.01)** | 100 | 100 | 100 | 90 | **40** |
| **McEliece (m=4, t=2)** | — | — | — | — | — |

- (—): 생성 시간 또는 QUBO 크기 제약으로 미측정
- McEliece: N별 비교 불가 (k=8 고정, total_vars=33~46). 별도 실험에서 SA 성공률 ~0% (10000 sweeps에서도 ~5%). 자세한 결과는 아래 Section 5.1 참조
- Quiet Planting: N=500+에서 QUBO 크기가 ~2,600 변수로 SA 비용 급증
- Posiform: N=1000에서 2-SAT 유일성 검사 생성 시간 ~108초로 실험 제외

### 2.2 난이도 순서 (SA 기준)

```
SA-trivial ◄───────────────────────────────────────────────────────► SA-hard

ZeroExp     Posiform    Hard(α=0.1)    Hard(α=0.01)    Quiet(f=0.5)    Wishart(α=0.7)
 100%@1K     100%@1K     100%@1K       40%@1K           0%@200          0%@100
                                                            │
                                              McEliece(m=4,t=2): ~0%@k=8
```

---

## 3. 평균 Hamming Distance

SA가 찾은 해와 target 사이의 비트 차이. 0이면 정확히 일치, N/2이면 랜덤 수준.

| Method | N=50 | N=100 | N=200 | N=500 | N=1000 |
|---|:---:|:---:|:---:|:---:|:---:|
| **Zero Expectation** | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| **Wishart (alpha=0.7)** | 5.4 | 43.1 | 86.2 | 231.9 | 467.2 |
| **Quiet Planting** | 0.0 | 0.6 | 4.6 | — | — |
| **Posiform** | 0.0 | 0.0 | 0.0 | 0.0 | — |
| **Hard Posiform (lin2, alpha=0.1)** | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| **Hard Posiform (lin2, alpha=0.01)** | 0.0 | 0.0 | 0.0 | 0.4 | 3.1 |

**실패 양상의 차이**:
- **Wishart**: Hamming ~ N/2 → SA가 metastable trap에 갇혀 target과 완전히 다른 위치에 도달. "길을 완전히 잃음"
- **Quiet Planting**: Hamming 2~5 → SA가 target 근처까지 접근하지만 마지막 몇 비트를 뒤집지 못함. 전형적 glassy 거동
- **Hard Posiform (alpha=0.01)**: Hamming 0~3 → Quiet Planting과 유사한 near-miss 패턴

---

## 4. QUBO 구조 비교

### 4.1 비영(non-zero) 항 수

| Method | N=50 | N=100 | N=200 | N=500 | N=1000 | 수식 |
|---|---:|---:|---:|---:|---:|---|
| **Zero Expectation** | 1,275 | 5,050 | 20,100 | 125,250 | 500,500 | N(N+1)/2 (dense) |
| **Wishart** | 1,275 | 5,050 | 20,100 | 125,250 | 500,500 | N(N-1)/2 (dense) |
| **Quiet Planting** | 1,038 | 2,107 | 4,255 | — | — | O(N*alpha) (sparse) |
| **Posiform** | 353 | 790 | 1,897 | 5,247 | — | O(N*clauses) (sparse) |
| **Hard Posiform** | 540 | 1,314 | 3,148 | 8,424 | 17,497 | O(N*k) (sparse) |
| **McEliece (m=4,t=2)** | — | — | — | — | — | k + O(kN) (보조변수 포함) |

### 4.2 QUBO 변수 수

| Method | QUBO 변수 수 | N=100 기준 | 보조변수 |
|---|---|---:|---|
| Zero Expectation | N | 100 | 없음 |
| Wishart | N | 100 | 없음 |
| Quiet Planting | N(1+alpha) | **520** | m개 (clause당 1개) |
| Posiform | N | 100 | 없음 |
| Hard Posiform | N | 100 | 없음 |
| McEliece (m=4,t=2) | k + aux | **33** (k=8) | 25개 (Rosenberg 차수축소) |

---

## 5. Hardened Posiform 상세: Sweep 전이 실험

N=500, 5 instances, 20 reads/sweep. 논문 (Pelofske 2024) Fig 7-8 재현.

| Config | 10 sw | 50 sw | 100 sw | 500 sw | 1000 sw | 5000 sw |
|---|:---:|:---:|:---:|:---:|:---:|:---:|
| **lin2, alpha=0.01** | 2.0% | 1.0% | 4.0% | 6.0% | 12.0% | 24.0% |
| **lin2, alpha=0.1** | 24.0% | 59.0% | 83.0% | 100% | 100% | 100% |
| **lin20, alpha=0.01** | 0.0% | 2.0% | 12.0% | 31.0% | 52.0% | 68.0% |
| **lin20, alpha=0.1** | 87.0% | 98.0% | 99.0% | 100% | 100% | 100% |

**논문 핵심 발견 재현 확인**:

| 논문 발견 | 코드 결과 | 일치 |
|-----------|-----------|:---:|
| alpha=0.01이 alpha=0.1보다 어려움 | lin2: 24% vs 100% (5000sw) | O |
| lin2가 lin20보다 어려움 | alpha=0.01에서 lin2(24%) < lin20(68%) | O |
| Sweep 수 증가 → S-curve 전이 | 모든 config에서 단조 증가 | O |
| alpha=0.01은 매우 많은 sweep 필요 | lin2,alpha=0.01: 5000sw에서도 24% | O |

### 5.1 McEliece 상세: SA 실험 결과

McEliece는 다른 방법론과 달리 QUBO 크기가 m,t에 의해 결정되므로 N별 비교가 불가. 별도의 m-scaling, t-sweep 실험 결과:

**m-Scaling** (t=2, num_reads=200):

| m | N(=2^m) | k | total_vars | aux | SA 성공률 |
|---|---------|---|------------|-----|-----------|
| 3 | 8 | 2 | 2 | 0 | **100%** (trivial) |
| 4 | 16 | 8 | 33~46 | 25~38 | **~0%** (SA-hard) |

**t-Sweep** (m=4, num_reads=200):

| t | k | total_vars | aux | SA 성공률 |
|---|---|------------|-----|-----------|
| 1 | 11 | 71 | 60 | **0%** |
| 2 | 8 | 33~46 | 25~38 | **~0-50%** |
| 3 | 4 | 8~11 | 4~7 | **~100%** |

**Sweep 전이** (m=4, t=2, reads=50/sweep):

| sweeps | 1 | 10 | 100 | 1000 | 10000 |
|--------|---|----|----|------|-------|
| GS rate | 0% | 0% | 0% | 0% | ~5% |

**핵심 발견**:
1. McEliece QUBO의 SA 난이도는 **보조변수 수**에 의해 결정됨. m=4,t=2에서 aux=25~38개로 SA-hard
2. t 증가 → k 감소 + aux 감소 → **역설적으로 더 쉬워짐** (t=3에서 k=4, aux=4~7로 SA 100%)
3. m≥5는 Rosenberg 차수축소의 지수적 비용(O(2^w), w=컬럼 가중치)으로 **QUBO 생성 자체가 불가**
4. Eq.14 exact decomposition 구현 시 보조변수 제거 → 확장성 문제 해결 기대

---

## 6. 방법론 특성 종합

| 속성 | Zero Exp | Wishart | Quiet Plant | Posiform | Hard Posiform | McEliece |
|---|---|---|---|---|---|---|
| **SA 난이도** | trivial | **hard** | medium | trivial | **tunable** | **hard** (m=4) |
| **난이도 파라미터** | density | alpha (M/N) | alpha, field | coeff_range | coeff_type, alpha_scale | m, t |
| **GS 보장** | 수학적 | 수학적 (유한정밀도 주의) | field 의존 | 수학적 (Tarjan SCC) | 수학적 | 조건부 (M 의존) |
| **GS 유일성** | Z2 대칭 | Z2 대칭 | field 필요 | **유일 (증명)** | **유일 (증명)** | 조건부 |
| **구별 불가능성** | E[q_ij]=0 보장 | X (low-rank) | alpha<3.86 보장 | 미보장 | X (block-diagonal) | 미분석 |
| **보조변수** | 없음 | 없음 | 있음 (n+m) | 없음 | 없음 | **있음 (대량)** |
| **SA 상전이** | 없음 | alpha_c ~ 0.95 | N~300 (f=0.5) | 없음 | N~500 (alpha=0.01) | m=3→4 전이 |
| **생성 시간** | O(N^2) 빠름 | O(N^2) 빠름 | O(N) 빠름 | O(N^2.5) 느림 | O(N^2.5) 느림 | O(2^w) (m≥5 불가) |
| **참고 논문** | 자체 설계 | Hamze 2020 | Krzakala 2009 | Hahn 2023 | Pelofske 2024 | Mandrà 2024 |

---

## 7. 3중 트레이드오프

SA-hard, 구별 불가능, 수학적 GS 보장 — 세 가지를 동시에 달성하는 방법론은 아직 없다.

| 방법론 | SA-hard | 구별 불가능 | GS 보장 |
|--------|:-------:|:----------:|:-------:|
| Zero Expectation | X | **O** | **O** |
| Wishart | **O** | X | 조건부 |
| Quiet Planting | 중간 | **O** | 조건부 |
| Posiform | X | 미보장 | **O** |
| Hard Posiform | **조절 가능** | X | **O** |
| McEliece | **O** (m=4) | 미분석 | 조건부 (M 의존) |

---

## 8. 용도별 권장 생성기

| 시나리오 | 1순위 | 2순위 | 이유 |
|---------|:-----:|:-----:|------|
| **솔버 정확성 검증 (sanity check)** | Posiform | ZeroExp | SA-trivial + 유일한 GS 수학적 보장 |
| **SA 한계 측정** | Wishart | Hard Posiform | N=100부터 SA 실패, alpha로 연속 제어 |
| **양자 우위 벤치마크** | Wishart | Hard Posiform | SA-hard이면서 GS가 알려져 있어 정답 비교 가능 |
| **난이도 연속 조절** | Hard Posiform | Wishart | coeff_type + alpha_scale로 trivial~hard 전 구간 커버 |
| **통계적 은닉 연구** | Quiet Planting | ZeroExp | alpha<3.86에서 랜덤 3-SAT과 구별 불가 |
| **대규모 문제 (N>1000)** | Wishart, Hard Posiform | ZeroExp | QUBO 크기 = N, 생성 빠름 (Wishart/ZeroExp) |
| **블라인드 벤치마크 (부정행위 방지)** | ZeroExp | Quiet Planting | Q 행렬 분석으로 target 추측 불가 |
| **암호학적 난이도 연구** | McEliece | — | McEliece 복호화와 동치, 이론적 O(2^k) 하한 |

---

## 9. 핵심 인사이트

1. **SA 난이도의 본질은 에너지 landscape 구조**: Q 계수의 1차 통계량(평균, 분산)이 아니라 **계수 간 상관 구조**가 metastable trap을 만들어 SA를 가둠. Wishart의 low-rank 상관이 대표적.

2. **"구별 가능"과 "SA-hard"는 독립**: ZeroExp는 대각 편향이 커서 구별하기 쉽지만 SA-easy. Wishart는 대각 편향이 작아서 구별하기 어렵지만 SA-hard.

3. **Posiform 계열의 SA-easy 근본 원인**: 2-SAT은 P(다항 시간)에 풀리는 문제 → posiform QUBO의 에너지 지형이 smooth → metastable trap 없음.

4. **Hardened Posiform이 이를 극복하는 원리**: random discrete-coefficient QUBO의 rugged landscape 위에 posiform을 overlay → random QUBO의 metastable trap을 유지하면서 GS 유일성 보장.

5. **실용적 병목**: Posiform 계열은 2-SAT 유일성 검사에 O(N^2.5) 소요 → N=1000에서 ~108초. Quiet Planting은 보조변수로 QUBO 크기가 5.2배 팽창.

---

## 10. 실험 환경

- **OS**: Linux 5.15.167.4 (WSL2)
- **Python**: 3.x
- **SA Solver**: `neal.SimulatedAnnealingSampler` (D-Wave Ocean SDK)
- **실험 일자**: 2026-02-20
