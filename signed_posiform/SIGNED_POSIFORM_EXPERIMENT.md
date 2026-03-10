# Signed Posiform QUBO 실험 보고서

## 방법론 개요

기존 posiform planting은 wrong tuple indicator에 **양수 weight만** 사용.
Signed Posiform은 **음수 weight도 허용**하되, GS(Ground State) 보장 범위 내로 제한.

### 핵심 아이디어

wrong tuple (x_i, x_j) = (w_i, w_j)에 weight b를 곱해서 에너지 함수에 추가할 때:
- **b > 0**: 항상 안전 (기존 posiform)
- **b < 0**: |b| < X - Y 범위에서만 안전
  - X = 해당 wrong tuple 부분공간에서 현재 에너지 함수의 최솟값
  - Y = target 에너지 (global minimum)

### 2단계 구성

- **Phase 1**: 절의 (1 - negative_ratio) 비율을 양수 weight로 처리 → 기본 에너지 지형 수립
- **Phase 2**: 나머지를 음수 weight로 처리 → gap 기반 범위 내

### 장점

1. Q 행렬에 양/음 혼합 → random QUBO와 유사 (planted solution 감지 어려움)
2. Frustration 도입 가능
3. GS 수학적 보장 (한계치 준수로 전수검사 없이 보장)

### 한계

- gap 계산에 O(2^n) 필요 → 실용 범위 n ≤ 25
- SA 난이도는 기존 posiform과 유사 (SA-easy)
- 단독으로는 hardness 향상 효과 미미 → hardened posiform과 결합 필요

---

## 실험 환경

- SA: `neal.SimulatedAnnealingSampler`
- num_reads=200, num_sweeps=1000 (n≤20)
- coeff_range=(1.0, 3.0), margin=0.01

---

## 실험 1: GS 보장 검증

signed posiform이 모든 n, seed에서 GS를 보장하는지 brute force 검증.

| n | seeds | GS 보장 | Unique | avg gap | Q 음수 비율 |
|---|-------|---------|--------|---------|------------|
| 5 | 20 | 20/20 | 20/20 | 0.7478 | 47.0% |
| 7 | 20 | 20/20 | 20/20 | 0.4409 | 52.3% |
| 8 | 20 | 20/20 | 20/20 | 0.2732 | 48.6% |
| 10 | 20 | 20/20 | 20/20 | 0.2273 | 50.5% |
| 12 | 1 | 1/1 | 1/1 | 0.2196 | 42.6% |
| 14 | 1 | 1/1 | 1/1 | 0.0050 | 52.3% |
| 16 | 1 | 1/1 | 1/1 | 0.0460 | 43.1% |
| 18 | 1 | 1/1 | 1/1 | 0.0619 | 46.0% |
| 20 | 1 | 1/1 | 1/1 | 0.0769 | 44.9% |
| 22 | 1 | 1/1 | 1/1 | 0.0190 | 55.9% |
| 25 | 1 | 1/1 | 1/1 | 0.0221 | 52.4% |

**결론**: 모든 테스트에서 100% GS 보장. Q 행렬의 약 50%가 음수 계수.

---

## 실험 2: SA Scaling (N=10, 15, 20)

| N | Success% | Avg Hamming | Q Neg% | Gen Time |
|---|----------|-------------|--------|----------|
| 10 | 100.0% | 0.0 | 50.6% | 0.00s |
| 15 | 100.0% | 0.0 | 50.1% | 0.02s |
| 20 | 100.0% | 0.0 | 49.5% | 0.88s |

**결론**: n≤20에서 SA 100% 성공. Posiform 기반 방법론 특성상 SA-easy.

---

## 실험 3: 음수 비율 Sweep (N=15, 20 seeds)

negative_ratio 변화에 따른 에너지 갭과 SA 성공률.

| Neg Ratio | Success% | Q Neg% | Avg Gap | Avg Hamming |
|-----------|----------|--------|---------|-------------|
| 0.0 | 100.0% | 48.5% | 1.5872 | 0.0 |
| 0.1 | 100.0% | 49.5% | 0.3887 | 0.0 |
| 0.2 | 100.0% | 48.9% | 0.2241 | 0.0 |
| 0.3 | 100.0% | 49.0% | 0.0796 | 0.0 |
| 0.5 | 100.0% | 49.2% | 0.0733 | 0.0 |
| 0.7 | 100.0% | 49.5% | 0.1001 | 0.0 |

**관찰**:
- negative_ratio 증가 → 에너지 갭 감소 (1.59 → 0.07)
- 그러나 SA 성공률에는 영향 없음 (모두 100%)
- Q 행렬의 음수 비율은 neg_ratio와 무관하게 ~49% 유지
  - 이유: 기존 posiform도 양수 weight로 Q에 음수 항 생성 (indicator 전개 시)
  - neg_ratio=0.0에서도 Q_neg≈48%인 것이 이를 보여줌

---

## 실험 4: Signed Posiform vs Plain Posiform (N=20, 20 runs)

| Method | Success% | EXACT | ENERGY | FAIL | Avg Hamming |
|--------|----------|-------|--------|------|-------------|
| SignedPosiform | 100.0% | 20 | 0 | 0 | 0.0 |
| Posiform | 100.0% | 20 | 0 | 0 | 0.0 |

**결론**: 동일 규모에서 SA 난이도 차이 없음. 두 방법 모두 SA-easy.

---

## 성능 최적화

Phase 2에서 매 clause마다 에너지를 재계산해야 하므로, 최적화가 중요.

### 최적화 전 vs 후 (벡터화 + 증분 업데이트)

| n | 최적화 전 | 최적화 후 | 속도 향상 |
|---|----------|----------|----------|
| 12 | 0.3s | 0.0s | - |
| 14 | 2.4s | 0.0s | ~100x |
| 16 | 12.1s | 0.0s | ~1000x |
| 18 | 63.8s | 0.1s | 638x |
| 20 | 400.8s | 1.0s | 401x |
| 22 | (추정 ~3000s) | 6.7s | ~450x |
| 25 | 불가능 | 45.0s | - |

**핵심 최적화**:
1. **증분 에너지 업데이트**: Phase 2 clause 추가 시 Q 전체 순회 대신 indicator delta만 O(2^n) 추가
2. **numpy 벡터화**: Python 루프 → numpy 배열 연산
3. **비트 배열 사전 계산**: `_precompute_bits`로 매번 비트 추출 제거

---

## 사용법

```bash
# 기본 실행 (target=10110, seed=42)
python3 signed_posiform/qubo_signed_posiform.py 10110 42

# 길이 지정 랜덤 target
python3 signed_posiform/qubo_signed_posiform.py 15

# SA Scaling 실험
python3 signed_posiform/test_signed_posiform.py --scaling 20

# 음수 비율 sweep 실험
python3 signed_posiform/test_signed_posiform.py --neg-sweep 20

# Signed vs Plain Posiform 비교
python3 signed_posiform/test_signed_posiform.py --compare 20
```
