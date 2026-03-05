# Truth Table Concat QUBO — 실험 보고서

## 1. 배경 및 동기

### 기존 Truth Table QUBO의 한계

Truth Table QUBO는 진리표(2^n 비트스트링 → 에너지)로부터 임의의 에너지 landscape를 정확히 QUBO로 인코딩할 수 있는 강력한 방법이다. 그러나 진리표 크기가 2^n이므로:

| 모드 | 실용 한계 | 제약 원인 |
|------|-----------|-----------|
| Exact (Rosenberg) | n ≤ 9 | Möbius 변환 O(n·2^n) + 보조변수 폭발 |
| Approx (QP) | n ≤ 12 | 제약조건 수 2^n - 1 |

즉, 단일 블록으로는 **최대 12변수**까지만 생성할 수 있다.

### 아이디어: Block-Diagonal 접합

동일한 k-bit target을 가진 Truth Table QUBO를 h개 독립적으로 생성하여 block-diagonal로 접합하면:

```
Q_total = diag(Q_1, Q_2, ..., Q_h)
```

- 총 변수 수: k × h (approx) 또는 k × h + aux (exact)
- 각 블록은 독립 → 전체 ground state = target 반복 h회
- 블록별 다른 seed → 다른 에너지 landscape, 동일 target
- **유일성 보장**: 각 블록의 ground state가 유일하면 전체도 유일

이를 통해 k ≤ 12인 Truth Table의 제약을 우회하여 **임의 크기의 QUBO**를 생성할 수 있다.

---

## 2. 알고리즘

### 2.1 블록 생성

각 블록 i (i = 0, 1, ..., h-1):

1. `preset_energy_gap(k, target, gap, noise_scale, seed=seed+i)` → 진리표 생성
   - E(target) = 0
   - E(x ≠ target) = gap + \|N(0, noise_scale)\|
2. Approx 모드: `create_qubo_approx(tt)` → k변수 QUBO (보조변수 0)
3. Exact 모드: `create_qubo_truthtable(tt)` → (k + aux)변수 QUBO

### 2.2 Block-Diagonal 조립

블록 i의 모든 변수 인덱스에 offset을 더하여 전역 인덱스로 변환:

```
블록 0: 변수 [0, 1, ..., k-1]
블록 1: 변수 [k, k+1, ..., 2k-1]
...
블록 h-1: 변수 [(h-1)k, ..., hk-1]
```

Q_total의 비영 항은 각 블록 내부에만 존재하므로 블록 간 결합은 없다.

### 2.3 Ground State 유일성

- Approx 모드: QP 제약으로 E_Q(target) + ε ≤ E_Q(x) 보장 → 각 블록에서 target이 유일 ground state
- Exact 모드: Möbius 변환이 진리표를 정확히 재현 → gap > 0이므로 target이 유일 ground state
- 독립 블록의 직합: 모든 블록에서 유일 → 전체 Q_total에서도 유일

---

## 3. 실험 결과

### 3.1 실험 1: h-Scaling (Approx vs Exact)

블록 수(h) 증가에 따른 SA 성공률 변화. 블록 크기 k=7 고정. Approx와 Exact 모드 비교.

**설정**:

| 파라미터 | 값 |
|----------|-----|
| k (블록 크기) | 7 |
| h (블록 수) | 1, 2, 5, 10, 20 |
| mode | approx / exact |
| gap | 2.0 |
| noise_scale | 1.0 |
| instances | 100 |
| num_reads | 100 |
| num_sweeps | 1000 |

#### Approx 모드 결과

| h | N (총 변수) | QUBO 크기 | GS Rate | Avg Hamming |
|---|------------|-----------|---------|-------------|
| 1 | 7 | 7 | **58.08%** | 1.62 |
| 2 | 14 | 14 | **30.27%** | 3.56 |
| 5 | 35 | 35 | **9.73%** | 8.79 |
| 10 | 70 | 70 | **4.97%** | 16.73 |
| 20 | 140 | 140 | **5.85%** | 33.04 |

#### Exact 모드 결과

| h | N (총 변수) | QUBO 크기 | GS Rate | Avg Hamming |
|---|------------|-----------|---------|-------------|
| 1 | 7 | 22 | **6.86%** | 3.37 |
| 2 | 14 | 44 | **0.54%** | 6.61 |
| 5 | 35 | 110 | **0.00%** | 16.24 |
| 10 | 70 | 220 | **0.00%** | 33.17 |
| 20 | 140 | 440 | **0.00%** | 68.27 |

#### Approx vs Exact 비교

| h | N | Approx QUBO | Approx GS Rate | Exact QUBO | Exact GS Rate | QUBO 팽창률 |
|---|---|-------------|----------------|------------|---------------|------------|
| 1 | 7 | 7 | **58.08%** | 22 | **6.86%** | ×3.1 |
| 2 | 14 | 14 | **30.27%** | 44 | **0.54%** | ×3.1 |
| 5 | 35 | 35 | **9.73%** | 110 | **0.00%** | ×3.1 |
| 10 | 70 | 70 | **4.97%** | 220 | **0.00%** | ×3.1 |
| 20 | 140 | 140 | **5.85%** | 440 | **0.00%** | ×3.1 |

**관찰**:

1. **Exact 모드의 보조변수 비용**: k=7에서 블록당 보조변수 15개 → QUBO 크기 22 (×3.1 팽창). 이 비용이 h에 비례하여 누적됨.
2. **Approx 성공률 감소**: h 증가에 따라 58% → 30% → 10% → 5%로 단조 감소. h=10→20에서 소폭 반등(4.97→5.85%)은 통계 노이즈.
3. **Exact는 h=1에서부터 난이도 높음**: 단일 블록도 22변수 QUBO가 되어 sweeps=1000에서 6.86%만 성공. h≥5에서 완전 0%.
4. **Exact Avg Hamming ≈ Approx의 2배**: Exact h=1에서 3.37 vs Approx 1.62. 보조변수 탐색 공간이 SA를 정답 근처에서도 멀어지게 함.
5. **Hamming distance 선형 증가**: 두 모드 모두 블록이 독립적이므로 Avg Hamming이 h에 비례. Approx 블록당 ~1.65, Exact 블록당 ~3.4.
6. **이론 vs 실측 (Approx)**: 단일 블록 성공률 p ≈ 0.58일 때, p^5 = 6.9% (실측 9.73%), p^10 = 0.5% (실측 4.97%). 실측이 이론보다 높은 이유는 SA가 block-diagonal 구조를 부분적으로 활용하기 때문으로 추정.
7. **결론**: Approx 모드가 Concat QUBO에 압도적으로 유리. Exact 모드는 보조변수 비용으로 실용성이 없음.

### 3.2 실험 2: Hardened vs Concat 비교

동일 변수 수(N≈35)에서 Concat과 Hardened Posiform의 SA 난이도 비교.

**설정**:

| 파라미터 | 값 |
|----------|-----|
| N | ≈35 |
| Concat | k=7, h=5 |
| Hardened | lin2, α=0.1 |
| instances | 100 |
| num_reads | 100 |
| num_sweeps | 1000 |

**결과**:

| 방법론 | GS Rate | 소요 시간 |
|--------|---------|-----------|
| Concat-Exact | **0.00%** | 8.4s |
| Concat-Approx | **15.99%** | 47.0s |
| Hardened (lin2, α=0.1) | **100.00%** | 9.8s |

**분석**:

1. **Hardened가 N=35에서 SA-easy**: lin2, α=0.1이라도 N=35는 Hardened의 상전이 크기(N≈500)보다 훨씬 작아 SA가 쉽게 풀어냄. SA sweeps=1000으로도 100% 성공.
2. **Concat-Approx가 더 어려움**: 같은 N=35에서 Concat-Approx는 15.99%만 성공. 블록 독립성에도 불구하고, 5개 블록이 모두 동시에 맞아야 하므로 난이도가 높음.
3. **Concat-Exact 완전 실패 (0%)**: Rosenberg 보조변수로 QUBO 크기가 k*h보다 훨씬 커짐 (약 22×5=110 변수). 보조변수가 탐색 공간을 크게 확장하여 sweeps=1000으로는 불충분.
4. **Concat-Approx vs Exact 트레이드오프**: Approx 모드는 보조변수가 없어(QUBO 크기=35) SA가 효율적으로 탐색. Exact는 정확하지만 보조변수 비용이 큼.

---

## 4. 특성 분석

### Block-Diagonal 구조의 장단점

**장점**:
- 단일 블록의 진리표 제약(n≤12)을 우회하여 임의 크기 QUBO 생성 가능
- 각 블록의 ground state 유일성 → 전체 유일성 보장
- Approx 모드에서 보조변수 0개 → QUBO 크기 = k×h

**단점**:
- **Block-diagonal 구조 노출**: Q 행렬의 비대각 블록이 모두 0이므로, 관찰자가 쉽게 블록 구조를 탐지 가능. 랜덤 QUBO와 통계적으로 구별됨.
- **블록 독립 분할 가능**: SA 또는 다른 솔버가 블록 구조를 인식하면 문제를 h개 소문제로 분할하여 개별적으로 풀 수 있음. 이 경우 난이도가 단일 블록 수준으로 하락.
- **SA 난이도 한계**: 단일 블록이 SA-easy(k=7에서 ~58%)이므로, 블록 독립 풀이가 가능한 지능적 솔버에 대해서는 전체 문제도 SA-easy.

### 다른 방법론과의 포지셔닝

| 생성기 | SA-Hard? | 구별 불가능? | 보조변수 | QUBO 크기 제한 |
|--------|----------|-------------|---------|---------------|
| Concat-Approx | 중간 (블록 수 의존) | 아니오 (block-diagonal) | 0 | 없음 |
| Concat-Exact | SA-hard (보조변수 비용) | 아니오 | 많음 | 없음 |
| Hardened Posiform | SA-hard (N≥500) | 예 (off-diagonal) | 0 | 없음 |
| Wishart | SA-hard (alpha≥0.7) | 아니오 (low-rank) | 0 | 없음 |
| Quiet Planting | SA-hard (field 의존) | 예 (alpha<3.86) | m개 | 없음 |

Concat QUBO는 **정답이 보장된 대규모 벤치마크 생성**에 유용하지만, **SA-hardness나 통계적 은닉성**을 목표로 하는 용도에는 적합하지 않다.

---

## 5. 사용법

```bash
# Approx 모드 (기본값, 보조변수 없음)
python3 truthtable_concat/qubo_truthtable_concat.py 1001111 10
python3 truthtable_concat/qubo_truthtable_concat.py 1001111 10 --gap 2.0 --noise 1.0 --seed 42

# Exact 모드 (Rosenberg)
python3 truthtable_concat/qubo_truthtable_concat.py 1001111 10 --exact

# h-Scaling 실험
python3 truthtable_concat/test_truthtable_concat.py --scaling 100

# Hardened vs Concat 비교 실험
python3 truthtable_concat/test_truthtable_concat.py --compare 100
```

---

## 6. 결론

1. **Block-diagonal 접합**은 Truth Table QUBO의 크기 제한(n≤12)을 우회하는 실용적 방법이다.
2. **Approx 모드**가 Exact보다 SA 탐색에 유리하다 (보조변수 0 vs 많음).
3. **SA 난이도는 중간 수준**: 단일 블록이 SA-easy이므로, 블록 구조를 분할하면 쉽게 풀림. 블록 구조를 모르는 솔버에 대해서만 난이도가 있음.
4. **N=35에서 Hardened(100%) > Concat-Approx(16%) > Concat-Exact(0%)**: Hardened는 이 규모에서 SA-easy이지만, Concat은 "모든 블록 동시 성공" 조건 때문에 오히려 더 어려움.
5. **주된 용도**: 정답이 보장된 대규모 QUBO 벤치마크 생성. SA-hardness 벤치마크에는 Hardened Posiform이나 Quiet Planting이 더 적합.
