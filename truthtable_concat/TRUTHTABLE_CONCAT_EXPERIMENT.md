# Truth Table Concat QUBO — 실험 보고서

## 1. 배경 및 동기

### 기존 Truth Table QUBO의 한계

Truth Table QUBO는 진리표(2^n 비트스트링 → 에너지)로부터 임의의 에너지 landscape를 정확히 QUBO로 인코딩할 수 있는 강력한 방법이다. 그러나 진리표 크기가 2^n이므로:

| 모드 | 실용 한계 | 제약 원인 |
|------|-----------|-----------|
| Exact (Rosenberg) | n ≤ 9 | Möbius 변환 O(n·2^n) + 보조변수 폭발 |
| Approx (QP) | n ≤ 12 | 제약조건 수 2^n - 1 |

즉, 단일 블록으로는 **최대 12변수**까지만 생성할 수 있다.

### 아이디어: Block-Diagonal 접합 + Posiform Hardening

k-bit Truth Table QUBO를 h개 독립적으로 생성하여 block-diagonal로 접합:

```
Q_concat = diag(Q_1, Q_2, ..., Q_h)
```

선택적으로 posiform hardening을 적용하여 cross-block coupling 추가:

```
Q_final = Q_concat + α × Q_posiform
```

- 총 변수 수: k × h (approx/random) 또는 k × h + aux (exact)
- **α=0 (posiform 없음)**: 모든 블록이 동일 target → 전체 target = target 반복 h회
- **α>0 (posiform hardening)**: 블록별 랜덤 target → 반복 패턴 없음, posiform이 전체 ground state 보장
- 블록별 다른 seed → 다른 에너지 landscape
- **유일성 보장**: 각 블록의 ground state가 유일 + posiform이 target에서 유일 최소 → 전체도 유일
- **Posiform hardening**: α > 0이면 블록 간 결합이 생겨 block-diagonal 분해 불가능

---

## 2. 알고리즘

### 2.1 Landscape 모드

#### Planted (기본)
- 블록별 target을 심은 진리표 → QUBO
- E(target) = 0, E(x ≠ target) = uniform(0.1, 5.0)
- 에너지 bowl 구조 → SA-easy per block

#### Random
- 블록별 무작위 진리표 → unconstrained QP fit → QUBO
- 에너지가 0 중심 uniform(-5, 5)으로 모든 상태에 대해 랜덤
- Target은 brute force로 발견 (심지 않음)
- Genuine local minima (QP fit 잔차로 인한 복합 landscape)
- **Posiform hardening 필수** (블록이 target을 보장하지 않으므로)

### 2.2 블록 Target 결정

- **planted, α=0**: 모든 블록이 사용자 지정 target → `full_target = target * h` (반복)
- **planted, α>0**: 블록별 랜덤 k-bit target 생성 (seed+7777 기반 RNG)
- **random**: 블록별 QP fit 후 brute force GS 탐색 → target 자동 결정

### 2.3 블록 생성

**Planted**:
1. `preset_random_landscape(k, block_target_i, seed=seed+i)` → 진리표 생성
2. Approx: `create_qubo_approx(tt)` → k변수 QUBO (보조변수 0)
3. Exact: `create_qubo_truthtable(tt, reduce_strategy='greedy')` → (k + aux)변수 QUBO

**Random**:
1. `create_random_block_qubo(k, energy_range, seed=seed+i)` → QUBO + GS
2. Random truth table → unconstrained lstsq (상수항 포함) → QUBO
3. 보조변수 0개

### 2.4 Block-Diagonal 조립

블록 i의 모든 변수 인덱스에 offset을 더하여 전역 인덱스로 변환.

### 2.5 Posiform Hardening (planted: 선택, random: 필수)

α > 0이면 전체 target(k×h bits)으로 posiform planted QUBO를 생성하여 결합:
- `Q_final = Q_concat + α × Q_posiform`
- Random landscape: `targeted=True` (강한 유일성 필요)
- Planted landscape: `targeted=False` (블록이 이미 target 보장)

---

## 3. 실험 결과

### 3.0 실험 0: Ground State 검증 (Random Landscape)

Random landscape + posiform hardening에서 brute force로 GS 유일성 검증 (n ≤ 21).

**설정**: 20 instances per config

| n | k | h | α | GS 일치 | GS 유일 | avg gap | avg LM |
|---|---|---|---|---------|---------|---------|--------|
| 7 | 7 | 1 | 0.01 | **20/20** | **20/20** | 0.3453 | 3.5 |
| 7 | 7 | 1 | 0.1 | **20/20** | **20/20** | 0.6251 | 3.5 |
| 14 | 7 | 2 | 0.01 | **20/20** | **20/20** | 0.2000 | 6.9 |
| 14 | 7 | 2 | 0.1 | **20/20** | **20/20** | 0.5264 | 6.9 |
| 21 | 7 | 3 | 0.01 | **20/20** | **20/20** | 0.1768 | 10.4 |
| 21 | 7 | 3 | 0.1 | **20/20** | **20/20** | 0.5010 | 10.4 |

**관찰**: 모든 config에서 100% GS 일치 + 유일성. α 증가 시 energy gap 증가 (posiform 기여). 블록당 평균 ~3.5개 local minima (k=7 → 2^7=128 상태 중).

### 3.1 실험 1: h-Scaling (6 Configs)

블록 수(h) 증가에 따른 SA 성공률. k=7 고정, 6 configs.

**설정**: 10 instances, 100 reads, 1000 sweeps

| h | N | Approx | Approx(α=0.01) | Random(α=0.01) | Random(α=0.1) | Greedy | Greedy(α=0.01) |
|---|---|--------|----------------|----------------|---------------|--------|----------------|
| 1 | 7 | **63.00%** | **74.00%** | **87.60%** | **96.60%** | **2.50%** | **5.20%** |
| 2 | 14 | **41.00%** | **53.40%** | **73.90%** | **95.20%** | **0.20%** | **0.40%** |
| 5 | 35 | **21.50%** | **24.30%** | **55.90%** | **99.30%** | **0.00%** | **0.00%** |
| 10 | 70 | **14.60%** | **8.90%** | **44.90%** | **100.00%** | **0.00%** | **0.00%** |
| 20 | 140 | **6.90%** | **1.40%** | **17.60%** | **98.10%** | **0.00%** | **0.00%** |

**관찰**:

1. **Random(α=0.1) SA-easy**: h=10에서 100%, h=20에서도 98.1%. Posiform 신호가 강하면 random landscape의 얕은 local minima를 쉽게 극복.
2. **Random(α=0.01) > Planted(α=0.01)**: 모든 h에서 Random이 우수 (h=5: 55.9% vs 24.3%). Random landscape의 QP fit 잔차가 planted의 sharp bowl보다 SA-friendly.
3. **Planted Approx(α=0.01) 역전**: h≤5에서는 bare보다 우수, h≥10에서 역전 (h=10: 8.9% vs 14.6%). 랜덤 target 복잡도.
4. **Greedy 전멸**: h≥2에서 0%. 보조변수 ×3.1 팽창이 지배적.

### 3.2 실험 2: 14-way 비교

동일 변수 수(N=35)에서 Concat-Approx / Concat-Greedy / Random / Hardened Posiform 비교.

**설정**: k=7, h=5, 10 instances, 100 reads, 1000 sweeps

| Method | GS Rate | Count | QUBO | Time |
|--------|---------|-------|------|------|
| Concat-Greedy | **0.00%** | 0/1000 | 110 | 0.8s |
| Concat-Greedy(α=0.001) | **0.00%** | 0/1000 | 110 | 1.2s |
| Concat-Greedy(α=0.01) | **0.00%** | 0/1000 | 110 | 1.3s |
| Concat-Greedy(α=0.1) | **0.00%** | 0/1000 | 110 | 1.3s |
| Concat-Approx(α=0.001) | **10.90%** | 109/1000 | 35 | 13.4s |
| Concat-Approx | **20.90%** | 209/1000 | 35 | 13.1s |
| Concat-Approx(α=0.01) | **23.30%** | 233/1000 | 35 | 12.5s |
| Random(α=0.001) | **39.20%** | 392/1000 | 35 | 0.7s |
| Random(α=0.01) | **56.60%** | 566/1000 | 35 | 0.7s |
| Hardened(α=0.001) | **81.80%** | 818/1000 | 35 | 0.8s |
| Hardened(α=0.01) | **84.70%** | 847/1000 | 35 | 0.9s |
| Concat-Approx(α=0.1) | **98.10%** | 981/1000 | 35 | 13.3s |
| Random(α=0.1) | **98.70%** | 987/1000 | 35 | 0.7s |
| Hardened(α=0.1) | **99.70%** | 997/1000 | 35 | 0.9s |

**분석**:

1. **Concat-Greedy 전멸**: 모든 α에서 0%. 보조변수 비용 극복 불가.
2. **Planted < Random < Hardened** (동일 α): Random landscape는 planted보다 쉽지만 hardened보다 어려움.
   - α=0.01: Planted 23.3% < Random 56.6% < Hardened 84.7%
   - α=0.001: Planted 10.9% < Random 39.2% < Hardened 81.8%
3. **α=0.1 수렴**: 세 방법 모두 ~98-99%. 강한 posiform에서 base QUBO 차이 무시.
4. **생성 속도**: Random은 ~0.7s (planted ~13s). QP fit이 constrained QP보다 빠름.

---

## 4. 특성 분석

### Base QUBO 비교 (동일 α에서)

| Base QUBO | 생성 방법 | QUBO 크기 | α=0.01 GS Rate | α=0.1 GS Rate |
|-----------|-----------|-----------|----------------|---------------|
| Planted TT QP (랜덤 target) | 진리표 planted → QP 최적화 | N | 23.30% | 98.10% |
| Random TT QP | random 진리표 → unconstrained lstsq | N | 56.60% | 98.70% |
| Random discrete-coeff (Hardened) | 랜덤 {-1,+1} 계수 | N | 84.70% | 99.70% |
| Rosenberg (Greedy) | 진리표 → Möbius → Rosenberg | ~3.1×N | 0.00% | 0.00% |

### Landscape 특성 비교

| 특성 | Planted | Random |
|------|---------|--------|
| Target 결정 | 사용자 지정 / 랜덤 | Brute force 발견 |
| E(target) | 0 (고정) | Landscape에서 자연 최소 |
| Local minima 수 | 적음 (bowl) | ~3.5/block (k=7) |
| SA 난이도 | 중간 | 쉬움 (QP fit residual 얕음) |
| Posiform 필수? | 아니오 | 예 |
| 생성 속도 | 느림 (constrained QP) | 빠름 (lstsq) |

### 다른 방법론과의 포지셔닝

| 생성기 | SA-Hard? | 구별 불가능? | 보조변수 | QUBO 크기 제한 |
|--------|----------|-------------|---------|---------------|
| Concat-Planted (α=0) | 중간 | 아니오 (block-diagonal) | 0/많음 | 없음 |
| Concat-Planted (α>0) | 중간~쉬움 | 부분적 | 0/많음 | 없음 |
| Concat-Random (α>0) | 쉬움~중간 | 부분적 | 0 | 없음 |
| Hardened Posiform | SA-hard (N≥500) | 예 | 0 | 없음 |
| Wishart | SA-hard (α≥0.7) | 아니오 (low-rank) | 0 | 없음 |
| Quiet Planting | SA-hard (field 의존) | 예 (α<3.86) | m개 | 없음 |

---

## 5. 사용법

```bash
# Planted Approx 모드 (기본값, 보조변수 없음, 동일 target)
python3 truthtable_concat/qubo_truthtable_concat.py 1001111 10
python3 truthtable_concat/qubo_truthtable_concat.py 1001111 10 --seed 42

# Planted Exact 모드 (Rosenberg greedy, 동일 target)
python3 truthtable_concat/qubo_truthtable_concat.py 1001111 10 --exact

# Planted + Posiform hardening (블록별 랜덤 target 자동 적용)
python3 truthtable_concat/qubo_truthtable_concat.py 1001111 10 --harden 0.1

# Random landscape + Posiform hardening
python3 truthtable_concat/qubo_truthtable_concat.py 1001111 10 --random --harden 0.01
python3 truthtable_concat/qubo_truthtable_concat.py 1001111 10 --random --harden 0.01 --seed 42

# GS brute force 검증 (random landscape, n≤21)
python3 truthtable_concat/test_truthtable_concat.py --verify 20

# h-Scaling 실험 (6 configs)
python3 truthtable_concat/test_truthtable_concat.py --scaling 10

# 14-way 비교 실험
python3 truthtable_concat/test_truthtable_concat.py --compare 10
```

---

## 6. 결론

1. **Block-diagonal 접합**은 Truth Table QUBO의 크기 제한(n≤12)을 우회하는 실용적 방법이다.
2. **Approx 모드**가 Greedy(exact)보다 SA 탐색에 압도적으로 유리하다 (보조변수 0 vs ×3.1 팽창). Greedy는 모든 α에서 0% — 실용성 없음.
3. **Random landscape**: planted보다 SA-easy하지만 생성이 빠르고 genuine local minima를 제공. 동일 α에서 planted보다 항상 높은 성공률.
4. **난이도 순서 (동일 α)**: Concat-Planted < Concat-Random < Hardened Posiform (SA-easy 순).
5. **α=0.1에서 수렴**: 세 방법 모두 ~98-99%. 강한 posiform 신호에서 base QUBO 차이 소멸.
6. **SA-hardness 벤치마크에는 Hardened Posiform이 가장 적합**. Concat은 큰 N에서도 생성 가능하다는 장점이 있으나, 작은 α에서도 SA가 쉽게 풀 수 있음.
7. **Random landscape는 기존 truthtable_posiform 방법론을 통합**한 것으로, 동일 framework에서 planted/random 전환이 가능.
