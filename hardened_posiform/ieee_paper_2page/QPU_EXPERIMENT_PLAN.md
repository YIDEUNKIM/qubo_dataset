# QPU 실험 계획 — QPU 계수 정밀도 바닥 측정

**목표 (재설정)**: SA에서 본 난이도 곡선을 QPU로 *재현*하는 것이 아니라, **QPU가 planted 항 αP를 α 얼마부터 분해하는지 = QPU의 유효 계수 정밀도 바닥을 α-sweep으로 측정**한다. SA는 무한 정밀도 산수라 αP가 아무리 작아도 보지만, QPU는 아날로그 소자라 αP 커플러 값이 노이즈(ICE) 위로 올라와야 본다. 그 전이점 α\*가 QPU 정밀도의 직접 지표다. 부차적으로 분해 가능한 α(≳α\*)에서 SA와 같은 인스턴스를 비교한다.

## 왜 "재현"이 아니라 "정밀도 측정"인가 (실측 근거)

인스턴스 실측: R 계수는 전부 ±1(lin2), P 계수는 |값| 1~28(중앙값 2). QPU 결합계수 노이즈 ICE ≈ 정규화 J의 1~3%(약 0.02~0.03). αP 개별 커플러 크기 = α·|P|:

| α | median \|αP_ij\| | max \|αP_ij\| | ICE(~0.02–0.03) 대비 |
|---|---|---|---|
| 1e-3 | 0.002 | 0.028 | 한참 아래 (안 보임) |
| 3e-3 | 0.006 | 0.084 | 아래 (거의 안 보임) |
| 5e-3 | 0.010 | 0.14 | 경계 아래 |
| **1e-2** (SA peak) | **0.020** | 0.28 | **딱 ICE 바닥** |
| 2e-2 | 0.040 | 0.56 | 위 (분해됨) |

→ α ≲ 5e-3에서는 planted 신호가 QPU 노이즈에 묻혀 QPU가 사실상 ΣR만 본다. SA 곡선의 핵심(Δ_R 바닥 → 상승 → α=0.01 peak)이 대부분 QPU 정밀도 아래에 있어 **"재현"은 실패가 예정**돼 있다. (총 ΔP가 ~220으로 커도 소용없다 — QPU는 커플러 하나하나(각 ~0.002~0.02)에 노이즈와 함께 작용한다.) 그래서 이 한계를 *결과*로 바꾼다.

## 측정 대상 (산출물)

1. **분해 threshold α\***: QPU 난이도 곡선이 α=0 baseline(ΣR만 보는 행동)에서 통계적으로 벗어나기 시작하는 α.
2. **SA−QPU 격차**: SA는 ~1e-3부터 ΔE가 상승하지만 QPU는 α\*부터 — 그 차이가 "QPU 정밀도 비용".
3. 분해 가능 α(≳α\*)에서 **같은 인스턴스 per-instance 성공/순위 비교** (양자 우위 탐색).

## 설정 (SA와 동일)

| 항목 | 값 | 비고 |
|---|---|---|
| 인스턴스 | `instances/instances_pegasus16_lin2_hard300.pkl` (300개) | SA 실패빈도로 선정. ⚠ **QPU에선 안 갇힐 수 있음** → per-α QPU 실패율 필수 보고 |
| 솔버 | **Advantage_system6.4** 이름 고정 | 인스턴스가 박힌 칩 |
| α (10점) | 0, 1e-9, 1e-5, 1e-4, 1e-3, 3e-3, 5e-3, 1e-2, 1.5e-2, 2e-2 | SA와 동일 그리드(threshold 비교용) |
| num_reads | 100 (게이지 4 × 25) | |
| 지표 | 실패(x≠target) 조건부 ΔE·HD, **원래 Q_α로 비트열에서 재계산** | QPU 스케일·ICE 무관 |

## 프로토콜 — QPU 필수 추가 (검증의 핵심 3가지)

1. **고정 스케일 (auto_scale 끔)**: α마다 max\|계수\|가 달라(1 → α=0.02에서 ~1.56) auto_scale이 문제마다 다르게 재정규화하면 물리 에너지 스케일(유효 온도)이 α마다 변해 sweep이 오염된다. **전 α·전 인스턴스 공통 상수 1개로 스케일**해 αP만 변하게 한다.
2. **spin-reversal(게이지) 평균**: ICE·누설의 *계통* 편향을 상쇄(D-Wave 벤치마킹 표준). 없으면 노이즈가 비대칭 계통오차로 남아 ΔE/HD를 왜곡.
3. **QPU 실패율(trapped fraction) 보고**: 집단이 SA-기준이라 QPU에선 일부 풀릴 수 있다. α마다 QPU가 몇 %나 실패하는지 보고하고, 안 갇히면 wrong-only가 생존편향됨을 명시.

```python
import pickle
from dwave.system import DWaveSampler

data = pickle.load(open('instances/instances_pegasus16_lin2_hard300.pkl','rb'))
insts, meta = data['instances'], data['meta']

qpu = DWaveSampler(solver={'name': meta['qpu_solver']})   # 6.4 이름 고정
assert qpu.solver.name == meta['qpu_solver']

# 전 α 공통 고정 스케일 (auto_scale=False) — α마다 재정규화 방지
cmax = max(abs(v) for inst in insts for a in ALPHAS
           for v in build_Q(inst, a).values())
SCALE = 1.0 / cmax

for inst in insts:
    for alpha in ALPHAS:
        Q = build_Q(inst, alpha)                 # 물리 큐빗 인덱스 키
        # active-set 검증: Q의 노드/엣지가 현재 칩에 다 있는지 (없으면 즉시 에러)
        Qs = {k: SCALE * v for k, v in Q.items()}
        resp = qpu.sample_qubo(Qs, num_reads=100, annealing_time=20,
                               auto_scale=False, num_spin_reversal_transforms=4,
                               answer_mode='raw')
        t_energy = inst['t_energy_r'] + alpha * inst['t_energy_p']
        for s, _e, occ in resp.data(['sample','energy','num_occurrences']):
            x = ''.join(str(s[nd]) for nd in inst['sorted_nodes'])
            hd = sum(a != b for a, b in zip(inst['target_str'], x))
            dE = qubo_energy(x, Q) - t_energy     # ⚠ 원래 Q_α 로 재계산 (스케일 안 함)
```

## 시간 예산 (6.4 타이밍 기준)

per-problem access ≈ 15.93 ms(programming) + 100 × (20 + 173.3 + 20.58) µs ≈ **37.3 ms**. 게이지 평균은 reads를 나눠 쓰므로 총량 동일. 300 inst × 10 α = 3000 problems → **≈ 1.9 분 access**(과금), wall-clock ≈ 19 분. 단일 세션.

## 분석

- QPU 곡선을 α=0 baseline과 비교 → **α\* 추출** (baseline에서 유의하게 벗어나는 첫 α).
- SA와 같은 그리드에 overlay하되 **비교는 (a) lift-off 위치(α\*), (b) 분해 α에서 per-instance 성공/순위**만. **절대 ΔE overlay는 하지 말 것** — T=20µs와 1000 sweeps는 대응이 아니라 절대값이 솔버·스케줄 의존이다(사과-오렌지).

## 한계 / 주의

- **T=20µs ↔ 1000 sweeps 대응 아님**: 정성적 모양·순위만 비교.
- **dual-budget(20 vs 200µs)**: QPU는 어닐링 길수록 좋아지지 않을 수 있음(디코히런스·열화로 비단조) → "예산↑→ΔE↓" 기대 금지. 보조 실험으로만.
- **집단 SA-선정**: QPU 실패율 보고 필수. 필요하면 QPU 실패빈도로 재선정한 별도 집단도 분석.
- **hardware-native 검증**: posiform 절이 KL 블록을 가로지르는 커플러를 쓰므로, Q_α의 모든 노드/엣지가 현재 active P16 set에 있는지 active-set 검증(없으면 silent 0 금지, 즉시 RuntimeError).

## 기대 (가설)

- QPU valley(α ≲ 5e-3) ≈ α=0 평평 (정밀도 바닥에 묻힘).
- α\* ~ 0.005–0.01 부근에서 lift-off.
- SA의 peak는 QPU에서 흐릿하거나 오른쪽으로 밀림.
- **α\* 자체가 핵심 결과**: "이 QPU·이 P-스케일에서 planted 신호를 분해하는 한계는 α ≈ X" — SA(무한 정밀도) 대비 아날로그 QPU의 정밀도 비용을 정량화.

## 실행

`hardened_posiform/hardware_native_qubo.ipynb` QPU 셀에서 인스턴스 pkl·solver·α·고정스케일·게이지를 위 조건으로 교체. QPU 샘플을 SA와 같은 samples.csv 스키마(alpha, instance_id, hd, delta_energy, is_target)로 저장 → `analyze_hard300.py` 로직 적용해 곡선 생성·overlay.
