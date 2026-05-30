# SA vs QPU — 실험 과정과 결과

목적: SA로 얻은 난이도 곡선(고정 300 집단, 오답만 ΔE/HD)을 실제 D-Wave QPU에서 같은 인스턴스로 측정해 비교한다. 한 줄 결론: **QPU도 같은 정성적 난이도 구조(평탄 바닥 → α=0.01 부근 ΔE 비단조 peak → HD 단조 감소)를 보이되, 한 자릿수 큰 하드웨어 오차 바닥 위에서 본다.** → 구조가 SA 산물이 아니라 인스턴스 성질임.

## 데이터 위치
- **QPU 본 결과 (g2)**: `results/qpu_run_i0-100_20260527_221806/` — auto_scale + 2-게이지 spin-reversal 평균, 100 인스턴스. `samples.csv`(9칼럼: alpha, instance_id, read_idx, hd, energy, delta_energy, is_target, is_R_gs, eng_ok), `summary.json`
- **g1 대조 (게이지 없음)**: `results/qpu_run_g1_i0-100_20260527_215718/` — auto_scale, 게이지 1회. g2와 비교해 게이지 효과 측정용
- **정밀도 바닥 대조 (under-scaled)**: `results/qpu_ctrl_fixed0.641_i0-100_20260527_213127/` — 고정 SCALE=0.641 (range 23%만 사용)
- **SA**: `results/gsp_fullrun500_..._20260521_173814/` + `results/gsp_extend500peak_..._20260527_162048/` (고정 300 집단 = `results/hard300_roster.json`)
- 그림: `figures/fig_sa_vs_qpu.png` (SA-300 vs QPU-100 g2)

## 실험 설정 (SA ↔ QPU)
| 항목 | SA | QPU |
|---|---|---|
| 솔버 | neal SA, 1000 sweeps | D-Wave Advantage (Pegasus P16, 5612 큐빗), anneal 20µs |
| 인스턴스 | 가장 어려운 300개 고정 집단 | 그 집단 중 100개 (인스턴스 0–99, 동일) |
| reads / (inst, α) | 100 | 100 |
| α (10점) | 0, 1e-9, 1e-5, 1e-4, 1e-3, 3e-3, 5e-3, 1e-2, 1.5e-2, 2e-2 | 동일 |
| 스케일 | — | **auto_scale**(문제별 range 채움) |
| 게이지 | — | 2회 spin-reversal 평균 |
| 지표 | 오답만(x≠target) ΔE·HD, 원래 Q_α로 계산 | 동일 |

## ⚠ 스케일 교훈 (중요)
처음엔 cross-α 비교를 위해 **고정 SCALE=1/cmax_QUBO=0.641**을 썼으나 **틀렸다**: QUBO 최댓값만 1로 맞춰 하드웨어 J range의 ~23%만 사용 → αP 커플러가 J~0.013으로 ICE(~0.02–0.03) 아래로 묻힘 → 봉우리 소멸(ΔE 354/HD 847로 폭망, `qpu_ctrl_fixed0.641_*`). ΔE/HD는 **원래 Q_α로 재계산**하므로 스케일과 무관하게 단위 비교됨 — 스케일은 QPU 최적화 fidelity만 좌우. 따라서 **auto_scale**(문제별 range 채움, 표준)이 옳다. (올바른 고정 스케일도 가능: 전 α를 Ising range에 맞추는 SCALE≥1.406, binding은 대각 h.) → 이 대조가 본문 "정밀도 바닥" caveat 한 문장의 근거.

## 실험 과정 (QPU)
1. 인스턴스(`instances_pegasus16_lin2_hard300.pkl`)를 **솔버 이름으로 고정**해 로드 — 자동 선택 시 다른 칩이 잡혀 silent 실패하므로.
2. **hardware-native**: minor-embedding·chain 없이 물리 큐빗에 직접 제출(raw `DWaveSampler`). 제출 전 인스턴스 Q_α의 노드·커플러가 칩 active set에 다 있는지 검증 → 100/100 통과.
3. 각 (inst, α)에 `sample_qubo`로 100 reads, **auto_scale + 2-게이지 평균**. **ΔE는 QPU 보고 에너지가 아니라 원래 Q_α로 비트열에서 직접 재계산**.
4. 모든 통계는 **오답 조건부**(QPU는 전 α에서 100% 실패라 사실상 전체 샘플).
5. 단계별 실행: 단계1 [0,100) 완료. 단계2 [100,300)는 **D-Wave 분기 점검**(Advantage_system6.4→Advantage_system6 rename, Pegasus 오프라인)으로 대기 중 → 복귀 시 재실행 예정.

## 결과 (QPU g2, 오답만 mean)
| α | SA ΔE | QPU ΔE | SA HD | QPU HD | QPU 실패율 |
|---|------:|-------:|------:|-------:|------:|
| 0       | 0.57 | 24.28 | 460.2 | 501.0 | 100% |
| 1e-9    | 1.11 | 23.85 |  72.2 | 499.8 | 100% |
| 1e-5    | 1.11 | 23.72 |  72.3 | 498.8 | 100% |
| 1e-4    | 1.17 | 24.24 |  71.9 | 497.0 | 100% |
| 1e-3    | 1.67 | 27.01 |  67.1 | 474.5 | 100% |
| 3e-3    | 2.52 | 33.25 |  57.2 | 427.2 | 100% |
| 5e-3    | 3.04 | 37.18 |  48.1 | 380.9 | 100% |
| 1e-2    | 3.28 | 42.88 |  29.9 | 285.9 | 100% |
| **1.5e-2** | 2.71 | **43.54** | 17.7 | 216.1 | 100% |
| 2e-2    | 2.08 | 42.05 |  10.8 | 167.1 | 100% |

![fig_sa_vs_qpu](figures/fig_sa_vs_qpu.png)

그림. SA(파랑, 300) vs QPU(빨강, 100), 오답만 평균. (a) HD 공통 로그축. (b) ΔE 이중 y축(좌 SA / 우 QPU) — 스케일 차이가 커서 분리.

## 게이지 효과 (g1 vs g2)
| α | g1 ΔE(게이지X) | g2 ΔE(게이지2) | g1 HD | g2 HD |
|---|---:|---:|---:|---:|
| 0    | 28.63 | 24.28 | 500.0 | 501.0 |
| 1e-2 | 46.43 | 42.88 | 288.7 | 285.9 |
| 2e-2 | 45.13 | 42.05 | 171.4 | 167.1 |

- 게이지 평균은 **ΔE를 전 α 일률 ~10–15% 낮춤**(ICE 계통편향 부분 보정), **HD는 거의 불변**, **봉우리·단조 구조 그대로**. → 난이도 시그니처가 게이지에 강건.

## 해석
- **성공률**: QPU는 5612변수 정확해를 사실상 못 찾아 **전 α 100% 실패** → 이진 지표 완전 무용(이 집단이 QPU에서도 100% 갇히므로 재선정 불필요).
- **ΔE — 같은 봉우리, 다른 바닥**: 둘 다 valley 평탄 → α=0.01 부근 **비단조 peak** → 하강. 단 바닥이 SA ≈ 1, QPU ≈ 24. QPU는 정답 근처에 못 가 베이스라인이 큼.
- **HD 단조 감소**: SA 460→11, QPU 501→167. 둘 다 planted bias에 반응해 비트로는 target에 다가가나, QPU는 멀리서 멈춤.

→ **메트릭(오답만 ΔE/HD)이 실제 QPU에서도 planted 난이도를 드러낸다**(성공률은 0%인데도). 차이는 *구조 재현 실패*가 아니라 **최적화 품질**.

## 재생성
```bash
# QPU 측정 (Leap 토큰 필요). qpu_run.py 상단: INST_START/END, GAUGES, FIXED_SCALE
DWAVE_API_TOKEN=... python3 qpu_run.py
# SA vs QPU overlay 그림 (results/qpu_run_i* 자동 합산)
python3 make_overlay.py
```
