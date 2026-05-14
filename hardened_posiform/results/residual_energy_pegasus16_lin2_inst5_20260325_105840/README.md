# Residual Energy + 에너지 갭 실험 결과

**날짜**: 2026-03-25 11:14

## 실험 설정

| 항목 | 값 |
|---|---|
| Topology | Pegasus P16 (n=5640) |
| coeff_type | lin2 |
| α | [0, 0.001, 0.01, 0.1] |
| sweeps | [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000] |
| num_instances | 5 |
| 총 소요 시간 | 910.2s |

## 에너지 갭 통계

| 지표 | Mean | Median | Min | Max |
|---|---|---|---|---|
| Δ_P | -32879.200 | -28776.000 | -46396.000 | -27568.000 |
| Δ_R | 1.000 | 1.000 | 1.000 | 1.000 |

## 장벽-갭 비율 ρ = Δ_R / (α·Δ_P)

| α | ρ mean | ρ median |
|---|---|---|
| 0.001 | -0 | -0 |
| 0.01 | -0 | -0 |
| 0.1 | -0 | -0 |

## 파일 목록

- `energy_gsp.png` — Energy GSP vs Sweeps
- `residual_energy.png` — Residual Energy vs Sweeps
- `delta_distribution.png` — Δ_P, Δ_R 분포
- `data.json` — 전체 데이터
