# GSP vs α 실험 결과 + HD CSV

**날짜**: 2026-05-19 21:23

## 실험 설정

| 항목 | 값 |
|---|---|
| Topology | Pegasus P16 (n=5627) |
| coeff_type | lin2 |
| α | [0, 1e-09, 1e-08, 1e-07, 1e-06, 1e-05, 0.0001, 0.001, 0.01, 0.1] |
| sweep | 1000 |
| num_reads | 100 |
| num_instances | 200 |
| 총 소요 시간 | 2333.5s |

**Main metric**: `gsp_gs`. **NEW**: `samples.csv`, `instances.csv`, HD 6-plot.

## 결과 요약 (GS / Bit / Energy / HD per-read)

| α | GS 성공률 | Bit 성공률 | Energy 성공률 | Avg HD |
|---|---|---|---|---|
| 0 | 11076/20000 (55.4%) | 0/20000 (0.0%) | 11076/20000 (55.4%) | 458.62 |
| 1e-09 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 69.70 |
| 1e-08 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 69.57 |
| 1e-07 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 69.37 |
| 1e-06 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 69.30 |
| 1e-05 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 69.36 |
| 0.0001 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 68.81 |
| 0.001 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 63.72 |
| 0.01 | 221/20000 (1.1%) | 221/20000 (1.1%) | 75/20000 (0.4%) | 25.48 |
| 0.1 | 19980/20000 (99.9%) | 19980/20000 (99.9%) | 0/20000 (0.0%) | 0.00 |

## 파일

- `gsp_vs_alpha.png/pdf` — 기존 GS 곡선 (2-panel)
- `samples.csv` — per-read raw (200,000 rows)
- `instances.csv` — per-instance summary (2,000 rows)
- `hd_A_violin.{png,pdf}` — HD 분포 violin per α
- `hd_B_cdf.{png,pdf}` — HD CDF per α
- `hd_C_heatmap.{png,pdf}` — instance × α heatmap (mean HD)
- `hd_D_trajectory.{png,pdf}` — α 따라가는 인스턴스 trajectory
- `hd_E_minmean.{png,pdf}` — per-inst min HD vs mean HD
- `hd_F_hexbin.{png,pdf}` — per-sample HD vs Δenergy hexbin
- `data.json` — params + α별 집계
