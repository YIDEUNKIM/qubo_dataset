# GSP vs α 실험 결과 + HD CSV

**날짜**: 2026-05-20 18:14

## 실험 설정

| 항목 | 값 |
|---|---|
| Topology | Pegasus P16 (n=5612) |
| coeff_type | lin20 |
| α | [0, 1e-09, 1e-08, 1e-07, 1e-06, 1e-05, 0.0001, 0.001, 0.01, 0.1] |
| sweep | 1000 |
| num_reads | 100 |
| num_instances | 200 |
| 총 소요 시간 | 1366.2s |

**Main metric**: `gsp_gs`. **NEW**: `samples.csv`, `instances.csv`, HD 6-plot.

## 결과 요약 (GS / Bit / Energy / HD per-read)

| α | GS 성공률 | Bit 성공률 | Energy 성공률 | Avg HD |
|---|---|---|---|---|
| 0 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 188.05 |
| 1e-09 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 137.48 |
| 1e-08 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 136.23 |
| 1e-07 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 134.44 |
| 1e-06 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 133.64 |
| 1e-05 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 132.04 |
| 0.0001 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 127.61 |
| 0.001 | 0/20000 (0.0%) | 0/20000 (0.0%) | 0/20000 (0.0%) | 97.29 |
| 0.01 | 2991/20000 (15.0%) | 2991/20000 (15.0%) | 2991/20000 (15.0%) | 9.68 |
| 0.1 | 19970/20000 (99.8%) | 19970/20000 (99.8%) | 19970/20000 (99.8%) | 0.00 |

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
