# QPU GSP vs α × anneal time 실험 결과

**날짜**: 2026-04-30 20:39

## 실험 설정

| 항목 | 값 |
|---|---|
| QPU | Advantage_system4.1 |
| Topology | Pegasus P16 (5627 active qubits) |
| coeff_type | lin2 |
| α | [0, 0.001, 0.01, 0.05, 0.1] |
| anneal_times (μs) | [0.5, 2000.0] |
| num_reads | 100 |
| num_instances | 0 (of 100) |
| 총 소요 시간 | 0.0s |

**Main metric**: `gsp_gs` (α=0 R-GS / α>0 target exact)

## 저장된 데이터

- `data.json` — 집계 + per-instance metrics + timing summary
- `raw/raw_alpha{α}_T{T}.pkl` — combo별 raw bits/energies/occurrences/timing (post-hoc 분석용)
- `qpu_gsp_vs_alpha_anneal.png/pdf` — main figure

## 결과 (GS / Bit / Energy / HD)

| α | T(μs) | GS | Bit | Energy | Avg HD |
|---|---|---|---|---|---|
| 0 | 0.5 | 0.000 | 0.000 | 0.000 | 0.0 |
| 0 | 2000.0 | 0.000 | 0.000 | 0.000 | 0.0 |
| 0.001 | 0.5 | 0.000 | 0.000 | 0.000 | 0.0 |
| 0.001 | 2000.0 | 0.000 | 0.000 | 0.000 | 0.0 |
| 0.01 | 0.5 | 0.000 | 0.000 | 0.000 | 0.0 |
| 0.01 | 2000.0 | 0.000 | 0.000 | 0.000 | 0.0 |
| 0.05 | 0.5 | 0.000 | 0.000 | 0.000 | 0.0 |
| 0.05 | 2000.0 | 0.000 | 0.000 | 0.000 | 0.0 |
| 0.1 | 0.5 | 0.000 | 0.000 | 0.000 | 0.0 |
| 0.1 | 2000.0 | 0.000 | 0.000 | 0.000 | 0.0 |

Δ_P mean = 1.7700, Δ_R mean = 1.00
