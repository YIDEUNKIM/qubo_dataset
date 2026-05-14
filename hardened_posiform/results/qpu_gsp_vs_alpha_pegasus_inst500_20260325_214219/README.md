# QPU GSP vs α 실험 결과

**날짜**: 2026-03-25 23:35

## 실험 설정

| 항목 | 값 |
|---|---|
| QPU | Advantage_system4.1 |
| Topology | pegasus (5627 qubits) |
| coeff_type | lin2 |
| α | [0, 0.001, 0.003, 0.005, 0.01, 0.03, 0.05, 0.1] |
| annealing_time | 20μs |
| num_reads | 100 |
| num_instances | 500 |
| energy_tol | 0.5 |
| 총 소요 시간 | 3787.1s |

## 결과

| α | QPU GSP |
|---|---|
| 0 | 0.004 |
| 0.001 | 0.000 |
| 0.003 | 0.000 |
| 0.005 | 0.000 |
| 0.01 | 0.000 |
| 0.03 | 0.000 |
| 0.05 | 0.002 |
| 0.1 | 0.040 |

## 파일 목록

- `qpu_vs_sa.png/pdf` — QPU vs SA 비교
- `qpu_gsp_vs_alpha.png/pdf` — QPU 단독
- `data.json` — 전체 데이터
