"""기존 lin2 200-inst 10α 결과 + 신규 3α 결과 → 13α 통합 결과 dir 생성.

생성물: results/gsp_vs_alpha_pegasus16_lin2_inst200_merged_<stamp>/
  samples.csv  — 26만 rows (200×13×100)
  instances.csv — 2600 rows (200×13)
  data.json    — 13 α 통합 메타
"""
import csv
import json
import os
import sys
from datetime import datetime

ORIG = "/home/yideun/qubo_dataset/hardened_posiform/results/gsp_vs_alpha_pegasus16_lin2_inst200_20260519_204139"
# 신규 dir은 인자로 받음
if len(sys.argv) < 2:
    print("usage: merge_alphas.py <new_results_dir>")
    sys.exit(1)
NEW = sys.argv[1]

stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
OUT = f"/home/yideun/qubo_dataset/hardened_posiform/results/gsp_vs_alpha_pegasus16_lin2_inst200_merged_{stamp}"
os.makedirs(OUT, exist_ok=True)
print(f"merging:\n  orig: {ORIG}\n  new:  {NEW}\n  out:  {OUT}")

# 1) samples.csv 병합
out_samples = os.path.join(OUT, 'samples.csv')
n_orig = n_new = 0
with open(out_samples, 'w', newline='') as fo:
    w = csv.writer(fo)
    with open(os.path.join(ORIG, 'samples.csv')) as f:
        r = csv.reader(f)
        header = next(r)
        w.writerow(header)
        for row in r:
            w.writerow(row); n_orig += 1
    with open(os.path.join(NEW, 'samples.csv')) as f:
        r = csv.reader(f); next(r)
        for row in r:
            w.writerow(row); n_new += 1
print(f"  samples.csv: orig={n_orig}, new={n_new}, total={n_orig+n_new}")

# 2) instances.csv 병합
out_inst = os.path.join(OUT, 'instances.csv')
i_orig = i_new = 0
with open(out_inst, 'w', newline='') as fo:
    w = csv.writer(fo)
    with open(os.path.join(ORIG, 'instances.csv')) as f:
        r = csv.reader(f)
        header = next(r)
        w.writerow(header)
        for row in r:
            w.writerow(row); i_orig += 1
    with open(os.path.join(NEW, 'instances.csv')) as f:
        r = csv.reader(f); next(r)
        for row in r:
            w.writerow(row); i_new += 1
print(f"  instances.csv: orig={i_orig}, new={i_new}, total={i_orig+i_new}")

# 3) data.json 병합
with open(os.path.join(ORIG, 'data.json')) as f:
    orig_data = json.load(f)
with open(os.path.join(NEW, 'data.json')) as f:
    new_data = json.load(f)

# 기존 data.json에 avg_dE 없음 → samples.csv에서 계산하여 보강
print("  computing avg_dE for orig αs from samples.csv...")
from collections import defaultdict
dE_sum_orig = defaultdict(float)
hd_check = defaultdict(int)
count_orig = defaultdict(int)
with open(os.path.join(ORIG, 'samples.csv')) as f:
    r = csv.reader(f); next(r)
    for row in r:
        a = row[0]
        # row[3]=hd, row[5]=delta_energy
        dE_sum_orig[a] += float(row[5])
        count_orig[a] += 1
        hd_check[a] += int(row[3])

# orig data.json results에 avg_dE/dE_sum 보강
for a_key, r in orig_data['results'].items():
    if a_key in dE_sum_orig:
        r['dE_sum'] = dE_sum_orig[a_key]
        r['avg_dE'] = dE_sum_orig[a_key] / count_orig[a_key] if count_orig[a_key] else 0.0
    # sanity: hd_sum 일치 확인
    if a_key in hd_check and abs(hd_check[a_key] - r.get('hd_sum', -1)) > 0:
        print(f"  [warn] α={a_key}: hd_sum mismatch ({hd_check[a_key]} vs {r.get('hd_sum')})")

# 신규 results에 이미 avg_dE 있음 — 그대로 사용
merged_results = dict(orig_data['results'])
for a_key, r in new_data['results'].items():
    merged_results[a_key] = r

# alpha key를 float으로 sort
def parse_alpha(k):
    try: return float(k)
    except: return 0.0
sorted_alphas = sorted(merged_results.keys(), key=parse_alpha)

# alphas list (params용)
all_alphas = sorted(set(orig_data['params']['alphas']) | set(new_data['params']['alphas']))

merged = {
    'params': {
        'topology': orig_data['params'].get('topology', new_data['params'].get('topology')),
        'topo_size': orig_data['params'].get('topo_size', new_data['params'].get('topo_size')),
        'coeff': orig_data['params'].get('coeff', new_data['params'].get('coeff')),
        'num_instances': orig_data['params'].get('num_instances'),
        'num_reads': orig_data['params'].get('num_reads'),
        'fixed_sweep': orig_data['params'].get('fixed_sweep'),
        'alphas': all_alphas,
        'elapsed_s_orig': orig_data['params'].get('elapsed_s'),
        'elapsed_s_new': new_data['params'].get('elapsed_s'),
        'n': new_data['params'].get('n'),  # 신규엔 n 있음 (5612)
        'merge_stamp': stamp,
    },
    'results': {a: merged_results[a] for a in sorted_alphas},
}
with open(os.path.join(OUT, 'data.json'), 'w') as f:
    json.dump(merged, f, indent=2)
print(f"  data.json: {len(merged['results'])} α entries → {os.path.join(OUT, 'data.json')}")
print(f"\n  ✓ merged dir: {OUT}")
