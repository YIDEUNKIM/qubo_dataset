"""
Quiet Planting field_strength 실험: Hardened Posiform α 실험과의 일반화 비교

Hardened Posiform에서 확인한 "Signal-Hardness Dilemma"가
Quiet Planting에서도 동일하게 나타나는지 검증.

대응 관계:
  Hardened Posiform     Quiet Planting
  ─────────────────     ──────────────
  α × P (posiform)  ↔  field_strength (planted field)
  Σ R_i (random)    ↔  SAT penalty QUBO (Rosenberg)
  α = 0 → 축퇴      ↔  field = 0 → 모든 SAT 해가 동일 에너지

측정:
  1. Target sampling rate (sweep 전이)
  2. GS 에너지 일치 vs target 일치 (에너지 분석)
  3. Hamming distance (signal 방향성)
  4. 음수 field (anti-signal) 검증
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from quiet_planting.qubo_quiet_planted import (
    create_qubo_quiet_planted,
    compute_auxiliary_values,
    extract_original_solution,
)
from qubo_utils import calculate_energy


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def add_planted_field(Q, target, field_strength, seed=None):
    """
    QUBO에 planted field를 수동으로 추가.
    field_strength > 0: target 방향 signal
    field_strength < 0: target 반대 방향 anti-signal

    Q를 복사하여 반환 (원본 수정 안 함).
    """
    rng = random.Random(seed)
    Q_new = dict(Q)
    n = len(target)

    for i in range(n):
        eps_i = rng.uniform(abs(field_strength) * 0.5, abs(field_strength) * 1.5)
        if field_strength > 0:
            # target 방향: target[i]=1이면 x_i=1 선호, target[i]=0이면 x_i=0 선호
            if target[i] == '1':
                Q_new[(i, i)] = Q_new.get((i, i), 0) - eps_i
            else:
                Q_new[(i, i)] = Q_new.get((i, i), 0) + eps_i
        else:
            # anti-signal: target 반대 방향
            if target[i] == '1':
                Q_new[(i, i)] = Q_new.get((i, i), 0) + eps_i
            else:
                Q_new[(i, i)] = Q_new.get((i, i), 0) - eps_i

    return Q_new


def run_field_experiment(n_bits=100, alpha=4.2, num_instances=10, reads_per_sweep=50):
    """
    field_strength를 음수~양수로 변화시키며 Hamming, 에너지 분석.
    """
    sweep_values = [1, 5, 10, 50, 100, 500, 1000, 5000]
    field_values = [-1.0, -0.1, 0.0, 0.1, 0.5, 1.0]

    print("=" * 120)
    print("Quiet Planting — Field Signal 실험")
    print(f"N={n_bits}, alpha={alpha}, QUBO size={n_bits + int(alpha * n_bits)}")
    print(f"num_instances={num_instances}, reads_per_sweep={reads_per_sweep}")
    print(f"field values: {field_values}")
    print("=" * 120)

    sampler = neal.SimulatedAnnealingSampler()
    all_results = []

    # base QUBO 생성 (field=0, 공유)
    print("\n[Base QUBO 인스턴스 생성 (field=0)]")
    base_instances = []
    t0 = time.time()
    for run in range(num_instances):
        target = ''.join(str(random.randint(0, 1)) for _ in range(n_bits))
        Q_base, clauses, info = create_qubo_quiet_planted(
            target, alpha=alpha, seed=run * 77, field_strength=0.0
        )
        base_instances.append((Q_base, clauses, info, target))
    gen_time = time.time() - t0
    print(f"  {len(base_instances)}개 생성 완료 ({gen_time:.1f}s)")

    for field in field_values:
        print(f"\n{'━' * 100}")
        print(f"  field_strength = {field}")
        print(f"{'━' * 100}")

        # field 적용
        instances = []
        for Q_base, clauses, info, target in base_instances:
            if abs(field) < 1e-15:
                Q = dict(Q_base)
            else:
                Q = add_planted_field(Q_base, target, field, seed=hash(target) & 0xFFFFFFFF)

            n = len(target)
            aux_str = compute_auxiliary_values(clauses, target, n)
            full_target = target + aux_str
            target_energy = calculate_energy(full_target, Q)
            instances.append((Q, clauses, info, target, full_target, target_energy))

        # local min 검사 (원래 변수만)
        print(f"\n  [Target Local Min 검사]")
        for idx, (Q, clauses, info, target, full_target, target_energy) in enumerate(instances):
            n = len(target)
            total_n = info['total_vars']
            min_flip_delta = float('inf')
            for i in range(total_n):
                flipped = list(full_target)
                flipped[i] = '0' if flipped[i] == '1' else '1'
                delta = calculate_energy(''.join(flipped), Q) - target_energy
                min_flip_delta = min(min_flip_delta, delta)
            is_local_min = min_flip_delta > -1e-10
            print(f"    inst {idx}: target_E={target_energy:>10.4f} | "
                  f"min_flip_delta={min_flip_delta:>+8.4f} | "
                  f"local_min={'O' if is_local_min else 'X'}")

        # sweep 전이
        for sweeps in sweep_values:
            total_samples = 0
            target_found = 0
            hamming_list = []

            t0 = time.time()
            for Q, clauses, info, target, full_target, target_energy in instances:
                n = len(target)
                ss = sampler.sample_qubo(Q, num_reads=reads_per_sweep,
                                         num_sweeps=sweeps)
                for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                    total_samples += 1
                    found_orig = extract_original_solution(sample, n)
                    if found_orig == target:
                        target_found += 1
                    hamming_list.append(hamming_distance(target, found_orig))

            solve_time = time.time() - t0
            rate = 100.0 * target_found / total_samples
            h_arr = np.array(hamming_list)

            print(f"  sw={sweeps:<6} | target rate: {target_found:>4}/{total_samples} "
                  f"({rate:>6.2f}%) | Hamming avg={h_arr.mean():>6.1f} "
                  f"med={np.median(h_arr):>5.0f} "
                  f"min={h_arr.min():>3} max={h_arr.max():>3} | {solve_time:.1f}s")

            all_results.append({
                'field': field, 'sweeps': sweeps,
                'target_rate': rate, 'avg_hamming': float(h_arr.mean()),
                'median_hamming': float(np.median(h_arr)),
            })

    # 에너지 분석 (sweep=1000 고정)
    print("\n" + "=" * 120)
    print("에너지 분석 (sweep=1000)")
    print("=" * 120)

    energy_summary = []
    for field in field_values:
        gs_energy_match = 0
        target_match = 0
        total_samples = 0
        hamming_all = []

        for Q_base, clauses, info, target in base_instances:
            n = len(target)
            if abs(field) < 1e-15:
                Q = dict(Q_base)
            else:
                Q = add_planted_field(Q_base, target, field, seed=hash(target) & 0xFFFFFFFF)

            aux_str = compute_auxiliary_values(clauses, target, n)
            full_target = target + aux_str
            target_energy = calculate_energy(full_target, Q)

            ss = sampler.sample_qubo(Q, num_reads=reads_per_sweep, num_sweeps=1000)
            for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                total_samples += 1
                found_orig = extract_original_solution(sample, n)

                if abs(energy - target_energy) < 1e-4:
                    gs_energy_match += 1
                if found_orig == target:
                    target_match += 1
                hamming_all.append(hamming_distance(target, found_orig))

        gs_pct = 100.0 * gs_energy_match / total_samples
        tgt_pct = 100.0 * target_match / total_samples
        h_arr = np.array(hamming_all)

        print(f"  field={field:<6} | GS에너지: {gs_energy_match}/{total_samples} ({gs_pct:.1f}%) | "
              f"Target: {target_match}/{total_samples} ({tgt_pct:.1f}%) | "
              f"차이: {gs_energy_match - target_match} | "
              f"Hamming med={np.median(h_arr):.0f}")

        energy_summary.append({
            'field': field, 'gs_pct': gs_pct, 'tgt_pct': tgt_pct,
            'hamming_med': float(np.median(h_arr)),
            'gs_not_target': gs_energy_match - target_match,
        })

    # ==========================================
    # 요약 테이블
    # ==========================================
    print("\n" + "=" * 120)
    print("요약 1: Target Sampling Rate (%)")
    print("=" * 120)
    header = f"  {'field':<8} |"
    for s in sweep_values:
        header += f" sw={s:<5}"
    print(header)
    print(f"  " + "-" * (10 + 8 * len(sweep_values)))
    for field in field_values:
        row = f"  {field:<8} |"
        for s in sweep_values:
            matching = [r for r in all_results
                        if r['field'] == field and r['sweeps'] == s]
            if matching:
                row += f" {matching[0]['target_rate']:>5.1f}% "
            else:
                row += f" {'N/A':>6} "
        print(row)

    print("\n" + "=" * 120)
    print("요약 2: Median Hamming Distance (sw=5000)")
    print("=" * 120)
    for field in field_values:
        matching = [r for r in all_results
                    if r['field'] == field and r['sweeps'] == 5000]
        if matching:
            print(f"  field={field:<6} | Hamming med={matching[0]['median_hamming']:>6.0f} "
                  f"avg={matching[0]['avg_hamming']:>6.1f}")

    print("\n" + "=" * 120)
    print("요약 3: 에너지 분석 (sweep=1000) — Degeneracy vs Hardness")
    print("=" * 120)
    print(f"  {'field':<8} | {'GS에너지%':>9} | {'Target%':>8} | {'Ham med':>7} | "
          f"{'GS찾았지만 target아닌':>20} | 해석")
    print(f"  " + "-" * 90)
    for s in energy_summary:
        if s['field'] < 0:
            interp = "anti-signal"
        elif s['field'] == 0:
            interp = "signal 없음 (축퇴)"
        else:
            interp = "signal"
        print(f"  {s['field']:<8} | {s['gs_pct']:>8.1f}% | {s['tgt_pct']:>7.1f}% | "
              f"{s['hamming_med']:>7.0f} | {s['gs_not_target']:>20} | {interp}")

    # Hardened Posiform과의 대응 표
    print("\n" + "=" * 120)
    print("요약 4: Hardened Posiform과의 대응")
    print("=" * 120)
    print("  Hardened Posiform          Quiet Planting")
    print("  ─────────────────          ──────────────")
    print("  α < 0  → anti-signal   ↔  field < 0  → anti-signal")
    print("  α = 0  → 축퇴           ↔  field = 0  → 축퇴 (모든 SAT해 동일 에너지)")
    print("  α > 0  → signal         ↔  field > 0  → signal")
    print()
    print("  동일 패턴 확인 여부:")
    for s in energy_summary:
        marker = ""
        if s['field'] == 0 and s['gs_not_target'] > 0:
            marker = "← GS 찾았지만 target 아님 = Hardened α=0과 동일!"
        elif s['field'] < 0 and s['gs_pct'] == 0 and s['tgt_pct'] == 0:
            marker = "← anti-signal로 target에서 밀려남"
        elif s['field'] > 0 and s['gs_pct'] > 0 and s['gs_not_target'] == 0:
            marker = "← GS = target (축퇴 제거됨)"
        print(f"    field={s['field']:<6} | GS={s['gs_pct']:>5.1f}% Target={s['tgt_pct']:>5.1f}% "
              f"Ham={s['hamming_med']:>3.0f} {marker}")

    return all_results, energy_summary


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    n_bits = 100
    num_instances = 10

    if len(sys.argv) > 1:
        n_bits = int(sys.argv[1])
    if len(sys.argv) > 2:
        num_instances = int(sys.argv[2])

    run_field_experiment(
        n_bits=n_bits, alpha=4.2,
        num_instances=num_instances, reads_per_sweep=50
    )
