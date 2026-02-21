"""
McEliece Cryptographic QUBO SA 실험 프레임워크

실험 1: m-Scaling — GF(2^m) 확장 차수 증가에 따른 SA 난이도 변화
실험 2: t-Sweep — 에러 정정 능력 t 증가에 따른 난이도 변화
실험 3: Sweep 전이 — SA sweep 수 vs 성공률 S-curve
실험 4: 6-Way 비교 — McEliece vs Hardened vs Posiform vs Quiet vs Wishart vs ZeroExp
"""

import random
import time
import sys
import os
import numpy as np
import neal

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mceliece.qubo_mceliece import (
    create_qubo_mceliece,
    extract_original_solution,
)
from qubo_utils import calculate_energy


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def classify_result(target, found_solution, target_energy, found_energy):
    """결과 분류: EXACT / ENERGY_MATCH / FAIL"""
    if found_solution == target:
        return "EXACT"
    if abs(found_energy - target_energy) < 1e-4:
        return "ENERGY_MATCH"
    return "FAIL"


# ============================================================
# 실험 1: m-Scaling
# ============================================================

def run_m_scaling(m_values=None, t=2, num_runs=10,
                  num_reads=200, seed_base=42):
    """
    실험 1: GF(2^m) 확장 차수 증가에 따른 SA 난이도 변화 측정.

    m이 커지면 N=2^m 증가 → 보조변수 증가 → QUBO 크기 증가 → SA 어려워짐.
    주의: m≥5는 Rosenberg 차수축소의 지수적 비용으로 QUBO 생성이 매우 느림.
    """
    if m_values is None:
        m_values = [3, 4]

    print("=" * 100)
    print(f"실험 1: m-Scaling (GF(2^m) 확장 차수 vs SA 성공률)")
    print(f"m_values={m_values}, t={t}")
    print(f"num_runs={num_runs}, num_reads={num_reads}")
    print("=" * 100)

    sampler = neal.SimulatedAnnealingSampler()
    results = []

    for m in m_values:
        # m과 t에 따라 가능한 target_k 자동 계산
        # t=1이면 Goppa 근 1개 제거로 N이 줄어남
        N_est = (1 << m) - (1 if t <= 1 else 0)
        target_k = max(2, min(5, N_est - m * t))

        print(f"\n{'─' * 80}")
        print(f"  m={m} (N=2^{m}={1 << m}, target_k={target_k})")
        print(f"{'─' * 80}")

        success_count = 0
        total_valid = 0
        hamming_dists = []
        gen_times = []
        solve_times = []
        infos = []

        for run in range(num_runs):
            rng = random.Random(seed_base + run * 7)
            target = ''.join(str(rng.randint(0, 1)) for _ in range(target_k))

            t0 = time.time()
            try:
                Q, info = create_qubo_mceliece(target, m=m, t=t, seed=seed_base + run)
            except (ValueError, RuntimeError) as e:
                print(f"  Run {run+1}: SKIP ({e})")
                continue
            gen_time = time.time() - t0
            gen_times.append(gen_time)

            num_sweeps = max(1000, 10 * info['total_vars'])

            t1 = time.time()
            ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
            solve_time = time.time() - t1
            solve_times.append(solve_time)

            # 원래 k 변수 추출
            best_sample = ss.first.sample
            found_orig = extract_original_solution(best_sample, info['n'])
            expected_orig = info['full_target'][:info['n']]

            total_valid += 1
            hdist = hamming_distance(expected_orig, found_orig)
            hamming_dists.append(hdist)

            if found_orig == expected_orig:
                success_count += 1
                result_str = "EXACT"
            else:
                result_str = "FAIL"

            if not infos:
                infos.append(info)

            print(f"  Run {run+1}: {result_str} | k={info['n']}, total_vars={info['total_vars']}, "
                  f"aux={info['num_aux']} | Hamming={hdist} | "
                  f"gen={gen_time:.2f}s, solve={solve_time:.2f}s")

        rate = 100.0 * success_count / total_valid if total_valid > 0 else 0
        avg_h = np.mean(hamming_dists) if hamming_dists else float('nan')
        avg_gen = np.mean(gen_times) if gen_times else float('nan')
        avg_solve = np.mean(solve_times) if solve_times else float('nan')

        info_sample = infos[0] if infos else {}
        print(f"\n  [m={m}] 성공률: {success_count}/{total_valid} ({rate:.1f}%)")
        print(f"    N={info_sample.get('N', '?')}, k={info_sample.get('n', '?')}, "
              f"total_vars={info_sample.get('total_vars', '?')}, "
              f"aux={info_sample.get('num_aux', '?')}")
        print(f"    평균 Hamming: {avg_h:.1f}, 평균 gen: {avg_gen:.2f}s, "
              f"평균 solve: {avg_solve:.2f}s")

        results.append({
            'm': m, 'N': info_sample.get('N', 0),
            'k': info_sample.get('n', 0),
            'total_vars': info_sample.get('total_vars', 0),
            'num_aux': info_sample.get('num_aux', 0),
            'success_rate': rate, 'success_count': success_count,
            'total_valid': total_valid,
            'avg_hamming': avg_h, 'avg_gen_time': avg_gen,
            'avg_solve_time': avg_solve,
        })

    # 요약 테이블
    print("\n" + "=" * 100)
    print("m-Scaling 요약")
    print("=" * 100)
    print(f"{'m':<4} | {'N':<6} | {'k':<6} | {'Vars':<6} | {'Aux':<6} | "
          f"{'Success%':<10} | {'Count':<8} | {'Avg Ham':<8} | {'Gen(s)':<8} | {'Solve(s)':<8}")
    print("-" * 90)
    for r in results:
        print(f"{r['m']:<4} | {r['N']:<6} | {r['k']:<6} | {r['total_vars']:<6} | "
              f"{r['num_aux']:<6} | {r['success_rate']:<10.1f} | "
              f"{r['success_count']}/{r['total_valid']:<5} | {r['avg_hamming']:<8.1f} | "
              f"{r['avg_gen_time']:<8.2f} | {r['avg_solve_time']:<8.2f}")

    return results


# ============================================================
# 실험 2: t-Parameter Sweep
# ============================================================

def run_t_sweep(t_values=None, m=4, num_runs=10,
                num_reads=200, seed_base=42):
    """
    실험 2: 에러 정정 능력 t 증가에 따른 난이도 변화.

    t 증가 → k(=N-m*t) 감소 + 보조변수 증가 (고차 상호작용 증가).
    m=4 사용 (m≥5는 차수축소 비용으로 생성이 느림).
    """
    if t_values is None:
        t_values = [1, 2, 3]

    N = 1 << m  # 2^m

    print("=" * 100)
    print(f"실험 2: t-Parameter Sweep (에러 정정 능력 vs SA 성공률)")
    print(f"m={m} (N={N}), t_values={t_values}")
    print(f"num_runs={num_runs}, num_reads={num_reads}")
    print("=" * 100)

    sampler = neal.SimulatedAnnealingSampler()
    results = []

    for t in t_values:
        expected_k = N - m * t
        if expected_k <= 0:
            print(f"\n  t={t}: SKIP (k={expected_k} ≤ 0)")
            continue

        # target 길이: 실제 코드 차원보다 작게 설정
        # t=1이면 N이 1 줄어남
        N_est = (1 << m) - (1 if t <= 1 else 0)
        actual_expected_k = N_est - m * t
        target_len = max(2, min(actual_expected_k, 5))

        print(f"\n{'─' * 80}")
        print(f"  t={t} (예상 k≈{expected_k}, target_len={target_len})")
        print(f"{'─' * 80}")

        success_count = 0
        total_valid = 0
        hamming_dists = []
        gen_times = []
        solve_times = []
        infos = []

        for run in range(num_runs):
            rng = random.Random(seed_base + run * 13 + t * 1000)
            target = ''.join(str(rng.randint(0, 1)) for _ in range(target_len))

            t0 = time.time()
            try:
                Q, info = create_qubo_mceliece(target, m=m, t=t,
                                                seed=seed_base + run + t * 100)
            except (ValueError, RuntimeError) as e:
                print(f"  Run {run+1}: SKIP ({e})")
                continue
            gen_time = time.time() - t0
            gen_times.append(gen_time)

            num_sweeps = max(1000, 10 * info['total_vars'])

            t1 = time.time()
            ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)
            solve_time = time.time() - t1
            solve_times.append(solve_time)

            best_sample = ss.first.sample
            found_orig = extract_original_solution(best_sample, info['n'])
            expected_orig = info['full_target'][:info['n']]

            total_valid += 1
            hdist = hamming_distance(expected_orig, found_orig)
            hamming_dists.append(hdist)

            if found_orig == expected_orig:
                success_count += 1
                result_str = "EXACT"
            else:
                result_str = "FAIL"

            if not infos:
                infos.append(info)

            print(f"  Run {run+1}: {result_str} | k={info['n']}, total_vars={info['total_vars']}, "
                  f"aux={info['num_aux']} | Hamming={hdist} | "
                  f"gen={gen_time:.2f}s, solve={solve_time:.2f}s")

        rate = 100.0 * success_count / total_valid if total_valid > 0 else 0
        avg_h = np.mean(hamming_dists) if hamming_dists else float('nan')
        avg_gen = np.mean(gen_times) if gen_times else float('nan')
        avg_solve = np.mean(solve_times) if solve_times else float('nan')

        info_sample = infos[0] if infos else {}
        print(f"\n  [t={t}] 성공률: {success_count}/{total_valid} ({rate:.1f}%)")
        print(f"    k={info_sample.get('n', '?')}, total_vars={info_sample.get('total_vars', '?')}, "
              f"aux={info_sample.get('num_aux', '?')}")
        print(f"    평균 Hamming: {avg_h:.1f}, 평균 gen: {avg_gen:.2f}s, "
              f"평균 solve: {avg_solve:.2f}s")

        results.append({
            't': t, 'k': info_sample.get('n', 0),
            'total_vars': info_sample.get('total_vars', 0),
            'num_aux': info_sample.get('num_aux', 0),
            'success_rate': rate, 'success_count': success_count,
            'total_valid': total_valid,
            'avg_hamming': avg_h, 'avg_gen_time': avg_gen,
            'avg_solve_time': avg_solve,
        })

    # 요약 테이블
    print("\n" + "=" * 100)
    print("t-Sweep 요약")
    print("=" * 100)
    print(f"{'t':<4} | {'k':<6} | {'Vars':<6} | {'Aux':<6} | "
          f"{'Success%':<10} | {'Count':<8} | {'Avg Ham':<8} | {'Gen(s)':<8} | {'Solve(s)':<8}")
    print("-" * 80)
    for r in results:
        print(f"{r['t']:<4} | {r['k']:<6} | {r['total_vars']:<6} | "
              f"{r['num_aux']:<6} | {r['success_rate']:<10.1f} | "
              f"{r['success_count']}/{r['total_valid']:<5} | {r['avg_hamming']:<8.1f} | "
              f"{r['avg_gen_time']:<8.2f} | {r['avg_solve_time']:<8.2f}")

    return results


# ============================================================
# 실험 3: Sweep 전이
# ============================================================

def run_sweep_transition(num_instances=10, reads_per_sweep=50, seed_base=42):
    """
    실험 3: SA sweep 수에 따른 ground-state 발견 확률의 S-curve 측정.

    인스턴스를 사전 생성 후 sweep 값별 반복 (hardened_posiform 패턴).
    """
    sweep_values = [1, 5, 10, 50, 100, 500, 1000, 5000, 10000]
    # m≥5는 차수축소 비용으로 QUBO 생성이 매우 느려 제외
    configs = [
        (3, 1),   # 작은 문제 (k≈4, total_vars≈7)
        (3, 2),   # 최소 문제 (k≈2, total_vars≈2)
        (4, 2),   # 중간 문제 (k≈8, total_vars≈33)
        (4, 3),   # 큰 문제 (k≈4, total_vars≈8)
    ]

    print("=" * 110)
    print(f"실험 3: SA Sweep Count 전이 실험")
    print(f"configs (m,t)={configs}, num_instances={num_instances}, "
          f"reads_per_sweep={reads_per_sweep}")
    print(f"Sweep values: {sweep_values}")
    print("=" * 110)

    sampler = neal.SimulatedAnnealingSampler()
    all_results = []

    for m, t in configs:
        label = f"m={m},t={t}"
        print(f"\n{'─' * 90}")
        print(f"  {label}")
        print(f"{'─' * 90}")

        # QUBO 인스턴스 미리 생성
        t0 = time.time()
        instances = []
        # 보수적 추정: t=1이면 Goppa 근 1개 제거로 N이 1 줄어남
        N_est = (1 << m) - (1 if t <= 1 else 0)
        target_k = max(2, min(5, N_est - m * t))

        for run in range(num_instances):
            rng = random.Random(seed_base + run * 11 + m * 100)
            target = ''.join(str(rng.randint(0, 1)) for _ in range(target_k))

            try:
                Q, info = create_qubo_mceliece(target, m=m, t=t,
                                                seed=seed_base + run + m * 50)
                instances.append((Q, info))
            except (ValueError, RuntimeError) as e:
                print(f"    인스턴스 {run+1}: SKIP ({e})")

        gen_time = time.time() - t0
        actual = len(instances)
        if actual > 0:
            sample_info = instances[0][1]
            print(f"  인스턴스 생성: {actual}/{num_instances} "
                  f"(k={sample_info['n']}, total_vars={sample_info['total_vars']}, "
                  f"aux={sample_info['num_aux']}, 생성 시간: {gen_time:.1f}s)")
        else:
            print(f"  인스턴스 생성 실패")
            continue

        for sweeps in sweep_values:
            total_samples = 0
            ground_state_found = 0
            total_hamming = 0

            t1 = time.time()
            for Q, info in instances:
                expected_orig = info['full_target'][:info['n']]
                n_vars = info['n']

                ss = sampler.sample_qubo(Q, num_reads=reads_per_sweep,
                                          num_sweeps=sweeps)

                for sample, energy, _ in ss.data(['sample', 'energy',
                                                   'num_occurrences']):
                    total_samples += 1
                    found = extract_original_solution(sample, n_vars)
                    if found == expected_orig:
                        ground_state_found += 1
                    total_hamming += hamming_distance(expected_orig, found)

            solve_time = time.time() - t1
            rate = (100.0 * ground_state_found / total_samples
                    if total_samples > 0 else 0)
            avg_h = (total_hamming / total_samples
                     if total_samples > 0 else float('nan'))

            print(f"  sweeps={sweeps:<6} | GS rate: {ground_state_found:>4}/{total_samples} "
                  f"({rate:>6.2f}%) | avg Hamming: {avg_h:>6.1f} | time: {solve_time:.1f}s")

            all_results.append({
                'm': m, 't': t, 'sweeps': sweeps,
                'gs_rate': rate, 'avg_hamming': avg_h,
                'gs_found': ground_state_found, 'total_samples': total_samples,
            })

    # 요약 테이블
    print("\n" + "=" * 110)
    print("Sweep 전이 요약 (Ground-State Sampling Rate %)")
    print("=" * 110)
    header = f"{'Config':<12} |"
    for s in sweep_values:
        header += f" {s:>7}"
    print(header)
    print("-" * (14 + 8 * len(sweep_values)))

    for m, t in configs:
        label = f"m={m},t={t}"
        row = f"{label:<12} |"
        for s in sweep_values:
            matching = [r for r in all_results
                        if r['m'] == m and r['t'] == t and r['sweeps'] == s]
            if matching:
                row += f" {matching[0]['gs_rate']:>6.1f}%"
            else:
                row += f" {'N/A':>7}"
        print(row)

    return all_results


# ============================================================
# 실험 4: 6-Way 비교
# ============================================================

def run_comparison(n_bits=8, num_runs=10, num_reads=200, num_sweeps=1000,
                   seed_base=42):
    """
    실험 4: McEliece를 포함한 6-way 방법론 비교.

    McEliece, Hardened Posiform, Posiform, Quiet Planting, Wishart, ZeroExp.
    n_bits=8: McEliece m=4,t=2에서 k=8에 맞춤.
    """
    from hardened_posiform.qubo_posiform_hardened import create_qubo_hardened_posiform
    from posiform.qubo_posiform import create_qubo_posiform
    from quiet_planting.qubo_quiet_planted import (
        create_qubo_quiet_planted,
        extract_original_solution as quiet_extract,
    )
    from wishart.qubo_wishart import create_qubo_wishart
    from zero_expectation.qubo_zero_expectation import create_qubo_precise

    quiet_alpha = 4.2
    wishart_alpha = 0.7

    print("=" * 120)
    print(f"실험 4: 6-Way 방법론 비교")
    print(f"n_bits={n_bits}, num_runs={num_runs}, num_reads={num_reads}, "
          f"num_sweeps={num_sweeps}")
    print(f"  McEliece:   m=auto, t=2")
    print(f"  Hardened:   lin2, α=0.1")
    print(f"  Posiform:   coeff=(1.0, 3.0)")
    print(f"  Quiet:      alpha={quiet_alpha}")
    print(f"  Wishart:    alpha={wishart_alpha}")
    print(f"  ZeroExp:    density=1.0")
    print("=" * 120)

    sampler = neal.SimulatedAnnealingSampler()

    methods = ['McEliece', 'Hardened', 'Posiform', 'Quiet', 'Wishart', 'ZeroExp']
    counts = {m: {"EXACT": 0, "ENERGY_MATCH": 0, "FAIL": 0} for m in methods}
    hammings = {m: [] for m in methods}

    # 헤더
    header = f"{'Run':<4}"
    for m_name in methods:
        header += f" | {m_name:<12} {'Ham':<4}"
    print(f"\n{header}")
    print("-" * 120)

    quiet_m = int(quiet_alpha * n_bits)
    quiet_total = n_bits + quiet_m
    quiet_sweeps = max(1000, 10 * quiet_total)

    for run in range(num_runs):
        target = ''.join(str(random.Random(seed_base + run * 7).randint(0, 1))
                         for _ in range(n_bits))

        row_parts = [f"{run+1:<4}"]

        # --- McEliece ---
        try:
            Q_mc, info_mc = create_qubo_mceliece(target, t=2,
                                                   seed=seed_base + run)
            mc_k = info_mc['n']
            mc_sweeps = max(num_sweeps, 10 * info_mc['total_vars'])
            ss_mc = sampler.sample_qubo(Q_mc, num_reads=num_reads,
                                         num_sweeps=mc_sweeps)
            found_mc_full = extract_original_solution(ss_mc.first.sample, mc_k)
            expected_mc = info_mc['full_target'][:mc_k]
            te_mc = info_mc['target_energy']
            # 성공 여부: 원래 k 변수가 모두 일치
            result_mc = classify_result(expected_mc, found_mc_full,
                                        te_mc, ss_mc.first.energy)
            # Hamming은 원래 n_bits에 대해 계산 (패딩 부분 제외)
            pad_len = mc_k - n_bits
            found_mc_orig = found_mc_full[pad_len:]
            hdist_mc = hamming_distance(target, found_mc_orig)
        except (ValueError, RuntimeError):
            result_mc = "FAIL"
            hdist_mc = n_bits

        # --- Hardened Posiform ---
        Q_hp, info_hp = create_qubo_hardened_posiform(
            n_bits, max_subgraph_size=15, coeff_type='lin2',
            posiform_scale=0.1, seed=seed_base + run * 53
        )
        te_hp = info_hp['target_energy']
        target_hp = info_hp['target']
        ss_hp = sampler.sample_qubo(Q_hp, num_reads=num_reads,
                                     num_sweeps=num_sweeps)
        found_hp = ''.join(str(ss_hp.first.sample[k]) for k in range(n_bits))
        result_hp = classify_result(target_hp, found_hp, te_hp,
                                     ss_hp.first.energy)
        hdist_hp = hamming_distance(target_hp, found_hp)

        # --- Posiform ---
        Q_p, info_p = create_qubo_posiform(target, coeff_range=(1.0, 3.0))
        te_p = info_p['target_energy']
        ss_p = sampler.sample_qubo(Q_p, num_reads=num_reads,
                                    num_sweeps=num_sweeps)
        found_p = ''.join(str(ss_p.first.sample[k]) for k in range(n_bits))
        result_p = classify_result(target, found_p, te_p, ss_p.first.energy)
        hdist_p = hamming_distance(target, found_p)

        # --- Quiet Planting ---
        Q_q, clauses_q, info_q = create_qubo_quiet_planted(
            target, alpha=quiet_alpha
        )
        te_q = info_q['target_energy']
        ss_q = sampler.sample_qubo(Q_q, num_reads=num_reads,
                                    num_sweeps=quiet_sweeps)
        found_q = quiet_extract(ss_q.first.sample, n_bits)
        result_q = classify_result(target, found_q, te_q, ss_q.first.energy)
        hdist_q = hamming_distance(target, found_q)

        # --- Wishart ---
        Q_w = create_qubo_wishart(target, alpha=wishart_alpha)
        te_w = calculate_energy(target, Q_w)
        ss_w = sampler.sample_qubo(Q_w, num_reads=num_reads,
                                    num_sweeps=num_sweeps)
        found_w = ''.join(str(ss_w.first.sample[k]) for k in range(n_bits))
        result_w = classify_result(target, found_w, te_w, ss_w.first.energy)
        # Wishart 대칭성 검사
        inverse_w = ''.join('0' if b == '1' else '1' for b in target)
        if found_w == inverse_w:
            result_w = "SYM_MATCH"
        hdist_w = hamming_distance(target, found_w)

        # --- Zero Expectation ---
        Q_z = create_qubo_precise(target, density=1.0)
        te_z = calculate_energy(target, Q_z)
        ss_z = sampler.sample_qubo(Q_z, num_reads=num_reads,
                                    num_sweeps=num_sweeps)
        found_z = ''.join(str(ss_z.first.sample[k]) for k in range(n_bits))
        result_z = classify_result(target, found_z, te_z, ss_z.first.energy)
        hdist_z = hamming_distance(target, found_z)

        # 집계
        for name, result, hdist in [
            ('McEliece', result_mc, hdist_mc),
            ('Hardened', result_hp, hdist_hp),
            ('Posiform', result_p, hdist_p),
            ('Quiet', result_q, hdist_q),
            ('Wishart', result_w, hdist_w),
            ('ZeroExp', result_z, hdist_z),
        ]:
            r_key = "EXACT" if result in ("EXACT", "SYM_MATCH") else result
            if r_key not in counts[name]:
                r_key = "FAIL"
            counts[name][r_key] += 1
            hammings[name].append(hdist)
            row_parts.append(f"{result:<12} {hdist:<4}")

        print(" | ".join(row_parts))

    # 요약
    print("\n" + "=" * 100)
    print("6-Way 비교 요약")
    print("=" * 100)
    print(f"{'Method':<14} | {'Success%':<10} | {'EXACT':<6} | {'ENERGY':<6} | "
          f"{'FAIL':<6} | {'Avg Hamming':<12}")
    print("-" * 70)

    for name in methods:
        c = counts[name]
        success = c.get('EXACT', 0)
        success_rate = 100.0 * success / num_runs
        avg_h = np.mean(hammings[name])
        print(f"{name:<14} | {success_rate:<10.1f} | {success:<6} | "
              f"{c.get('ENERGY_MATCH', 0):<6} | {c.get('FAIL', 0):<6} | "
              f"{avg_h:<12.1f}")

    return counts, hammings


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    mode = "m-scaling"
    num_runs = 10

    if len(sys.argv) > 1:
        if sys.argv[1] == "--m-scaling":
            mode = "m-scaling"
        elif sys.argv[1] == "--t-sweep":
            mode = "t-sweep"
        elif sys.argv[1] == "--sweep":
            mode = "sweep"
        elif sys.argv[1] == "--compare":
            mode = "compare"
        else:
            print(f"사용법: python3 {sys.argv[0]} [--m-scaling|--t-sweep|--sweep|--compare] [num_runs]")
            sys.exit(1)

    if len(sys.argv) > 2:
        num_runs = int(sys.argv[2])

    if mode == "m-scaling":
        run_m_scaling(num_runs=num_runs)
    elif mode == "t-sweep":
        run_t_sweep(num_runs=num_runs)
    elif mode == "sweep":
        run_sweep_transition(num_instances=num_runs)
    elif mode == "compare":
        run_comparison(num_runs=num_runs)
