"""
순수 랜덤 QUBO (n=500)에서 SA 실험.
planted solution 없이 랜덤 Q 행렬 생성 → SA 실행 → 해의 일관성 분석.
"""

import numpy as np
import neal
import time


def hamming_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def create_random_qubo(n, seed=None):
    """순수 랜덤 QUBO 생성. Q_ij ~ N(0, 1)."""
    rng = np.random.RandomState(seed)
    Q = {}
    for i in range(n):
        Q[(i, i)] = rng.normal(0, 1)
        for j in range(i + 1, n):
            Q[(i, j)] = rng.normal(0, 1)
    return Q


def run_experiment(n=500, num_instances=5, num_reads=50, sweep_list=None):
    if sweep_list is None:
        sweep_list = [100, 500, 1000, 5000]

    sampler = neal.SimulatedAnnealingSampler()

    print("=" * 100)
    print(f"순수 랜덤 QUBO SA 실험")
    print(f"N={n}, instances={num_instances}, reads={num_reads}")
    print("=" * 100)

    for num_sweeps in sweep_list:
        print(f"\n{'─' * 80}")
        print(f"  sweeps = {num_sweeps}")
        print(f"{'─' * 80}")

        all_best_energies = []
        all_hamming_between_reads = []
        all_energy_spreads = []

        t0 = time.time()
        for inst in range(num_instances):
            Q = create_random_qubo(n, seed=inst * 97)

            ss = sampler.sample_qubo(Q, num_reads=num_reads, num_sweeps=num_sweeps)

            # 에너지 수집
            energies = []
            solutions = []
            for sample, energy, _ in ss.data(['sample', 'energy', 'num_occurrences']):
                energies.append(energy)
                sol = ''.join(str(sample[k]) for k in range(n))
                solutions.append(sol)

            energies = np.array(energies)
            best_idx = np.argmin(energies)
            best_energy = energies[best_idx]
            best_sol = solutions[best_idx]

            all_best_energies.append(best_energy)
            all_energy_spreads.append(energies.max() - energies.min())

            # best solution과 다른 read들 사이 Hamming distance
            hammings = [hamming_distance(best_sol, s) for s in solutions]
            all_hamming_between_reads.extend(hammings)

            # 같은 최저 에너지를 찾은 read 수
            gs_count = np.sum(np.abs(energies - best_energy) < 1e-6)

            print(f"  inst {inst}: best_E={best_energy:>10.2f} | "
                  f"best 도달: {gs_count}/{num_reads} | "
                  f"E범위: [{energies.min():.2f}, {energies.max():.2f}] | "
                  f"Hamming(best vs others): med={np.median(hammings):.0f}")

        elapsed = time.time() - t0
        h_arr = np.array(all_hamming_between_reads)

        print(f"\n  *** 종합 (sweeps={num_sweeps}) ***")
        print(f"  Best 에너지: mean={np.mean(all_best_energies):.2f}, "
              f"std={np.std(all_best_energies):.2f}")
        print(f"  에너지 범위 (max-min): mean={np.mean(all_energy_spreads):.2f}")
        print(f"  Hamming (best vs all reads): "
              f"mean={h_arr.mean():.1f}, median={np.median(h_arr):.0f}, "
              f"min={h_arr.min()}, max={h_arr.max()}")
        print(f"  Hamming=0 (best와 동일): {np.sum(h_arr == 0)}/"
              f"{len(h_arr)} ({100*np.sum(h_arr==0)/len(h_arr):.1f}%)")
        print(f"  시간: {elapsed:.1f}s")


if __name__ == "__main__":
    run_experiment(n=500, num_instances=5, num_reads=50,
                   sweep_list=[100, 500, 1000, 5000])
