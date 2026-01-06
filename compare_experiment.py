
import random
import neal
import time
from qubo_zero_expectation import create_qubo_precise, calculate_energy

def run_experiment(num_runs=4, n_bits=500, balance=False):
    mode_str = " (Balanced)" if balance else ""
    print(f"Starting comparison experiment: {num_runs} runs of {n_bits}-bit problems{mode_str}")
    print("-" * 80)
    print(f"{'Run':<5} | {'Target (First 10)':<15} | {'Target Energy':<15} | {'Found Energy':<15} | {'Time(s)':<8} | {'Result':<10}")
    print("-" * 80)

    sampler = neal.SimulatedAnnealingSampler()
    
    for i in range(num_runs):
        # 1. Generate Random Target
        target = "".join(str(random.randint(0, 1)) for _ in range(n_bits))
        
        # 2. Create QUBO (Measure generation time)
        start_time = time.time()
        Q = create_qubo_precise(target, density=1.0, balance_rows=balance)
        gen_time = time.time() - start_time
        
        # 3. Calculate expected target energy
        target_energy = calculate_energy(target, Q)
        
        # 4. Solve using Simulated Annealing
        solve_start_time = time.time()
        sampleset = sampler.sample_qubo(Q, num_reads=50) # Fast check
        solve_time = time.time() - solve_start_time
        
        best_sample = sampleset.first.sample
        best_energy = sampleset.first.energy
        # reconstruct bitstring from sample dict (keys are indices 0..n-1)
        found_solution = "".join(str(best_sample[k]) for k in range(n_bits))
        
        # 5. Compare
        energy_match = abs(best_energy - target_energy) < 1e-4
        exact_match = (found_solution == target)
        
        # Check for global spin flip (Ising symmetry)
        inverse_target = "".join('0' if b == '1' else '1' for b in target)
        symmetry_match = (found_solution == inverse_target)
        
        if exact_match:
            result_str = "EXACT_MATCH"
        elif symmetry_match:
            result_str = "SYM_MATCH"
        elif energy_match:
            result_str = "ENERGY_MATCH"
        else:
            result_str = "FAIL"
        
        # Calculate Hamming Distance for info
        diff_bits = sum(1 for a, b in zip(target, found_solution) if a != b)
        
        print(f"{i+1:<5} | {target[:10]}...    | {target_energy:<15.4f} | {best_energy:<15.4f} | {solve_time:<8.2f} | {result_str:<12} (Diff: {diff_bits})")
        print(f"{target[:10]}, {found_solution[:10]}")
    print("-" * 80)

if __name__ == "__main__":
    import sys
    n_bits = 500
    num_runs = 4
    
    if len(sys.argv) > 1:
        n_bits = int(sys.argv[1])
    if len(sys.argv) > 2:
        num_runs = int(sys.argv[2])
        
    balance_mode = "balance" in sys.argv
        
    run_experiment(num_runs=num_runs, n_bits=n_bits, balance=balance_mode)
