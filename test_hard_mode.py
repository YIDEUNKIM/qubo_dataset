
import random
import neal
import time
import sys
from qubo_hard_mode import create_qubo_hard
from qubo_zero_expectation import calculate_energy

def run_hard_experiment(num_runs=4, n_bits=300, noise_ratio=0.1, density=1.0):
    print(f"Starting HARD MODE experiment: {num_runs} runs, N={n_bits}, Noise={noise_ratio}, Density={density}")
    print("-" * 110)
    print(f"{'Run':<4} | {'Target (Prefix)':<12} | {'Target E':<12} | {'Found E':<12} | {'Diff':<4} | {'Result':<12} | {'Bit Diff':<8}")
    print("-" * 110)

    sampler = neal.SimulatedAnnealingSampler()

    for i in range(num_runs):
        # 1. Generate Random Target
        target = "".join(str(random.randint(0, 1)) for _ in range(n_bits))
        
        # 2. Create Hard QUBO
        Q = create_qubo_hard(target, density=density, noise_ratio=noise_ratio)
        
        # 3. Calculate expected target energy
        # Note: In Hard Mode, Target might NOT be the Ground State anymore!
        target_energy = calculate_energy(target, Q)
        
        # 4. Solve using Simulated Annealing
        # Increase num_reads to give solver a fair chance
        sampleset = sampler.sample_qubo(Q, num_reads=100) 
        
        best_sample = sampleset.first.sample
        best_energy = sampleset.first.energy
        found_solution = "".join(str(best_sample[k]) for k in range(n_bits))
        
        # 5. Compare
        energy_diff = best_energy - target_energy
        
        # Check for Ising symmetry
        inverse_target = "".join('0' if b == '1' else '1' for b in target)
        
        exact_match = (found_solution == target)
        sym_match = (found_solution == inverse_target)
        
        if exact_match:
            result_str = "EXACT"
        elif sym_match:
            result_str = "SYM_MATCH"
        elif abs(energy_diff) < 1e-4:
            result_str = "ENERGY_TIE"
        elif energy_diff < 0:
            result_str = "LOWER_FOUND" # Found a state deeper than Target!
        else:
            result_str = "HIGHER_FOUND" # Failed to reach Target energy
            
        bit_diff = sum(1 for a, b in zip(target, found_solution) if a != b)
        
        print(f"{i+1:<4} | {target[:10]}.. | {target_energy:<12.2f} | {best_energy:<12.2f} | {energy_diff:<4.2f} | {result_str:<12} | {bit_diff:<8}")

    print("-" * 100)

if __name__ == "__main__":
    n = 600
    noise = 0.1
    density = 1.0
    
    if len(sys.argv) > 1:
        noise = float(sys.argv[1])
    if len(sys.argv) > 2:
        density = float(sys.argv[2])
        
    run_hard_experiment(num_runs=4, n_bits=n, noise_ratio=noise, density=density)
