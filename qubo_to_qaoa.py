# 수정 필요


import random
import sys
import os

# Import the QUBO generator from the existing file
try:
    from qubo_zero_expectation import create_qubo_precise
except ImportError:
    # Fallback if running in a different directory structure
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from qubo_zero_expectation import create_qubo_precise

def qubo_to_ising(Q, n):
    """
    Converts a QUBO matrix Q (for binary x_i in {0, 1}) to an Ising Hamiltonian 
    (for spin z_i in {+1, -1}).
    
    Relation: x_i = (1 - z_i) / 2
    
    Returns:
        h (dict): Linear coefficients for Z_i (term: h[i] * Z_i)
        J (dict): Quadratic coefficients for Z_i Z_j (term: J[(i,j)] * Z_i * Z_j)
        offset (float): CONSTANT energy offset
    """
    h = {}
    J = {}
    offset = 0.0
    
    # Initialize h for all variables
    for i in range(n):
        h[i] = 0.0
    
    for (i, j), weight in Q.items():
        if i == j:
            # Diagonal term Q_ii * x_i
            # x_i = (1 - z_i)/2
            # term = Q_ii/2 - Q_ii/2 * z_i
            offset += weight / 2.0
            h[i] -= weight / 2.0
        else:
            # Off-diagonal term Q_ij * x_i * x_j
            # x_i * x_j = (1 - z_i - z_j + z_i*z_j) / 4
            # term = Q_ij/4 - Q_ij/4 * z_i - Q_ij/4 * z_j + Q_ij/4 * z_i * z_j
            w_4 = weight / 4.0
            offset += w_4
            h[i] -= w_4
            h[j] -= w_4
            
            # Use tuple key (min, max) for J to handle symmetry if needed
            # Assuming input Q always has i < j for off-diagonals
            if i < j:
                key = (i, j)
            else:
                key = (j, i)
            
            J[key] = J.get(key, 0.0) + w_4
            
    return h, J, offset

def formula_string(h, J, offset):
    """Returns a string representation of the Hamiltonian."""
    terms = []
    
    # Constant
    if abs(offset) > 1e-6:
        terms.append(f"{offset:.4f}")
    
    # Linear Terms
    for i in sorted(h.keys()):
        val = h[i]
        if abs(val) > 1e-6:
            sign = "+" if val > 0 else "-"
            terms.append(f"{sign} {abs(val):.4f}*Z_{i}")
            
    # Quadratic Terms
    for (i, j) in sorted(J.keys()):
        val = J[(i, j)]
        if abs(val) > 1e-6:
            sign = "+" if val > 0 else "-"
            terms.append(f"{sign} {abs(val):.4f}*Z_{i}Z_{j}")
            
    full_str = " ".join(terms)
    # Clean up leading "+" if present
    if full_str.startswith("+ "):
        full_str = full_str[2:]
        
    return f"H = {full_str}"

def generate_qiskit_script(h, J, offset, n, target_str):
    """
    Generates the content of a Python script that uses Qiskit to run QAOA for this Hamiltonian.
    """
    # Prepare the Pauli list for SparsePauliOp
    # Format: ("ZZ...Z", coeffs) where Z is at specific positions
    # Qiskit uses little-endian (rightmost bit is 0), but often it's easier to map index i to position i directly 
    # and just be consistent.
    # Standard Qiskit convention: "II...Z...I" where Z is at index k means qubit k.
    # Note: Qiskit's from_list takes [("Label", coeff), ...]
    # Label is length N. "ZIII" -> Z on qubit 3 (if N=4). 
    # WAIT: Qiskit Operator labeling is usually "q_{n-1} ... q_0".
    # So index 0 corresponds to the RIGHTMOST character.
    
    pauli_lines = []
    
    # Constant term is usually Identity "II...I"
    if abs(offset) > 1e-9:
        pauli_lines.append(f'    ("I" * {n}, {offset}),')

    # Linear terms h_i * Z_i
    # To place Z at index i (0-indexed from right):
    # Label construction helper needed in the script?
    # Actually, let's write a helper function IN the generated script to keep the string cleaner.
    
    script_content = f"""
import numpy as np
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms import QAOA
from qiskit_algorithms.optimizers import COBYLA
from qiskit.primitives import Estimator
# For newer Qiskit versions, we might use different imports, but this is standard for recent stable.
from qiskit.circuit.library import QAOAAnsatz

def run_qaoa_solver():
    print("="*60)
    print("QAOA Solver for Target: {target_str} (N={n})")
    print("="*60)

    num_qubits = {n}
    
    # Define Hamiltonian terms (Pauli strings, Coefficient)
    # Note: Qiskit string format is q_{{n-1}} ... q_1 q_0
    # So index 0 is the LAST character.
    
    pauli_list = []
    
    # 1. Constant offset (Identity)
    if {abs(offset)} > 1e-9:
        pauli_list.append(("I" * num_qubits, {offset}))
        
    # 2. Linear terms (h_i * Z_i)
    # We constructed h dictionary in generation step
    h_coeffs = {h} # Dictionary {{index: coeff}}
    
    for idx, coeff in h_coeffs.items():
        if abs(coeff) < 1e-9: continue
        # Create string with Z at pos idx (from right)
        # e.g. for N=4, idx=0 => "IIIZ"
        op_list = ["I"] * num_qubits
        op_list[num_qubits - 1 - idx] = "Z" # Reversed index
        label = "".join(op_list)
        pauli_list.append((label, coeff))

    # 3. Quadratic terms (J_ij * Z_i * Z_j)
    J_coeffs = {J} # Dictionary {{(i,j): coeff}}
    
    for (i, j), coeff in J_coeffs.items():
        if abs(coeff) < 1e-9: continue
        op_list = ["I"] * num_qubits
        op_list[num_qubits - 1 - i] = "Z"
        op_list[num_qubits - 1 - j] = "Z"
        label = "".join(op_list)
        pauli_list.append((label, coeff))

    # Create Operator
    hamiltonian = SparsePauliOp.from_list(pauli_list)
    
    print(f"Hamiltonian encoded with {{len(pauli_list)}} terms.")
    
    # Setup QAOA
    optimizer = COBYLA(maxiter=200)
    estimator = Estimator()
    
    # In newer Qiskit Algorithms, we often use SamplingVQE or QAOA directly
    # Let's use the standard QAOA class
    qaoa = QAOA(estimator, optimizer, reps=1)
    
    # Run optimization
    # Note: QAOA finds the MINIMUM eigenvalue
    print("Running QAOA optimization...")
    result = qaoa.compute_minimum_eigenvalue(hamiltonian)
    
    print("\\nOptimization Complete!")
    print(f"Minimum Energy Found: {{result.eigenvalue.real:.4f}}")
    print(f"Optimal Parameters: {{result.optimal_point}}")
    
    # To get the best bitstring, we ideally interpret the circuit or sample it.
    # But compute_minimum_eigenvalue focuses on energy.
    # Let's simple-print the result for now.
    return result

if __name__ == "__main__":
    try:
        run_qaoa_solver()
    except Exception as e:
        print(f"Error running Qiskit job: {{e}}")
        print("Please ensure you have installed: pip install qiskit qiskit-algorithms numpy")
"""
    return script_content.strip()

def main():
    print("=" * 60)
    print("QUBO -> QAOA (Ising Hamiltonian) Converter")
    print("=" * 60)
    
    # 1. Create QUBO
    target = "10110"  # Default 5-bit example
    # Check if user provided an argument
    if len(sys.argv) > 1:
        target = sys.argv[1]
        
    print(f"Target Solution: {target}")
    
    # Generate precise QUBO
    Q = create_qubo_precise(target, density=1.0)
    n = len(target)
    
    print(f"QUBO generated with {len(Q)} terms.")
    
    # 2. Convert to Ising
    h, J, offset = qubo_to_ising(Q, n)
    
    # 3. Print Formula
    print("\nIsing Hamiltonian Form:")
    print(formula_string(h, J, offset))
    
    # 4. Generate Qiskit Script
    script_filename = f"run_qaoa_n{n}.py"
    qiskit_code = generate_qiskit_script(h, J, offset, n, target)
    
    with open(script_filename, "w") as f:
        f.write(qiskit_code)
        
    print("\n" + "-"*60)
    print(f"✓ GENERATED QISKIT SCRIPT: {script_filename}")
    print("-"*60)
    print(f"To run the QAOA simulation, install qiskit and run:")
    print(f"  pip install qiskit qiskit-algorithms")
    print(f"  python3 {script_filename}")

if __name__ == "__main__":
    main()
