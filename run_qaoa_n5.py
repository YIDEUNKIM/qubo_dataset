import numpy as np
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms import QAOA
from qiskit_algorithms.optimizers import COBYLA
from qiskit.primitives import Estimator
# For newer Qiskit versions, we might use different imports, but this is standard for recent stable.
from qiskit.circuit.library import QAOAAnsatz

def run_qaoa_solver():
    print("="*60)
    print("QAOA Solver for Target: 10110 (N=5)")
    print("="*60)

    num_qubits = 5
    
    # Define Hamiltonian terms (Pauli strings, Coefficient)
    # Note: Qiskit string format is q_{n-1} ... q_1 q_0
    # So index 0 is the LAST character.
    
    pauli_list = []
    
    # 1. Constant offset (Identity)
    if 1.0876903045175745 > 1e-9:
        pauli_list.append(("I" * num_qubits, -1.0876903045175745))
        
    # 2. Linear terms (h_i * Z_i)
    # We constructed h dictionary in generation step
    h_coeffs = {0: 0.18767595902098, 1: -1.1225604522606052, 2: 1.3041055878992567, 3: 4.045402231266415, 4: -2.9267597856144807} # Dictionary {index: coeff}
    
    for idx, coeff in h_coeffs.items():
        if abs(coeff) < 1e-9: continue
        # Create string with Z at pos idx (from right)
        # e.g. for N=4, idx=0 => "IIIZ"
        op_list = ["I"] * num_qubits
        op_list[num_qubits - 1 - idx] = "Z" # Reversed index
        label = "".join(op_list)
        pauli_list.append((label, coeff))

    # 3. Quadratic terms (J_ij * Z_i * Z_j)
    J_coeffs = {(0, 1): 0.733318846038385, (0, 2): -2.1310284782891324, (0, 3): -2.5726042204327118, (0, 4): 1.3560150483548563, (1, 2): 1.712593608436046, (1, 3): 0.6727621178371507, (1, 4): -0.23395739262428644, (2, 3): -2.3745578577111686, (2, 4): 1.219229043948784, (3, 4): 1.2180560486480854} # Dictionary {(i,j): coeff}
    
    for (i, j), coeff in J_coeffs.items():
        if abs(coeff) < 1e-9: continue
        op_list = ["I"] * num_qubits
        op_list[num_qubits - 1 - i] = "Z"
        op_list[num_qubits - 1 - j] = "Z"
        label = "".join(op_list)
        pauli_list.append((label, coeff))

    # Create Operator
    hamiltonian = SparsePauliOp.from_list(pauli_list)
    
    print(f"Hamiltonian encoded with {len(pauli_list)} terms.")
    
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
    
    print("\nOptimization Complete!")
    print(f"Minimum Energy Found: {result.eigenvalue.real:.4f}")
    print(f"Optimal Parameters: {result.optimal_point}")
    
    # To get the best bitstring, we ideally interpret the circuit or sample it.
    # But compute_minimum_eigenvalue focuses on energy.
    # Let's simple-print the result for now.
    return result

if __name__ == "__main__":
    try:
        run_qaoa_solver()
    except Exception as e:
        print(f"Error running Qiskit job: {e}")
        print("Please ensure you have installed: pip install qiskit qiskit-algorithms numpy")