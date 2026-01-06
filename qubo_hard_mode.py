
import random
import sys
import os
import csv
from qubo_zero_expectation import print_q_matrix, print_qubo_formula

def create_qubo_hard(target_binary_string, density=1.0, base_range=(1, 3), noise_ratio=0.1):
    """
    Hard Mode QUBO 생성기
    
    기반: Ising-Derived Balanced Model
    추가: Noise Injection (Frustration)
    
    Args:
        noise_ratio (float): Frustration probability for non-backbone edges (0.0 ~ 1.0).
                             NOTE: High noise refers to the probability of flipping weak edges.
    """
    n = len(target_binary_string)
    Q = {}
    
    # Ising Spin 변환 (-1, 1)
    spins = [1 if b == '1' else -1 for b in target_binary_string]
    
    # --- Strategy: Backbone (Strong Signal) + Noise (Weak Frustration) ---
    # Strong Weight (Backbone): 20.0 (Safety Margin)
    # Weak Weight (Noise): 0.2 (Local perturbation)
    # Ratio 100:1 absolutely guarantees ground state uniqueness against any cluster flip.
    W_STRONG = 20.0
    W_WEAK = 0.2
    
    # 1. Create Backbone (Star Graph centered at 0)
    # Why Star? In a Star, every node i is connected to center 0. 
    # To flip node i, you must overcome the bond (0, i).
    # If W_STRONG > Sum(all noise weights on i), flipping is impossible.
    # Max expected noise degree ~ N * density. W_STRONG must be > N * density * W_WEAK.
    # 300 * 0.1 * 0.2 = 6.0. W_STRONG=20.0 is safe.
    connected_edges = set()
    center_node = 0
    
    for i in range(1, n):
        J_ij = W_STRONG * spins[center_node] * spins[i] # Aligned with target
        
        qt_i, qt_j, qt_w = center_node, i, -4 * J_ij 
        
        Q[(qt_i, qt_j)] = Q.get((qt_i, qt_j), 0) + qt_w
        Q[(qt_i, qt_i)] = Q.get((qt_i, qt_i), 0) + 2 * J_ij
        Q[(qt_j, qt_j)] = Q.get((qt_j, qt_j), 0) + 2 * J_ij
        
        connected_edges.add((min(center_node, i), max(center_node, i)))
        
    # 2. Add Noise Edges (Frustration)
    # Randomly add extra edges. If they flip (frustration), they are weak.
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in connected_edges:
                continue
                
            if random.random() > density:
                continue
            
            # Decide if this edge supports Target or fights it
            is_noise = random.random() < noise_ratio
            
            # Base interaction: Aligned with target
            interaction_sign = 1.0
            if is_noise:
                interaction_sign = -1.0 # Flip! (Frustration)
                
            J_ij = W_WEAK * interaction_sign * spins[i] * spins[j]
            
            # Convert to QUBO
            q_ij = -4 * J_ij
            q_i = 2 * J_ij
            q_j = 2 * J_ij
            
            Q[(i, j)] = Q.get((i, j), 0) + q_ij
            Q[(i, i)] = Q.get((i, i), 0) + q_i
            Q[(j, j)] = Q.get((j, j), 0) + q_j
            
    return Q

if __name__ == "__main__":
    random.seed(42)
    
    target = "11001"
    noise = 0.1
    
    if len(sys.argv) > 1:
        if sys.argv[1].isdigit():
            length = int(sys.argv[1])
            target = "".join(str(random.randint(0, 1)) for _ in range(length))
        else:
            target = sys.argv[1]
            
    if len(sys.argv) > 2:
        noise = float(sys.argv[2])
        
    print(f"Hard Mode QUBO 생성 (Target Len: {len(target)}, Noise: {noise})")
    print(f"Target: {target}")
    
    Q = create_qubo_hard(target, noise_ratio=noise)
    
    if len(target) <= 20:
        print_q_matrix(Q, len(target))
        print_qubo_formula(Q)
