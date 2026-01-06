

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
import random
from qubo_hard_mode import create_qubo_hard

def visualize_qubo_surface(Q, n_bits, filename="qubo_smooth_viz.png"):
    """
    QUBO 행렬을 MATLAB 스타일의 부드러운 3D Surface Plot으로 시각화합니다.
    """
    print(f"Generating Smooth 3D Surface Plot for {n_bits}x{n_bits} QUBO matrix...")
    
    # 1. Prepare Data
    # For surface plots, showing the whole structure is often better, 
    # but strictly discrete data being smoothed might be misleading.
    # We will use the full range but can subsample if huge.
    
    display_n = n_bits
    if n_bits > 100:
        print(f"Warning: N={n_bits} is large. Visualization might be crowded.")
        
    matrix = np.zeros((display_n, display_n))
    
    # Populate matrix (Symmetric for better visual symmetry, or Upper Triangular)
    # Surface plots look best when symmetric or fully filled.
    # QUBO is upper triangular, so let's mirror it to make the surface look like a "Mountain/Valley".
    for (i, j), value in Q.items():
        if i < display_n and j < display_n:
            matrix[i, j] = value
            matrix[j, i] = value # Mirror for symmetry
            
    # 2. Setup 3D Plot
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Coordinates
    X = np.arange(display_n)
    Y = np.arange(display_n)
    X, Y = np.meshgrid(X, Y)
    Z = matrix
    
    # 3. Plot Surface
    # cmap='coolwarm' (Blue-Red) or 'viridis' or 'jet'
    # antialiased=True makes it smooth
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, 
                           linewidth=0, antialiased=True, alpha=0.9)
    
    # Customize aesthetics
    ax.set_title(f"QUBO Energy Landscape (Matrix Form) - {n_bits}x{n_bits}", fontsize=15)
    ax.set_xlabel('Variable i')
    ax.set_ylabel('Variable j')
    ax.set_zlabel('Interaction Strength')
    
    # Add color bar
    fig.colorbar(surf, shrink=0.5, aspect=10)
    
    # Adjust view angle for dramatic effect
    ax.view_init(elev=30, azim=45)
    
    # Save
    plt.savefig(filename, dpi=300)
    print(f"Visualization saved to {filename}")

if __name__ == "__main__":
    # Generate a sample QUBO
    target_len = 50 # Default to 50 for good demo
    if len(sys.argv) > 1:
        target_len = int(sys.argv[1])
        
    target = "".join(str(random.randint(0, 1)) for _ in range(target_len))
    
    # Use Hard Mode with Star Backbone for interesting visualization
    # Density 0.5 makes it look like a rugged terrain
    Q = create_qubo_hard(target, density=0.5, noise_ratio=0.2)
    
    output_file = f"qubo_smooth_viz_N{target_len}.png"
    visualize_qubo_surface(Q, target_len, output_file)
