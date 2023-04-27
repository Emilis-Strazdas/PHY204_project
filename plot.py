"""
This file contains the plot function, which plots the potential grid.
"""

# ---------- Imports ---------- #
import numpy as np
import matplotlib.pyplot as plt

# ---------- Plot ---------- #
def plot_efficient(Square_grid, Diagonal_grid):
    """
    Plots the potential grid.
    
    Args:
        Potential_grid (numpy.ndarray): The potential grid.
        grid_parameters (dict): A dictionary containing the grid parameters.
        
    Returns:
        None
    """

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.imshow(Square_grid, cmap='inferno')
    ax2.imshow(Diagonal_grid, cmap='inferno')
    plt.show()    

def plot_finished(Square_grid, Diagonal_grid, grid_parameters, plot_type='2D'):
    """
    Plots the potential grid.

    The plot type can be either 2D or 3D.
    
    Args:
        Potential_grid (numpy.ndarray): The potential grid.
        grid_parameters (dict): A dictionary containing the grid parameters.
        
    Returns:
        None
        
    """

    # Extract the size of the grid
    L = grid_parameters['size']

    # Combine the square and diagonal grids
    Potential_grid = np.zeros((2*L, 2*L))

    # Add the sqaure parts
    Potential_grid[0:L, 0:L] = Square_grid
    Potential_grid[L:2*L, L:2*L] = np.flip(Square_grid).T       # Anti-transpose the square grid

    # Complete the diagonal matrix
    anti_diagonal = np.flip(np.diag(np.flip(Diagonal_grid, 1).diagonal()), 1)
    Potential_grid[L:2*L, 0:L] = Diagonal_grid + np.flip(Diagonal_grid).T - anti_diagonal
    
    if plot_type == '2D':
        # Plot the potential grid
        plt.imshow(Potential_grid, cmap='inferno')
        plt.show()
    elif plot_type == '3D':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(np.arange(0, 2*L, 2*L), np.arange(0, 2*L, 2*L), Potential_grid)
        plt.show()
