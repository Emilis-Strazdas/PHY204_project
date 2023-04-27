"""
This file contains the iteration_step function, which performs one iteration
step.
"""

# ---------- Imports ---------- #
import numpy as np

# ---------- Iteration Step ---------- #
def iteration_step_efficient(Square_grid, Diagonal_grid, grid_parameters):
    """
    Performs one iteration step.
    
    Args:
        Square_grid (numpy.ndarray): The square part of the potential grid.
        Diagonal_grid (numpy.ndarray): The diagonal part of the potential grid.
        grid_parameters (dict): A dictionary containing the grid parameters.
        
    Returns:
        Potential_grid_new (numpy.ndarray): The updated Square_grid grid.
        Diagonal_grid_new (numpy.ndarray): The updated Diagonal_grid grid.
    """

    # Get the size of the grid
    L = grid_parameters['size']

    # Initialise the new potential grid
    Square_grid_new = np.copy(Square_grid)
    Diagonal_grid_new = np.copy(Diagonal_grid)

    # Update the square part
    for i in range(1, L-1):
        for j in range(L-2, 0, -1):
            Square_grid_new[i, j]     = (Square_grid[i+1, j]     + Square_grid[i-1, j]   + Square_grid[i, j+1]   + Square_grid[i, j-1]) * 0.25
    # Update the boundary of the square part
    for j in range(L-2, 0, -1):
            Square_grid_new[L-1, j]   = (Square_grid[L-2, j]     + Diagonal_grid[0, j]   + Square_grid[L-1, j-1] + Square_grid[L-1, j+1]) * 0.25

    # # Update the boundary of the diagonal part
    for j in range(L-2, 0, -1):
            Diagonal_grid_new[0, j]   = (Diagonal_grid[1, j]     + Square_grid[L-1, j]   + Diagonal_grid[0, j-1] + Diagonal_grid[0, j+1]) * 0.25

    # Update the diagonal part
    for i in range(1, L-2):
        for j in range(L-2-i, 1, -1):
            Diagonal_grid_new[i, j]   = (Diagonal_grid[i+1, j]   + Diagonal_grid[i-1, j] + Diagonal_grid[i, j+1] + Diagonal_grid[i, j-1]) * 0.25

    # Update the anti-diagonal
    for i in range(1, L-1):
            Diagonal_grid_new[i, L-1-i] = (Diagonal_grid[i-1, L-1-i] + Diagonal_grid[i, L-2-i]) * 0.5

    return Square_grid_new, Diagonal_grid_new
