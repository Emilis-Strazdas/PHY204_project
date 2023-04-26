"""
This file contains the iteration_step function, which performs one iteration
step.
"""

# ---------- Imports ---------- #
import numpy as np
from initialise import initialise
from plot import plot

# ---------- Iteration Step ---------- #
def iteration_step(Potential_grid, grid_parameters):
    """
    Performs one iteration step.
    
    Args:
        Potential_grid (numpy.ndarray): The potential grid.
        grid_parameters (dict): A dictionary containing the grid parameters.
        
    Returns:
        Potential_grid_new (numpy.ndarray): The updated potential grid.
    """

    # Get the size of the grid
    L = grid_parameters['size']

    # Initialise the new potential grid
    Potential_grid_new = np.copy(Potential_grid)

    # Update the upper left square
    for i in range(1, L+1):
        for j in range(1, L):
            Potential_grid_new[i, j] = (Potential_grid[i+1, j] + Potential_grid[i-1, j] + Potential_grid[i, j+1] + Potential_grid[i, j-1]) * 0.25
    
    # Update the lower two squares
    for i in range(L+1, 2*L-1):
        for j in range(1, 2*L-1):
            Potential_grid_new[i, j] = (Potential_grid[i+1, j] + Potential_grid[i-1, j] + Potential_grid[i, j+1] + Potential_grid[i, j-1]) * 0.25

    return Potential_grid_new