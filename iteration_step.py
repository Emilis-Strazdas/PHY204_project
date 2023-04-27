"""
This file contains the iteration_step function, which performs one iteration
step.
"""

# ---------- Imports ---------- #
import numpy as np
from initialise import initialise
from plot import plot
from time_function import time_function

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

def iteration_step_efficient(Square_grid, Diagonal_grid, grid_parameters):
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

def iteration_step_efficient_omega(Square_grid, Diagonal_grid, grid_parameters):
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
    w = grid_parameters['w']

    # Initialise the new potential grid
    Square_grid_new = np.copy(Square_grid)
    Diagonal_grid_new = np.copy(Diagonal_grid)

    # Update the square part
    for i in range(1, L-1):
        for j in range(L-2, 0, -1):
            # Update using Jacobi
            Square_grid_new[i, j]     = (Square_grid[i+1, j]     + Square_grid[i-1, j]   + Square_grid[i, j+1]   + Square_grid[i, j-1]) * 0.25
            # Update using SOR
            Square_grid_new[i, j]     = w * Square_grid_new[i, j] + (1-w) * Square_grid[i, j]
    
    # Update the boundary of the square part
    for j in range(L-2, 0, -1):
            # Update using Jacobi
            Square_grid_new[L-1, j]   = (Square_grid[L-2, j]     + Diagonal_grid[0, j]   + Square_grid[L-1, j-1] + Square_grid[L-1, j+1]) * 0.25
            # Update using SOR
            Square_grid_new[L-1, j]   = w * Square_grid_new[L-1, j] + (1-w) * Square_grid[L-1, j]

    # # Update the boundary of the diagonal part
    for j in range(L-2, 0, -1):
            # Update using Jacobi
            Diagonal_grid_new[0, j]   = (Diagonal_grid[1, j]     + Square_grid[L-1, j]   + Diagonal_grid[0, j-1] + Diagonal_grid[0, j+1]) * 0.25
            # Update using SOR
            Diagonal_grid_new[0, j]   = w * Diagonal_grid_new[0, j] + (1-w) * Diagonal_grid[0, j]

    # Update the diagonal part
    for i in range(1, L-2):
        for j in range(L-2-i, 1, -1):
            # Update using Jacobi
            Diagonal_grid_new[i, j]   = (Diagonal_grid[i+1, j]   + Diagonal_grid[i-1, j] + Diagonal_grid[i, j+1] + Diagonal_grid[i, j-1]) * 0.25
            # Update using SOR
            Diagonal_grid_new[i, j]   = w * Diagonal_grid_new[i, j] + (1-w) * Diagonal_grid[i, j]

    # Update the anti-diagonal
    for i in range(1, L-1):
            # Update using Jacobi
            Diagonal_grid_new[i, L-1-i] = (Diagonal_grid[i-1, L-1-i] + Diagonal_grid[i, L-2-i]) * 0.5
            # Update using SOR
            Diagonal_grid_new[i, L-1-i] = w * Diagonal_grid_new[i, L-1-i] + (1-w) * Diagonal_grid[i, L-1-i]

    return Square_grid_new, Diagonal_grid_new

def iteration_step_square(Potential_grid, Error_grid, grid_parameters):
     # Get the size of the grid
    L = grid_parameters['size']

    # Initialise the new potential grid
    Potential_grid_new = np.copy(Potential_grid)

    # Update the upper left square
    for j in range(1, L-1):
            for i in range(1, L-1, -1):
                Potential_grid_new[i, j] = (Potential_grid[i+1, j]     + Potential_grid[i-1, j]   + Potential_grid[i, j+1]   + Potential_grid[i, j-1]) * 0.25

def iteration_step_square_error(Potential_grid, Error_grid, grid_parameters):
     # Get the size of the grid
    L = grid_parameters['size']

    # Initialise the new potential grid
    Potential_grid_new = np.copy(Potential_grid)

    # Update the upper left square
    for j in range(1, L-1):
        if Error_grid[0, j] > 1e-5:
            for i in range(1, L-1, -1):
                Potential_grid_new[i, j] = (Potential_grid[i+1, j]     + Potential_grid[i-1, j]   + Potential_grid[i, j+1]   + Potential_grid[i, j-1]) * 0.25
                Error_grid[i, j] = np.abs(Potential_grid_new[i, j] - Potential_grid[i, j])
                if Error_grid[i, j] > Error_grid[0, j]:
                    Error_grid[0, j] = Error_grid[i, j]
                    