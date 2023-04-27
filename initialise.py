"""
This file contains the initialisation function, which initialises the potential grid.
"""

# ---------- Imports ---------- #
import numpy as np

# ---------- Initialisation ---------- #
def initialise_efficient(L):
    """
    Initialises the potential grid.

    The potential grid is initialised to zero everywhere except for the two lines of
    length L, where the potential is linearly increasing from 0 to V0 from the boundary
    to the centre of the grid.

    Args:
        grid_parameters (dict): A dictionary containing the grid parameters.

    Returns:
        Square_grid (numpy.ndarray): The initialised square part of the potential grid.
        Diagonal_grid (numpy.ndarray): The initialised diagonal part of the potential grid.
    """
    # Calculate the potential step
    dV = 1 / L

    # Initialise the potential grids
    Square_grid = np.zeros((L, L))
    Diagonal_grid = np.copy(Square_grid)

    # Update the square grid
    for i in range(1, L):
        Square_grid[i, L-1] = i * dV

    # Update the diagonal grid
    Diagonal_grid[0, L-1] = 1

    return Square_grid, Diagonal_grid
