"""
This file contains the initialisation function, which initialises the potential grid.
"""

# ---------- Imports ---------- #
import numpy as np

# ---------- Initialisation ---------- #
def initialise(grid_parameters):
    """
    Initialises the potential grid.

    The potential grid is initialised to zero everywhere except for the two lines of
    length L, where the potential is linearly increasing from 0 to V0 from the boundary
    to the centre of the grid.

    Args:
        grid_parameters (dict): A dictionary containing the grid parameters.

    Returns:
        Potential_grid (numpy.ndarray): The initialised potential grid.
    """

    # Extract grid parameters
    L = grid_parameters['size']
    V0 = grid_parameters['V0']

    # Calculate the potential step
    dV = V0 / L

    # Initialise the potential grid
    Potential_grid = np.zeros((2*L, 2*L))

    # Update the upper left square
    for i in range(1, L):
        Potential_grid[i, L] = i * dV
        Potential_grid[L, L+i] = V0 - i * dV

    Potential_grid[L, L] = V0

    return Potential_grid