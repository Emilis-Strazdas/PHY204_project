"""
This file contains the build function, which computes the potential grid.
"""

# ---------- Imports ---------- #
from initialise import initialise_efficient
from iteration import iteration_efficient
from construct import construct

# ---------- Build ---------- #
def build(L):
    """
    Computes the potential grid.
    
    Args:
        L (int): The number of grid points in one direction.

    Returns:
        Potential_grid (numpy.ndarray): The potential grid.
    """
    Square_grid, Diagonal_grid = initialise_efficient(L)
    Square_grid, Diagonal_grid = iteration_efficient(Square_grid, Diagonal_grid, L)
    return construct(Square_grid, Diagonal_grid)
