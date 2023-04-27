"""
This file contains the function that constructs the full potential grid from the
square and diagonal parts.
"""

# ---------- Imports ----------#
import numpy as np
from matrix_operations import anti_transpose
from matrix_operations import complete_anti_diagonal_symmetric

# ---------- Construct ----------#
def construct(Square_grid, Diagonal_grid):
    """
    Constructs the full potential grid from the square and diagonal parts.
    
    Args:
        Square_grid (numpy.ndarray): The square part of the potential grid.
        Diagonal_grid (numpy.ndarray): The diagonal part of the potential grid.
        
    Returns:
        Potential_grid (numpy.ndarray): The full potential grid.
    """
    # Extract the size of the grid
    L = len(Square_grid)

    # Initialize full grid
    Potential_grid = np.zeros((2*L, 2*L))

    # Add the sqaure parts
    Potential_grid[0:L, 0:L] = Square_grid
    Potential_grid[L:2*L, L:2*L] = anti_transpose(Square_grid)

    # Complete the diagonal matrix
    Potential_grid[L:2*L, 0:L] = complete_anti_diagonal_symmetric(Diagonal_grid)

    return Potential_grid
