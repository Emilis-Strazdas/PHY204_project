"""
This file contains functions for matrix operations. These functions are used in
the iteration_step.py file.
"""

# ---------- Imports ---------- #
import numpy as np

# ---------- Matrix Operations ---------- #
def anti_transpose(A):
    """
    Returns the anti-transpose of a matrix.
    
    Args:
        A (numpy.ndarray): The matrix to be anti-transposed.
        
    Returns:
        A (numpy.ndarray): The anti-transposed matrix.
    """
    return np.flip(A).T

def complete_anti_diagonal_symmetric(A):
    """
    Returns the fully computed matrix given a symmetric matrix' upper triagular
    part.
    
    Args:
        A (numpy.ndarray): The initial triangular matrix.
        
    Returns:
        A (numpy.ndarray): The anti-diagonal symmetric matrix.
    """
    anti_diagonal = np.flip(np.diag(np.flip(A, 1).diagonal()), 1)
    A = A + anti_transpose(A) - anti_diagonal
    return A