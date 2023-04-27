"""
"""

# ---------- Imports ---------- #
import numpy as np

# ---------- Matrix Operations ---------- #
def anti_transpose(A):
    return np.flip(A).T

def complete_anti_diagonal_symmetric(A):
    anti_diagonal = np.flip(np.diag(np.flip(A, 1).diagonal()), 1)
    A = A + anti_transpose(A) - anti_diagonal
    return A