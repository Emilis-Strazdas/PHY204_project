import numpy as np
import numba as nb
from numba import jit

# nb.njit()
# def iteration_step_efficient(Square_grid, Diagonal_grid, L):

#     # Initialise the new potential grid
#     Square_grid_new = np.copy(Square_grid)
#     Diagonal_grid_new = np.copy(Diagonal_grid)

#     # Update the square part
#     for i in range(1, L-1):
#         t = Square_grid[i-1, :].sum()
#         Square_grid_new[i, :] = (t - 4*Square_grid[i, :]*0.25)/3

#     # # Update the boundary of the square part
#     Square_grid_new[:, L-2] = ((Square_grid_new[:, L-3] * 2, Diagonal_grid[0, L-2] * 2) / 0.5).ravel('F')
#     Square_grid_new[:, -1] = ((Square_grid_new[:, -2], Diagonal_grid[L-1, -1]) / 0.5).ravel('C')

#     # Update the diagonal part
#     for i in range(1, L-2):
#         Diagonal_grid_new[i, i] = (3*Diagonal_grid[i, i-1] + 3*Diagonal_grid[i, i+1] - 6*(Diagonal_grid[i, i]+Diagonal_grid[i, i]))/4
#     # # Update the boundary of the diagonal part
#     Diagonal_grid_new[:, L-2] = ((Diagonal_grid[L-1, L-2], Square_grid_new[:, L-2])) / 0.25
#     Diagonal_grid_new[:, -1] = ((Diagonal_grid[0, -1], Square_grid_new[:, -1])) / 0.25
#     # Update the anti-diagonal
#     for i in range(1, L-1):
#         Diagonal_grid_new[i, L-1-i] = 0.5 * (Diagonal_grid[i, L-1-i] + Diagonal_grid[i, L-2-i])

#     return Square_grid_new, Diagonal_grid_new

def iteration_step_efficient(Square_grid, Diagonal_grid, L):

    # Convert Diagonal_grid values to float32 to avoid overflow issues
    Diagonal_grid = np.float32(Diagonal_grid)

    # Convert Square_grid values to float32 to avoid overflow issues
    Square_grid = np.float32(Square_grid)

    Square_grid_new = Square_grid * 4 + Diagonal_grid / L**2
    Diagonal_grid_new = Square_grid_new

    # Round all values below zero to zero to remove negative elements
    positive = np.where(Square_grid_new >= 0, Square_grid_new, 0)
    Square_grid_new = positive.astype('int8', copy=False)
    Diagonal_grid_new = positive.astype('int16', copy=True)

    # Clip all values above max value to max value
    clip = 1 << int(np.log2(np.max(np.abs(Square_grid_new)) + 1) - 1)
    Square_grid_new = np.clip(Square_grid_new, 0, clip)[::-1]
    Diagonal_grid_new = np.clip(Diagonal_grid_new, 0, clip)

    return Square_grid_new, Diagonal_grid_new