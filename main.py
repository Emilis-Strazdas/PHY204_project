#==========   PHY205 Project      ==========
# Date Created: 2023-04-26  
# Authors:      E. Strazdas
#               L. Rousseau
#=========================================

"""
    This code is used to solve the Laplace equation with fixed boundary
    conditions using the Jacobi method. The code is split into three files:
    initialise.py, iteration.py and iteration_step.py. The initialise.py
    file contains the function that initialises the grid and the boundary
    conditions. The iteration.py file contains the function that iterates
    over the grid until the error is below a certain tolerance. The
    iteration_step.py file contains the function that performs one iteration
    step. The main.py file contains the main function that runs the code.
    
    TODO:
        * Add comments.
        * Add docstrings to the import functions.
"""

# ---------- Imports ---------- #
from initialise import initialise_efficient
from iteration import iteration_efficient
from plot import plot_efficient
from plot import plot_finished
from construct import construct
from time_function import time_function


# ---------- Main ---------- #

def main_efficient():
    # Define grid parameters:
    L = 100
    # Define iteration parameters ('max_iterations' is optional):
    tolerance = 1e-5

    # Compile half of the grid
    # Square_grid, Diagonal_grid = build(grid_parameters, iteration_parameters)

    print(f'Without numba: {time_function(lambda: build(L, tolerance))}')

    # Plot the results
    # plot_finished(Square_grid, Diagonal_grid, grid_parameters, plot_type='3D')

def build(L, tolerance):
    Square_grid, Diagonal_grid = initialise_efficient(L)
    Square_grid, Diagonal_grid = iteration_efficient(Square_grid, Diagonal_grid, L, tolerance)
    return construct(Square_grid, Diagonal_grid)
# -------------------------- #

if __name__ == "__main__":
    main_efficient()
