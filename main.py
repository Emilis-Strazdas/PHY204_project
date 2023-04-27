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
from time_function import time_function

# ---------- Main ---------- #

def main_efficient():
    # Define grid parameters:
    grid_parameters = {
        "size": 20,
        "V0": 1,
        "w": 0.5
    }
    # Define iteration parameters ('max_iterations' is optional):
    iteration_parameters = {
        "max_iterations": 1000,
        "tolerance": 1e-5
    }

    # Compile half of the grid
    Square_grid, Diagonal_grid = build(grid_parameters, iteration_parameters)

    # Construct full grid

    # Plot the results
    plot_finished(Square_grid, Diagonal_grid, grid_parameters, plot_type='3D')

def build(grid_parameters, iteration_parameters):
    Square_grid, Diagonal_grid = initialise_efficient(grid_parameters)
    return iteration_efficient(Square_grid, Diagonal_grid, grid_parameters, iteration_parameters)

# -------------------------- #

if __name__ == "__main__":
    main_efficient()
