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
from initialise import initialise
from iteration import iteration
from plot import plot

from initialise import initialise_efficient
from iteration import iteration_efficient
from plot import plot_efficient
from plot import plot_finished

from time_function import time_function

import numpy as np
import matplotlib.pyplot as plt

# ---------- Main ---------- #

def main():
    # Define grid parameters:
    grid_parameters = {
        "size": 10,
        "V0": 1
    }
    # Define iteration parameters ('max_iterations' is optional):
    iteration_parameters = {
        "max_iterations": 1000,
        "tolerance": 1e-5
    }

    # Initialise the grid
    Potential_grid = initialise(grid_parameters)
    
    # Iterate over the grid until the error is below the tolerance
    Potential_grid = iteration(Potential_grid, grid_parameters, iteration_parameters)

    # Plot the results
    plot(Potential_grid, grid_parameters)

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

    Square_grid, Diagonal_grid, max_iter = build(grid_parameters, iteration_parameters)

    # Plot the results
    plot_finished(Square_grid, Diagonal_grid, grid_parameters)

def main_efficient_omega():
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

    Square_grid, Diagonal_grid, max_iter = build(grid_parameters, iteration_parameters, omega=True)

    # Plot the results
    # plot_finished(Square_grid, Diagonal_grid, grid_parameters)

def build(grid_parameters, iteration_parameters, omega=False):
    Square_grid, Diagonal_grid = initialise_efficient(grid_parameters)
    return iteration_efficient(Square_grid, Diagonal_grid, grid_parameters, iteration_parameters, omega=omega)

# -------------------------- #

if __name__ == "__main__":
    main_efficient()
    # main_efficient_omega()
    # Define grid parameters:
    # grid_parameters = {
    #     "size": 10,
    #     "V0": 1,
    #     "w": 1
    # }
    # # Define iteration parameters ('max_iterations' is optional):
    # iteration_parameters = {
    #     "max_iterations": 1000,
    #     "tolerance": 1e-5
    # }

    # for i in np.linspace(0.05, 2, 40):
    #     grid_parameters['w'] = i
    #     Square_grid, Diagonal_grid, max_iter = build(grid_parameters, iteration_parameters, omega=True)
    #     print(f'for {np.round(i, 2)} we have {10000 - max_iter}')
