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
import numpy as np
import matplotlib.pyplot as plt

# ---------- Main ---------- #

def main():
    # Define grid parameters:
    grid_parameters = {
        "size": 30,
        "V0": 1
    }
    # Define iteration parameters ('max_iterations' is optional):
    iteration_parameters = {
        "max_iterations": 1000000,
        "tolerance": 1e-5
    }

    # Initialise the grid
    Potential_grid = initialise(grid_parameters)
    
    # Iterate over the grid until the error is below the tolerance
    Potential_grid = iteration(Potential_grid, grid_parameters, iteration_parameters)

    # Plot the results
    plot(Potential_grid, grid_parameters, plot_type='3D')

# -------------------------- #

if __name__ == "__main__":
    main()