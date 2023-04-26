"""
This file contains the iteration function, which iterates the potential grid
until the error is below a certain tolerance.
"""

# ---------- Imports ---------- #
import numpy as np
from iteration_step import iteration_step

# ---------- Iteration ---------- #
def iteration(Potential_grid, grid_parameters, iteration_parameters):
    """
    Iterates the potential grid until the error is below a certain tolerance.

    Max_iter (maximum number of iterations) is only here as a safety measure
    in case the error is never below the tolerance. It is not necessary for
    the code to run and it will converge without it.

    Args:
        Potential_grid (numpy.ndarray): The potential grid.
        grid_parameters (dict): A dictionary containing the grid parameters.
        iteration_parameters (dict): A dictionary containing the iteration
            parameters.
    
    Returns:
        Potential_grid (numpy.ndarray): The updated potential grid.
    """

    # Extracting the iteration parameters
    max_iter = iteration_parameters['max_iterations']
    tolerance = iteration_parameters['tolerance']

    # Iterating the potential grid
    while max_iter > 0:
        # Counting the number of iterations
        max_iter -= 1
        # Performing one iteration step
        Potential_grid_new = iteration_step(Potential_grid, grid_parameters)
        # Calculating the error
        error = np.max(np.abs(Potential_grid - Potential_grid_new))
        # Checking if the error is below the tolerance
        if error < tolerance:
            break
        # Updating the potential grid
        Potential_grid = Potential_grid_new

    return Potential_grid