"""
This file contains the iteration function, which iterates the potential grid
until the error is below a certain tolerance.
"""

# ---------- Imports ---------- #
import numpy as np
import matplotlib.pyplot as plt
from iteration_step import iteration_step

from iteration_step import iteration_step_efficient
from iteration_step import iteration_step_efficient_omega
from iteration_step import *

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

def iteration_efficient(Square_grid, Diagonal_grid, grid_parameters, iteration_parameters, omega=False):
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

    if omega:
        # Iterating the potential grid
        while max_iter > 0:
            # Counting the number of iterations
            max_iter -= 1
            # Performing one iteration step
            Square_grid_new, Diagonal_grid_new = iteration_step_efficient_omega(Square_grid, Diagonal_grid, grid_parameters)
            # Calculating the error
            if np.max(Square_grid_new - Square_grid) < tolerance:
                break
            # Updating the potential grid
            Square_grid, Diagonal_grid = Square_grid_new, Diagonal_grid_new  

    else:
        # Iterating the potential grid
        while max_iter > 0:
            # Counting the number of iterations
            max_iter -= 1
            # Performing one iteration step
            Square_grid_new, Diagonal_grid_new = iteration_step_efficient(Square_grid, Diagonal_grid, grid_parameters)
            # Calculating the error
            if np.max(Square_grid_new - Square_grid) < tolerance:
                break
            # Updating the potential grid
            Square_grid, Diagonal_grid = Square_grid_new, Diagonal_grid_new

    return Square_grid_new, Diagonal_grid_new, max_iter

def iteration_square(Square_grid, grid_parameters, iteration_parameters, Error_grid=False):
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

    if isinstance(Error_grid, np.ndarray):
        while max_iter > 0:
            # Counting the number of iterations
            max_iter -= 1
            # Performing one iteration step
            Square_grid_new = iteration_step_square_error(Square_grid, Error_grid, grid_parameters)
            # Calculating the error
            error = np.max(np.abs(Square_grid - Square_grid_new))
            # Checking if the error is below the tolerance
            if error < tolerance:
                break
            # Updating the potential grid
            Square_grid = Square_grid_new

    # Iterating the potential grid
    while max_iter > 0:
        # Counting the number of iterations
        max_iter -= 1
        # Performing one iteration step
        Square_grid_new = iteration_step_square(Square_grid, Error_grid, grid_parameters)
        # Calculating the error
        error = np.max(np.abs(Square_grid - Square_grid_new))
        # Checking if the error is below the tolerance
        if error < tolerance:
            break
        # Updating the potential grid
        Square_grid = Square_grid_new

    return Square_grid

def initialise_square():
    L = 100
    V0 = 1
    dV = V0/L
    Potential_grid = np.zeros((L, L))
    Error_grid = np.zeros((L, L))
    for i in range(L):
        Potential_grid[i, L-1] = dV*i
    for i in range(L):
        for j in range(L):
            Error_grid[i, j] = 1
    return Potential_grid, Error_grid

# Potential_grid, Error_grid = initialise_square()
# grid_parameters = {'size': 100, 'w': 1.9}
# iteration_parameters = {'max_iterations': 10000, 'tolerance': 1e-5}
# print(time_function(lambda: iteration_square(Potential_grid, grid_parameters, iteration_parameters)))
# print(time_function(lambda: iteration_square(Potential_grid, grid_parameters, iteration_parameters, Error_grid=Error_grid)))