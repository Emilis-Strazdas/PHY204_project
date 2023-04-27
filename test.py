import numpy as np
from time_function import *
import matplotlib.pyplot as plt

def initialise_square(grid_parameters):
    L = grid_parameters['size']
    V0 = grid_parameters['V0']
    dV = V0/L
    Square_grid = np.zeros((L, L))
    Error_grid = np.zeros(L)
    for i in range(L):
        Square_grid[i, L-1] = dV*i
        Error_grid[i] = 1
    return Square_grid, Error_grid

def iteration_step_square(Potential_grid, grid_parameters, Error_grid):
     # Get the size of the grid
    L = grid_parameters['size']

    # Initialise the new potential grid
    Potential_grid_new = np.copy(Potential_grid)

    if isinstance(Error_grid, np.ndarray):
        # Update the upper left square
        for j in range(L-2, 0, -1):
            # if Error_grid[j] > 1e-5 and Error_grid[j] != -1:
            
                Error_grid[j] = -1
                for i in range(1, L-1):
                    Potential_grid_new[i, j] = (Potential_grid[i+1, j] + Potential_grid[i-1, j] + Potential_grid[i, j+1] + Potential_grid[i, j-1]) * 0.25
                    Error = np.abs(Potential_grid_new[i, j] - Potential_grid[i, j])
                    if Error > Error_grid[j] and Potential_grid_new[i, j] != 0:
                        Error_grid[j] = Error
    else:
        # Update the upper left square
        for j in range(L-2, 0, -1):
                for i in range(1, L-1):
                    Potential_grid_new[i, j] = (Potential_grid[i+1, j] + Potential_grid[i-1, j] + Potential_grid[i, j+1] + Potential_grid[i, j-1]) * 0.25

    return Potential_grid_new, Error_grid

def iteration_square(Square_grid, grid_parameters, iteration_parameters, Error_grid=None):
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
        Square_grid_new, Error_grid_new = iteration_step_square(Square_grid, grid_parameters, Error_grid=Error_grid)
        # Calculating the error
        error = np.max(Square_grid_new - Square_grid)
        # Checking if the error is below the tolerance
        if error < tolerance:
            print('Converged after', iteration_parameters['max_iterations'] - max_iter, 'iterations.')
            break
        # Updating the potential grid
        Square_grid = Square_grid_new
        Error_grid = Error_grid_new
        print(Error_grid)

    # plt.imshow(Square_grid)
    # plt.show()  

    return Square_grid

def main():
    grid_parameters = {'size': 10, 'V0': 1}
    iteration_parameters = {'max_iterations': 10000, 'tolerance': 1e-5}
    Square_grid, Error_grid = initialise_square(grid_parameters)

    # Without error grid
    time_simple = time_function(lambda: iteration_square(Square_grid, grid_parameters, iteration_parameters))
    
    # With error grid
    time_error  = time_function(lambda: iteration_square(Square_grid, grid_parameters, iteration_parameters, Error_grid=Error_grid))

    print('Time without error grid: ', np.round(time_simple, 5), 's')
    print('Time with error grid: ', np.round(time_error, 5), 's')

if __name__ == '__main__':
    main()
