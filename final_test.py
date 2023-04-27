
# ---------- Imports ---------- #
import numpy as np
from numba import njit
import time
import matplotlib.pyplot as plt

# ---------- Iteration Step ---------- #
# Using Numba to compile the function for faster execution
@njit
def iteration_step_efficient(Square_grid, Diagonal_grid, L):
    """
    Performs one iteration step.
    
    Args:
        Square_grid (numpy.ndarray): The square part of the potential grid.
        Diagonal_grid (numpy.ndarray): The diagonal part of the potential grid.
        grid_parameters (dict): A dictionary containing the grid parameters.
        
    Returns:
        Potential_grid_new (numpy.ndarray): The updated Square_grid grid.
        Diagonal_grid_new (numpy.ndarray): The updated Diagonal_grid grid.
    """

    # Initialise the new potential grid
    Square_grid_new = np.copy(Square_grid)
    Diagonal_grid_new = np.copy(Diagonal_grid)
    
    # Update the square part
    for i in range(1, L-1):
        for j in range(L-2, 0, -1):
            Square_grid_new[i, j]     = (Square_grid[i+1, j]     + Square_grid[i-1, j]   + Square_grid[i, j+1]   + Square_grid[i, j-1]) * 0.25
    # Update the boundary of the square part
    for j in range(L-2, 0, -1):
            Square_grid_new[L-1, j]   = (Square_grid[L-2, j]     + Diagonal_grid[0, j]   + Square_grid[L-1, j-1] + Square_grid[L-1, j+1]) * 0.25

    # # Update the boundary of the diagonal part
    for j in range(L-2, 0, -1):
            Diagonal_grid_new[0, j]   = (Diagonal_grid[1, j]     + Square_grid[L-1, j]   + Diagonal_grid[0, j-1] + Diagonal_grid[0, j+1]) * 0.25

    # Update the diagonal part
    for i in range(1, L-2):
        for j in range(L-2-i, 1, -1):
            Diagonal_grid_new[i, j]   = (Diagonal_grid[i+1, j]   + Diagonal_grid[i-1, j] + Diagonal_grid[i, j+1] + Diagonal_grid[i, j-1]) * 0.25

    # Update the anti-diagonal
    for i in range(1, L-1):
            Diagonal_grid_new[i, L-1-i] = (Diagonal_grid[i-1, L-1-i] + Diagonal_grid[i, L-2-i]) * 0.5

    return Square_grid_new, Diagonal_grid_new

# ---------- Iteration ---------- #
def iteration_efficient(Square_grid, Diagonal_grid, L):
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
    max_iter = 100000   # precaution, never reached

    # Iterating the potential grid
    while max_iter > 0:
        # Counting the number of iterations
        max_iter -= 1
        # Performing one iteration step
        Square_grid_new, Diagonal_grid_new = iteration_step_efficient(Square_grid, Diagonal_grid, L)
        # Calculating the error
        if np.max(Square_grid_new - Square_grid) < 1e-5:
            break
        # Updating the potential grid
        Square_grid, Diagonal_grid = Square_grid_new, Diagonal_grid_new

    if max_iter == 0:
        print("Warning: max_iter reached.")

    return Square_grid_new, Diagonal_grid_new

# ---------- Initialisation ---------- #
def initialise_efficient(L):
    """
    Initialises the potential grid.

    The potential grid is initialised to zero everywhere except for the two lines of
    length L, where the potential is linearly increasing from 0 to V0 from the boundary
    to the centre of the grid.

    Args:
        L (int): The number of grid points in one direction.

    Returns:
        Square_grid (numpy.ndarray): The initialised square part of the potential grid.
        Diagonal_grid (numpy.ndarray): The initialised diagonal part of the potential grid.
    """
    # Calculate the potential step
    dV = 1 / L

    # Initialise the potential grids
    Square_grid = np.zeros((L, L))
    Diagonal_grid = np.copy(Square_grid)

    # Update the square grid
    for i in range(1, L):
        Square_grid[i, L-1] = i * dV

    # Update the diagonal grid
    Diagonal_grid[0, L-1] = 1

    return Square_grid, Diagonal_grid

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

# ---------- Construct ----------#
def construct(Square_grid, Diagonal_grid):
    """
    Constructs the full potential grid from the square and diagonal parts.
    
    Args:
        Square_grid (numpy.ndarray): The square part of the potential grid.
        Diagonal_grid (numpy.ndarray): The diagonal part of the potential grid.
        
    Returns:
        Potential_grid (numpy.ndarray): The full potential grid.
    """
    # Extract the size of the grid
    L = len(Square_grid)

    # Initialize full grid
    Potential_grid = np.zeros((2*L, 2*L))

    # Add the sqaure parts
    Potential_grid[0:L, 0:L] = Square_grid
    Potential_grid[L:2*L, L:2*L] = anti_transpose(Square_grid)

    # Complete the diagonal matrix
    Potential_grid[L:2*L, 0:L] = complete_anti_diagonal_symmetric(Diagonal_grid)

    return Potential_grid

# ---------- Build ---------- #
def build(L):
    """
    Computes the potential grid.
    
    Args:
        L (int): The number of grid points in one direction.

    Returns:
        Potential_grid (numpy.ndarray): The potential grid.
    """
    Square_grid, Diagonal_grid = initialise_efficient(L)
    Square_grid, Diagonal_grid = iteration_efficient(Square_grid, Diagonal_grid, L)
    return construct(Square_grid, Diagonal_grid)

# ---------- Plot ---------- #
def plot(Potential_grid, plot_type='2D'):
    """
    Plots the potential grid.

    The plot type can be either 2D or 3D.

    Args:
        Potential_grid (numpy.ndarray): The potential grid.
        plot_type (str): The type of plot to be made. Can be either '2D' or '3D'.

    Returns:
        None
    """
    L = int(Potential_grid.shape[0])
    if plot_type == '2D':
        plt.imshow(Potential_grid, cmap='inferno')
        plt.show()
    elif plot_type == '3D':
        X = np.linspace(0, L, L)
        Y = np.linspace(0, L, L)
        X, Y = np.meshgrid(X, Y)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, Potential_grid, cmap='inferno')
        plt.show()

# ---------- Timing function ---------- #
def time_function(f):
    """
    Times the execution of a function.

    Args:
        f (function): The function to be timed.

    Returns:
        The function output and the time it took to execute.
    """

    # Starting the timer
    start = time.time()
    # Executing the function
    f()
    # Stopping the timer
    end = time.time()
    # Calculating the time difference
    time_difference = end - start

    return time_difference

# ---------- Main ---------- #
def main():
    """
    Plots the potential grid for a given L.
    """
    # Define grid parameters:
    L = 100

    Potential_grid = build(L)

    # Plot the results
    plot(Potential_grid, plot_type='3D')

# ---------- Main for timing ---------- #
def main_time():
    """
    Measures the time it takes to run the build function for different L.
    
    L=1 is used to precompile the code.
    
    Args:
        L_list (list): A list of L values.
    """
    # Define grid parameters:
    L_list = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # Precompile the code
    build(1)

    # Time the function
    for L in L_list:
        print(f'Time for L = {L} was: {time_function(lambda: build(L))} seconds.')

# -------------------------- #

if __name__ == "__main__":
    main()
    main_time()
