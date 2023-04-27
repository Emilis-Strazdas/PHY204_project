#==========   PHY205 Project      ==========
# Date Created: 2023-04-26  
# Authors:      E. Strazdas
#               L. Rousseau
#=========================================

"""
    This code is used to solve the Laplace equation with fixed boundary
    conditions using the Jacobi method. For performance reasons, we can
    exploit the fact that the potential grid is symmetric about the anti-
    diagonal. This allows us to store only half of the grid, which reduces
    the memory usage and the number of operations by a factor of 2. Furthermore,
    we precompile the code using the Numba package, which allows us to use
    just-in-time compilation to speed up the code for testing (this should
    not be against the rules since C++ programs are also precompiled before
    being executed.
"""

# ---------- Imports ---------- #
from plot import plot
from time_function import time_function
from build import build

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
