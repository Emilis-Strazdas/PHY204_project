"""
This file contains the plot function, which plots the potential grid.
"""
import numpy as np
import matplotlib.pyplot as plt

def plot(Potential_grid, grid_parameters, plot_type='2D'):
    """
    Plots the potential grid.

    The plot type can be either 2D or 3D.
    
    Args:
        Potential_grid (numpy.ndarray): The potential grid.
        grid_parameters (dict): A dictionary containing the grid parameters.
        
    Returns:
        None
        
    """

    # Extract the size of the grid
    L = grid_parameters['size']

    if plot_type == '2D':
        plt.imshow(Potential_grid, cmap='inferno')
        plt.show()

    elif plot_type == '3D':
        fig = plt.figure()
        
        plt.plot_surface(np.arange(0, 2*L), np.arange(0, 2*L), Potential_grid, cmap='inferno')
        plt.show()

    else:
        print('#ERROR: Invalid plot type.')