"""
This file contains the plot function, which plots the potential grid.
"""

# ---------- Imports ---------- #
import numpy as np
import matplotlib.pyplot as plt

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
