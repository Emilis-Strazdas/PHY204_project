# PHY204_project

 Iterative way to solve the Laplasian equation for $\Delta V = 0$, knowing the boundary conditions. The code utilises the Mean-Value Theorem for the electric potential and a Jaccobi iteration method. It is split into three main parts:

* iteration.py file contains the function that iterates over the grid until the error is below a certain tolerance;
* iteration_step.py file contains the function that performs one iteration step;
* main.py file contains the main function that runs the code.

There is also plot.py for plotting the function in either 2D or 3D, but has nothing to do with calculating the grid.

The code has been written by **Emilis Strazdas** and **Leopold Rousseau**.
