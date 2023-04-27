"""
This file contains the time_function function, which times the execution of a
function.
"""

# ---------- Imports ---------- #
import time

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
