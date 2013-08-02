import numpy as np


def gauss(x, SD):
    """A simple gaussian function with a mean of 0
    """
    return 1.0 * np.exp(-(x) ** 2 / (2 * SD ** 2))

