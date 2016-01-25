########################################################
# Data for type A~_n Coxeter groups

from .. import utils as utils
from . import typ1A

import numpy as np

def cartanmat(n):
    """Takes a zero matrix and returns the type A~_n Cartan matrix as a
    NumPy array."""
    if n == 2:
        return np.array([[2, -2], [-2, 2]], dtype = 'int8')
    elif n <= 1:
        raise ValueError("The rank of A~_n must be at least 2.")
    else:
        C = typ1A.cartanmat(n)
        C[n-1][0], C[0][n-1] = -1, -1
        return C

