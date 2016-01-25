########################################################
# Data for type B~_n Coxeter groups

from .. import utils as utils
from . import typ1A

import numpy as np

def cartanmat(n):
    """Takes a zero matrix and returns the type B~_n Cartan matrix as a
    NumPy array."""
    if n < 3:
        raise ValueError("The rank of B~_n must be at least 3.")
    if n == 3:
        return np.array([[2, -1, 0], [-2, 2, -2], [0, -1, 2]], dtype = 'int8')
    else:
        C = np.zeros((n,n), dtype='int8')
        np.fill_diagonal(C, 2)
        C[range(1,n-1), range(2,n)] = [-1]*(n-2)
        C[range(2,n), range(1,n-1)] = [-1]*(n-2)
        C[1][2] = -2
        C[0][n-2], C[n-2][0] = -1, -1
        return C
