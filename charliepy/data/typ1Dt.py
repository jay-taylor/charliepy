########################################################
# Data for type D~_n Coxeter groups

from .. import utils as utils

import numpy as np

def cartanmat(n):
    """Takes a zero matrix and returns the type A_n Cartan matrix as a
    NumPy array."""
    if n < 5:
        raise ValueError("The rank of D~_n must be at least 2.")
    else:
        C = np.zeros((n,n), dtype='int8')
        np.fill_diagonal(C, 2)
        C[range(1,n-1), range(2,n)] = [-1]*(n-2)
        C[range(2,n), range(1,n-1)] = [-1]*(n-2)
        C[1:4,1:4] = [[2, 0, -1], [0, 2, -1], [-1, -1, 2]]
        C[0][n-2], C[n-2][0] = -1, -1
        return C

