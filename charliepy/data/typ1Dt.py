########################################################
# Data for type D~_n Coxeter groups
# The Dynkin diagram here follows Bourbaki so we have.
#
#           1 o                         o n-1
#               \     3        n-3    /
#               2 o---o-- ... --o---o n-2
#               /                     \
#           0 o                         o n

from .. import utils as utils

import numpy as np

def cartanmat(n):
    """Takes a zero matrix and returns the type D~n Cartan matrix as a
    NumPy array."""
    if n < 5:
        raise ValueError("The rank of D~_n must be at least 5.")
    else:
        C = np.zeros((n,n), dtype='int8')
        np.fill_diagonal(C, 2)
        C[range(1, n-2), range(2, n-1)] = [-1]*(n-3)
        C[range(2, n-1), range(1, n-2)] = [-1]*(n-3)
        C[n-1][n-3], C[n-3][n-1] = -1, -1
        #C[1:4,1:4] = [[2, 0, -1], [0, 2, -1], [-1, -1, 2]]
        C[0][2], C[2][0] = -1, -1
        return C

