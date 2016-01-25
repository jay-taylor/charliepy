########################################################
# Data for type C~_n Coxeter groups
# The Dynkin diagram here follows Bourbaki so we have
#
#          0   1   2        n-2 n-1  n
#          o=>=o---o-- ... --o---o=<=o
#

from .. import utils as utils

import numpy as np

def cartanmat(n):
    """Takes a zero matrix and returns the type C_n Cartan matrix as a
    NumPy array."""
    if n < 3:
        raise ValueError("The rank of C~_n must be at least 3.")
    #if n == 3:
    #    return np.array([[2, -2, -2], [-1, 2, 0], [-1, 0, 2]], dtype = 'int8')
    #else:
    C = np.zeros((n,n), dtype='int8')
    np.fill_diagonal(C, 2)
    C[range(1,n-1), range(2,n)] = [-1]*(n-2)
    C[range(2,n), range(1,n-1)] = [-1]*(n-2)
    C[n-1][n-2] = -2
    C[0][1], C[1][0] = -2, -1
    return C

