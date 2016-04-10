########################################################
# Data for type B~_n Coxeter groups.
# The Dynkin diagram here follows Bourbaki so we have
#
#             0   2   1
#       B3 :  o=>=o=<=o
#
# and for n >= 4 we have the diagram is given by
#
#          1 o
#              \          n-1  n
#             2  o-- ... --o=>=o
#              /
#          0 o
#

from .. import utils
from . import typ1A

import numpy as np

def cartanmat(n):
    """Takes a zero matrix and returns the type B~_n Cartan matrix as a
    NumPy array."""
    if n < 3:
        raise ValueError("The rank of B~_n must be at least 3.")
    if n == 3:
        return np.array([[2, 0, -2], [0, 2, -2], [-1, -1, 2]], dtype = 'int8')
        #return np.array([[2, -1, 0], [-2, 2, -2], [0, -1, 2]], dtype = 'int8')
    else:
        C = np.zeros((n,n), dtype='int8')
        np.fill_diagonal(C, 2)
        C[range(1,n-1), range(2,n)] = [-1]*(n-2)
        C[range(2,n), range(1,n-1)] = [-1]*(n-2)
        C[n-2][n-1] = -2
        C[0][2], C[2][0] = -1, -1
        return C
