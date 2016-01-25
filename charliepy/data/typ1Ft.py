########################################################
# Data for type F~_4 Coxeter groups
#

import numpy as np

def cartanmat(n):
    if n == 5:
        return np.array([[ 2, -1,  0,  0,  0],
                         [-1,  2, -1,  0,  0],
                         [ 0, -1,  2, -1,  0],
                         [ 0,  0, -2,  2, -1],
                         [ 0,  0,  0, -1,  2]],
                         dtype = 'int8')
    else:
        raise ValueError("Rank of F~_4 must be 5!")

