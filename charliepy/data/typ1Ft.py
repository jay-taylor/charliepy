########################################################
# Data for type F~_4 Coxeter groups
#
#              0   1   2   3   4
#              o---o---o=>=o---o
#

import numpy as np

def cartanmat(n):
    if n == 5:
        return np.array([[ 2, -1,  0,  0,  0],
                         [-1,  2, -1,  0,  0],
                         [ 0, -1,  2, -2,  0],
                         [ 0,  0, -1,  2, -1],
                         [ 0,  0,  0, -1,  2]],
                         dtype = 'int8')
    else:
        raise ValueError("Rank of F~_4 must be 5!")

