########################################################
# Data for type G~_2 Coxeter groups
#

import numpy as np

def cartanmat(n):
    if n == 3:
        return np.array([[ 2, -1,  0],
                         [-1,  2, -1],
                         [ 0, -3,  2]],
                         dtype = 'int8') 
    else:
        raise ValueError("Rank of G~_2 must be 3!")

