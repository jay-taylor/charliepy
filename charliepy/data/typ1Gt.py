########################################################
# Data for type G~_2 Coxeter groups
# The Dynkin diagram here follows Bourbaki so we have
#
#                       1   2   0
#                       o-<-o---o
#                         6

import numpy as np

def cartanmat(n):
    if n == 3:
        return np.array([[ 2,  0, -1],
                         [ 0,  2, -1],
                         [-1, -3,  2]],
                         dtype = 'int8') 
    else:
        raise ValueError("Rank of G~_2 must be 3!")

