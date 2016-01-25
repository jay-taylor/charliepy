########################################################
# Data for type G_2 Coxeter groups and reductive groups.
# The Dynkin diagram here follows Bourbaki so we have
#
#                       0   1
#                       o-<-o
#                         6
#

import numpy as np
import math

def cartanmat(n):
    if n == 2:
        return np.array([[ 2, -1],
                        [-3,  2]],
                        dtype = 'int8') 
    else:
        raise ValueError("Rank of G_2 must be 2!")

def diagram(inds):
    """Prints the Dynkin diagram."""
    print("G2 {} < {}".format(*inds))
    return None

def degrees(n):
    return [2, 6]

def conjclassdata(ind, **kwargs):
    # stores the data: representatives, centraliser orders, names
    repcentnam = list() 

    # Stored as an iterable because we're about to loop through it.
    reps = iter([[], [0], [1], [0, 1], [0, 1, 0, 1], [0, 1, 0, 1, 0, 1]])

    # create words based on indices of the irreducible component
    repcentnam.append([[ind[j] for j in rep] for rep in reps])

    repcentnam.append([12, 4, 4, 6, 6, 12])

    repcentnam.append(['A_0', '~A_1', 'A_1', 'G_2', 'A_2', 'A_1+~A_1'])

    return repcentnam

def irrchardata(n, **kwargs):
    binv = [0, 6, 3, 3, 1, 2]
    ainv = [0, 6, 1, 1, 1, 1]
    nam = ["phi_{1,0}", "phi_{1,6}", "phi_{1,3}'", "phi_{1,3}''",
           "phi_{2,1}", "phi_{2,2}"]

    return [nam, ainv, binv]

def chartable(n, **kwargs):
    return np.array([[1,  1,  1,  1,  1,  1], 
                     [1, -1, -1,  1,  1,  1], 
                     [1,  1, -1, -1,  1, -1], 
                     [1, -1,  1, -1,  1, -1], 
                     [2,  0,  0,  1, -1, -2], 
                     [2,  0,  0, -1, -1,  2]],
                     dtype='int')

