########################################################
# Data for type 3D_4 Coxeter groups and reductive groups
#
# As in the main part of PyCox we assume that the Dynkin diagram is
# labelled as follows:
#
#               1
#               |
#           0 - 2 - 3
#
# Recall that we have the following two possibilities for the
# automorphism phi of order 3:
#
#           phi:  0 -> 1 -> 3 -> 0
#           phi': 0 -> 3 -> 1 -> 0
#
# The algorithms below adjust the data based on both the labelling of
# the Dynkin diagram and whether we're dealing with the automorphism phi
# or phi'.

import numpy as np

def conjclassdata(ind, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type D_4 and phi is a graph automorphism of order
    3. The data is adapted from the corresponding files in GAP-Chevie.

    """
    phi = kwargs['phi']
    
    # Stores the data: representatives, centraliser orders, names.
    repcentnam = list() 

    # The first instance assumes we're dealing with the automorphism phi
    # and the second case assumes we're dealing with the automorphism
    # phi'. Passing from one case to the other is simply a matter of
    # interchanging the roles of 1 and 3.
    if phi[ind[0]] == ind[1]:
        reps = iter([[0], [], [0, 1, 2, 0, 1, 2], [2], [0, 2],
                     [0, 1, 2, 0, 1, 3, 2, 1], [0, 1, 2, 1]])
    elif phi[ind[0]] == ind[3]:
        reps = iter([[0], [], [0, 3, 2, 0, 3, 2], [2], [0, 2],
                     [0, 3, 2, 0, 3, 1, 2, 3], [0, 3, 2, 3]])

    # Relabel the reps based on the indices passed
    repcentnam.append([[ind[j] for j in rep] for rep in reps])

    repcentnam.append([4, 12, 12, 4, 4, 24, 24])

    repcentnam.append(["C_3", "\\tilde A_2", "C_3+A_1", "\\tilde A_2+A_1",
                       "F_4", "\\tilde A_2+A_2","F_4(a_1)"])

    return repcentnam

def irrchardata(n, **kwargs):
    return [[".4", ".1111", ".22", "11.2", "1.3", "1.111", "1.21"]]

def chartable(n, **kwargs):
    # Here we use Lusztig's preferred extension. This is the unique
    # extension which is defined over the rationals.
    return np.array([[  1, 1,  1,  1,  1,  1,  1 ],
                     [ -1, 1,  1, -1,  1,  1,  1 ],
                     [  0, 2,  2,  0, -1, -1, -1 ],
                     [  0, 0,  0,  0, -1,  3,  3 ],
                     [  1, 1, -1, -1,  0, -2,  2 ],
                     [ -1, 1, -1,  1,  0, -2,  2 ],
                     [  0, 2, -2,  0,  0,  2, -2 ]],
                    dtype='int')

