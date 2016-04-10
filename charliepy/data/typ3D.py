########################################################
# Data for type 3D_4 Coxeter groups and reductive groups
#
# As in the main part of PyCox we assume that the Dynkin diagram is
# labelled as follows:
#
#                               2
#                             /
#                       0 - 1
#                             \
#                               3
#
# Recall that we have the following two possibilities for the
# automorphism phi of order 3:
#
#           phi:  0 -> 2 -> 3 -> 0
#           phi': 0 -> 3 -> 2 -> 0
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
    
    # The first instance assumes we're dealing with the automorphism phi
    # and the second case assumes we're dealing with the automorphism
    # phi'. Passing from one case to the other is simply a matter of
    # interchanging the roles of 1 and 3.
    if phi[ind[0]] == ind[2]:
        classes = (
            ([0], 4, 'C3'), ([], 12, '~A2'), ([0, 2, 1, 0, 2, 1], 12, 'C3+A1'),
            ([1], 4, '~A2+A1'), ([0, 1], 4, 'F4'),
            ([0, 2, 1, 0, 2, 3, 1, 2], 24, '~A2+A2'),
            ([0, 2, 1, 2], 24, 'F4(a1)')
        )
    else:
        classes = (
            ([0], 4, 'C3'), ([], 12, '~A2'), ([0, 3, 1, 0, 3, 1], 12, 'C3+A1'),
            ([1], 4, '~A2+A1'), ([0, 1], 4, 'F4'),
            ([0, 3, 1, 0, 3, 2, 1, 3], 24, '~A2+A2'),
            ([0, 3, 1, 3], 24, 'F4(a1)')
        )

    return (([ind[j] for j in cls[0]], cls[1], cls[2]) for cls in classes)

def conjclassdata_min(ind, **kwargs):
    classes = (
        (4, 'C3'), (12, '~A2'), (12, 'C3+A1'), (4, '~A2+A1'), (4, 'F4'),
        (24, '~A2+A2'), (24, 'F4(a1)')
    )

    return (cls for cls in classes)

def irrchardata(n, **kwargs):
    charnames = (".4", ".1111", ".22", "11.2", "1.3", "1.111", "1.21")
    return ((nam, None, None) for nam in charnames)

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

