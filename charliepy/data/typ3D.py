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

from . import rawdata as _rd
from .. import permutat
import numpy as np

def conjclasses(inds, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type D_4 and phi is a graph automorphism of order
    3. The data is adapted from the corresponding files in GAP-Chevie.

    """
    phi = kwargs['phi']

    if inds[0]^phi == inds[2]:
        p = permutat.Perm([0, 1, 2, 3])
    else:
        p = permutat.Perm([0, 1, 3, 2])

    return (
        (nam, cent, [inds[i^p] for i in rep])
        for (nam, cent, rep) in _rd.D4_3conjclasses
    )


def conjclasses_min(inds, **kwargs):
    return (cls[:2] for cls in _rd.D4_3conjclasses)

def irrchars(n, labels=(), **kwargs):
    return _rd.D4_3irrchardata

def chartable(n, **kwargs):
    # Here we use Lusztig's preferred extension. This is the unique
    # extension which is defined over the rationals.
    return _rd.D4_3chartable

