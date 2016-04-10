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
from . import rawdata as _rd

def cartanmat(n):
    if n == 2:
        return np.array([[ 2, -1],
                         [-3,  2]],
                        dtype = 'int8') 
    else:
        raise ValueError("Rank of G_2 must be 2!")

def rootlengths(n, **kwargs):
    """
    Returns a tuple giving the relative root lengths.

    """
    return (1, 3)

def diagram(inds):
    """Prints the Dynkin diagram."""
    print("G2 :  {} < {}".format(*inds))
    return None

def degrees(n):
    return [2, 6]

def longestword(inds):
    return [inds[0], inds[1]]*3

def conjclasses(inds, **kwargs):
    return (
        (nam, cent, [inds[i] for i in rep])
        for (nam, cent, rep) in _rd.G2conjclasses
    )

def conjclasses_min(inds, **kwargs):
    return (cls[:2] for cls in _rd.G2conjclasses)

def irrchars(n, **kwargs):
    return (
        ('phi{{{d},{b}}}{p}'.format(d=deg, b=b, p=prim), a, b)
        for (deg, a, b, prim) in _rd.G2irrchardata
    )

def chartable(n, **kwargs):
    return _rd.G2chartable

def coxeterclasses(inds):
    return (
        (nam, orblen, [inds[i] for i in rep])
        for (nam, orblen, rep) in _rd.G2coxeterclasses
    )










