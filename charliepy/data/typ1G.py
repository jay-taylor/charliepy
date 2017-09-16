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
from .. import permutat

# We use the standard permutation generators below in class recognition.
_G2permgens = [permutat.cyclestoperm(*_rd.G2permgens[0]),
               permutat.cyclestoperm(*_rd.G2permgens[1])]

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

def maxparachain(inds):
    return inds

def conjclasses(inds, **kwargs):
    return (
        (nam, cent, [inds[i] for i in rep])
        for (nam, cent, rep) in _rd.G2conjclasses
    )

def conjclasses_min(inds, **kwargs):
    return (cls[:2] for cls in _rd.G2conjclasses)

def wordtoclass(n, w):
    """
    Returns the name of the conjugacy class in the Weyl group of type G_2
    containing the element w given as a word in the standard generators.

    """
    # First convert the word to a permutation.
    perm = permutat.Perm(range(12))
    for i in w:
        perm *= _G2permgens[i]

    # Each class is uniquely determined by the pair of cycletypes of the
    # restricted permutation on each orbit.
    param = (permutat.restrictedperm(perm, _rd.G2rootorbits[0]).cyclestructure,
             permutat.restrictedperm(perm, _rd.G2rootorbits[1]).cyclestructure)

    return _rd.G2classparams[param]

def irrchars(n, **kwargs):
    return (
        ('phi{{{d},{b}}}{p}'.format(d=deg, b=b, p=prim), a, b)
        for (deg, a, b, prim) in _rd.G2irrchardata
    )

def chartable(n, **kwargs):
    return _rd.G2chartable

def coxeterclasses(inds):
    return (
        (nam, orblen, norm, [inds[i] for i in rep])
        for (nam, orblen, norm, rep) in _rd.G2coxeterclasses
    )










