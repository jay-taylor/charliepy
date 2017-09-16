########################################################
# Data for type C_n Coxeter groups and reductive groups.
# The Dynkin diagram here follows Bourbaki so we have
#
#            0   1                n-1
#            o---o-- ... --o---o=<=o

from .. import utils
from . import typ1A as typA
from . import typ1B as typB

def cartanmat(n):
    if n < 1:
        raise ValueError("Rank of C_n must be at least 1!")
    elif n == 1:
        return np.array([2], dtype = 'int8')

    cmat = typA.cartanmat(n)
    cmat[n-1][n-2] = -2
    return cmat

def diagram(inds):
    """Prints the Dynkin diagram."""
    out = "C{} :  ".format(len(inds))
    out += " - ".join(str(x) for x in inds[:-1])
    out += " < {}".format(inds[-1])
    print(out)
    return None

def rootlengths(n, **kwargs):
    """
    Returns a generator giving the relative root lengths.

    """
    if n == 1:
        yield 1
    else:
        for i in range(n-1):
            yield 1
        yield 2

def degrees(n):
    return list(range(2, 2*n+1, 2))

def longestword(inds):
    return (inds[0::2] + inds[1::2])*len(inds)

def maxparachain(inds):
    return inds

def conjclasses(ind, **kwargs):
    n = len(ind)

    # If n == 2 we permute the labellings of the conjugacy classes. This means
    # that the computation of the character table stays the same for both B and
    # C.
    if n == 2:
        return typB.conjclasses(ind[::-1])
    else:
        return typB.conjclasses(ind)

# These are all identical to the type B case.
wordtoclass = typB.wordtoclass
conjclasses_min = typB.conjclasses_min
irrchars = typB.irrchars
chartable = typB.chartable
