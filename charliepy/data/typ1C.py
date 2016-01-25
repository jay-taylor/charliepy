########################################################
# Data for type C_n Coxeter groups and reductive groups.
# The Dynkin diagram here follows Bourbaki so we have
#
#            0   1                n-1
#            o---o-- ... --o---o=<=o

from .. import utils as utils
from . import typ1A as typA
from . import typ1B as typB

def cartanmat(n):
    if n < 1:
        raise ValueError("Rank of C_n must be at least 1!")

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

def degrees(n):
    return list(range(2, 2*n+1, 2))

def conjclassdata(ind, **kwargs):
    # stores the data: representatives, centraliser orders, names
    n = len(ind)

    # If n == 2 we permute the labellings of the conjugacy classes. This means
    # that the computation of the character table stays the same for both B and
    # C.
    if n == 2:
        return typB.conjclassdata(ind[::-1])
    else:
        return typB.conjclassdata(ind)

def irrchardata(n, **kwargs):
    return typB.irrchardata(n)

def chartable(n, **kwargs):
    return typB.chartable(n)
