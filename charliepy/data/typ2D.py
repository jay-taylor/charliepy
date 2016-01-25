########################################################
# Data for type 2D_n Coxeter groups and reductive groups
#
# As in main.chartableD we label the generators of B_n
# as 0, ..., n-1 and the generators of D_n as
# 0', ..., (n-1)' such that we have i' = i if
# i > 1 and 0' = 1, 1' = 010.

import itertools

from .. import utils as utils
from . import typ1B

def conjclassdata(ind, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type D_n and phi is the unique graph automorphism
    of order 2. The data is adapted from the corresponding files in GAP-Chevie.
    """
    # Stores the data: representatives, centraliser orders, names
    repcentnam = [[],[],[]]
    n = len(ind)

    # We construct the reps w(alpha,beta).0 as in [GKP, 4.2] and [GP, 3.4.7].
    for mu in utils.partitiontuplesgen(n,2):
        m = 0
        word = list()
        if len(mu[1])%2:
            # mu[1] <-> neg blocks and mu[0] <-> pos blocks
            # so need mu[1] increasing and mu[0] decreasing
            for d in mu[1][::-1]:
                if m == 0:
                    # b^-(0,d).0 = 0...(d-1)0 = 1'...(d-1)' 
                    word += ind[1:m+d]
                else:
                    # b^-(m,d).0 = m...101...(m+d-1)0 = m'...2'0'1'...(m+d-1)'
                    # kill all 0's as even no. and 0's commute with blocks
                    word += ind[m:1:-1] + ind[:m+d]
                m += d

            for d in mu[0]:
                # b^+(m,d) = (m+1)...(m+d-1)
                word += ind[m+1:m+d]
                m += d

            repcentnam[0].append(word)
            # phi-cent order is order of cent of w(alpha,beta) in W_n' 
            repcentnam[1].append(utils.centralisertuple(n,2,mu)/2)
            repcentnam[2].append(utils.parttupletostring(mu))
        else:
            continue

    return repcentnam

def irrchardata(n, **kwargs):
    return [[utils.parttupletostring(x)
             for x in filter(testprefext, utils.partitiontuplesgen(n, 2))]]

def chartable(n, **kwargs):
    # Get the indices of the characters labelling the preferred extensions and
    # the non-degenerate conjugacy classes of type B_n.
    exts = [x[0] for x in itertools.ifilter(testprefexttup,
                            enumerate(utils.partitiontuplesgen(n, 2)))]
    classes = [item[0] for item in itertools.ifilter(
            lambda x: len(x[1][1])%2, enumerate(utils.partitiontuplesgen(n,2)))]
    
    # Return the appropriate slice of the character table of type B_n.
    return typ1B.chartable(n)[exts,:][:,classes]

#########################################
## Utility Functions
##

def testprefexttup(tup, good=False):
    return testprefext(tup[1], good)

def testprefext(bipart, good=False):
    """Tests if a bipartition represents the preferred extension of an
    irreducible character of type D_n or not. If the optional flag good is set
    to True then tests if the bipartition is Shoji's good extension.
    """
    # Don't want to shred the input because we're going to print it.
    dupe = [bipart[0][:], bipart[1][:]]
    dlen = len(dupe[0]) - len(dupe[1])
    dupe[1] += [0]*dlen
    dupe[0] += [0]*(-dlen)

    symb = [[val+i for i, val in enumerate(reversed(dupe[0]))],
            [val+i for i, val in enumerate(reversed(dupe[1]))]]
    
    if good:
        return symb[0] < symb[1]
    else:
        return symb[1] < symb[0]

