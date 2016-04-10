########################################################
# Data for type 2D_n Coxeter groups and reductive groups
#
# As in main.chartableD we label the generators of B_n
# as 0, ..., n-1 and the generators of D_n as
# 0', ..., (n-1)' such that we have i' = i if
# i > 1 and 0' = 1, 1' = 010.


from .. import utils
from . import typ1B

import itertools

def conjclassdata(ind, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type D_n and phi is the unique graph automorphism
    of order 2. The data is adapted from the corresponding files in GAP-Chevie.
    """
    n = len(ind)

    # We construct the reps w(alpha,beta).0 as in [GKP, 4.2] and [GP, 3.4.7].
    for mu in typ1B._conjlabels(n):
        m = 0
        word = list()
        if len(mu[1]) % 2:
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

            yield (
                word,
                utils.centralisertuple(n, 2, mu)//2,
                utils.parttupletostring(mu)
            )

        else:
            continue

def conjclassdata_min(ind, **kwargs):
    n = len(ind)

    for mu in typ1B._conjlabels(n):
        if len(mu[1]) % 2:
            yield (typ1B.centraliser(mu)//2, utils.parttupletostring(mu))
        else:
            continue

def irrchardata(n, **kwargs):
    return ((utils.parttupletostring(x), None, None)
            for x in filter(testprefext, typ1B._charlabels(n)))

def chartable(n, **kwargs):
    # Get the indices of the characters labelling the preferred extensions and
    # the non-degenerate conjugacy classes of type B_n.
    chars = [i for i, x in enumerate(typ1B._charlabels(n)) if testprefext(x)]
    classes = [i for i, x in enumerate(typ1B._conjlabels(n))
        if len(x[1]) % 2]
    
    # Return the appropriate slice of the character table of type B_n.
    ctB = typ1B.chartable(n)
    return [[ctB[i][j] for j in classes] for i in chars]

#########################################
## Utility Functions
##

def testprefext(bipart):
    """
    Returns True if bipart is a bipartition representing the preferred
    extension of an irreducible character of type D_n and False
    otherwise.

    """
    # Assume (alpha; beta) is the bipartition we are considering with
    #
    #       alpha = (alpha_1, ..., alpha_s)
    #       beta = (beta_1, ..., beta_t)
    #
    # If t < s then the reduced symbol associated to (alpha; beta) will
    # have a 0 in the lower row but no 0 in the upper row. Thus this
    # will label the preferred extension. Conversely if s < t then this
    # will not label the preferred extension. Finally if s == t then the
    # bipartition labels the preferred extension if and only if beta <
    # alpha in the lexicographic ordering.
    if len(bipart[0]) != len(bipart[1]):
        return len(bipart[1]) < len(bipart[0])
    else:
        return bipart[1] < bipart[0]

