########################################################
# Data for type 2A_n Coxeter groups and reductive groups

from .. import permutat
from .. import utils
from . import typ1A

import numpy as np

def conjclasses(ind, **kwargs):
    """ Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type A_n and phi is the unique graph automorphism
    of order 2. The data is adapted from the corresponding files in GAP-Chevie.
    """
    n = len(ind)

    # This is an auxiliary function to determine a reduced word for a
    # permutation in the standard Coxeter generators of S_(n+1). It is a
    # classic bubble sort algorithm and is an adaptation of [GP00, pg.
    # 9, Algorithm A]. The difference is our algorithm multiplies the
    # word by generators on the right and not on the left. This has the
    # side effect that the word is reduced but not lexicographically
    # smallest.
    def redword(w):
        s = 0
        lim = w.degree - 1
        word = []
        while s < lim:
            if s^w > (s+1)^w:
                word.append(ind[s])
                w = (s, s+1)*w
                if s != 0:
                    s -= 1
                else:
                    s += 1
            else: 
                s += 1
  
        return word
  
    # The longest word in S_(n+1) as a permutation 
    w0 = permutat.Perm(range(n,-1,-1))
  
    # Construct helper list for constructing elements in stair case
    # form. Note that the indexing must be changed to reflect the fact
    # that we index from 0. In particular we have helper[i] is
    # a_(i+1) - 1 in the notation of [GKP] (for all i in range(n+1)),
    helper = [0]*(n+1)
    helper[::2] = range((n+1)//2 + (n+1) % 2)
    helper[1::2] = [n - i for i in range((n+1)//2)]
  
    for mu in typ1A._conjlabels(n+1):
        # Make mu a maximal composition of n+1 following [GKP] then
        # construct an element in stair case form as in [GKP, 3.1]. We
        # assume that the partitions provided are weakly decreasing
        # (which they are).
        mu = ([x for x in mu if x%2 == 0] + [x for x in mu if x%2 == 1])

        lmu = len(mu)
  
        # Create the permutation in stair form.
        w = [None]*lmu 
        j = 0
        for i in range(lmu):
            w[i] = helper[j:j+mu[i]]
            j += mu[i]

        # We work with w0*w to compensate for the fact that our reduced
        # word algorithm multiplies on the right. This means we get a
        # lexicographically smallest reduced word. Note that w0*w and
        # w*w0 are conjuagte by w0, hence they have the same cycle type.
        w = w0*permutat.cyclestoperm(*w)

        # Note that cycletype comes back as a sorted partition.
        yield (
            utils.intlisttostring(mu),
            typ1A.centraliser(w.cycletype),
            redword(w)
        )

def conjclasses_min(ind, **kwargs):
    """ Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type A_n and phi is the unique graph automorphism
    of order 2. The data is adapted from the corresponding files in GAP-Chevie.
    """
    n = len(ind)
  
    # The longest word in S_(n+1) as a permutation 
    w0 = permutat.Perm(range(n,-1,-1))
  
    # Construct helper list for constructing elements in stair case
    # form. Note that the indexing must be changed to reflect the fact
    # that we index from 0. In particular we have helper[i] is
    # a_(i+1) - 1 in the notation of [GKP] (for all i in range(n+1)),
    helper = [0]*(n+1)
    helper[::2] = range((n+1)//2 + (n+1) % 2)
    helper[1::2] = [n - i for i in range((n+1)//2)]
  
    for mu in typ1A._conjlabels(n+1):
        # Make mu a maximal composition of n+1 following [GKP] then
        # construct an element in stair case form as in [GKP, 3.1]. We
        # assume that the partitions provided are weakly decreasing
        # (which they are).
        mu = ([x for x in mu if x%2 == 0] + [x for x in mu if x%2 == 1])

        lmu = len(mu)
  
        # Create the permutation in stair form.
        w = [None]*lmu 
        j = 0
        for i in range(lmu):
            w[i] = helper[j:j+mu[i]]
            j += mu[i]

        # We work with w0*w to compensate for the fact that our reduced
        # word algorithm multiplies on the right. This means we get a
        # lexicographically smallest reduced word. Note that w0*w and
        # w*w0 are conjuagte by w0, hence they have the same cycle type.
        w = w0*permutat.cyclestoperm(*w)

        # Note that cycletype comes back as a sorted partition.
        yield (
            utils.intlisttostring(mu),
            typ1A.centraliser(w.cycletype)
        )


def irrchars(n, **kwargs):
    return ((utils.intlisttostring(mu), None, None)
            for mu in typ1A._charlabels(n+1))

def chartable(n, **kwargs):
    # Here we use Lusztig's preferred extension defined in [L85, IV,
    # 17.2]. In particular w.phi acts on the representation E as
    # w(-1)^a_E. This involves multiplying all rows of the character
    # table by appropriate signs.
    chartabsplit = typ1A.chartable(n)
    signs = [(-1)**sum(i*val for i, val in enumerate(mu)) 
             for mu in typ1A._charlabels(n + 1)]
    return [[s*x for x in row] for s, row in zip(signs, chartabsplit)]
