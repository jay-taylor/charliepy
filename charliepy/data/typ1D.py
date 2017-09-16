########################################################
# Data for type D_n Coxeter groups and reductive groups
# The Dynkin diagram here follows Bourbaki so we have.
#
#                                  o n-2
#            0   1        n-4    /
#            o---o-- ... --o---o n-3
#                                \
#                                  o n-1
#
# Note that the labelling is consistent when read from the right. This
# means the digram of D_3 is given by
#
#                    1   0   2
#                    o---o---o
#
# Now assume W is a Weyl group of type D_n and W' is a 
# Weyl group of type B_n with Dynkin diagram
#
#            0'  1'              (n-1)'
#            o---o-- ... --o---o=>=o
#
# Then we embed W -> W' such that i = i' for any
# 0 <= i < n-1 and (n-1) = (n-1)'(n-2)'(n-1)'.

from .. import utils
from .. import permutat
from . import typ1A as typA
from . import typ1B as typB
from . import _cdata

import numpy as np
import itertools

########################################################################
########################################################################
##                                                                    ##
##                          Cartan Matrix                             ##
##                                                                    ##

def cartanmat(n):
    if n < 2:
        raise ValueError("Rank of D_n must be at least 2!")

    if n == 2:
        return np.array([[2, 0],
                         [0, 2]], dtype = 'int8')
    
    C = typA.cartanmat(n)
    C[-3:,-3:] = [[2, -1, -1], [-1, 2, 0], [-1, 0, 2]]
    return C

def rootlengths(n, **kwargs):
    """
    Returns a generator giving the relative root lengths.

    """
    return (1 for _ in range(n))

def diagram(inds):
    """Prints the Dynkin diagram."""
    n = len(inds)
    if n == 2:
        print("D2 :  {}   {}".format(*inds))
        return None
    elif n == 3:
        print("D3 :  {} - {} - {}".format(inds[1], inds[0], inds[2]))
        return None

    out3 = "D{} :  ".format(len(inds))
    out3 += " - ".join(str(x) for x in inds[:-2])
    out2 = " "*len(out3) + " /"
    out4 = " "*len(out3) + " \\"
    out1 = " "*len(out2) + " {}".format(inds[-2])
    out5 = " "*len(out2) + " {}".format(inds[-1])
    print(out1, out2, out3, out4, out5, sep="\n")
    return None

def degrees(n):
    return list(range(2, 2*n-1, 2)) + [n]

def longestword(inds):
    n = len(inds)
    L = [inds[0:n-1:2], inds[1:n-1:2]]
    L[n % 2].append(inds[-1])
    return (L[0] + L[1])*(n - 1)

def maxparachain(inds):
    n = len(inds)
    
    if n > 4:
        return inds
    elif n == 4:
        return [inds[i] for i in (0, 2, 1, 3)]
    elif n == 3:
        return [inds[i] for i in (1, 0, 2)]
    else:
        return inds

########################################################################
########################################################################
##                                                                    ##
##                        Conjugacy Classes                           ##
##                                                                    ##

def _evenpartitions(n):
    # For an even number n this produces all the even partitions of n,
    # i.e., all partitions where every part is divisible by 2. The
    # ordering of the partitions is the same as that produced by
    # typA._conjlabels(n); the algorithm is essentially the same.
    # However now we start with a list containing only 2s instead of
    # only 1s. We also have:
    #
    #   i - index of smallest entry != 2
    #   m - length of list.

    part = [2 for j in range(n//2)]
    i, m = -1, n//2
    
    while True:
        yield part[:m]

        # Easy case [2, 2] -> [4].
        if 2 <= m - i - 1:
            i += 1
            m -= 1
            part[i] = 4
            continue

        if m == 1:
            return
        
        # Do the replacement [x, ..., x, r] -> [(x+2), 2, ..., 2].
        # First find the index of the first copy of x; note this will
        # become the index of the largest entry not equal to 2.
        i = m-3
        while 0 <= i:
            if part[i] == part[i+1]:
                i -= 1
            else:
                break
        i += 1
        
        # Now set the right number of 2s. First correcting the length.
        m = i + (m-2-i)*part[m-2]//2 + part[m-1]//2
        j = i + 1
        while j < m:
            part[j] = 2
            j += 1
        part[i] += 2

def _conjlabels(n):
    # We use here the fact that the even partitions of n produced by
    # typB._conjlabels are in the same order as the following list
    #
    #    [[2*i for i in part] for part in typA._conjlabels(n//2)]
    #
    # This is the list produced by _evenpartitions(n) which uses the
    # same algorithm as for producing partitions.

    # This produces all bipartitions. If L is the list of all
    # bipartitions returned then the sublist L[::2] is ordered in graded
    # reverse lexicographic ordering and the sublist L[1::2] is obtained
    # from L[::2] by swapping the entries in the bipartition.
    if not n % 2:
        for bpart in typB._conjlabels(n):
            if not len(bpart[1]) % 2 and (bpart[1]
                                          or any(x % 2 for x in bpart[0])):
                yield bpart
        for part in _evenpartitions(n):
            yield [part, '+']
            yield [part, '-']
    else:
        for bpart in typB._conjlabels(n):
            if not len(bpart[1]) % 2:
                yield bpart

# As in [GP00] we get representatives of the classes in type Dn from
# those of type Bn. Assume (alpha;beta) is a bipartition labelling a
# class then the representative is given by
#
#   b^+(m_1,alpha_1)...b^+(m_k,alpha_k)b^-(n_1,beta_1)...b^-(n_l,beta_l)
#
# as before. Just as before we define the positive block to be
#
#   b^+(m, d) = m'...(m+d-2)' = m....(m+d-2)
#
# Here the ' notation denotes the embedding Dn -> Bn defined at the
# beginning of this file. To get a word in the generators of Dn we write
# each negative block in the form (n-1)'.w, with w a reduced expression
# in the Dn generators. It's clear that the generator (n-1)' commutes
# with all negative blocks except possibly the far most right negative
# block. Hence, we obtain a reduced expression in the generators of Dn
# by simply striking off the (n-1)' terms.
#
# The desired expression for the negative blocks is given as follows
#
#   b^-(m, d) = m'...(n-2)'(n-1)'(n-2)'...(m+d-1)'
#             = (n-1)'m...(n-3)(n-1)(n-2)...(m+d-1)
#
# Note: One might be tempted to beautify the output by swapping the
# terms (n-1) and (n-2) but this can only be done assuming we're not in
# the last block. To avoid always having to treat the last block
# differently we don't do this.
def conjclasses(inds, **kwargs):
    n = len(inds)

    for mu in _conjlabels(n):
        alpha, beta = mu
        w = []
        m = 0

        if beta == '+':
            # First build the + representative.
            for d in alpha:
                # Very good in the symmetric group (see typ1A).
                w += inds[m:m+d-1:2] + inds[m+1:m+d-1:2]
                m += d

            c = typB.centraliser([mu[0], []])
            spart = utils.intlisttostring(mu[0])

            yield spart + '.+', c, w

            # Now the - representative. Note we must make a copy of the
            # word otherwise we'll change what we've just yielded.
            w = w[:]

            # The - representative is obtained by flipping the generator
            # n-2 in the word to n-1. This must occurr in the last
            # segment of the word. Either it's at the end of the odd
            # list (hence the end of the word) or the even list.
            t = d - 2
            if t % 2:
                w[-1] = inds[-1]
            else:
                w[-t//2-1] = inds[-1]

            yield spart + '.-', c, w

        elif isinstance(beta, list):
            for d in alpha:
                # Very good in the symmetric group (see typ1A).
                w += inds[m:m+d-1:2] + inds[m+1:m+d-1:2]
                m += d
            for d in beta:
                # Note the block b^-(n-1, 1) = (n-1)' so we just delete
                # this in this case.
                if d == 1:
                    if m == 0:
                        w += inds[m:-2] + inds[-1::-1]
                    elif m != n-1:
                        w += inds[m:-2] + inds[-1:m+d-2:-1]
                else:
                    w += inds[m:-2] + inds[-1:m+d-2:-1]
                m += d

            yield utils.parttupletostring(mu), typB.centraliser(mu)//2, w

def conjclasses_min(inds, **kwargs):
    n = len(inds)

    for mu in _conjlabels(n):
        if mu[1] == '+':
            c = typB.centraliser([mu[0], []])
            spart = utils.intlisttostring(mu[0])

            yield spart + '.+', c
            yield spart + '.-', c

        elif isinstance(mu[1], list):
            yield utils.parttupletostring(mu), typB.centraliser(mu)//2

def wordtoclass(n, w):
    """
    Returns the name of the conjugacy class in the Weyl group of type D_n
    containing the element w given as a word in the standard generators.

    """
    lab = _cdata.wordtoclassD(n, w)
    
    if isinstance(lab[1], tuple):
        return utils.parttupletostring(lab)
    else:
        return "{0}.{1}".format(utils.intlisttostring(lab[0]), lab[1])

#def wordtoclass(n, elm):
#    # We use the following representation of the Weyl group of type D_n
#    # into S_{2n}. If s_0, ..., s_{n-1} are the generators then we have
#    #
#    #       s_i -> (i, i+1)(N-i, N-i-1) if i = 0, ..., n-2
#    #   s_{n-1} -> (n-2, n)(n-1, n+1)
#    #
#    # Here N = 2n-1. This is compatible with our embedding of the type
#    # B_n Weyl group.
#    N = 2*n - 1
#    p = permutat.Perm(range(2*n))
#    for i in elm:
#        if i == n-1:
#            p *= (n-2, n)
#            p *= (n-1, n+1)
#        else:
#            p *= (i, i+1)
#            p *= (N-i, N-i-1)
#
#    # Get the B_n label first.
#    seen = list(range(n))
#    pos, neg = [], []
#    for i in range(n):
#        if seen[i] is not None:
#            cyc = p.cycle(i)
#            if N-i in cyc:
#                neg.append(len(cyc)//2)
#            else:
#                pos.append(len(cyc))
#        for j in cyc:
#            if j >= n:
#                seen[N-j] = None
#            else:
#                seen[j] = None
#
#    if neg or any(i % 2 for i in pos):
#        return (tuple(sorted(pos, reverse=True)),
#                tuple(sorted(neg, reverse=True)))
#    else:
#        # We now need to work out whether with have a + class or a -
#        # class. We start by finding an element x of the Weyl group of
#        # type B_n which conjugates p to an element of the + copy of
#        # S_n. Such an element must exist. If x lies in D_n then we have
#        # a + element. If x doesn't lie in D_n then s_{n-1}*x does but
#        # s_{n-1} sends a + element to a - element, so in that case it's
#        # a minus element. As the D_n subgroup is obtained by
#        # intersecting the B_n subgroup with the alternating group
#        # A_{2n} < S_{2n} we need only look at the parity of the
#        # permutation to determine whether it's in D_n or not.
#        sgn = 1
#        flag = True
#        while flag:
#            flag = False
#            for i in range(n):
#                j = i^p
#                if j >= n:
#                    p = (j, N-j)*p*(j,N-j)
#                    sgn *= -1
#                    flag = True
#                    break
#
#        if sgn == 1:
#            return (tuple(sorted(pos, reverse=True)), "+")
#        else:
#            return (tuple(sorted(pos, reverse=True)), "-")


########################################################################
########################################################################
##                                                                    ##
##                    Irreducible Characters                          ##
##                                                                    ##

def _charlabels(n):
    # The ordering here is such that the sum of the entries of the first
    # partition is greater than or equal to the sum of the entries in the
    # second partition.
    end = n//2
    for i in range(n, end, -1):
        for pair in itertools.product(utils.partitions(i),
                                      utils.partitions(n-i)):
            yield list(pair)
    if not n % 2:
        for pair in itertools.combinations(utils.partitions(end), 2):
            yield list(pair)
        for part in utils.partitions(end):
            yield [part, '+']
            yield [part, '-']

def irrchars(n, **kwargs):
    for mu in _charlabels(n):
        # We always know that sum(mu[0]) >= sum(mu[1]) from _charlabels.
        # Moreover we always know that the degenerate labels come in pairs
        # ordered by + then -.
        if mu[1] == '+':
            # label, a-value, b-value.
            a = utils.ainvsymbol(utils.symbol([mu[0], mu[0]], 0), n)
            b = 4*sum(i*val for i, val in enumerate(mu[0])) + n//2
            spart = utils.intlisttostring(mu[0])
            yield (spart + '.+', a, b)
            yield (spart + '.-', a, b)
        elif isinstance(mu[1], list):
            # label, a-value, b-value.
            yield (
                utils.parttupletostring(mu),
                utils.ainvsymbol(utils.symbol(mu, 0), n),
                2*sum(i*val for i, val in enumerate(mu[0]))
                    + 2*sum(i*val for i, val in enumerate(mu[1]))
                    + sum(mu[1])
            )

def chartable(n, **kwargs):
    ctBhalf, d = typB.chartablehalf(n)

    if n % 2:
        # In this case we simply need to select the appropriate columns
        # corresponding to the classes in the Weyl group of type D.
        inds = [i for i, mu in enumerate(typB._conjlabels(n))
                if not len(mu[1]) % 2]
        return [[row[i] for i in inds] for row in ctBhalf]

    else:
        # Need the degenerate and non-degenerate conjugacy class indices.
        conjinds = {'deg' : [], 'nondeg' : []}
        for ind, mu in enumerate(typB._conjlabels(n)):
            if not len(mu[1]) % 2:
                if not mu[1] and all(not x % 2 for x in mu[0]):
                    conjinds['deg'].append(ind)
                else:
                    conjinds['nondeg'].append(ind)

        e = n//2
        m = 2*d
        chartabdim = len(ctBhalf) + d
        out = [None]*chartabdim

        # First the non-degenerate characters. Note that on the
        # degenerate classes the value is the same as the corresponding
        # value in typeB.
        out[:-m] = [[row[i] for i in conjinds['nondeg']]
                    + list(itertools.chain.from_iterable(
                        zip(*([row[i] for i in conjinds['deg']],)*2)))
                    for row in ctBhalf[:-d]]

        # Degenerate characters on non-degenerate classes. We pad here
        # with None type so that we can slice the lists below.
        out[-m::2] = [[row[i]//2 for i in conjinds['nondeg']]
                      + [None]*m
                      for row in ctBhalf[-d:]]
        out[-m+1::2] = [[row[i]//2 for i in conjinds['nondeg']]
                        + [None]*m
                        for row in ctBhalf[-d:]]

        # Degenerate characters on degenerate classes. Let S = x_+ + x_-
        # be the sum of 2 degenerate characters, i.e., the restriction
        # of a character from B_n. We denote by D = x_+ - x_- the
        # difference character. Clearly we have
        #
        #   x_+(w) = (S(w) + D(w))//2       x_-(w) = (S(w) - D(w))//2
        #
        # The difference character is described in [Gec15]. Below D
        # describes the values D(w_+) on the positive classes. We then
        # have D(w_-) == -D(w_+).
        scal = [(-1)**e*2**(len(part)-1) for part in typA._conjlabels(e)]
        D = [[a*b for a, b in zip(row, scal)] for row in typA.chartable(e-1)]
        S = [[row[i]//2 for i in conjinds['deg']] for row in ctBhalf[-d:]]

        plusD = [[a + b for a, b in zip(*row)] for row in zip(S, D)]
        minusD = [[a - b for a, b in zip(*row)] for row in zip(S, D)]

        # Positive chars on pos then neg classes.
        for i, row in enumerate(out[-m::2]):
            row[-m::2], row[-m+1::2] = plusD[i], minusD[i]

        # Negative chars on pos then neg classes.
        for i, row in enumerate(out[-m+1::2]):
            row[-m::2], row[-m+1::2] = minusD[i], plusD[i]
        return out



