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
# Note that the labelling is consistent so that we have
# the diagram of D_3 is given by
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
from . import typ1B as typB
from . import typ1A as typA

import numpy as np
import itertools

def cartanmat(n):
    if n < 2:
        raise ValueError("Rank of D_n must be at least 2!")

    if n == 2:
        return np.array([[2, 0],
                         [0, 2]], dtype = 'int8')
    
    C = typA.cartanmat(n)
    C[-3:,-3:] = [[2, -1, -1], [-1, 2, 0], [-1, 0, 2]]
    return C

def diagram(inds):
    """Prints the Dynkin diagram."""
    n = len(inds)
    if n == 2:
        print("D2 {} - {}".format(inds))
        return None
    elif n == 3:
        print("D3 {} - {} -- {}".format(inds[1], inds[0], inds[2]))
        return None

    out3 = "D{} ".format(len(inds))
    out3 += " - ".join(str(x) for x in inds[:-2])
    out2 = " "*len(out3) + " /"
    out4 = " "*len(out3) + " \\"
    out1 = " "*len(out2) + " {}".format(inds[-2])
    out5 = " "*len(out2) + " {}".format(inds[-1])
    print(out1, out2, out3, out4, out5, sep="\n")
    return None

def degrees(n):
    return list(range(2, 2*n+1, 2)) + [n]

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

def _evenpartitions(n):
    # For an even number n this produces all the even partitions of n, i.e., all
    # partitions where every part is divisible by 2. The ordering of the
    # partitions is the same as that produced by typA._conjlabels(n); the
    # algorithm is essentially the same. However now we start with a list
    # containing only 2s instead of only 1s. We also have:
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
    # This is the list produced by _evenpartitions(n) which uses the same
    # algorithm as for producing partitions.

    # This produces all bipartitions. If L is the list of all bipartitions
    # returned then the sublist L[::2] is ordered in graded reverse
    # lexicographic ordering and the sublist L[1::2] is obtained from L[::2] by
    # swapping the entries in the bipartition. The optional flag partinds gives
    # the bipartitions by listing their indices in the sequence of partitions.
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

def conjclassdata(ind, **kwargs):
    # stores the data: representatives, centraliser orders, names
    repcentnam = [[], [], []]
    n = len(ind)

    for mu in _conjlabels(n):
        # We construct the reps w(alpha,beta) as in [GP00, 3.4.12]. This
        # means mu[1] <-> neg blocks and mu[0] <-> pos blocks so need
        # mu[1] increasing and mu[0] decreasing. Note partitions from
        # _conjlabels(n) are decreasing. Also by [GP00, 3.4.12] a
        # bipartition labels a class only if len(mu[1]) is even.
        w = []
        m = 0

        if mu[1] == '+':
            # First build the + representative.
            for d in mu[0]:
                # We make the representative very good by having it be a very
                # good representative in the symmetric group (see typ1A). 
                if not m and d != 1:
                    # b^+(0,d) = 1...(d-1) = 0'2'...(d-1)'
                    # Need d != 1 test otherwise ind[0] is added not [].
                    w += ind[0:1] + ind[3:d:2] + ind[2:d:2]
                else:
                    # b^+(m,d) = (m+1)...(m+d-1) = (m+1)'...(m+d-1)'
                    # Note this does nothing when m = 0 and d = 1.
                    w += ind[m+1:m+d:2] + ind[m+2:m+d:2]
                m += d

            repcentnam[0].append(w)

            # Now the - representative. Note we must make a copy of the
            # word otherwise we'll change what we've just appended.
            w = w[:]
            w[0] = ind[1]
            repcentnam[0].append(w)

            c = utils.centralisertuple(n, 2, [mu[0], []])
            spart = utils.intlisttostring(mu[0])
            
            repcentnam[1] += [c, c]
            repcentnam[2] += [spart + '.+', spart + '.-']

        elif isinstance(mu[1], list):
            for d in reversed(mu[1]):
                # As there is an even no. of neg. blocks, and 0's
                # commute with blocks, we must simply write each block
                # as b.0 and then strike off the 0's.
                if not m and d != 1:
                    # b^-(0,d) = 0...(d-1) = 1'...(d-1)'0
                    # However, for cosmetics, we conjugate w, hence
                    # b^-(0,d), by 0 so that the word starts with 0'
                    # instead of 1'.
                    w += ind[0:1] + ind[2:d]
                elif m:
                    # b^-(m,d) = m...101...(m+d-1)
                    #          = m'...2'0'1'...(m+d-1)'0
                    w += ind[m:1:-1] + ind[:m+d]
                m += d

            for d in mu[0]:
                if not m and d != 1:
                    # b^+(0,d) = 1...(d-1) = 0'2'...(d-1)'
                    # Need d != 1 otherwise ind[0] is added not [].
                    w += ind[0:1] + ind[3:d:2] + ind[2:d:2]
                else:
                    # b^+(m,d) = (m+1)...(m+d-1) = (m+1)'...(m+d-1)'
                    # Note this does nothing when m = 0 and d = 1.
                    w += ind[m+1:m+d:2] + ind[m+2:m+d:2]
                m += d

            repcentnam[0].append(w)
            repcentnam[1].append(utils.centralisertuple(n, 2, mu)//2)
            repcentnam[2].append(utils.parttupletostring(mu))

    return repcentnam

def irrchardata(n, **kwargs):
    ainv, binv, nam = [], [], []

    for mu in _charlabels(n):
        # We always know that sum(mu[0]) >= sum(mu[1]) from _charlabels.
        # Moreover we always know that the degenerate labels come in pairs
        # ordered by + then -.
        if mu[1] == '+':
            a = utils.ainvsymbol(utils.symbol([mu[0], mu[0]], 0), n)
            b = 4*sum(i*val for i, val in enumerate(mu[0])) + n//2
            spart = utils.intlisttostring(mu[0])
            binv += [b, b]
            ainv += [a, a]
            nam += [spart + '.+', spart + '.-']
        elif isinstance(mu[1], list):
            b = (2*sum(i*val for i, val in enumerate(mu[0]))
                    + 2*sum(i*val for i, val in enumerate(mu[1]))
                    + sum(mu[1]))
            binv.append(b)
            ainv.append(utils.ainvsymbol(utils.symbol(mu, 0), n))
            nam.append(utils.parttupletostring(mu))

    return [nam, ainv, binv]

def chartable(n, **kwargs):
    ctBhalf, d = typB.chartablehalf(n)

    if n % 2:
        return ctBhalf[:, [i for i, mu in enumerate(typB._conjlabels(n))
                              if not len(mu[1]) % 2]]
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
        chartabdim = ctBhalf.shape[0] + d
        out = np.empty((chartabdim, chartabdim), dtype='int')

        # First the totally nondegenerate part.
        out[:-m, :-m] = ctBhalf[:-d, conjinds['nondeg']]

        # Non-degenerate characters on degenerate classes.
        out[:-m, -m::2] = ctBhalf[:-d, conjinds['deg']]
        out[:-m, -m+1::2] = ctBhalf[:-d, conjinds['deg']]

        # Degenerate characters on non-degenerate classes.
        rows = ctBhalf[-d:, conjinds['nondeg']]//2
        out[-m::2, :-m], out[-m+1::2, :-m] = rows, rows

        # Degenerate characters on degenerate classes. Let S = x_+ + x_- be the
        # sum of 2 degenerate characters, i.e., the restriction of a character
        # from B_n. We denote by D = x_+ - x_- the difference character. Clearly
        # we have
        #
        #   x_+(w) = (S(w) + D(w))//2       x_-(w) = (S(w) - D(w))//2
        #
        # The difference character is described in [Gec15]. Below D describes
        # the values D(w_+) on the positive classes. We then have
        # D(w_-) == -D(w_+).
        scal = np.array([(-1)**e*2**(len(part)-1) for part in
                                                        typA._conjlabels(e)])
        D = typA.chartable(e-1)*scal
        S = ctBhalf[-d:, conjinds['deg']]//2

        plusD = S + D
        minusD = S - D

        # Positive chars on pos then neg classes.
        out[-m::2, -m::2], out[-m::2, -m+1::2] = plusD, minusD
        
        # Negative chars on pos then neg classes.
        out[-m+1::2, -m::2], out[-m+1::2, -m+1::2] = minusD, plusD
        return out



