########################################################
# Data for type B_n Coxeter groups and reductive groups.
# The Dynkin diagram here follows Bourbaki so we have
#
#            0   1                n-1
#            o---o-- ... --o---o=>=o

from .. import utils
from .. import permutat
from . import typ1A
from . import _cdata

import numpy as np
import itertools
#import operator
import collections
#import multiprocessing as mp

########################################################################
########################################################################
##                                                                    ##
##                          Cartan Matrix                             ##
##                                                                    ##

def cartanmat(n):
    if n < 1:
        raise ValueError("Rank of B_n must be at least 1!")
    elif n == 1:
        return np.array([2], dtype = 'int8')

    C = typ1A.cartanmat(n)
    C[n-2][n-1] = -2
    return C

def rootlengths(n, **kwargs):
    """
    Returns a generator giving the relative root lengths.

    """
    if n == 1:
        yield 1
    else:
        for i in range(n-1):
            yield 2
        yield 1

def diagram(inds):
    """Prints the Dynkin diagram."""
    out = "B{} :  ".format(len(inds))
    out += " - ".join(str(x) for x in inds[:-1])
    out += " > {}".format(inds[-1])
    print(out)
    return None

def degrees(n):
    return list(range(2, 2*n+1, 2))

def longestword(inds):
    return (inds[0::2] + inds[1::2])*len(inds)

def maxparachain(inds):
    return inds

########################################################################
########################################################################
##                                                                    ##
##                        Conjugacy Classes                           ##
##                                                                    ##

# We wish to split each partition of n into a bipartition. This is done
# as follows. Let P = [k_1*1, k_2*2, ..., k_n*n] be the part count form
# of a partition of n. Now assume [L, R] is a bipartition of n such that
# the multiset of entries of L cup R is the same as that of P. We may
# uniquely describe the bipartition [L, R] using the partcount form as
#
#       L = [l_1*1, l_2*2, ..., l_n*n]
#       R = [r_1*1, r_2*2, ..., r_n*n]
#
# with l_i + r_i = k_i.
#
# All possible lists [l_1, ..., l_n] and [r_1, ..., r_n] are then
# obtained by taking all the compositions of each number occuring in the
# list [k_1, ..., k_n] into 2 parts.
#
# Pfeiffer's algorithm does this like an odometer. He starts with the
# lists
#
#       L = [k_1, ..., k_n]
#       R = [0, ..., 0]
#
# and then increases the rightmost term on each iteration. This is
# exactly the behaviour of itertools.product, which we use here.  Note
# however that in the end we must return our lists in reverse because we
# want our partitions to be decreasing.
def _conjlabels(n):
    for part in typ1A._conjlabels(n):
        pcount = collections.Counter(part)
        # We aren't gauranteed that the elements of Counter are sorted.
        entries = sorted(pcount)

        def _bipartint(k):
            for i in range(k + 1):
                yield (k-i, i)

        gen = (_bipartint(pcount[x]) for x in entries)
        for item in itertools.product(*gen):
            # Note that the generator produces items of the form
            # ((l_1,r_1), ..., (l_n, r_n)). So we have to reverse the
            # lists at the end.
            lpart, rpart = [], []
            for i, inds in enumerate(item):
                lpart += (entries[i],)*inds[0]
                rpart += (entries[i],)*inds[1]
            yield [lpart[::-1], rpart[::-1]]

def centraliser(mu):
    """
    returns the order of the centraliser of an element of a given type
    (specified by an r-tuple of partitions mu) in the wreath product of a cyclic
    group of order r with the full symmetric group of degree n. (The program is
    taken from the GAP library and re-written almost 1-1 in python.)

    """
    cent = 1
    for i in range(2):
        last, k = 0, 1

        for x in mu[i]:
            cent *= 2*x
            if x == last:
                k += 1
                cent *= k
            else:
                k = 1
            last = x

    return cent

# For the representatives of the conjugacy classes we largely follow
# section 3.4 of [GP00] but we have to make some slight modifications
# because we have a different labelling of the Dynkin diagram. For any
# integers m >= 0, d >= 1 we define
#
#   b^+(m, d) = m....(m+d-2)
#
# where the indices denote simple reflections. This is called a positive
# block. We then define a negative block by setting
#
#   b^-(m, d) = m...(m+d-2)(m+d-1)...(n-2)(n-1)(n-2)...(m+d-1)
#
# Now assume (alpha;beta) is a bipartition labelling a class then we
# have a representative for this class is given by the product of blocks
#
#   b^+(m_1,alpha_1)...b^+(m_k,alpha_k)b^-(n_1,beta_1)...b^-(n_l,beta_l)
#
# Here the m_i and the n_j simply sum all the proceeding parts of the
# bipartitions that have come before.
def conjclasses(inds, **kwargs):
    n = len(inds)

    for mu in _conjlabels(n):
        alpha, beta = mu
        w = []    
        m = 0

        for d in alpha:
            # We make the representative very good by having it be a
            # very good representative in the symmetric group (see
            # typ1A).
            w += inds[m:m+d-1:2] + inds[m+1:m+d-1:2]
            m += d
        for d in beta:
            # Have to put a special case here because of the way Python
            # interprets -1 in its slice notation.
            if d == 1 and m == 0:
                w += inds[m:] + inds[-2::-1]
                m += d
            else:
                w += inds[m:] + inds[-2:m+d-2:-1]
                m += d

        yield utils.parttupletostring(mu), centraliser(mu), w

def conjclasses_min(inds, **kwargs):
    n = len(inds)

    return ((utils.parttupletostring(mu), centraliser(mu))
            for mu in _conjlabels(n))

def wordtoclass(n, w):
    """
    Returns the name of the conjugacy class in the Weyl group of type B_n
    containing the element w given as a word in the standard generators.

    """
    return utils.parttupletostring(_cdata.wordtoclassB(n, w))



########################################################################
########################################################################
##                                                                    ##
##                    Irreducible Characters                          ##
##                                                                    ##

def _charlabels(n):
    # This produces all bipartitions. If L is the list of all
    # bipartitions returned then the sublist L[::2] is ordered in graded
    # reverse lexicographic ordering and the sublist L[1::2] is obtained
    # from L[::2] by swapping the entries in the bipartition. The
    # optional flag partinds gives the bipartitions by listing their
    # indices in the sequence of partitions.
    end = n//2
    for i in range(n, end, -1):
        for pair in itertools.product(utils.partitions(i),
                                      utils.partitions(n-i)):
            yield [pair[0], pair[1]]
            yield [pair[1], pair[0]]
    if not n % 2:
        for pair in itertools.combinations(utils.partitions(end), 2):
            yield [pair[0], pair[1]]
            yield [pair[1], pair[0]]
        for part in utils.partitions(end):
            yield [part, part]

def _charlabelshalf(n):
    # This produces precisely half the bipartitions of n. If [L, R] is a
    # bipartition of n then we represent this bipartition as a list of
    # lists [[a, x], [b, y]] where a is the sum of the entries of L and
    # x is the index of L in list(typ1A._charlabels(a)) and similarly
    # for R. The bipartitions are produced in graded reverse
    # lexicographic order. In particular, they are ordered in the
    # following way:
    #
    #   - a > b,
    #   - a == b and x < y,
    #   - a == b and x == y.
    #
    # If n is odd then exactly half the bipartitions of n satisfy these
    # conditions. If n is even then this covers slightly more than half
    # as we also get the so-called degenerate bipartitions.
    end = n//2
    allnrs = utils.nrpartitions(n, lessthan=True)

    for i in range(n, end, -1):
        for pair in itertools.product(range(allnrs[i]),
                                      range(allnrs[n-i])):
            yield [[i, pair[0]], [n-i, pair[1]]]
    if not n % 2:
        for pair in itertools.combinations(range(allnrs[end]), 2):
            yield [[end, pair[0]], [end, pair[1]]]
        for i in range(allnrs[end]):
            yield [[end, i], [end, i]]

def irrchars(n, **kwargs):
    # For the a-value and b-value see 6.4.3 and 5.5.3 of [GP00].
    for mu in _charlabels(n):
        # label, a-value, b-value.
        yield (
            utils.parttupletostring(mu),
            utils.ainvsymbol(utils.symbol(mu, d=1), n),
            2*sum(i*val for i, val in enumerate(mu[0]))
                + 2*sum(i*val for i, val in enumerate(mu[1]))
                + sum(mu[1])
        )

# Introduce a named tuple instead of a dict to make the code below
# cleaner.
_Charcol = collections.namedtuple('Charcol', ['col', 'pos', 'nlen'])

def chartablehalf(n, **kwargs):
    # Get the induction scheme and construct the character table to be
    # filled. d is the number of partitions of n//2. This is simply
    # being passed through from inductionscheme so it doesn't have to be
    # recomputed.
    scheme, d = inductionscheme(n)
    ncols = utils.nrbipartitions(n) 
    cthalf = [None]*ncols

    charcol = _cdata.charcolB
    
    ## We now specify how to construct a new column of the character
    ## table of W_(m+1) from the column t of the character table of
    ## W_(k+1). Note the parameter p keeps track of whether we're adding
    ## to the column an unsigned, p = 0, or signed, p = 1, block. The
    ## parameter q denotes whether there are an even, q = 0, or odd, q =
    ## 1, number of signed blocks in the element of W_(k+1).
    #def charcol(schm, t, k, p, q):
    #    col = [None]*len(schm)
    #    for ind, pi in enumerate(schm):
    #        val = 0
    #        for j, swp in pi[0][k]:
    #            if j < 0:
    #                if swp and q:
    #                    val += t[-j-1]
    #                else:
    #                    val -= t[-j-1]
    #            else:
    #                if swp and q:
    #                    val -= t[j-1]
    #                else:
    #                    val += t[j-1]
    #        for j in pi[1][k]:
    #            if j < 0:
    #                if p:
    #                    val += t[-j-1]
    #                else:
    #                    val -= t[-j-1]
    #            else:
    #                if p:
    #                    val -= t[j-1]
    #                else:
    #                    val += t[j-1]
    #        col[ind] = val
    #    return col
    
    # Now construct all the columns of the character table.
    cols = [[] for i in range(n-1)]

    # First the m-cycle in all possible places, keeping track of where
    # you put it. The set contains all parts which have appeared on the
    # right hand side.  Once something has appeared on the right hand
    # side it can then only be put on the right hand side, not the left.
    # The integer gives the parity of the number of signed blocks in the
    # representative.
    for k in range(n//2):
        cols[k].append(_Charcol(charcol(scheme[k], [1], k, 0, 0), set(), 0))
        cols[k].append(_Charcol(charcol(scheme[k], [1], k, 1, 1), set([k]), 1))

        # Add the lists constructed above to everything you know.
        for m in range(k+1, n-k-1):
            for t in cols[m-k-1]:
                # If you haven't ever added a cycle of the length k to
                # the right hand side then it's OK to add it to the
                # left. Otherwise only consider the version with the
                # k-cycle added on the right.
                if not k in t.pos:
                    cols[m].append(
                        _Charcol(charcol(scheme[m], t.col, k, 0, t.nlen),
                                 t.pos, t.nlen)
                    )
                    cols[m].append(
                        _Charcol(charcol(scheme[m], t.col, k, 1, t.nlen),
                                 t.pos | {k}, (t.nlen + 1) % 2)
                    )
                else:
                    cols[m].append(
                        _Charcol(charcol(scheme[m], t.col, k, 1, t.nlen),
                                 t.pos, (t.nlen + 1) % 2)
                    )
    
    # Fill the first n-2 columns of the array and then the last two.
    c = 0
    for k in range(n-1):
        for t in cols[n-k-2]:
            if not k in t.pos:
                cthalf[c] = charcol(scheme[n-1], t.col, k, 0, t.nlen)
                cthalf[c+1] = charcol(scheme[n-1], t.col, k, 1, t.nlen)
                c += 2
            else:
                cthalf[c] = charcol(scheme[n-1], t.col, k, 1, t.nlen)
                c += 1
    cthalf[ncols-2] = charcol(scheme[n-1], [1], n-1, 0, 0)
    cthalf[ncols-1] = charcol(scheme[n-1], [1], n-1, 1, 1)

    # d is the number of partitions of n//2
    return list(map(list, zip(*cthalf))), d


def chartable(n):
    cthalf, d = chartablehalf(n)

    out = [None]*len(cthalf[0])
    signchar = [(-1)**len(bpair[1]) for bpair in _conjlabels(n)]

    if not n % 2:
        out[:-d:2] = cthalf[:-d]
        out[1:-d:2] = [[a*b for a, b in zip(row, signchar)]
                       for row in cthalf[:-d]]
        out[-d:] = cthalf[-d:]
    else:
        out[::2] = cthalf
        out[1::2] = [[a*b for a, b in zip(row, signchar)]
                     for row in cthalf]
     
    return out

# Be careful if changing the inductionscheme function below. The C
# version of charcol relies on the specific format of the ouput from
# inductionscheme.
def inductionscheme(n):
    # Get the type A induction scheme for the symmetric group S_n.
    schemesym = typ1A.inductionscheme(n)
    nrparts = utils.nrpartitions(n, lessthan=True)

    # We store these values here as they will be used frequently by
    # bipartind. See the comments below for their meaning.
    cummul = [[None]*(m+1) for m in range(n+1)]
    for m in range(n+1):
        for a in range(m, m//2-1, -1):
            cummul[m][a] = sum([nrparts[j]*nrparts[m-j]
                                    for j in range(m, a, -1)])

    def bipartind(bpart, m):
        """
        Assume m is an integer and bpart = [[a, x], [b, y]] is a bipartition
        of m where [a, x] records a partition by letting a be the sum of the
        entries and x be the index of the partition in the list
        typ1A._charlabels(m).

        """
        # Assume bpart = [[a, x], [b, y]] is a bipartition of m where
        # each partition is represented by a pair [s, i] where s is the
        # sum of the entries in the partition and i is the index of the
        # partition in list(typ1A._charlabels(s)). Note that any such
        # pair generated by _charlabelshalf(m) satisfies a >= b.
        #
        # We will need to count the number of bipartitions that have
        # come before bpart in the list _charlabelshalf(m). This
        # certainly includes everything of the form:
        #
        #   - [[a', x'], [m-a', y']] with a+1 <= a' <= m,
        #                                 0 <= x' < p(a'),
        #                                 0 <= y' < p(m-a').
        #
        # As the number of such elements is called frequently it is
        # stored in the list cummul and is accessed via cummul[m][a].
        #
        # Now assume for the moment that a != b. Then we also have
        # everything of the form:
        #
        #   - [[a, x'], [b, y']] with 0 <= x' < x and 0 <= y' < p(b).
        #   - [[a, x], [b, y']] with 0 <= y' < y.
        #
        # From this one easily computes the index of bpart from x and y.
        #
        # Now assume that a == b. If x == y then we also have everything
        # of the form:
        #
        #   - [[a, x'], [a, y']] with 0 <= x' < p(m//2) and x <= y' < p(m//2).
        #   - [[a, y'], [a, y']] with 0 <= x' < y.
        #
        # The contrbution of the first part is just p(m//2)Cr2, where Cr
        # denotes the number of combinations. If x < y then we have
        # everything of the form:
        #
        #   - [[a, x'], [a, y']] with 0 <= x' < x and x' < y' < p(m//2)
        #   - [[a, x], [a, y']] with x < y' < y
        #
        # An expression for the number of elements in the first row is
        # given by
        #
        #       sum_{i=1}^x (k-i) = T_{k-1} - T_{k-x-1}
        #                         = x*(2*k - x - 1)//2
        #
        # where T_i = sum_{j=0}^i j denotes the ith triangular number.
        # The contributions from the last row is clearly y-x-1. This
        # easily treats the case x <= y.
        a, x = bpart[0]
        b, y = bpart[1]

        if a == b:
            if x == y:
                return cummul[m][a] + nrparts[a]*(nrparts[a]-1)//2 + y
            else:
                return cummul[m][a] + x*(2*nrparts[a] - x - 3)//2 + y - 1
        else:
            return cummul[m][a] + x*nrparts[m-a] + y

    # If bpair is an entry in allbiparts[m] then this auxilliary
    # function computes all partitions obtained from bpair[0] and
    # bpair[1] by removing a hook from the partition. It also keeps
    # track of the leg parity corresponding to that hook.
    def hooks(bpair, m):
        hks = [[[] for i in range(m)], [[] for i in range(m)]]

        # We now check to see if any (j+1)-hooks can be removed from
        # bpair[0].  Note that bpair consists of a partition of a and a
        # partition of m-a.  When we remove the (j+1)-hook we obtain a
        # bipartition of m-j-1 consisting of a partition of a-j-1 and a
        # partition of m-a. We need to keep track of whether or not this
        # bipartition is contained in the list
        # list(_charlabelshalf(m-j-1)). This will be the case unless:
        #       - a-j-1 < m-a  iff  2a-m-1 < j.
        #       - a-j-1 == m-a and bpair[0][1] > bpair[1][1].
        #
        # In this case we not only record the index of the resulting
        # hook but also an additional integer which is 0 or 1. A value
        # of 1 indicates that we had to swap the entries in the
        # bipartition to be in the list produced by _charlabelshalf and
        # a value of 0 indicates we didn't.
        a, ind = bpair[0]
        for j in range(a):
            # If a (j+1)-hook exists then we construct the bipartition
            # obtained from bpair by removing this (j+1)-hook. We have
            # to first make a duplicate of bpair because we can't change
            # what we're looping over.
            tmp = [None]*len(schemesym[a-1][ind][j])
            diff = 2*a - m - 1

            for c, k in enumerate(schemesym[a-1][ind][j]):
                if k < 0:
                    dupe = [a-j-1, -k-1]
                    sgn = -1
                else:
                    dupe = [a-j-1, k-1]
                    sgn = 1

                if j == m:
                    tmp[c] = [sgn, 1]
                else:
                    if diff < j or (diff == j and dupe[1] > bpair[1][1]):
                        tmp[c] = [sgn*(bipartind([bpair[1], dupe], m-j-1) + 1),
                                  1]
                    else:
                        tmp[c] = [sgn*(bipartind([dupe, bpair[1]], m-j-1) + 1),
                                  0]

            hks[0][j] = tmp

        # We now check to see if (j+1)-hooks can be removed from
        # bpair[1]. This time we don't have to be careful. Because the
        # bipartition we started with is assumed to consist of a pair of
        # partitions whose entries sum to a and b respectively with a >=
        # b. If we remove a (j+1)-hook from the second partition then we
        # clearly have a > b-j-1 so the resulting bipartition is
        # produced by _charlabelshalf.
        a, ind = bpair[1]
        for j in range(a):
            # If a (j+1)-hook exists then we construct the bipartition
            # obtained from bpair by removing this (j+1)-hook. We have
            # to first make a duplicate of bpair because we can't change
            # what we're looping over.
            tmp = [None]*len(schemesym[a-1][ind][j])
            for c, k in enumerate(schemesym[a-1][ind][j]):
                if k < 0:
                    dupe = [a-j-1, -k-1]
                    sgn = -1
                else:
                    dupe = [a-j-1, k-1]
                    sgn = 1

                if j == m:
                    tmp[c] = sgn
                else:
                    tmp[c] = sgn*(bipartind([bpair[0], dupe], m-j-1) + 1)

            hks[1][j] = tmp
            
        return hks

    return ([[hooks(bpair, m) for bpair in _charlabelshalf(m)] 
             for m in range(1, n+1)], len(schemesym[n//2-1]))

