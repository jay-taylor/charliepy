########################################################
# Data for type A_n Coxeter groups and reductive groups.
# The Dynkin diagram here follows Bourbaki so we have
#
#            0   1   2        n-1
#            o---o---o-- ... --o

from .. import utils
from .. import permutat
from . import __USE_C__
if __USE_C__:
    from . import _chartabs

import numpy as np
from collections import OrderedDict

########################################################################
########################################################################
##                                                                    ##
##                          Cartan Matrix                             ##
##                                                                    ##

def cartanmat(n, **kwargs):
    """Takes a zero matrix and returns the type A_n Cartan matrix as a
    NumPy array."""
    if n < 0:
        raise ValueError("Rank must be positive!")
    elif n == 0:
        return np.array([], dtype = 'int8')
    elif n == 1:
        return np.array([2], dtype = 'int8')
    else:
        C = np.zeros((n,n), dtype='int8')
        np.fill_diagonal(C, 2)
        C[range(n-1), range(1,n)] = [-1]*(n-1)
        C[range(1,n), range(n-1)] = [-1]*(n-1)
        return C


def rootlengths(n, **kwargs):
    """
    Returns a generator giving the relative root lengths.

    """
    return (1 for _ in range(n))


def diagram(inds):
    """Prints the Dynkin diagram."""
    print("A{} :  ".format(len(inds)) + " - ".join(str(x) for x in inds))
    return None


def degrees(n):
    return list(range(2, n+2))

def longestword(inds):
    n = len(inds) + 1
    if n % 2:
        return (inds[0::2] + inds[1::2])*(n//2) + inds[0::2]
    else:
        return (inds[0::2] + inds[1::2])*(n//2)


########################################################################
########################################################################
##                                                                    ##
##                        Conjugacy Classes                           ##
##                                                                    ##

def _conjlabels(n):
    # See clp.utils.partitions for more information concerning the
    # algorithm.
    part = [1 for j in range(n)]
    i, m = -1, n
    
    while True:
        yield part[:m]

        if 2 <= m - i - 1:
            i += 1
            m -= 1
            part[i] = 2
            continue

        if m == 1:
            return
        
        # Do the replacement [x, ..., x, r] -> [(x+1), 1, ..., 1].
        # First find the index of the first copy of x; note this will
        # become the index of the largest entry not equal to 1.
        i = m-3
        while 0 <= i:
            if part[i] == part[i+1]:
                i -= 1
            else:
                break
        i += 1

        # Now set the right number of 1s. First correcting the length.
        m = i + (m-2-i)*part[m-2] + part[m-1]
        j = i + 1
        while j < m:
            part[j] = 1
            j += 1
        part[i] += 1


def centraliser(mu):
    """
    Returns the order of the centraliser of an element of cycle type of a
    given partition in the symmetric group. See [Pfe94, 1.2] for more details.

    """
    cent, last, k = 1, 0, 1
    for x in mu:
        cent *= x
        if x == last:
            k += 1
            cent *= k
        else:
            k = 1
        last = x

    return cent

def conjclasses(inds, **kwargs):
    n = len(inds)

    if not n:
        return ([], 1, '()')

    for mu in _conjlabels(n+1):
        w, i = [], -1
        c, last, k = 1, 0, 1
        # The representative is simply a Coxter element in a parabolic subgroup
        # of W labelled by the partition mu. As in GAP-CHEVIE this is a 'very
        # good' representative as in [GM, Theorem 1.1].
        for x in mu:
            # This creates the representative.
            w += inds[i+1:i+x:2] + inds[i+2:i+x:2]
            i += x
            
            # This determines the centraliser order. See [Pfe94, 1.2].
            c *= x
            if x == last:
                k += 1
                c *= k
            else:
                k = 1
            last = x

        yield utils.intlisttostring(mu), c, w

def conjclasses_min(ind, **kwargs):
    n = len(ind)

    if not n:
        return ((1, '()'))

    return ((utils.intlisttostring(mu), centraliser(mu))
            for mu in _conjlabels(n+1))


def wordtoclass(elm, inds, n):
    """
    Takes a word in the standard generators and returns the parameter of its
    conjugacy class.

    """
    p = permutat.Perm(range(n + 1))
    for i in elm:
        p *= (i, i + 1)
    return utils.intlisttostring(p.cycletype(True))



########################################################################
########################################################################
##                                                                    ##
##                    Irreducible Characters                          ##
##                                                                    ##

def _charlabels(n):
    # See ..utils.partitions for more information concerning the
    # algorithm.
    part = [1 for j in range(n)]
    part[0] = n
    i, m = 0, 1

    while m < n:
        yield part[:m]

        # Easy case [2] -> [1, 1]; this case is the most common. Note we
        # back up the list because we've made the current entry a 1.
        if part[i] == 2:
            part[i] -= 1
            m += 1
            i -= 1
            continue

        # Do the replacement [(x+1), 1, ..., 1] -> [x, ..., x, r].
        part[i] -= 1
        m += 1
        while m-i-1 > part[i]:
            part[i+1] = part[i]
            m -= part[i]-1
            i += 1

        # If there is more than 1 trailing 1 then add r.
        if m-i-1 > 1:
            part[i+1] = m-i-1
            m -= m-i-2
            i += 1

    yield part


def irrchars(n, **kwargs):
    # For the a-value and b-value see 5.4.2, 5.4.4 of [GP00]. For the
    # A-value we use 5.11.5 of [Lus84] to get that
    #
    #       A(mu) = N - a(mu*) = N - a*(mu)
    #
    # where N is the number of positive roots, mu* is the dual
    # partition, and a*(mu) is as in 5.4.1 of [GP00].
    if n:
        N = n*(n+1)//2
        for mu in _charlabels(n+1):
            a = sum(i*val for i, val in enumerate(mu))
            yield (utils.intlisttostring(mu), a, a)
    else:
        yield ('()', 0, 0)
     
    # A-value : N - sum(val*(val-1) for val in mu)//2


def chartable(r, **kwargs):
    """
    returns the character table of the irreducible Coxeter group of
    type A_r, i.e., the symmetric group S_(r+1). The rows and columns
    are indexed by the partitons of r+1, as ordered in partitions(r+1).
    The function uses the fast algorithm constructed by Goetz Pfeiffer in
    [Pfe94].

    >>> partitions(4)
    [[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]]
    >>> chartable(3)
    array([[ 1, -1,  1,  1, -1],                # sign character
           [ 3, -1, -1,  0,  1],                # reflection character
           [ 2,  0,  2, -1,  0],
           [ 3,  1, -1,  0, -1],
           [ 1,  1,  1,  1,  1]])               # trivial character

    """
    n = r+1
    scheme = inductionscheme(n)

    chartabdim = len(scheme[r])

    # Construct the character table to be filled.
    out = [None]*chartabdim

    # Define how to construct a column of the character table of S_(m+1) from
    # the column t of the character table of S_(k+1).
    if __USE_C__:
        charcol = _chartabs.charcolA
    else:
        def charcol(schm, t, k):
            col = [None]*len(schm)
            for ind, pi in enumerate(schm):
                val = 0
                for i in pi[k]:
                    if i < 0:
                        val -= t[-i-1]
                    else:
                        val += t[i-1]
                col[ind] = val

            return col

    # Construct all necessary columns.
    colsm = [[] for i in range(n)]
    for k in range(n//2):
        colsm[k].append(charcol(scheme[k], [1], k))

        for m in range(k+1, n-k-1):
            for t in colsm[m-k-1]:
                colsm[m].append(charcol(scheme[m], t, k))

    # Fill the first n-1 columns of the array and then the last one. The
    # variable c keeps track of which column we're filling.
    c = 0
    for k in range(r):
        for t in colsm[r-k-1]:
            out[c] = charcol(scheme[r], t, k)
            c += 1
    out[c] = charcol(scheme[r], [1], r)

    return list(map(list, zip(*out)))


def inductionscheme(n):
    """
    This constructs an induction scheme following section 3 of
    [Pfe94]. It returns a list of length n. Suppose alpha is a partition
    of 1 <= m <= n and i is the index of alpha in partitions(m). For all
    1 <= k <= m we have

        inductionscheme(n)[m-1][i][k-1]
    
    is a list of integers. For each integer p in the list we have
    abs(p)-1 gives the index of a partition of m-k which occurs from
    alpha by removing a k-hook. The sign cmp(p,0) of p gives the leg
    parity of that hook.  Note that we must record the positions as
    non-zero integers to keep track of the leg parity.

    """
    # allparts stores all partitions. Note that allparts[m] stores the
    # partitions of m+1. Similarly for scheme.
    allparts = [None]*n
    scheme = [None]*n

    # This auxilliary function computes all partitions obtained from beta by
    # removing a hook. It also keeps track of the leg parity corresponding
    # to that hook. Note that beta is a betaset corresponding to a partition of
    # m+1.
    def hooks(beta, m):
        # First list contains indices and second contains signs.
        hks = [[] for i in range(m+1)]

        for i in beta:
            sgn = 1

            for j in range(1, i+1):
                if i-j in beta:
                    sgn *= -1
                elif j == m+1:
                    hks[m].append(sgn)
                else:
                    # If we add a 0 to gamma then we have to obtain a reduced
                    # beta set to find it in the list. Here we assume that beta
                    # contains no zeroes to start with.
                    if i == j:
                        gamma = beta - {i} | {i-j}

                        # Find the largest subset {0,1,...,i} contained in
                        # gamma, remove it and then subtract i+1 from all
                        # remaining entries.
                        k = 0
                        while k in gamma:
                            k += 1
                        gamma = frozenset([x - k for x in gamma -
                                                set(range(k))])
                    else:
                        gamma = frozenset(beta - {i} | {i-j})

                    # We have to store the position of gamma in the list of
                    # partitions as a number >0 otherwise we can't keep track
                    # of the parity.
                    hks[j-1].append( sgn*(allparts[m-j][gamma]+1) )
        return hks

    for i in range(n):
        allparts[i] = OrderedDict([utils.betaset(part, frozen=True), i]
                                   for i, part in enumerate(_charlabels(i+1)))
        scheme[i] = [hooks(beta, i) for beta in allparts[i].keys()]
    
    return scheme









