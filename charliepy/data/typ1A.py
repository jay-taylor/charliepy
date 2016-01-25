########################################################
# Data for type A_n Coxeter groups and reductive groups.
# The Dynkin diagram here follows Bourbaki so we have
#
#            0   1   2        n-1
#            o---o---o-- ... --o

from .. import utils as utils
from . import __USE_C__
if __USE_C__:
    from . import _chartabs

import numpy as np
from collections import OrderedDict

def cartanmat(n):
    """Takes a zero matrix and returns the type A_n Cartan matrix as a
    NumPy array."""
    if n < 0:
        raise ValueError("Rank of A_n must be positive!")
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

def diagram(inds):
    """Prints the Dynkin diagram."""
    print("A{} ".format(len(inds)) + " - ".join(str(x) for x in inds))
    return None

def degrees(n):
    return list(range(2, n+2))

def _conjlabels(n):
    # See clp.utils.combinat.partitions for more information concerning the
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

def _charlabels(n):
    # See clp.utils.combinat.partitions for more information concerning the
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

def conjclassdata(ind, **kwargs):
    # stores the data: representatives, centraliser orders, names
    repcentnam = [[], [], []]
    n = len(ind)

    for mu in _conjlabels(n+1):
        w = []
        i = -1
        # The representative is simply a Coxter element in a parabolic subgroup
        # of W labelled by the partition mu. As in GAP-CHEVIE this is a 'very
        # good' representative as in [GM, Theorem 1.1].
        for x in mu:
            w += ind[i+1:i+x:2] + ind[i+2:i+x:2]
            i += x
        repcentnam[0].append(w)

        repcentnam[1].append(utils.centraliserpartition(n+1, mu))

        if n > 0:
            repcentnam[2].append(utils.intlisttostring(mu))
        else:
            repcentnam[2].append(' ')

    return repcentnam

def irrchardata(n, **kwargs):
    if n == 0: 
        nam = ['1']
        ainv = [0]
        binv = [0]

    else:
        # According to timeit this is faster than putting everything together
        # in a single for loop and using list.append.
        nam = list(map(utils.intlisttostring, _charlabels(n+1)))
        ainv = [sum(i*val for i, val in enumerate(mu)) for mu in
                _charlabels(n+1)]

        return [nam, ainv, ainv[:]]

def chartable(r, **kwargs):
    """returns the character table of the irreducible Coxeter group of
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
    out = np.empty((chartabdim, chartabdim), dtype='int')

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
            out[:,c] = charcol(scheme[r], t, k)
            c += 1
    out[:,c] = charcol(scheme[r], [1], r)

    return out


def inductionscheme(n):
    """This constructs an induction scheme following section 3 of
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
