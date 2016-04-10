import itertools
import operator
import sys
import collections
import numpy as np

# define what is imported under import *
__all__ = [
    'betaset',
    'betasets',
    'bipartitions',
    'dualpartition',
    'flatten',
    'matrixblocks',
    'nrbipartitions',
    'nrpartitions',
    'partitions',
    'partitionsfixedlength',
    'symbol',
    'symbolrank',
    'symbols'
]

##################################################
#                Partitions

def partitions(n, k=0, reverse=False):
    """
    This function returns a generator which returns all the integer
    partitions of an integer n in the reverse lexicographic ordering.

    If an integer k is specified then the generator will only give those
    partitions whose maximal entry is less than or equal to k.

    If reverse=True is set then the partitions will be generated in the
    reverse order.

    """
    # This an implementation of Algorithm P from Knuth's The Art of Computer
    # Programming. The algorithm produces all partitions in reverse
    # lexicographic ordering. The steps of the algorithm are as follows.
    #
    # Input: [n]
    # Loop : Let the tail of the current partition be [(x+1), 1, ..., 1] with
    #        possibly no 1s. Then the next lexicographically smallest partition
    #        is obtained by replacing the tail with [x, ..., x, r] for some
    #        possible remainder r.
    # Terminate: [1, ..., 1]
    #
    # So we don't have to keep appending 1s to a list we simply start with a
    # list of the form [n, 1, ..., 1], which has length n. We then modify this
    # list in place and slice the result to obtain the partition.
    #
    # The variables do the following:
    #   m - tracks the length of current partition
    #   i - index of largest entry not equal to 1
    

    # Trivial cases first.
    if n < 0:
        raise ValueError('The integer n must be positive.')
    elif not n:
        yield []
        return

    part = [1 for j in range(n)]

    if reverse:
        i, m = -1, n

        if k == 1:
            yield part
            return
        elif k < 0:
            raise ValueError('We must have k is a positive integer.')
        elif not k or k >= n:
            k = n
        
        while m:
            yield part[:m]

            # Easy case [1, 1] -> [2].
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

            if not i and part[i] > k:
                return

    else:
        i, m = -1, 0
       
        # For 0 < k <= n we simply start the process at the lexicographically
        # largest partition beginning with k. We have to deal with some weird
        # edge cases separately.
        if not k or k >= n:
            part[0] = n
            i, m = 0, 1
        elif k == 1:
            yield part
            return
        elif k < 0:
            raise ValueError('We must have k is a positive integer.')
        else:
            tot = n
            while 1 < k:
                while k <= tot:
                    i += 1
                    m += 1
                    part[i] = k
                    tot -= k
                k -= 1
            if tot:
                m += 1

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

            # Replace as many trailing 1s with copies of x as possible.
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


def dualpartition(mu):
    """
    Returns the dual (or conjugate) partition to mu.

    """
    if mu:
        return [len([l for l in mu if l > j]) for j in range(mu[0])]
    else:
        return list()


def nrpartitions(n, lessthan=False):
    # Uses the following recurrence relation for p(n) = the number of partitions
    # of n:
    #
    #   p(n) = sum_{k>0} (-1)^{k-1}*( p(n-f(k)) + p(n-(f(k)+k)) )
    #
    # where f(k) = (3*k^2-k)/2. Setting lessthan=True will return a list p such
    # that p[i] is the number of partitions of i.
    p = [0]*(n+1)
    p[0] = 1
    for m in range(1, n+1):
        k, f, val = -1, 0, 0
        # Note that f(k+1) = f(k) + 3k + 1.
        while True:
            k += 1
            f += 3*k + 1
            if not f+k < m:
                break
            val += (-1)**k * (p[m-f] + p[m-f-k-1])
        if f <= m:
            val += (-1)**k * p[m-f]
        p[m] = val
    if lessthan:
        return p
    else:
        return p[n]


def partitionsfixedlength(n, k):
    # This is Algorithm H in Knuth's The Art of Programming. It generates all
    # partitions of n of length k.
    if k == 1:
        yield [n]
        return
    part = [1 for i in range(k+1)]
    part[0] = n-k+1
    part[k] = -1

    while True:
        yield part[:k]

        if part[1] < part[0] - 1:
            part[0] -= 1
            part[1] += 1
            continue
      
        j = 2
        s = part[0] + part[1] - 1
        while part[j] >= part[0] - 1:
            s += part[j]
            j += 1

        if j+1 > k:
            break
        part[j] += 1
        j -= 1
        while j > 0:
            part[j] = part[j+1]
            s -= part[j]
            j -= 1
        part[0] = s


def bipartitions(n, k=0, reverse=False):
    # This produces all bipartitions. If L is the list of all bipartitions
    # returned then the sublist L[::2] is ordered in graded reverse
    # lexicographic ordering. The optional flag partinds gives the bipartitions
    # by listing their indices in the sequence of partitions.
    if reverse:
        for i in range(n+1):
            for pair in itertools.product(partitions(i, k, reverse=True),
                                          partitions(n-i, k, reverse=True)):
                yield list(pair)
    else:
        for i in range(n, -1, -1):
            for pair in itertools.product(partitions(i, k), partitions(n-i, k)):
                yield list(pair)

# This function was written together with Mike Davis.
def partitiontuples(n, r):
    """
    Returns a generator which gives all posible lists L of length r such
    that sum(sum(L[i]) for i in range(r)) is n. In particular when r=2 this
    gives all bipartitions of n. For example,

    """
    # Auxilliary function. It takes a list and returns all partitions of that
    # list into r parts. This is done by slicing the list at r-1 places. For
    # instance, if r=2 then this returns vlist[:i] + vlist[i:] for all i in
    # range(len(vlist)+1). Note that we have to go one more than the length of
    # the list to get two copies of the whole list.
    def expandLikeGroup(vlist, r):
        expandedGroup = list()
        vlen = len(vlist)
        for item in itertools.combinations_with_replacement(
                range(vlen+1), r-1):
            temp = list()
            # We start slicing from the end of the list so that the first result
            # is [vlist, [], ..., []], which is simply for cosmetics.
            i = vlen
            for j in item:
                temp.append(vlist[j:i])
                i = j
            temp.append(vlist[:j])
            expandedGroup.append(temp)

        return expandedGroup

    for part in partitionsgen(n):
        grplist = list()
        plen = len(part)

        # Keep track of the number of distinct integers in the partition.
        c = 0
        for val, vallist in itertools.groupby(part):
            # Group like terms in the partition together and then return all
            # possible partitions of those groups into r parts. For example,
            # if part=[2, 1] and r=2 then
            # --> [[2], [1]]
            # --> [[[[], [2]], [[2], []]], [[[], [1]], [[1], []]]]
            grplist.append(expandLikeGroup(list(vallist), r))
            c += 1

        if c > 1:
            # If there is more than one distinct integer in the partition then
            # we have to gather things together. For example,
            #
            # [[[[], [2]], [[2], []]], [[[], [1]], [[1], []]]]
            # --> [[[], [2, 1]], [[1], [2]], [[2], [1]], [[2, 1], []]]
            result = grplist[0]
            for i in range(1, len(grplist)):
                if i == len(grplist)-1:
                    for x in itertools.product(result, grplist[i]):
                        yield map(operator.add, *x)
                else:
                    result = [map(operator.add, *x)
                              for x in itertools.product(result, grplist[i])]
   
        else:
            # If there is only one distinct integer in the partition then we're
            # done with this partition.
            for x in grplist[0]:
                yield x

def nrbipartitions(n, lessthan=False):
    allnrs = nrpartitions(n, lessthan=True)

    if lessthan:
        p = [None]*(n+1)
        p[0] = 1
        for m in range(1, n+1):
            val = 0
            for i in range(m, m//2, -1):
                val += 2*allnrs[i]*allnrs[m-i]
            if not m % 2:
                val += allnrs[m//2]**2
            p[m] = val
        return p
    else:
        val = sum(2*allnrs[i]*allnrs[n-i] for i in range(n, n//2, -1))
        if n % 2:
            return val
        else:
            return val + allnrs[n//2]**2



##################################################
#                   Symbols                      #

def betaset(part, frozen=False):
    """
    Returns a betaset corresponding to the partition part. If part contains
    no 0 entries then the betaset will be reduced. Note this assumes that the
    partition is weakly decreasing.
    
    """
    if frozen:
        return frozenset([val + len(part) - i - 1 for i, val in
                                                        enumerate(part)])
    else:
        return {val + len(part) - i - 1 for i, val in enumerate(part)}


def betasets(n, frozen=False):
    if frozen:
        for part in partitions(n):
            yield frozenset([val + len(part) - i - 1 for i, val in
                                                        enumerate(part)])
    else:
        for part in partitions(n):
            yield {val + len(part) - i - 1 for i, val in enumerate(part)}


def symbol(bipart, d=0):
    """
    Returns a pair of beta sets corresponding to the bipartition bipart with
    defect d. If part contains no 0 entries then the betaset will be reduced.

    """
    x = len(bipart[0]) - len(bipart[1])
    if x <= d:
        return [betaset(bipart[0] + [0]*(d-x)), betaset(bipart[1])]
    else:
        return [betaset(bipart[0]), betaset(bipart[1] + [0]*(x-d))]

def symbols(n, d=0):
    """
    Returns a list of all reduced symbols of rank n.

    """
    for bipart in bipartitions(n):
        yield symbol(bipart, d=d)

def symbolrank(symb):
    """
    Get the rank of a symbol.

    """
    return (sum(symb[0]) + sum(symb[1])
            - (len(symb[0]) + len(symb[1]) - 1)**2//4)

def ainvsymbol(symb, n=False):
    """Returns the a-invariant of a symbol, if the rank is known then this may
    be optionally passed. Otherwise the rank is automatically computed.
    """
    # If we don't know the rank then compute this first.
    if n == False:
        n = symbolrank(symb)

    # Get all entries of the symbol in a sorted list.
    allentries = sorted(x for x in itertools.chain(symb[0], symb[1]))

    s = len(allentries)

    return sum((s - i - 1)*(x - i//2) for i, x in enumerate(allentries))


##################################################
#                   Formatting                   #

def intlisttostring(ilist, brackets='()'):
    """This turns a list of integers into a string for displaying in such
    things as character tables. When the list contains an integer >=10 then the
    elements of the list are seperated by commas and surrounded by brackets.
    One may optionally specify the brackets if required.

    >>> intlisttostring([1, 2, 2, 6])
    '1226'
    >>> intlisttostring([1, 2, 10, -4])
    '(1,2,10,-4)'
    >>> intlisttostring([1, 2, 10, -4], brackets='[]')
    '[1,2,10,-4]'
    """
    if all(0 <= i < 10 for i in ilist):
        return ''.join(str(j) for j in ilist)
    else:
        return '{}{}{}'.format(
            brackets[0],
            ','.join(str(j) for j in ilist),
            brackets[1]
        )


def parttupletostring(part, brackets='()'):
    """Takes as argument a list containing lists of integers. It then applies
    intlisttostring to each list and seperates them by a '.'. If one passes the
    optional argument brackets then this is passed to all calls of
    intlisttostring.

    >>> parttupletostring([[1, 2], [7, 7, 8]])
    '12.778'
    >>> parttupletostring([[10, 22], []])
    '(10,22).'
    >>> parttupletostring([[10, 22], [7]], brackets='{}')
    '{10,22}.7'

    """
    return '.'.join(intlisttostring(i, brackets) for i in part)

def arraytostring(array, rowlabels=(), collabels=(), maxwidth=78):
    """
    Creates a string representation of an array which is suitable for printing.
    Setting the maxwidth=0 will mean that the strings width is not capped.

    """
    # Maps elements of the array to strings, except for 0 it returns
    # a ".".  This relies on the fact that the matrix is integral.
    def zerostodots(elm):
        if elm:
            return str(elm)
        else:
            return "."

    mat = [[zerostodots(x) for x in row] for row in array]

    # Get column widths.
    if collabels:
        widths = [len(max(rowlabels, key=len))]
        widths += [max(len(max(row, key=len)), len(collabels[i]))
                    for i, row in enumerate(zip(*mat))]
    else:
        widths = [0]
        widths += [len(max(row, key=len)) for row in zip(*mat)]

    # Justify the columns.
    for i in range(len(array)):
        mat[i] = [s.rjust(j) for s, j in zip(mat[i], widths[1:])]

    # We break the columns of the matrix up into blocks so we don't exceed the
    # specified width.
    blocks = []

    if maxwidth:
        start, end = 0, 0
        blockwidth = widths[0]

        for i, width in enumerate(widths[1:]):
            blockwidth += width + 1

            if i == len(widths)-2:
                if blockwidth <= maxwidth:
                    blocks.append(slice(start, end + 1))
                else:
                    blocks.append(slice(start, end))
                    blocks.append(slice(end, end + 1))
            elif blockwidth <= maxwidth:
                end += 1
            else:
                blocks.append(slice(start, end))
                blockwidth = widths[0] + width + 1
                start = end
                end += 1
    else:
        blocks = [slice(0,len(mat))]

    # First justify all the characternames so we don't keep doing it
    # in the loop below.
    if rowlabels:
        rowlabelsjust = [nam.ljust(widths[0]) for nam in rowlabels]

    # Now sort out the rows.
    rowblocks = []
    for block in blocks:
        printrow = []

        # Header row first.
        if collabels:
            printrow.append(" "*(widths[0] + 1) + " ".join(
                            (s.rjust(widths[block.start + i + 1])
                            for i, s in enumerate(collabels[block]))
                           ))
        if rowlabels:
            printrow.append("-"*(widths[0]
                            + sum(widths[block.start + 1 : block.stop + 1])
                            + block.stop - block.start))

            for i, row in enumerate(mat):
                printrow += [rowlabelsjust[i] + " " + " ".join(row[block])]

        else:
            for i, row in enumerate(mat):
                printrow += [" ".join(row[block])]

        rowblocks.append(printrow)

    return "\n\n".join(("\n".join(row) for row in rowblocks))

##################################################
#                     Arrays                     #

def matrixblocks(X):
    """
    If the rows and columns of the matrix can be permuted so that X has
    a block diagonal form then the result is a list of lists describing
    the corresponding blocks. It is assumed that X[i][j] != 0 if and
    only if X[j][i] != 0.
    
    """
    L = range(len(X))
    orbs = []

    while L!=[]:
        orb = [L[0]]
        for x in orb:
            for i in L:
                if X[x][i] != 0 and not i in orb:
                    orb.append(i)
        for i in orb: 
            L.remove(i)
        orb.sort()
        orbs.append(orb)
    return orbs


def npindex(a, b, cols=False):
    """
    Let a be a numpy array of shape (n,n) for some n>=1 and let b be a
    numpy array of shape (n). This function finds the first occurance of
    b as a row of a. If the optional value cols is set to True then
    columns will be searched instead of rows. If b is not contained in a
    then None is returned.

    """
    n = a.shape[0]
    if cols:
        return next((i for i in range(n) if np.array_equal(a[:,i], b)), None)
    else:
        return next((i for i in range(n) if np.array_equal(a[i], b)), None)


def npperm(a, b, cols=False):
    """
    Let a and b be numpy arrays of shape (n,n) such that a is obtained
    from b by permuting the rows. This function returns a list such that
    the ith row of a is the same as the L[i]th row of b. If the optional
    flag cols is set to True then the same happens but for a permutation
    of columns.

    """
    n = a.shape[0]
    out = [None]*n
    S = set(range(n))
    if cols:
        for i in range(n):
            for j in S:
                if np.array_equal(a[:,i], b[:,j]):
                    out[i] = j
                    break
            S.remove(j)
    else:
        for i in range(n):
            for j in S:
                if np.array_equal(a[i], b[j]):
                    out[i] = j
                    break
            S.remove(j)

    return out

##################################################
#                      Lists                     #

def flatten(L):
    """
    Returns a generator whose list is all elements in L.

    """
    # The following is taken from http://stackoverflow.com/a/2158532.
    for el in L:
        if (isinstance(el, collections.Iterable)
                and not isinstance(el, (str, bytes))):
            yield from flatten(el)
        else:
            yield el

##################################################
#                     System                     #
def getmodule(typ, twist=1):
    return sys.modules['charliepy.data.typ{}{}'.format(twist, typ)]


