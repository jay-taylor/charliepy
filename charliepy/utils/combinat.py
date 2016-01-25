# JT 11/03/14 : Replaced partitions and partitiontuples with more efficient
# versions, which include generator versions. Partitions keeps the same
# ordering as the old function but partitiontuples order has changed

# JT 13/03/14: Added intlisttostring and parttupletostring to obtain names as
# in CHEVIE.

# JT 14/03/14: Removed for loops from ainvbipartition to make it more
# efficient.

import itertools
import operator
import functools

def partitions(n, k=0, reverse=False):
    """This function returns a generator which returns all the integer
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

def betaset(part, frozen=False):
    """Returns a betaset corresponding to the partition part. If part contains
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
    """Returns a pair of beta sets corresponding to the bipartition bipart with
    defect d. If part contains no 0 entries then the betaset will be reduced.
    """
    x = len(bipart[0]) - len(bipart[1])
    if x <= d:
        return [betaset(bipart[0] + [0 for i in range(d-x)]),
                betaset(bipart[1])]
    else:
        return [betaset(bipart[0]),
                betaset(bipart[1] + [0 for i in range(x-d)])]

def symbols(n, d=0):
    """Returns a list of all reduced symbols of rank n.
    """
    for bipart in bipartitions(n):
        yield symbol(bipart, d=d)

def symbolrank(symb):
    """Get the rank of a symbol.
    """
    return sum(symb[0]) + sum(symb[1]) - (len(symb[0])+len(symb[1])-1)**2//4

def partitionsgen(n):
    """Constructs a generator which produces all integer pattitions of n.
    Applying list to this generator gives a list of all integer partitions. This
    algorithm, which is free of recursion, is due to Pfeiffer and is taken 1-1
    from [P94, pg. 167]. The ordering of the partitions is the same as that in
    the GAP library. For example,

    >>> partitions(4)
    <generator object partitions at 0x10367eeb0>
    >>> list(partitions(4))
    [[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]]
    """
    # First deal with the trivial case.
    if n == 0:
        yield []
        return

    # Consider B(k, m) the set of all partitions of k whose maximal entry is at
    # most m. If m >= k then this set is simply the set P(k) of all partitions
    # of k. If we write B(k) for the set B(k, min(k, n-k)) then we have P(n)
    # is the union of all sets
    #
    #      [ [k] + t for t in B(n-k) for k in range(n+1) ]
    #
    # The set B(k, m) is the union of the set B(k, m-1) and
    # 
    #     [ [m] + t for t in B(k-m, m) ]
    #
    # As min(k, n-k) <= n/2 we will need to construct the sets B(k, m) for all
    # 0 <= m <= n/2. This is how we will construct each such set.

    partm = [[] for i in range(n)]

    # This constructs the sets B(m).
    for m in range(1, n//2+1):
        partm[m].append([m])
        for k in range(m+1, n-m+1):
            partm[k] += [ [m] + t for t in partm[k-m] ]

    # Now construct the partitions from the sets B(n-k).
    for k in range(1, n):
        for t in partm[n-k]:
            yield [k] + t

    yield [n]

def dualpartition(mu): #pycox
    """Returns the dual (or conjugate) partition to mu.
    """
    if mu:
        return [len([l for l in mu if l>j]) for j in range(mu[0])]
    else:
        return list()

def centraliserpartition(n, mu): #pycox
    """Returns the order of the centraliser of an element of cycle type of a
    given partition in the full symmetric group. (The program is taken from the
    GAP library and re-written almost 1-1 in python.)
    """
    res, last, k = 1, 0, 1
    for p in mu:
        res *= p
        if p == last:
            k += 1
            res *= k
        else:
            k = 1
        last = p
    return res

def centralisertuple(n, r, mu): #pycox
    """returns the order of the centraliser of an element of a given type
    (specified by an r-tuple of partitions mu) in the wreath product of a cyclic
    group of order r with the full symmetric group of degree n. (The program is
    taken from the GAP library and re-written almost 1-1 in python.)
    """
    res = 1
    for i in range(r):
        last, k = 0, 1

        for p in mu[i]:
            res *= r*p
            if p == last:
                k += 1
                res *= k
            else:
                k = 1
            last = p

    return res

#F ainvbipartition
def ainvsymbol(symb, n=False):
    """Returns the a-invariant of a symbol, if the rank is known then this may
    be optionally passed. Otherwise the rank is automatically computed.
    """
    # If we don't know the rank then compute this first.
    if n == False:
        n = symbolrank(symb)

    # Get all entries of the symbol in a sorted list.
    allentries = [x for x in itertools.chain(symb[0], symb[1])]
    allentries.sort()

    s = len(allentries)

    return sum((s - i - 1)*(x - i//2) for i, x in enumerate(allentries))

##F differencepartitions
#def differencepartitions(gamma,alpha):
#  """returns  a dictionary with information about the  difference of
#  two partitions (if it exists); this function is  needed for the
#  computation of character values in type B_n. It is taken almost 
#  1-1 from the gap-chevie library.
# 
#  See also 'heckevalueB'. 
#  """
#  dp={'cc':0, 'll':0}
#  if len(alpha)>len(gamma):
#    return False
#  old=[]
#  inhook=False
#  alpha=alpha[:]
#  for i in range(len(alpha),len(gamma)):
#    alpha.append(0)
#  for i in range(len(gamma)):
#    if alpha[i]>gamma[i]:
#      return False
#    new=list(range(alpha[i],gamma[i]))
#    intsec=[r for r in old if r in new]
#    if len(intsec)>1:
#      return False
#    elif len(intsec)==1:
#      dp['ll']+=1
#    else:
#      if inhook:
#        dp['cc']+=1
#        dp['d']=old[0]-i
#        inhook=False
#    if new!=[]:
#      inhook=True
#    old=new
#  if inhook:
#    dp['cc']+=1
#    dp['d']=old[0]-len(gamma)
#  return dp

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
    out = str()
    if all((i<10 and i>=0) for i in ilist):
        out += ''.join(str(j) for j in ilist)
    else:
        out += brackets[0] + ','.join(str(j) for j in ilist) + brackets[1]

    return out

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
    if brackets == '()':
        return '.'.join(map(intlisttostring, part))
    else:
        slist = [intlisttostring(i, brackets) for i in part]
        return '.'.join(slist)

##F bipartitions 
#def bipartitions1(n):
#  """returns the list of all bipartitions of n (as in gap). The
#  ordering is different from that of partitiontuples(n,2).
#  """
#  if n==0:
#    return [[[],[]]]
#  pm=[[] for i in range(n)]
#  for m in range(1,n+1):
#    pm[m-1].append([[],[m]])
#    for k in range(m+1,n+1):
#      for t in pm[k-m-1]:
#        s=[[],[m]]
#        s[1].extend(t[1])
#        pm[k-1].append(s)
#  for m in range(1,n//2+1):
#    pm[m-1].append([[m],[]])
#    for k in range(m+1,n-m+1):
#      for t in pm[k-m-1]:
#        s=[[m],t[1]]
#        s[0].extend(t[0])
#        pm[k-1].append(s)
#  res=[]
#  for k in range(1,n):
#    for t in pm[n-k-1]:
#      s=[[k],t[1]]
#      s[0].extend(t[0])
#      res.append(s)
#  res.append([[n],[]])
#  res.extend(pm[n-1])
#  return res;

#def partitiontuples(n, r):
#    return list(partitiontuplesgen(n, r))

#def partitiontuplesgen(n,r):
#  """returns the list of all r-tuples of partitions of n.
#
#  >>>partitiontuples(3,2)
#  [[[1,1,1],[]],
#   [[1,1],[1]],
#   [[1],[1,1]],
#   [[],[1,1,1]],
#   [[2,1],[]],
#   [[1],[2]],
#   [[2],[1]],
#   [[],[2,1]],
#   [[3],[]],
#   [[],[3]]]
#
#  The ordering is exactly the same as in gap. (The program is actually 
#  taken from the gap library and re-written almost 1-1 in python.)
#  
#  See also 'partitions'.
#  """
#  empty={'tup':[[] for x in range(r)],'pos':(n-1)*[1]}
#  if n==0:
#    yield [empty['tup']]
#  pm=[[] for x in range(1,n)]
#  for m in range(1,n//2+1):
#    for i in range(1,r+1):
#      s={'tup':[l[:] for l in empty['tup']],'pos':empty['pos'][:]}
#      s['tup'][i-1]=[m]
#      s['pos'][m-1]=i
#      pm[m-1].append(s)
#    for k in range(m+1,n-m+1):
#      for t in pm[k-m-1]:
#        for i in range(t['pos'][m-1],r+1):
#          t1={'tup':[l[:] for l in t['tup']],'pos':t['pos'][:]}
#          s=[m]
#          s.extend(t['tup'][i-1])
#          t1['tup'][i-1]=s
#          t1['pos'][m-1]=i
#          pm[k-1].append(t1)
#  for k in range(1,n):
#    for t in pm[n-k-1]:
#      for i in range(t['pos'][k-1],r+1):
#        t1=[l[:] for l in t['tup']]
#        s=[k]
#        s.extend(t['tup'][i-1])
#        t1[i-1]=s
#        yield t1
#  for i in range(1,r+1):
#    s=[l[:] for l in empty['tup']]
#    s[i-1]=[n]
#    yield s
  
# This function was written together with Mike Davis.
def fastpartitiontuples(n, r):
    """Returns a generator which gives all posible lists L of length r such
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

##F partitiontuples
#def partitiontuples(n,r):
#  """returns the list of all r-tuples of partitions of n.
#
#  >>>partitiontuples(3,2)
#  [[[1,1,1],[]],
#   [[1,1],[1]],
#   [[1],[1,1]],
#   [[],[1,1,1]],
#   [[2,1],[]],
#   [[1],[2]],
#   [[2],[1]],
#   [[],[2,1]],
#   [[3],[]],
#   [[],[3]]]
#
#  The ordering is exactly the same as in gap. (The program is actually 
#  taken from the gap library and re-written almost 1-1 in python.)
#  
#  See also 'partitions'.
#  """
#  empty={'tup':[[] for x in range(r)],'pos':(n-1)*[1]}
#  if n==0:
#    return [empty['tup']]
#  pm=[[] for x in range(1,n)]
#  for m in range(1,n//2+1):
#    for i in range(1,r+1):
#      s={'tup':[l[:] for l in empty['tup']],'pos':empty['pos'][:]}
#      s['tup'][i-1]=[m]
#      s['pos'][m-1]=i
#      pm[m-1].append(s)
#    for k in range(m+1,n-m+1):
#      for t in pm[k-m-1]:
#        for i in range(t['pos'][m-1],r+1):
#          t1={'tup':[l[:] for l in t['tup']],'pos':t['pos'][:]}
#          s=[m]
#          s.extend(t['tup'][i-1])
#          t1['tup'][i-1]=s
#          t1['pos'][m-1]=i
#          pm[k-1].append(t1)
#  res=[]
#  for k in range(1,n):
#    for t in pm[n-k-1]:
#      for i in range(t['pos'][k-1],r+1):
#        t1=[l[:] for l in t['tup']]
#        s=[k]
#        s.extend(t['tup'][i-1])
#        t1[i-1]=s
#        res.append(t1)
#  for i in range(1,r+1):
#    s=[l[:] for l in empty['tup']]
#    s[i-1]=[n]
#    res.append(s)
#  return res

##F lusztigsymbolB
#def lusztigsymbolB(n,vs,vt,dblpart):
#  """returns the symbol associated with a  bipartition, as defined
#  by Lusztig,  taking into account  weights.  In this form, the 
#  total number of entries in a symbol only depends on n and the 
#  parameters (but not on the bipartition).
#
#  >>> bipartitions(2)
#  [[[1], [1]], [[1, 1], []], [[2], []], [[], [1, 1]], [[], [2]]]
#  >>> [lusztigsymbolB(2,1,1,pi) for pi in bipartitions(2)]
#  [[[0, 1, 3], [0, 2]], 
#   [[0, 2, 3], [0, 1]], 
#   [[0, 1, 4], [0, 1]], 
#   [[0, 1, 2], [1, 2]], 
#   [[0, 1, 2], [0, 3]]]
#
#  See also 'redlusztigsymbolB'.
#  """
#  q,r=vt//vs,vt%vs
#  a,b=dblpart[0][:],dblpart[1][:]
#  if len(a)>len(b)+q:
#    b.extend((len(a)-len(b)-q)*[0])
#  elif len(b)+q>len(a):
#    a.extend((len(b)-len(a)+q)*[0])
#  while len(a)+len(b)<2*n+q:
#    a.append(0)
#    b.append(0)
#  la,mu=(n+q)*[0],n*[0]
#  for i in range(1,len(b)+1):
#    la[i-1]=vs*(a[len(a)-i]+i-1)+r
#    mu[i-1]=vs*(b[len(b)-i]+i-1)
#  for i in range(len(b)+1,len(a)+1):
#    la[i-1]=vs*(a[len(a)-i]+i-1)+r
#  return [la,mu]

##F redlusztigsymbolB
#def redlusztigsymbolB(vs,vt,dblpart):
#  """similar to 'lusztigsymbolB' but now the number of entries in a 
#  symbol is as small as possible (depending on the bipartition).
#
#  >>> bipartitions(2)
#  [[[1], [1]], [[1, 1], []], [[2], []], [[], [1, 1]], [[], [2]]]
#  >>> [redlusztigsymbolB(1,1,pi) for pi in bipartitions(2)]
#  [[[0, 2], [1]], 
#   [[1, 2], [0]], 
#   [[2], []], 
#   [[0, 1, 2], 
#   [1, 2]], [[0, 1], [2]]]
#
#  See also 'lusztigsymbolB'.
#  """
#  q,r=vt//vs,vt%vs
#  a,b=dblpart[0][:],dblpart[1][:]
#  if len(a)>len(b)+q:
#    b.extend((len(a)-len(b)-q)*[0])
#  elif len(b)+q>len(a):
#    a.extend((len(b)-len(a)+q)*[0])
#  n=(len(a)+len(b)-q)//2
#  la,mu=(n+q)*[0],n*[0]
#  for i in range(1,len(b)+1):
#    la[i-1]=vs*(a[len(a)-i]+i-1)+r
#    mu[i-1]=vs*(b[len(b)-i]+i-1)
#  for i in range(len(b)+1,len(a)+1):
#    la[i-1]=vs*(a[len(a)-i]+i-1)+r
#  return [la,mu]



##F ainvbipartition
#def ainvbipartition(n,vs,vt,bip):
#  """returns the a-invariant of a  bipartition,  computed from 
#  the associated Lusztig symbol. See also 'lusztigsymbolB'.
#  """
#  q,r=vt//vs,vt%vs
#  p=[]
#  s=lusztigsymbolB(n,vs,vt,bip)
#  for i in s[0]:
#    p.append(i)
#  for i in s[1]:
#    p.append(i)
#  p.sort()
#  N=(len(p)-q)//2
#  z0=[i*vs for i in range(N)]
#  z0.extend([i*vs+r for i in range(N+q)])
#  z0.sort()
#  return sum(i*(p[-i-1]-z0[-i-1]) for i in range(len(p)))

