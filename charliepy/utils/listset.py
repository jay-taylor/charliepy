import itertools
import collections

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

##F flatlist
#def flatlist(lst):
#  """returns the list of all elements that are contained in a given
#  list or in any of its sublists. (Taken from the gap library.)
#
#  >>> flatlist([1,[2,3],[[1,2],3]])
#  [1,2,3,1,2,3]
#  >>> flatlist([]);
#  []
#  """
#  flt=[]
#  for elm in lst:
#    if type(elm)!=type([]):
#      flt.append(elm)
#    else:
#      flt.extend(flatlist(elm))
#  return flt

##F transclos
#def transclos(n,r):
#  """returns the transitive closure of a relation on the integers
#  0,1,...,n-1 given by a list of pairs in r.
#  """
#  m=[]
#  for i in range(n):
#    l=n*[False]
#    l[i]=True
#    for p in r:
#      if p[0]==i:
#        l[p[1]]=True
#    m.append(l)
#  for i in range(n):
#    for j in range(n):
#      if m[j][i]==1:
#        m[j]=[m[i][k] or m[j][k] for k in range(n)]
#  return m
 
##F noduplicates
#def noduplicates(seq):
#  """returns a new list in which all duplicates in the original list
#  have been removed.  Also works for non-hashable types. (Learned  
#  this from stackoverflow.)
#
#  >>> noduplicates([[1, 1], [2, 1], [1, 2], [1, 1], [3, 1]])
#  [[1, 1], [2, 1], [1, 2], [3, 1]]
#  """
#  seen = set()
#  return [x for x in seq if str(x) not in seen and not seen.add(str(x))] 

##F partitioncap
#def partitioncap(p1,p2):
#  """returns the common refinement of two partitions of the same set.
#  """
#  s2=[set([str(x) for x in l]) for l in p2]
#  neu=[]
#  for p in p1:
#    s=[str(x) for x in p]
#    for j in range(len(s2)):
#      nn=[p[i] for i in range(len(s)) if s[i] in s2[j]]
#      if nn!=[]:
#        neu.append(nn)
#  return neu

##F permmult
#def permmult(p,q):
#  """returns the composition of two permutations (acting from the
#  right.)
#  """
#  return tuple([q[i] for i in p])

##F perminverse
#def perminverse(p):
#  """returns the inverse of a permutation.
#  """
#  np=len(p)*[0]
#  for i in range(len(p)):
#    np[p[i]]=i
#  return tuple(np)
 
##F cartesian2
#def cartesian2(liste,n,tup,i):
#  if i==n:
#    tups=[tup[:]]
#  else:
#    tups=[]
#    for l in liste[i]:
#      if i==len(tup):
#        tup.append(l)
#      else:
#        tup[i]=l
#      tups.extend(cartesian2(liste,n,tup,i+1))
#  return tups

# This is equivalent to itertools.product. In particular, the itertools
# cartesian product of two sets A and B is equivalent to the list comprehension
#
#   [(x, y) for x in A for y in B]
#
# Here the for loops cycle like an odometer, with the rightmost loop advancing
# on every iteration.

##F cartesian
#def cartesian(*arg):
#  """returns the cartesian product of lists. 
#
#  >>> cartesian([1,2],[3,4],[4,5])
#  [[1,3,4],[1,3,5],[1,4,4],[1,4,5],[2,3,4],[2,3,5],[2,4,4],[2,4,5]]
#  >>> cartesian([1,2,2],[1,1,2])
#  [[1,1],[1,1],[1,2],[2,1],[2,1],[2,2],[2,1],[2,1],[2,2]]
#
#  In the first form the argument is a comma-separated sequence l1, l2,
#  ..., and the function returns the cartesian product of l1, l2, ...
#    
#  In the second form the argument is a list of lists [l1,l2,,...], and 
#  and the function returns the cartesian product of those lists.
#
#  If more than two lists are given,  cartesian(l1,l2,...) is the  same 
#  (up to some nested bracketing) as cartesian(cartesian(l1,l2),...).
#
#  >>> cartesian(cartesian([1,2],[3,4]),[4,5])
#  [[1,3],4],[[1,3],5],[[1,4],4],[[1,4],5],[[2,3],4],
#                                  [[2,3],5],[[2,4],4],[[2,4],5]]
#
#  The ordering is exactly the same as in gap. (The program is actually
#  taken from the gap library and re-written almost 1-1 in python.)
#  """
#  if len(arg)==1:
#    return cartesian2(arg[0],len(arg[0]),[],0)
#  else:
#    return cartesian2(arg,len(arg),[],0)
#
