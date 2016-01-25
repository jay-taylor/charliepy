import numpy as np

#from misc import lprint

##F idmat
#def idmat(rng,scalar):
#  """returns the scalar matrix of size len(rng) with a given scalar 
#  on the diagonal.
#
#  >>> idmat([0,1,2],-3)
#  [[-3, 0, 0], [0, -3, 0], [0, 0, -3]]
#  """
#  m=[len(rng)*[0] for x in rng]
#  for x in range(len(rng)): m[x][x]=scalar
#  return m
#
##F transposemat
#def transposemat(mat):
#  """returns the transpose of a matrix.
#  
#  >>> transposemat([[1,2,3],[4,5,6]])
#  [[1,4], [2,5], [3,6]]
#  """
#  return list(map(lambda *row: list(row), *mat))
##  return zip(*mat)
#
##F flatblockmat
#def flatblockmat(blmat):
#  """flattens a block matrix.
#  """
#  a=[]
#  for b in range(len(blmat)):
#    for i in range(len(blmat[b][b])):
#      a.append(flatlist(l[i] for l in blmat[b]))
#  return a
#
##F matadd
#def matadd(a,b):
#  """returns the sum of two matrices.
#  """
#  return [[a[i][j]+b[i][j] for j in range(len(a[0]))] for i in range(len(a))]
#
##F matsub
#def matsub(a,b):
#  """returns the difference of two matrices.
#  """
#  return [[a[i][j]-b[i][j] for j in range(len(a[0]))] for i in range(len(a))]
#
##F matmult
#def matmult(a,b):
#  """returns the matrix product of the matrices a and b.
#
#  See also 'matadd' and 'scalmatmult'.
#  """
#  return [[sum(row[k]*b[k][j] for k in range(len(b)))
#                        for j in range(len(b[0]))] for row in a]
#
##F scalmatmult
#def scalmatmult(a,b):
#  """multiplies a matrix b with a scalar a.
#  """
#  return [[a*b[i][j] for j in range(len(b[0]))] for i in range(len(b))]
#
##F directsummat
#def directsummat(a,b):
#  """returns the matrix direct sum of the matrices a and b.
#
#  >>> c=directsummat(cartanmat("A",2),cartanmat("G",2))
#  [[2,-1,0,0],[-1,2,0,0],[0,0,2,-1],[0,0,-3,2]]
#  """
#  if a==[[]]: 
#    return b
#  elif b==[[]]: 
#    return a
#  else:
#    c=idmat(range(len(a[0])+len(b[0])),0)
#    for i in range(len(a[0])):
#      for j in range(len(a[0])):
#        c[i][j]=a[i][j]
#    for i in range(len(b[0])):
#       for j in range(len(b[0])):
#         c[len(a[0])+i][len(a[0])+j]=b[i][j]
#    return c
#
##F kroneckerproduct
#def kroneckerproduct(mat1,mat2):
#  """returns the  Kronecker product of the matrices mat1 and mat2. If 
#  mat1  has size m times n and mat2  has  size p times q, then the 
#  Kronecker product is a matrix of size m*p times n*q.
#
#  >>> mat1=[[ 0,-1, 1],
#            [-2, 0,-2]]
#  >>> mat2=[[1,1],
#            [0,1]]
#  >>> kroneckerproduct(mat1,mat2)
#  [[ 0, 0,-1,-1, 1, 1], 
#   [ 0, 0, 0,-1, 0, 1], 
#   [-2,-2, 0, 0,-2,-2], 
#   [ 0,-2, 0, 0, 0,-2]]
#
#  (The program is taken from the gap library and re-written almost 
#  1-1 in python.)
#  """
#  krpr=[]
#  for row1 in mat1:
#    for row2 in mat2:
#      row=[]
#      for i in row1:
#        row.extend([i*x for x in row2])
#      krpr.append(row)
#  return krpr
#
##F decomposemat
#def decomposemat(mat):
#  """tests if a matrix can be decomposed in block diagonal form; it 
#  is assumed that  mat[i][j]=0  if and only if  mat[j][i]=0. The 
#  result is range(len(mat[0]))  if  mat  can not be  decomposed. 
#  Otherwise, the lists of block indices are returned. 
#
#  >>>d=decomposemat([[ 2, 0,-1, 0, 0], 
#                     [ 0, 2, 0,-3, 0], 
#                     [-2, 0, 2, 0,-1], 
#                     [ 0,-1, 0, 2, 0], 
#                     [ 0, 0,-1, 0, 2]]
#  [[0, 2, 4], [1, 3]]
#
#  Thus, there are two blocks, obtained by taking the submatrices 
#  with row and colum indices (0,2,4) and (1,3), respectively.
#  """
#  l=list(range(len(mat[0])))
#  orbs=[]
#  while l!=[]:
#    orb=[l[0]]
#    for o in orb:
#      for i in l:
#        if mat[o][i]!=0 and not i in orb:
#          orb.append(i)
#    for i in orb: 
#      l.remove(i)
#    orb.sort()
#    orbs.append(orb)
#  return orbs
#
##F determinantmat1
#def determinantmat1(mat):
#  n=len(mat[0])
#  if n==1:
#    return mat[0][0]
#  else:
#    d=0
#    for k in range(n):
#      l=list(range(n))
#      l.remove(k)
#      if mat[k][0]!=0:
#        d+=(-1)**k*mat[k][0]*determinantmat1([[mat[i][j] 
#                           for j in range(1,n)] for i in l])
#    return d
#
##F determinantmat
#def determinantmat(mat):
#  """returns the determinant of a matrix. If all coefficients are integers,
#  the function uses a simplified version of the algorithm for  computing  
#  the elementary divisors of a matrix; in particular,  it does  not  use 
#  fractions. In general, it uses induction on the size of the matrix and 
#  expansion  along the first column.  (Thus,  it will  work for matrices 
#  of moderate size over any commutative ring).
#  """
#  a=[list(l) for l in mat]
#  if not all(type(x)==type(0) for x in flatlist(a)):
#    return determinantmat1(mat)
#  n=len(mat[0])
#  d=1
#  for p in range(n-1):
#    i=p
#    while i<n and a[i][p]==0:
#      i+=1
#    if i>=n:
#      return 0
#    if i!=p:
#      d=-d
#      for j in range(n):
#        x=a[p][j]
#        a[p][j]=a[i][j]
#        a[i][j]=x
#    fertig=False
#    i=p+1
#    while i<n and not fertig:
#      if a[i][p]!=0:
#        q=a[i][p]//a[p][p]
#        for j in range(n):
#          a[i][j]-=q*a[p][j]
#        if a[i][p]!=0:
#          for j in range(n):
#            x=a[p][j]
#            a[p][j]=a[i][j]
#            a[i][j]=x
#          d=-d
#        else:
#          i+=1
#      else:
#        i+=1
#    d*=a[p][p]
#  d*=a[n-1][n-1]
#  return d
#
##F inversematp
#def inversematp(mat,p):
#  """checks if an integer matrix is invertible; if this is the case, the 
#  function  returns the inverse of that matrix  modulo a prime number.
#  """
#  n=len(mat[0])
#  a=[]
#  for i in range(n):
#    l=[mat[i][j]%p for j in range(n)]
#    l.extend(n*[0])
#    l[n+i]=1
#    a.append(l)
#  for k in range(n):
#    k1=k
#    while k1<n and a[k1][k]%p==0:
#      k1+=1
#    if k1==n:
#      return False
#    if k!=k1:
#      for j in range(2*n):
#        a[k][j],a[k1][j]=a[k1][j],a[k][j]
#    p1=1
#    while (p1*a[k][k])%p!=1:
#      p1+=1
#    for j in range(2*n):
#      a[k][j]=(p1*a[k][j])%p
#    for i in range(n): 
#      if i!=k:
#        q=a[i][k]
#        for j in range(2*n):
#          a[i][j]=(a[i][j]-(q*a[k][j])%p)%p
#  return [l[n:] for l in a]
#
##F displaymat
#def displaymat(mat,rows=[],cols=[],width=78):
#  """displays a matrix, where the optional arguments 'rows' and 'cols'
#  can be used to specify labels for the rows and columns.  There is 
#  a further  optional  argument by which one can set the 'width' of 
#  the display, i.e., the maximum number of characters printed  (the 
#  default value is 78 characters per line).
#
#  >>> displaymat(chartable(coxeter("H",3))['irreducibles'])
#   1 -1     1  1  1     -1     1 -1     -1 -1
#   1  1     1  1  1      1     1  1      1  1
#   5 -1     .  1 -1      .     .  1      . -5
#   5  1     .  1 -1      .     . -1      .  5
#   3 -1   ir5 -1  .  1-ir5 1-ir5  .    ir5  3
#   3 -1 1-ir5 -1  .    ir5   ir5  .  1-ir5  3
#   3  1   ir5 -1  . -1+ir5 1-ir5  .   -ir5 -3
#   3  1 1-ir5 -1  .   -ir5   ir5  . -1+ir5 -3
#   4  .    -1  .  1      1    -1 -1      1 -4
#   4  .    -1  .  1     -1    -1  1     -1  4
#
#  (Values equal to 0 are printed as '.'.)
#
#  See also 'displaychartable'.
#  """
#  m=len(mat)
#  n=len(mat[0])
#  csp=[max([len(repr(mat[i][j])) for i in range(m)]) for j in range(n)]
#  if cols!=[]:
#    csp=[max(len(str(cols[j])),csp[j])+1 for j in range(n)]
#  else:
#    csp=[csp[j]+1 for j in range(n)]
#  if rows!=[]:
#    maxr=max([len(str(r)) for r in rows])
#  else:
#    maxr=0
#  co=0
#  cut=[0]
#  for j in range(len(csp)):
#    if co+csp[j]<=width-maxr:
#      co+=csp[j]
#    else:
#      cut.append(j)
#      co=csp[j]
#  if cut[:-1]!=n:
#    cut.append(n)
#  for k in range(len(cut)-1): 
#    if cols!=[]:
#      lprint(width*'-')
#      lprint('\n')
#      lprint(maxr*' ') 
#      for j in range(cut[k],cut[k+1]):
#        lprint((csp[j]-len(str(cols[j])))*' '+str(cols[j]))
#      lprint('\n')
#      lprint(width*'-')
#      lprint('\n')
#    for i in range(m):
#      if rows!=[]:
#        lprint((maxr-len(str(rows[i])))*' '+str(rows[i]))
#      for j in range(cut[k],cut[k+1]):
#        if mat[i][j]==0:
#          lprint((csp[j]-len(repr(mat[i][j])))*' '+'.')
#        elif type(mat[i][j])==type(0):
#          lprint((csp[j]-len(repr(mat[i][j])))*' '+str(mat[i][j]))
#        else:
#          lprint((csp[j]-len(repr(mat[i][j])))*' '+repr(mat[i][j]))
#      lprint('\n')
#    lprint('\n')
#  return None

def npindex(a, b, cols=False):
    """Let a be a numpy array of shape (n,n) for some n>=1 and let b be a
    numpy array of shape (n). This function finds the first occurance of b as a
    row of a. If the optional value cols is set to True then columns will be
    searched instead of rows. If b is not contained in a then None is returned.
    """
    n = a.shape[0]
    if cols:
        return next((i for i in range(n) if np.array_equal(a[:,i], b)), None)
    else:
        return next((i for i in range(n) if np.array_equal(a[i], b)), None)

def npperm(a, b, cols=False):
    """Let a and b be numpy arrays of shape (n,n) such that a is obtained from
    b by permuting the rows. This function returns a list such that the ith row
    of a is the same as the L[i]th row of b. If the optional flag cols is set
    to True then the same happens but for a permutation of columns.
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






