# import common functions into the global namespace
from functools import reduce

#F isprime
def isprime(n):
  if n==0 or n==1:
    return False
  if n<0:
    return isprime(-n)
  i=2
  while i*i<=n:
    if n%i==0:
      return False
    i+=1
  return True

#F nextprime
def nextprime(n):
  """returns the next prime number greater than a given integer. 
  (Taken from the gap-library.)

  >>> nextprime(10**8)
  100000007
  >>> isprime(100000007)
  True

  See also 'isprime'.
  """
  if -3==n:
    n=-2
  elif -3<n<2:
    n=2
  elif n%2==0:
    n+=1
  else:
    n+=2
  while not isprime(n):
    if n%6==1:
      n+=4
    else:
      n+=2
  return n

#F intlcm
def intlcm(a,b):
  """returns the least common multiple of two integers. 

  See also 'intlcmlist'.
  """
  gcd,tmp=a,b
  while tmp!=0:
    gcd,tmp=tmp,gcd%tmp
  return (a*b)//gcd

def intlcmlist(l):
  """returns the least common multiple of a list of integers.
  """
  return reduce(intlcm,l)

#F gcdex
def gcdex(m,n):
  """Extended gcd algorithm for integers (taken from the  gap-library). The 
  function returns a dictionary with entries gcd, coeff1, coeff2, coeff3 
  and coeff4. We have
 
      gcd = coeff1 * m + coeff2 * n  and  0 = coeff3 * m + coeff4 * n.

  >>> gcdex(4,15)
  {'coeff3': -15, 'coeff2': -1, 'coeff1': 4, 'gcd': 1, 'coeff4': 4}
 
  (Thus, gcd(4,15)=1; we have 1=4*4+(-1)*15 and 0=(-15)*4+4*15.)
  """
  if 0<=m:
    f,fm=m,1
  else:
    f,fm=-m,-1
  if 0<=n:
    g,gm=n,0
  else:
    g,gm=-n,0
  while g!=0:
    q,h,hm=f//g,g,gm
    g=f-q*g
    gm=fm-q*gm
    f,fm=h,hm
  if n==0:
    return {'gcd':f,'coeff1':fm,'coeff2':0,'coeff3':gm,'coeff4':1}
  else:
    return {'gcd':f,'coeff1':fm,'coeff2':(f-fm*m)//n,'coeff3':gm,
            'coeff4':(0-gm*m)//n}

