#import sys # stdout used in lprint
#import json # dump used in 'writeto'
#
## define our own print
#def lprint(x):
#  sys.stdout.write(x)
#  sys.stdout.flush()
#
#def writeto(fname,l):
##  import simplejson
#  f=open(fname, 'w')
#  json.dump(l,f)
#  f.close()
#  print("finished")
#
##F printfunction
#def printfunction(f):
#  """prints the source code of a function.
#
#  >>> printfunction(transposemat)
#  def transposemat(mat):
#    return list(map(lambda *row: list(row), *mat))
#  """
#  import inspect
#  print(''.join(inspect.getsourcelines(f)[0]))

def mkint(s):
    """Turns any string into an integer, including empty strings.
    """
    s = s.strip()
    return int(s) if s else 0
