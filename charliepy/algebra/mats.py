import numpy as np

class Matrix(np.ndarray):
    def __new__(cls, obj, dtype=None, copy=True, order=None, subok=False,
            ndmin=0):
        return np.array(obj, dtype=dtype, copy=copy, order=order, subok=subok,
                ndmin=ndmin).view(cls)

    def __repr__(self):
        return "Matrix(\n" + np.array2string(self, separator=', ') + "\n)"

    def __str__(self):
        return np.array2string(self, separator=', ')

    #def __eq__(self, other):
    #    try:
    #        return np.array_equal(self, other)
    #    except:
    #        return False

def equalmatrices(C, D):
    """Decides if two matrices are the same, i.e. same shape and entries.
    """
    return np.array_equal(C, D)

def zeromatrix(n, dtype='int'):
    """Return an nxn matrix filled with 0s. The default datatype is int.
    """
    return Matrix(np.zeros((n,n), dtype=dtype), copy=False)

def identitymatrix(n, dtype='int'):
    """Return an nxn matrix filled with 1s. The default datatype is int.
    """
    return Matrix(np.identity(n, dtype=dtype), copy=False)

def matrixblocks(X):
    """If the rows and columns of the matrix can be permuted so that X has a
    block diagonal form then the result is a list of lists describing the
    corresponding blocks. It is assumed that X[i][j] != 0 if and only if X[j][i]
    != 0.
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

