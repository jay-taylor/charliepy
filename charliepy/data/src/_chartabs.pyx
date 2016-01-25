cimport cython

def charcolA(list schm, list t, unsigned int k):
    cdef list col, pi
    cdef unsigned int ind
    cdef int val, i

    col = [None]*len(schm)
    for ind, pi in enumerate(schm):
        val = 0
        for i in pi[k]:
            if i < 0:
                val -= <int>(t[<unsigned int>(-i-1)])
            else:
                val += <int>(t[<unsigned int>(i-1)])
        col[ind] = val

    return col

def charcolB(list schm, list t, unsigned int k, unsigned int p, unsigned int q):
    cdef list col, pi
    cdef unsigned int ind, swp
    cdef int val, j

    col = [None]*len(schm)
    for ind, pi in enumerate(schm):
        val = 0
        for j, swp in pi[0][k]:
            if j < 0:
                if swp and q:
                    val += <int>(t[<unsigned int>(-j-1)])
                else:
                    val -= <int>(t[<unsigned int>(-j-1)])
            else:
                if swp and q:
                    val -= <int>(t[<unsigned int>(j-1)])
                else:
                    val += <int>(t[<unsigned int>(j-1)])
        for j in pi[1][k]:
            if j < 0:
                if p:
                    val += <int>(t[<unsigned int>(-j-1)])
                else:
                    val -= <int>(t[<unsigned int>(-j-1)])
            else:
                if p:
                    val -= <int>(t[<unsigned int>(j-1)])
                else:
                    val += <int>(t[<unsigned int>(j-1)])
        col[ind] = val
    return col


