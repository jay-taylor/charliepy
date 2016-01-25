########################################################
# Data for type I_2(m) Coxeter groups
#

import numpy as np

def degrees(n, **kwargs):
    return [2, kwargs['m']]

def conjclassdata(ind, **kwargs):
    """
    Returns the conjugacy class data for the ireducible Coxeter groups of type
    I_2(m). The function takes the value m as a keyword argument. For example,
    
            conjclassdata([0, 1], m=8)
    
    gives the data for the Coxeter group of type I_2(8). The data is adapted
    from the corresponding files in GAP-Chevie.
    """
    m = kwargs['m']
    mpar = m % 2

    # Stores the data: representatives, centraliser orders, names.
    if mpar == 0:
        repcentnam = [[[], ind[0:1], ind[1:2]], [2*m, 4, 4],
                      [' ', str(ind[0]), str(ind[1])]]
    elif mpar == 1:
        repcentnam = [[[], ind[0:1]], [2*m, 2], [' ', str(ind[0])]]

    # Create words and names based on indices of the irreducible component.
    for i in range(1, (m - mpar)//2 + 1):
        rep = ind[:]*i
        repcentnam[0] += [rep]
        repcentnam[1] += [m]
        repcentnam[2] += [''.join(str(x) for x in rep)]

    # If m is even then the last representative is in the centre of W.
    if mpar == 0:
        repcentnam[1][-1] = 2*m

    return repcentnam

def irrchardata(n, **kwargs):
    # Note that here n should in fact be such that the irreducible component is
    # of type I_2(n).
    m = kwargs['m']

    if m % 2:
        r = (m-1)//2 + 1
        binv = [0, m] + range(1, r)
        ainv = [0, m] + [1]*(r-1)
        nam = ["phi_{1,0}","phi_{1," + str(m) + "}"] + [
                "phi_{2," + str(j) + "}" for j in range(1, r)]
    else:
        r = (n-1)//2 + 1
        binv = [0, m//2, m//2, m] + range(1, r)
        ainv = [0, 1, 1, m] + [1]*(r-1)
        nam = ["phi_{1,0}", "phi_{1," + str(m//2) + "}'",
               "phi_{1," + str(m//2) + "}''", "phi_{1," + str(m) + "}"] + [
                "phi_{2," + str(j) + "}" for j in range(1, r)]

def chartable(n, **kwargs):
    pass
