########################################################
# Data for type I_2(m) Coxeter groups
#

import numpy as np

def degrees(n, **kwargs):
    return [2, kwargs['m']]

def conjclassdata(ind, **kwargs):
    m = kwargs['m']
    s, t = str(ind[0]), str(ind[1])

    if m % 2:
        yield [], 2*m, ' '
        yield ind[0:1], 2, s

        for i in range(1, (m + 1)//2):
            yield ind[:]*i, m, '.'.join((s, t)*i)
    else:
        yield [], 2*m, ' '
        yield ind[0:1], 4, s
        yield ind[1:2], 4, t

        for i in range(1, m//2):
            yield ind[:]*i, m, '.'.join((s, t)*i)

        yield ind[:]*(m//2), 2*m, '.'.join((s, t)*i)

def conjclassdata_min(ind, **kwargs):
    m = kwargs['m']
    s, t = str(ind[0]), str(ind[1])

    if m % 2:
        yield 2*m, ' '
        yield 2, s

        for i in range(1, (m + 1)//2):
            yield m, '.'.join((s, t)*i)

    else:
        yield (2*m, ' ')
        yield (4, s)
        yield (4, t)

        for i in range(1, m//2):
            yield (m, '.'.join((s, t)*i))

        yield (2*m, '.'.join((s, t)*i))

def irrchardata(n, **kwargs):
    # For the a-values and b-values see [GP00, 6.5.10].
    m = kwargs['m']

    if m % 2:
        # Linear characters first.
        yield ("phi_{1,0}", 0, 0)
        yield ("phi_{{1,{}}}".format(m), m, m)

        # Now the degree 2 characters.
        for j in range(1, (m-1)//2 + 1):
            yield ("phi_{{2,{}}}".format(j), 1, j)

    else:
        # Linear characters first.
        yield ("phi_{1,0}", 0, 0)
        yield ("phi_{{1,{}}}'".format(m//2), 1, m//2)
        yield ("phi_{{1,{}}}''".format(m//2), 1, m//2)
        yield ("phi_{{1,{}}}".format(m), m, m)

        # Now the degree 2 characters.
        for j in range(1, m//2):
            yield ("phi_{{2,{}}}".format(j), 1, j)

def chartable(n, **kwargs):
    pass
