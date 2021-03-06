########################################################
# Data for type F_4 Coxeter groups and reductive groups
# The Dynkin diagram here follows Bourbaki so we have
#
#                  0   1   2   3
#                  o---o=>=o---o

import numpy as np
from . import rawdata as _rd
from .. import permutat

# We use the standard permutation generators below in class recognition.
_F4permgens = [permutat.cyclestoperm(*_rd.F4permgens[0]),
               permutat.cyclestoperm(*_rd.F4permgens[1]),
               permutat.cyclestoperm(*_rd.F4permgens[2]),
               permutat.cyclestoperm(*_rd.F4permgens[3])]

def cartanmat(n):
    if n == 4:
        return np.array([[ 2, -1,  0,  0],
                         [-1,  2, -2,  0],
                         [ 0, -1,  2, -1],
                         [ 0,  0, -1,  2]],
                         dtype = 'int8')
    else:
        raise ValueError("Rank of F_4 must be 4!")

def rootlengths(n, **kwargs):
    """
    Returns a tuple giving the relative root lengths.

    """
    return (2, 2, 1, 1)

def diagram(inds):
    """Prints the Dynkin diagram."""
    print("F4 :  {} - {} > {} - {}".format(*inds))
    return None

def degrees(n):
    return [2, 6, 8, 12]

def longestword(inds):
    return [inds[k] for k in [0, 2, 1, 3]]*6

def maxparachain(inds):
    return [inds[i] for i in (1, 2, 0, 3)]

def coxeterclasses(inds):
    return (
        (nam, orblen, norm, [inds[i] for i in rep])
        for (nam, orblen, norm, rep) in _rd.F4coxeterclasses
    )

def conjclasses(inds, **kwargs):
    return (
        (nam, cent, [inds[i] for i in rep])
        for (nam, cent, rep) in _rd.F4conjclasses
    )

def conjclasses_min(inds, **kwargs):
    return (cls[:2] for cls in _rd.F4conjclasses)

def wordtoclass(n, w):
    """
    Returns the name of the conjugacy class in the Weyl group of type F_4
    containing the element w given as a word in the standard generators.

    """
    # First convert the word to a permutation.
    perm = permutat.Perm(range(48))
    for i in w:
        perm *= _F4permgens[i]

    # Each class is uniquely determined by the pair of cycletypes of the
    # restricted permutation on each orbit.
    param = (permutat.restrictedperm(perm, _rd.F4rootorbits[0]).cyclestructure,
             permutat.restrictedperm(perm, _rd.F4rootorbits[1]).cyclestructure)

    return _rd.F4classparams[param]

def irrchars(n, labels=(), **kwargs):
    if "kondo" in labels:
        return (
            ('{d}{p}'.format(d=deg, p=kondo), a, b)
            for (deg, a, b, prim, kondo) in _rd.F4irrchardata
        )
    else:
        return (
            ('phi{{{d},{b}}}{p}'.format(d=deg, b=b, p=prim), a, b)
            for (deg, a, b, prim, kondo) in _rd.F4irrchardata
        )

def chartable(n, **kwargs):
    return _rd.F4chartable

#def conjclassdata(ind, **kwargs):
#    classes = (
#        ([], 1152, 'A0'),
#        ([0, 1, 0, 2, 1, 0, 2, 1, 2, 3, 2, 1, 0, 2, 1, 2, 3, 2, 1, 0, 2, 1, 2,
#            3], 1152, '4A1'),
#        ([1, 2, 1, 2], 64, '2A1'),
#        ([1, 0], 36, 'A2'),
#        ([0, 1, 2, 1, 2, 3, 2, 1, 2, 3], 36, 'D4'),
#        ([0, 1, 0, 2, 1, 0, 2, 3, 2, 1, 2, 3], 96, 'D4(a1)'),
#        ([3, 2], 36, '~A2'),
#        ([0, 1, 0, 2, 1, 0, 2, 1, 2, 3], 36, 'C3+A1'),
#        ([0, 1, 0, 2, 1, 0, 2, 3, 2, 1, 0, 2, 1, 2, 3, 2], 72, 'A2+~A2'),
#        ([0, 1, 0, 2, 1, 2, 3, 2], 72, 'F4(a1)'),
#        ([0, 1, 2, 3], 12, 'F4'),
#        ([0], 96, 'A1'),
#        ([1, 2, 1, 2, 3, 2, 1, 2, 3], 96, '3A1'),
#        ([0, 3, 2], 12, '~A2+A1'),
#        ([3, 2, 1], 12, 'C3'),
#        ([1, 2, 1, 0, 2], 16, 'A3'),
#        ([2], 96, '~A1'),
#        ([0, 1, 0, 2, 1, 0, 2, 1, 2], 96, '2A1+~A1'),
#        ([1, 0, 3], 12, 'A2+~A1'),
#        ([2, 1, 0], 12, 'B3'),
#        ([1, 3, 2, 1, 2], 16, 'B2+A1'),
#        ([0, 2], 16, 'A1+~A1'),
#        ([2, 1], 32, 'B2'),
#        ([0, 1, 0, 2, 1, 0, 2, 1, 2, 3, 2, 1, 2, 3], 32, 'A3+~A1'),
#        ([0, 1, 2, 1, 3, 2], 8, 'B4')
#    )
#
#    return (([ind[j] for j in cls[0]], cls[1], cls[2]) for cls in classes)
#
#def conjclassdata_min(ind, **kwargs):
#    return (
#        (1152, 'A0'), (1152, '4A1'), (64, '2A1'), (36, 'A2'), (36, 'D4'),
#        (96, 'D4(a1)'), (36, '~A2'), (36, 'C3+A1'), (72, 'A2+~A2'),
#        (72, 'F4(a1)'), (12, 'F4'), (96, 'A1'), (96, '3A1'),
#        (12, '~A2+A1'), (12, 'C3'), (16, 'A3'), (96, '~A1'),
#        (96, '2A1+~A1'), (12, 'A2+~A1'), (12, 'B3'), (16, 'B2+A1'),
#        (16, 'A1+~A1'), (32, 'B2'), (32, 'A3+~A1'), (8, 'B4')
#    )
#    
#def irrchardata(n, **kwargs):
#    return (
#        ('1_1', 0, 0), ('1_2', 4, 12), ('1_3', 4, 12), ('1_4', 24, 24),
#        ('2_1', 1, 4), ('2_2', 13, 16), ('2_3', 1, 4), ('2_4', 13, 16),
#        ('4_1', 4, 8), ('9_1', 2, 2), ('9_2', 4, 6), ('9_3', 4, 6),
#        ('9_4', 10, 10), ('6_1', 4, 6), ('6_2', 4, 6), ('12', 4, 4),
#        ('4_2', 1, 1), ('4_3', 4, 7), ('4_4', 4, 7), ('4_5', 13, 13),
#        ('8_1', 3, 3), ('8_2', 9, 9), ('8_3', 3, 3), ('8_4', 9, 9),
#        ('16', 4, 5)
#    )
#    
#
#def chartable(n, **kwargs):
#    return [
#        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#            1], 
#        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1,
#            -1, -1, -1], 
#        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1,
#            -1, -1, -1], 
#        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
#            -1, 1, 1, 1, 1], 
#        [2, 2, 2, 2, 2, 2, -1, -1, -1, -1, -1, 2, 2, -1, -1, 2, 0, 0, 0, 0, 0,
#            0, 0, 0, 0], 
#        [2, 2, 2, 2, 2, 2, -1, -1, -1, -1, -1, -2, -2, 1, 1, -2, 0, 0, 0, 0, 0,
#            0, 0, 0, 0], 
#        [2, 2, 2, -1, -1, 2, 2, 2, -1, -1, -1, 0, 0, 0, 0, 0, 2, 2, -1, -1, 2,
#            0, 0, 0, 0], 
#        [2, 2, 2, -1, -1, 2, 2, 2, -1, -1, -1, 0, 0, 0, 0, 0, -2, -2, 1, 1, -2,
#            0, 0, 0, 0], 
#        [4, 4, 4, -2, -2, 4, -2, -2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#            0, 0, 0], 
#        [9, 9, 1, 0, 0, -3, 0, 0, 0, 0, 0, 3, 3, 0, 0, -1, 3, 3, 0, 0, -1, 1, 1,
#            1, -1], 
#        [9, 9, 1, 0, 0, -3, 0, 0, 0, 0, 0, 3, 3, 0, 0, -1, -3, -3, 0, 0, 1, -1,
#            -1, -1, 1], 
#        [9, 9, 1, 0, 0, -3, 0, 0, 0, 0, 0, -3, -3, 0, 0, 1, 3, 3, 0, 0, -1, -1,
#            -1, -1, 1], 
#        [9, 9, 1, 0, 0, -3, 0, 0, 0, 0, 0, -3, -3, 0, 0, 1, -3, -3, 0, 0, 1, 1,
#            1, 1, -1], 
#        [6, 6, -2, 0, 0, 2, 0, 0, 3, 3, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2,
#            -2, 0], 
#        [6, 6, -2, 0, 0, 2, 0, 0, 3, 3, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 2,
#            2, 0], 
#        [12, 12, -4, 0, 0, 4, 0, 0, -3, -3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#            0, 0, 0], 
#        [4, -4, 0, 1, -1, 0, 1, -1, -2, 2, 0, 2, -2, -1, 1, 0, 2, -2, -1, 1, 0,
#            0, 2, -2, 0], 
#        [4, -4, 0, 1, -1, 0, 1, -1, -2, 2, 0, 2, -2, -1, 1, 0, -2, 2, 1, -1, 0,
#            0, -2, 2, 0], 
#        [4, -4, 0, 1, -1, 0, 1, -1, -2, 2, 0, -2, 2, 1, -1, 0, 2, -2, -1, 1, 0,
#            0, -2, 2, 0], 
#        [4, -4, 0, 1, -1, 0, 1, -1, -2, 2, 0, -2, 2, 1, -1, 0, -2, 2, 1, -1, 0,
#            0, 2, -2, 0], 
#        [8, -8, 0, 2, -2, 0, -1, 1, 2, -2, 0, 4, -4, 1, -1, 0, 0, 0, 0, 0, 0, 0,
#            0, 0, 0], 
#        [8, -8, 0, 2, -2, 0, -1, 1, 2, -2, 0, -4, 4, -1, 1, 0, 0, 0, 0, 0, 0, 0,
#            0, 0, 0], 
#        [8, -8, 0, -1, 1, 0, 2, -2, 2, -2, 0, 0, 0, 0, 0, 0, 4, -4, 1, -1, 0, 0,
#            0, 0, 0], 
#        [8, -8, 0, -1, 1, 0, 2, -2, 2, -2, 0, 0, 0, 0, 0, 0, -4, 4, -1, 1, 0, 0,
#            0, 0, 0], 
#        [16, -16, 0, -2, 2, 0, -2, 2, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#            0, 0, 0]
#    ]
#
#
#def coxeterclasses(inds):
#    classes = (
#        ([], "A0", 1), ([0], "A1", 2), ([2], "~A1", 2), ([0, 2], "A1+~A1", 3), 
#        ([0, 1], "A2", 1), ([2, 3], "~A2", 1), ([1, 2], "B2", 1), 
#        ([0, 1, 3], "A2+~A1", 1), ([0, 2, 3], "~A2+A1", 1), 
#        ([0, 1, 2], "B3", 1), ([1, 2, 3], "C3", 1), ([0, 1, 2, 3], "F4", 1)
#    )
#
#    return (([inds[i] for i in x[0]], x[1], x[2]) for x in classes)
#
