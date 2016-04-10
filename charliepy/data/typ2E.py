########################################################
# Data for type 2E_6 Coxeter groups and reductive groups
#

from . import typ1E

import numpy as np

def conjclasses(inds, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type E_6 and phi is the unique graph automorphism
    of order 2. The data is adapted from the corresponding files in GAP-Chevie.

    """
    classes = (
        ([0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3,
            1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0], 51840, 'A0'),
        ([], 1152, '4A1'),
        ([2, 3, 2, 4, 3, 2], 192, '2A1'),
        ([0, 1, 3, 2, 0, 4, 3, 2, 5, 4, 3, 2], 648, '3A2'),
        ([0, 1, 2, 0, 3, 2, 0, 4, 3, 2, 0, 5, 4, 3, 2, 0], 216, 'A2'),
        ([1, 2, 3, 1, 2, 3, 5, 4, 3, 1, 2, 3, 4, 5], 108, '2A2'),
        ([0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 5, 4, 3, 2, 0], 96, 'D4(a1)'),
        ([0, 1], 16, 'A3+A1'),
        ([3, 4, 3, 1, 2, 0, 3, 4], 10, 'A4'),
        ([3, 1, 4, 3, 1, 2, 3, 4, 5, 4, 3, 1, 2, 3, 4, 5], 72, 'E6(a2)'),
        ([1, 3], 36, 'D4'),
        ([0, 4], 36, 'A5+A1'),
        ([4, 3], 24, 'A2+2A1'),
        ([0, 1, 4, 3], 9, 'E6(a1)'),
        ([0, 1, 2, 0, 3, 2], 12, 'E6'),
        ([0, 2, 0, 3, 2, 0, 4, 3, 2, 0, 5, 4, 3, 2, 0], 1440, 'A1'),
        ([1], 96, '3A1'),
        ([0], 96, 'A3+2A1'),
        ([1, 2, 3, 2, 4, 3, 2], 32, 'A3'),
        ([0, 2, 3, 2, 4, 3, 2], 36, 'A2+A1'),
        ([0, 2, 0, 3, 2], 36, '2A2+A1'),
        ([0, 1, 4], 12, 'A5'),
        ([1, 4, 3], 8, 'D5'),
        ([0, 4, 3], 10, 'A4+A1'),
        ([0, 1, 3], 12, 'D5(a1)')
    )

    return (([inds[j] for j in cls[0]], cls[1], cls[2]) for cls in classes)

def conjclasses_min(inds, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type E_6 and phi is the unique graph automorphism
    of order 2. The data is adapted from the corresponding files in GAP-Chevie.

    """
    classes = (
        (51840, 'A0'), (1152, '4A1'), (192, '2A1'), (648, '3A2'), (216, 'A2'),
        (108, '2A2'), (96, 'D4(a1)'), (16, 'A3+A1'), (10, 'A4'),
        (72, 'E6(a2)'), (36, 'D4'), (36, 'A5+A1'), (24, 'A2+2A1'),
        (9, 'E6(a1)'), (12, 'E6'), (1440, 'A1'), (96, '3A1'), (96, 'A3+2A1'),
        (32, 'A3'), (36, 'A2+A1'), (36, '2A2+A1'), (12, 'A5'), (8, 'D5'),
        (10, 'A4+A1'), (12, 'D5(a1)')
    )

    return (cls for cls in classes)

irrchars = typ1E.irrchars

def chartable(n, **kwargs):
    # Here we use Lusztig's preferred extension defined in [L85, IV,
    # 17.2]. In particular w.phi acts on the representation E as
    # w(-1)^a_E. This involves multiplying all rows of the character
    # table by appropriate signs.
    chartabsplit = typ1E.chartable(6)
    ainv = typ1E.irrchardata(6)[1]
    signs = [(-1)**a for a in ainv]
    return [[s*x for x in row] for s, row in zip(signs, chartabsplit)]
