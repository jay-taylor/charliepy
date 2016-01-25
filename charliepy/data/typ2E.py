########################################################
# Data for type 2E_6 Coxeter groups and reductive groups
#

from . import typ1E

import numpy as np

def conjclassdata(ind, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type E_6 and phi is the unique graph automorphism
    of order 2. The data is adapted from the corresponding files in GAP-Chevie.

    """
    # Stores the data: representatives, centraliser orders, names
    repcentnam = list() 

    repcentnam.append([[0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4,
                        3, 1, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0],
                       [],
                       [2, 3, 2, 4, 3, 2],
                       [0, 1, 3, 2, 0, 4, 3, 2, 5, 4, 3, 2],
                       [0, 1, 2, 0, 3, 2, 0, 4, 3, 2, 0, 5, 4, 3, 2, 0],
                       [1, 2, 3, 1, 2, 3, 5, 4, 3, 1, 2, 3, 4, 5],
                       [0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 5, 4, 3, 2, 0],
                       [0, 1],
                       [3, 4, 3, 1, 2, 0, 3, 4],
                       [3, 1, 4, 3, 1, 2, 3, 4, 5, 4, 3, 1, 2, 3, 4, 5],
                       [1, 3], [0, 4], [4, 3], [0, 1, 4, 3], [0, 1, 2, 0, 3, 2],
                       [0, 2, 0, 3, 2, 0, 4, 3, 2, 0, 5, 4, 3, 2, 0], [1], [0],
                       [1, 2, 3, 2, 4, 3, 2], [0, 2, 3, 2, 4, 3, 2],
                       [0, 2, 0, 3, 2], [0, 1, 4], [1, 4, 3], [0, 4, 3],
                       [0, 1, 3]])

    # Create words based on indices of the irreducible component
    repcentnam[0] = [[ind[j] for j in rep] for rep in repcentnam[0]]

    repcentnam.append([51840, 1152, 192, 648, 216, 108, 96, 16, 10, 72, 36,
                       36, 24, 9, 12, 1440, 96, 96, 32, 36, 36, 12, 8, 10, 12])

    repcentnam.append(['A_0', '4A_1', '2A_1', '3A_2', 'A_2',
                       '2A_2', 'D_4(a_1)', 'A_3+A_1', 'A_4', 'E_6(a_2)',
                       'D_4', 'A_5+A_1', 'A_2+2A_1', 'E_6(a_1)', 'E_6', 'A_1',
                       '3A_1', 'A_3+2A_1', 'A_3', 'A_2+A_1', '2A_2+A_1', 'A_5',
                       'D_5', 'A_4+A_1', 'D_5(a_1)'])

    return repcentnam

def chartable(n, **kwargs):
    # Here we use Lusztig's preferred extension defined in [L85, IV,
    # 17.2]. In particular w.phi acts on the representation E as
    # w(-1)^a_E. This involves multiplying all rows of the character
    # table by appropriate signs.
    chartabsplit = typ1E.chartable(6)
    ainv = typ1E.irrchardata(6)[1]
    sgn = np.array([[(-1)**a] for a in ainv], dtype='int')
    return sgn*chartabsplit
