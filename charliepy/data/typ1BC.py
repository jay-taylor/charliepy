########################################################
# Data for type BC_n root systems.
# The Dynkin diagram here is
#
#                0   1
#                o=>=o
#
# when n = 2 and when n >= 2 we have
#
#
#       0   1   2        n-1  n
#       o=>=o---o-- ... --o=>=o
#

from .. import utils as utils
from . import typ1A

import numpy as np

def cartanmat(n):
    if n < 2:
        raise ValueError("Type BC root system must have rank at least 2.")
    elif n == 2:
        C = np.array([[2, -4], [-1, 2]], dtype = 'int8')
    else:
        C = typ1A.cartanmat(n)
        C[0][1], C[n-2][n-1] = -2, -2
    return C
