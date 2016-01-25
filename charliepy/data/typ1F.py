########################################################
# Data for type F_4 Coxeter groups and reductive groups
# The Dynkin diagram here follows Bourbaki so we have
#
#                  0   1   2   3
#                  o---o=>=o---o

import numpy as np

def cartanmat(n):
    if n == 4:
        return np.array([[ 2, -1,  0,  0],
                         [-1,  2, -2,  0],
                         [ 0, -1,  2, -1],
                         [ 0,  0, -1,  2]],
                         dtype = 'int8')
    else:
        raise ValueError("Rank of F_4 must be 4!")

def diagram(inds):
    """Prints the Dynkin diagram."""
    print("F4 {} - {} > {} - {}".format(*inds))
    return None

def degrees(n):
    return [2, 6, 8, 12]

def conjclassdata(ind, **kwargs):
    # stores the data: representatives, centraliser orders, names
    repcentnam = list() 

    # Stored as an iterable because we're about to loop through it.
    reps = iter([[], [0, 1, 0, 2, 1, 0, 2, 1, 2, 3, 2, 1, 0, 2, 1, 2, 3, 2,
        1, 0, 2, 1, 2, 3],
    [1, 2, 1, 2], [1, 0], [0, 1, 2, 1, 2, 3, 2, 1, 2, 3],
    [0, 1, 0, 2, 1, 0, 2, 3, 2, 1, 2, 3], [3, 2],
    [0, 1, 0, 2, 1, 0, 2, 1, 2, 3],
    [0, 1, 0, 2, 1, 0, 2, 3, 2, 1, 0, 2, 1, 2, 3, 2],
    [0, 1, 0, 2, 1, 2, 3, 2], [0, 1, 2, 3], [0], [1, 2, 1, 2, 3, 2, 1, 2, 3],
    [0, 3, 2], [3, 2, 1], [1, 2, 1, 0, 2], [2], [0, 1, 0, 2, 1, 0, 2, 1, 2],
    [1, 0, 3], [2, 1, 0], [1, 3, 2, 1, 2], [0, 2], [2, 1],
    [0, 1, 0, 2, 1, 0, 2, 1, 2, 3, 2, 1, 2, 3], [0, 1, 2, 1, 3, 2]])

    # Create words based on indices of the irreducible component.
    repcentnam.append([[ind[j] for j in rep] for rep in reps])

    repcentnam.append([1152, 1152, 64, 36, 36, 96, 36, 36, 72, 72, 12, 96,
    96, 12, 12, 16, 96, 96, 12, 12, 16, 16, 32, 32, 8])

    repcentnam.append([' ', '4A_1', '2A_1', 'A_2', 'D_4', 'D_4(a_1)', '~A_2',
    'C_3+A_1', 'A_2+~A_2', 'F_4(a_1)', 'F_4', 'A_1', '3A_1', '~A_2+A_1',
    'C_3', 'A_3', '~A_1', '2A_1+~A_1', 'A_2+~A_1', 'B_3', 'B_2+A_1',
    'A_1+~A_1', 'B_2', 'A_3+~A_1', 'B_4'])

    return repcentnam

def irrchardata(n, **kwargs):
    ainv = [0, 4, 4, 24, 1, 13, 1, 13, 4, 2, 4, 4, 10, 4, 4, 4, 1, 4, 4, 13, 3,
            9, 3, 9, 4]
    binv = [0, 12, 12, 24, 4, 16, 4, 16, 8, 2, 6, 6, 10, 6, 6, 4, 1, 7, 7, 13,
            3, 9, 3, 9, 5]
    nam = ["1_1", "1_2", "1_3", "1_4", "2_1", "2_2", "2_3", "2_4", "4_1", "9_1",
           "9_2", "9_3", "9_4", "6_1", "6_2", "12", "4_2", "4_3", "4_4", "4_5",
           "8_1",  "8_2", "8_3", "8_4", "16"]
    
    return [nam, ainv, binv]

def chartable(n, **kwargs):
    return np.array([
        [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
        [1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1],
        [1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1],
        [2,2,2,2,2,2,-1,-1,-1,-1,-1,2,2,-1,-1,2,0,0,0,0,0,0,0,0,0],
        [2,2,2,2,2,2,-1,-1,-1,-1,-1,-2,-2,1,1,-2,0,0,0,0,0,0,0,0,0],
        [2,2,2,-1,-1,2,2,2,-1,-1,-1,0,0,0,0,0,2,2,-1,-1,2,0,0,0,0],
        [2,2,2,-1,-1,2,2,2,-1,-1,-1,0,0,0,0,0,-2,-2,1,1,-2,0,0,0,0],
        [4,4,4,-2,-2,4,-2,-2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [9,9,1,0,0,-3,0,0,0,0,0,3,3,0,0,-1,3,3,0,0,-1,1,1,1,-1],
        [9,9,1,0,0,-3,0,0,0,0,0,3,3,0,0,-1,-3,-3,0,0,1,-1,-1,-1,1],
        [9,9,1,0,0,-3,0,0,0,0,0,-3,-3,0,0,1,3,3,0,0,-1,-1,-1,-1,1],
        [9,9,1,0,0,-3,0,0,0,0,0,-3,-3,0,0,1,-3,-3,0,0,1,1,1,1,-1],
        [6,6,-2,0,0,2,0,0,3,3,-1,0,0,0,0,0,0,0,0,0,0,2,-2,-2,0],
        [6,6,-2,0,0,2,0,0,3,3,-1,0,0,0,0,0,0,0,0,0,0,-2,2,2,0],
        [12,12,-4,0,0,4,0,0,-3,-3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [4,-4,0,1,-1,0,1,-1,-2,2,0,2,-2,-1,1,0,2,-2,-1,1,0,0,2,-2,0],
        [4,-4,0,1,-1,0,1,-1,-2,2,0,2,-2,-1,1,0,-2,2,1,-1,0,0,-2,2,0],
        [4,-4,0,1,-1,0,1,-1,-2,2,0,-2,2,1,-1,0,2,-2,-1,1,0,0,-2,2,0],
        [4,-4,0,1,-1,0,1,-1,-2,2,0,-2,2,1,-1,0,-2,2,1,-1,0,0,2,-2,0],
        [8,-8,0,2,-2,0,-1,1,2,-2,0,4,-4,1,-1,0,0,0,0,0,0,0,0,0,0],
        [8,-8,0,2,-2,0,-1,1,2,-2,0,-4,4,-1,1,0,0,0,0,0,0,0,0,0,0],
        [8,-8,0,-1,1,0,2,-2,2,-2,0,0,0,0,0,0,4,-4,1,-1,0,0,0,0,0],
        [8,-8,0,-1,1,0,2,-2,2,-2,0,0,0,0,0,0,-4,4,-1,1,0,0,0,0,0],
        [16,-16,0,-2,2,0,-2,2,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]],
        dtype='int')

