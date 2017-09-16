########################################################
# Data for type H_3 and H_4 Coxeter groups
#

import numpy as np

#def cartanmat(n):
#    if n == 3:
#        return np.array([[           2, -utils.ir(5),  0],
#                         [-utils.ir(5),            2, -1],
#                         [           0,           -1,  2]])
#
#    if n == 4:
#        return np.array([[           2, -utils.ir(5),  0,  0],
#                         [-utils.ir(5),            2, -1,  0],
#                         [           0,           -1,  2, -1],
#                         [           0,            0, -1,  2]])

def degrees(n):
    if n == 3:
        return [2, 6, 10]
    if n == 4:
        return [2, 12, 20, 30]

def conjclassdata(ind, **kwargs):
    # stores the data: representatives, centraliser orders, names
    n = len(ind)

    # Stored as an iterable because we're about to loop through it.
    if n == 3:
        reps = [
            [], [0], [0, 1], [0, 2], [1, 2], [0, 1, 2], [0, 1, 0, 1],
            [0, 1, 0, 1, 2], [0, 1, 0, 1, 2, 1, 0, 1, 2],
            [0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2]
        ]

        cents = [120, 8, 10, 8, 6, 10, 10, 6, 10, 120]

    elif n == 4:
        reps = [
            [], [0], [0, 1], [0, 2], [1, 2], [0, 1, 2], [0, 1, 3],
            [0, 2, 3], [1, 3, 2], [0, 1, 0, 1], [0, 1, 2, 3], [0, 1, 0, 1, 2],
            [0, 1, 0, 1, 3], [0, 1, 0, 1, 2, 3], [0, 1, 2, 1, 0, 1, 2, 3],
            [0, 1, 0, 1, 2, 1, 0, 1, 2], [0, 1, 0, 1, 2, 1, 0, 1, 2, 3],
            [0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3],
            [0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3],
            [0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2],
            [0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3], 
            [0, 1, 0, 1, 2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 2, 3], 
            [0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 2, 3, 2, 1, 0, 1, 2, 3], 
            [0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 1, 2, 3], 
            [0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 2, 3], 
            [0, 1, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 2, 1, 0, 1,
                2, 3], 
            [1, 0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 2, 1,
                0, 1, 2, 3], 
            [0, 1, 0, 1, 2, 1, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 1, 0, 2, 1, 0, 3,
                2, 1, 0, 1, 2, 3], 
            [0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 3, 2, 1, 0, 1, 0, 2, 1,
                0, 3, 2, 1, 0, 1, 2, 3], 
            [0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 2, 1, 0,
                1, 2, 3, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3], 
            [0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 2,
                1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3], 
            [0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 3, 2, 1, 0, 1, 0, 2, 1,
                0, 1, 3, 2, 1, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 1, 2, 3], 
            [0, 1, 0, 1, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 2, 1,
                0, 1, 2, 3, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 2,
                1, 0, 1, 2, 3], 
            [0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 2,
                1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0,
                2, 1, 0, 1, 2, 3, 2, 1, 0, 1, 0, 2, 1, 0, 1, 2, 3]
        ]

        cents = [
            14400, 240, 100, 32, 36, 20, 20, 12, 8, 100, 30, 12, 20, 20, 30, 20,
            12, 600, 50, 240, 100, 30, 20, 360, 36, 600, 50, 30, 240, 600, 100,
            360, 600, 14400
        ]


    for i, elm in enumerate(reps):
        yield ([ind[j] for j in elm], cents[i],
               ''.join(str(ind[j]) for j in elm))

def conjclassdata_min(ind, **kwargs):
    for elm in conjclassdata(ind, **kwargs):
        yield (elm[1], elm[2])

def irrchardata(n, **kwargs):
    if n == 3:
        return (
            ("1_r'", 15, 15), ('1_r', 0, 0), ("5_r'", 5, 5), ('5_r', 2, 2),
            ('3_s', 6, 6), ('overline{3}_s', 6, 8), ("3_s'", 1, 1),
            ("overline{3}_s'", 1, 3), ("4_r'", 3, 3), ('4_r', 3, 4)
        )

    elif n == 4:
        return (
            ('1_r', 0, 0), ("1_r'", 60, 60), ('4_t', 1, 1), ("4_t'", 31, 31),
            ('overline{4}_t', 1, 7), ("overline{4}_t'", 31, 37),
            ('6_s', 6, 12), ('overline{6}_s', 6, 20), ('8_r', 6, 12),
            ('8_{rr}', 6, 13), ('9_s', 2, 2), ("9_s'", 22, 22),
            ('overline{9}_s', 2, 6), ("overline{9}_s'", 22, 26),
            ('10_r', 6, 12), ('16_t', 6, 11), ('overline{16}_t', 6, 13),
            ('16_{rr}', 3, 3), ("16_{rr}'", 18, 21), ('16_r', 3, 6),
            ("16_r'", 18, 18), ('18_r', 6, 10), ('24_t', 6, 11),
            ('overline{24}_t', 6, 7), ('24_s', 6, 12),
            ('overline{24}_s', 6, 6), ('25_r', 4, 4), ("25_r'", 16, 16),
            ('30_s', 6, 10), ('overline{30}_s', 6, 10), ('36_{rr}', 5, 5),
            ("36_{rr}'", 15, 15), ('40_r', 6, 8), ('48_{rr}', 6, 9)
        )

def chartable(n, **kwargs):
    pass
