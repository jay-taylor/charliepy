import charliepy
from charliepy.data import typ1A
import numpy as np
import math
import unittest

data = [{'ctab' : None, 'class' : None, 'chars' : None} for i in range(8)]

###############################################################################
#                             CHARACTER TABLES
#
# Tables of irreducible character values taken from CHEVIE. From ranks 1 to 6.
data[1]['ctab'] = np.array([
    [1, -1],
    [1,  1]], dtype='int')

data[2]['ctab'] = np.array([
    [1, -1,  1],
    [2,  0, -1],
    [1,  1,  1]], dtype='int')

data[3]['ctab'] = np.array([
    [1, -1,  1,  1, -1],
    [3, -1, -1,  0,  1],
    [2,  0,  2, -1,  0],
    [3,  1, -1,  0, -1],
    [1,  1,  1,  1,  1]], dtype='int')

data[4]['ctab'] = np.array([
    [1, -1,  1,  1, -1, -1,  1],
    [4, -2,  0,  1,  1,  0, -1],
    [5, -1,  1, -1, -1,  1,  0],
    [6,  0, -2,  0,  0,  0,  1],
    [5,  1,  1, -1,  1, -1,  0],
    [4,  2,  0,  1, -1,  0, -1],
    [1,  1,  1,  1,  1,  1,  1]], dtype='int')

data[5]['ctab'] = np.array([
    [1,  -1,  1, -1,  1, -1,  1, -1,  1,  1, -1],
    [5,  -3,  1,  1,  2,  0, -1, -1, -1,  0,  1],
    [9,  -3,  1, -3,  0,  0,  0,  1,  1, -1,  0],
    [5,  -1,  1,  3, -1, -1,  2,  1, -1,  0,  0],
    [10, -2, -2,  2,  1,  1,  1,  0,  0,  0, -1],
    [16,  0,  0,  0, -2,  0, -2,  0,  0,  1,  0],
    [5,   1,  1, -3, -1,  1,  2, -1, -1,  0,  0],
    [10,  2, -2, -2,  1, -1,  1,  0,  0,  0,  1],
    [9,   3,  1,  3,  0,  0,  0, -1,  1, -1,  0],
    [5,   3,  1, -1,  2,  0, -1,  1, -1,  0, -1],
    [1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1]], dtype='int')

data[6]['ctab'] = np.array([
    [1,  -1,  1, -1,  1, -1,  1,  1, -1,  1, -1,  1, -1, -1,  1],
    [6,  -4,  2,  0,  3, -1, -1,  0, -2,  0,  1,  1,  1,  0, -1],
    [14, -6,  2, -2,  2,  0,  2, -1,  0,  0,  0, -1, -1,  1,  0],
    [14, -4,  2,  0, -1, -1, -1,  2,  2,  0, -1, -1,  1,  0,  0],
    [15, -5, -1,  3,  3,  1, -1,  0, -1, -1, -1,  0,  0,  0,  1],
    [35, -5, -1, -1, -1,  1, -1, -1,  1,  1,  1,  0,  0, -1,  0],
    [21, -1,  1,  3, -3, -1,  1,  0,  1, -1,  1,  1, -1,  0,  0],
    [21,  1,  1, -3, -3,  1,  1,  0, -1, -1, -1,  1,  1,  0,  0],
    [20,  0, -4,  0,  2,  0,  2,  2,  0,  0,  0,  0,  0,  0, -1],
    [35,  5, -1,  1, -1, -1, -1, -1, -1,  1, -1,  0,  0,  1,  0],
    [14,  4,  2,  0, -1,  1, -1,  2, -2,  0,  1, -1, -1,  0,  0],
    [15,  5, -1, -3,  3, -1, -1,  0,  1, -1,  1,  0,  0,  0,  1],
    [14,  6,  2,  2,  2,  0,  2, -1,  0,  0,  0, -1,  1, -1,  0],
    [6,   4,  2,  0,  3,  1, -1,  0,  2,  0, -1,  1, -1,  0, -1],
    [1,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1]], dtype='int')

###############################################################################
#                             CONJUGACY CLASSES
#
# Tables of conjugacy classes taken from CHEVIE. From ranks 1 to 7.
data[1]['class'] = [
    [[1, 1], [2]],
    [[], [1]],
    [2, 2],
    ["11", "2"]]

data[2]['class'] = [
    [[1, 1, 1], [2, 1], [3]],
    [[], [1], [1, 2]],
    [6, 2, 3],
    ["111", "21", "3"]]

data[3]['class'] = [
    [[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]],
    [[], [1], [1, 3], [1, 2], [1, 3, 2]],
    [24, 4, 8, 3, 4],
    ["1111", "211", "22", "31", "4"]]

data[4]['class'] = [
    [[1, 1, 1, 1, 1], [2, 1, 1, 1], [2, 2, 1], [3, 1, 1], [3, 2], [4, 1], [5]],
    [[], [1], [1, 3], [1, 2], [1, 2, 4], [1, 3, 2], [1, 3, 2, 4]],
    [120, 12, 8, 6, 6, 4, 5],
    ["11111", "2111", "221", "311", "32", "41", "5"]]

data[5]['class'] = [
    [[1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1], [2, 2, 1, 1], [2, 2, 2],
        [3, 1, 1, 1 ], [3, 2, 1], [3, 3], [4, 1, 1], [4, 2], [5, 1], [6]],
    [[], [1], [1, 3], [1, 3, 5], [1, 2], [1, 2, 4], [1, 2, 4, 5],
         [1, 3, 2], [1, 3, 2, 5], [1, 3, 2, 4], [1, 3, 5, 2, 4]],
    [720, 48, 16, 48, 18, 6, 18, 8, 8, 5, 6],
    ["111111", "21111", "2211", "222", "3111", "321", "33", "411", "42", "51",
         "6"]]

data[6]['class'] = [
    [[1, 1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1, 1], [2, 2, 1, 1, 1], [2, 2, 2, 1],
        [3, 1, 1, 1, 1], [3, 2, 1, 1], [3, 2, 2], [3, 3, 1], [4, 1, 1, 1],
        [4, 2, 1], [4, 3], [5, 1, 1], [5, 2], [6, 1], [7]],
    [[], [1], [1, 3], [1, 3, 5], [1, 2], [1, 2, 4], [1, 2, 4, 6],
         [1, 2, 4, 5], [1, 3, 2], [1, 3, 2, 5], [1, 3, 2, 5, 6], [1, 3, 2, 4],
         [1, 3, 2, 4, 6], [1, 3, 5, 2, 4], [1, 3, 5, 2, 4, 6]],
    [5040, 240, 48, 48, 72, 12, 24, 18, 24, 8, 12, 10, 10, 6, 7],
    ["1111111", "211111", "22111", "2221", "31111", "3211", "322", "331",
         "4111", "421", "43", "511", "52", "61", "7" ]]

data[7]['class'] = [
    [[1, 1, 1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1, 1, 1], [2, 2, 1, 1, 1, 1],
        [2, 2, 2, 1, 1], [2, 2, 2, 2], [3, 1, 1, 1, 1, 1], [3, 2, 1, 1, 1],
        [3, 2, 2, 1], [3, 3, 1, 1], [3, 3, 2], [4, 1, 1, 1, 1], [4, 2, 1, 1],
        [4, 2, 2], [4, 3, 1], [4, 4], [5, 1, 1, 1], [5, 2, 1], [5, 3],
        [6, 1, 1], [6, 2], [7, 1], [8]],
    [[], [1], [1, 3], [1, 3, 5], [1, 3, 5, 7], [1, 2], [1, 2, 4],
        [1, 2, 4, 6], [1, 2, 4, 5], [1, 2, 4, 5, 7], [1, 3, 2],
        [1, 3, 2, 5], [1, 3, 2, 5, 7], [1, 3, 2, 5, 6],
        [1, 3, 2, 5, 7, 6], [1, 3, 2, 4], [1, 3, 2, 4, 6],
        [1, 3, 2, 4, 6, 7], [1, 3, 5, 2, 4], [1, 3, 5, 2, 4, 7],
        [1, 3, 5, 2, 4, 6], [1, 3, 5, 7, 2, 4, 6]],
    [40320, 1440, 192, 96, 384, 360, 36, 24, 36, 36, 96, 16, 32, 12, 32, 30,
        10, 15, 12, 12, 7, 8],
    ["11111111", "2111111", "221111", "22211", "2222", "311111", "32111",
        "3221", "3311", "332", "41111", "4211", "422", "431", "44", "5111",
        "521", "53", "611", "62", "71", "8"]]

###############################################################################
#                                CHARACTERS 
#
# Tables of irreducible characters taken from CHEVIE. From ranks 1 to 7.
data[1]['chars'] = [
    [[1, 1], [2]],
    ["11", "2"],
    [1, 0],
    [1, 0]]

data[2]['chars']= [
    [[1, 1, 1], [2, 1], [3]],
    ["111", "21", "3"],
    [3, 1, 0],
    [3, 1, 0]]

data[3]['chars']= [
    [[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]],
    ["1111", "211", "22", "31", "4"],
    [6, 3, 2, 1, 0],
    [6, 3, 2, 1, 0]]

data[4]['chars']= [
    [[1, 1, 1, 1, 1], [2, 1, 1, 1], [2, 2, 1], [3, 1, 1], [3, 2], [4, 1], [5]],
    ["11111", "2111", "221", "311", "32", "41", "5"],
    [10, 6, 4, 3, 2, 1, 0],
    [10, 6, 4, 3, 2, 1, 0]]

data[5]['chars']= [
    [[1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1], [2, 2, 1, 1], [2, 2, 2], [3, 1, 1, 1],
        [3, 2, 1], [3, 3], [4, 1, 1], [4, 2], [5, 1], [6]],
    ["111111", "21111", "2211", "222", "3111", "321", "33", "411", "42", "51",
        "6"],
    [15, 10, 7, 6, 6, 4, 3, 3, 2, 1, 0],
    [15, 10, 7, 6, 6, 4, 3, 3, 2, 1, 0]]

data[6]['chars']= [
    [[1, 1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1, 1], [2, 2, 1, 1, 1], [2, 2, 2, 1],
        [3, 1, 1, 1, 1], [3, 2, 1, 1], [3, 2, 2], [3, 3, 1], [4, 1, 1, 1],
        [4, 2, 1], [4, 3], [5, 1, 1], [5, 2], [6, 1], [7]],
    ["1111111", "211111", "22111", "2221", "31111", "3211", "322", "331",
        "4111", "421", "43", "511", "52", "61", "7"],
    [21, 15, 11, 9, 10, 7, 6, 5, 6, 4, 3, 3, 2, 1, 0],
    [21, 15, 11, 9, 10, 7, 6, 5, 6, 4, 3, 3, 2, 1, 0]]

data[7]['chars']= [
    [[1, 1, 1, 1, 1, 1, 1, 1], [2, 1, 1, 1, 1, 1, 1], [2, 2, 1, 1, 1, 1],
        [2, 2, 2, 1, 1], [2, 2, 2, 2], [3, 1, 1, 1, 1, 1], [3, 2, 1, 1, 1],
        [3, 2, 2, 1], [3, 3, 1, 1], [3, 3, 2], [4, 1, 1, 1, 1], [4, 2, 1, 1],
        [4, 2, 2], [4, 3, 1], [4, 4], [5, 1, 1, 1], [5, 2, 1], [5, 3],
        [6, 1, 1], [6, 2], [7, 1], [8]],
    ["11111111", "2111111", "221111", "22211", "2222", "311111", "32111",
        "3221", "3311", "332", "41111", "4211", "422", "431", "44", "5111",
        "521", "53", "611", "62", "71", "8"],
    [28, 21, 16, 13, 12, 15, 11, 9, 8, 7, 10, 7, 6, 5, 4, 6, 4, 3, 3, 2, 1, 0],
    [28, 21, 16, 13, 12, 15, 11, 9, 8, 7, 10, 7, 6, 5, 4, 6, 4, 3, 3, 2, 1, 0]]

class TestAData(unittest.TestCase):

    def test_chartabs(self):
        # Make sure the character tables agree with those from CHEVIE.
        for i in range(1, 7):
            ctab = typ1A.chartable(i)

            for j, mu in enumerate(typ1A._charlabels(i+1)):
                for k, nu in enumerate(typ1A._conjlabels(i+1)):
                    x = data[i]['chars'][0].index(mu)
                    y = data[i]['class'][0].index(nu)
                    self.assertTrue(ctab[j][k] == data[i]['ctab'][x][y],
                        'Rank ' + str(i+1) + ' character table not correct!')

        # Test that the character tables are internally consistent.
        for i in range(2, 16):
            ctab = typ1A.chartable(i)
            cents = typ1A.conjclassdata(list(range(1, i+1)))[1]
            size = math.factorial(i+1)
            classlens = [size//c for c in cents]

            # Row orthogonality relations.
            self.assertTrue(
                np.array_equal(
                    size*np.eye(ctab.shape[0], dtype='int'),
                    np.inner(ctab, ctab*classlens)
                ),
                'Rank ' + str(i) + ' character table not consistent!'
            )

            # Column orthogonality relations.
            self.assertTrue(
                np.array_equal(
                    cents*np.eye(ctab.shape[0], dtype='int'),
                    np.inner(ctab.T, ctab.T)
                ),
                'Rank ' + str(i) + ' character table not consistent!'
            )

    def test_conjclasses(self):
        # Make sure the conjugacy class information is being produced
        # correctly.
        for i in range(1, 8):
            classes = typ1A.conjclassdata(list(range(1, i+1)))

            for ind, part in enumerate(typ1A._conjlabels(i+1)):
                x = data[i]['class'][0].index(part)

                self.assertEqual(data[i]['class'][1][x], classes[0][ind],
                    'Rank ' + str(i) + ' conjugacy class representatives are '
                    'not correct!')
                self.assertEqual(data[i]['class'][2][x], classes[1][ind],
                    'Rank ' + str(i) + ' centraliser orders are '
                    'not correct!')
                self.assertEqual(data[i]['class'][3][x], classes[2][ind],
                    'Rank ' + str(i) + ' conjugacy class names are '
                    'not correct!')

    def test_characters(self):
        # Make sure the conjugacy class information is being produced
        # correctly.
        for i in range(1, 8):
            chars = typ1A.irrchardata(i)

            for ind, part in enumerate(typ1A._charlabels(i+1)):
                x = data[i]['class'][0].index(part)

                self.assertEqual(data[i]['chars'][1][x], chars[0][ind],
                    'Rank ' + str(i+1) + ' character names are not correct!')
                self.assertEqual(data[i]['chars'][2][x], chars[1][ind],
                    'Rank ' + str(i+1) + ' a-values are not correct!')
                self.assertEqual(data[i]['chars'][3][x], chars[2][ind],
                    'Rank ' + str(i+1) + ' b-values are not correct!')


def suite():
    tests = ['test_chartabs']#, 'test_conjclasses', 'test_characters']

    return unittest.TestLoader().loadTestsFromTestCase(TestAData)

# Run the test suite.
unittest.TextTestRunner(verbosity=2).run(suite())
