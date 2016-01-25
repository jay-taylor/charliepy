import charliepy
from charliepy.data import typ2A
import numpy as np
import unittest

###############################################################################
#                             CHARACTER TABLES
#
# Tables of irreducible character values taken from CHEVIE. From ranks 2 to 6.
ct2 = np.array([[-1, 1, -1],
                [-2, 0,  1],
                [ 1, 1,  1]], dtype='int')

ct3 = np.array([[ 1, -1,  1,  1, -1],
                [-3,  1,  1,  0, -1],
                [ 2,  0,  2, -1,  0],
                [-3, -1,  1,  0,  1],
                [ 1,  1,  1,  1,  1]], dtype='int')

ct4 = np.array([[ 1, -1, 1,  1, -1, -1,  1],
                [ 4, -2, 0,  1,  1,  0, -1],
                [ 5, -1, 1, -1, -1,  1,  0],
                [-6,  0, 2,  0,  0,  0, -1],
                [ 5,  1, 1, -1,  1, -1,  0],
                [-4, -2, 0, -1,  1,  0,  1],
                [ 1,  1, 1,  1,  1,  1,  1]], dtype='int')

ct5 = np.array([[ -1,  1, -1, 1, -1,  1, -1,  1, -1, -1,  1],
                [  5, -3,  1, 1,  2,  0, -1, -1, -1,  0,  1],
                [ -9,  3, -1, 3,  0,  0,  0, -1, -1,  1,  0],
                [  5, -1,  1, 3, -1, -1,  2,  1, -1,  0,  0],
                [ 10, -2, -2, 2,  1,  1,  1,  0,  0,  0, -1],
                [ 16,  0,  0, 0, -2,  0, -2,  0,  0,  1,  0],
                [ -5, -1, -1, 3,  1, -1, -2,  1,  1,  0,  0],
                [-10, -2,  2, 2, -1,  1, -1,  0,  0,  0, -1],
                [  9,  3,  1, 3,  0,  0,  0, -1,  1, -1,  0],
                [ -5, -3, -1, 1, -2,  0,  1, -1,  1,  0,  1],
                [  1,  1,  1, 1,  1,  1,  1,  1,  1,  1,  1]], dtype='int')

ct6 = np.array([[ -1,  1, -1, 1, -1,  1, -1, -1,  1, -1,  1, -1,  1,  1, -1],
                [ -6,  4, -2, 0, -3,  1,  1,  0,  2,  0, -1, -1, -1,  0,  1],
                [-14,  6, -2, 2, -2,  0, -2,  1,  0,  0,  0,  1,  1, -1,  0],
                [-14,  4, -2, 0,  1,  1,  1, -2, -2,  0,  1,  1, -1,  0,  0],
                [ 15, -5, -1, 3,  3,  1, -1,  0, -1, -1, -1,  0,  0,  0,  1],
                [-35,  5,  1, 1,  1, -1,  1,  1, -1, -1, -1,  0,  0,  1,  0],
                [ 21, -1,  1, 3, -3, -1,  1,  0,  1, -1,  1,  1, -1,  0,  0],
                [-21, -1, -1, 3,  3, -1, -1,  0,  1,  1,  1, -1, -1,  0,  0],
                [ 20,  0, -4, 0,  2,  0,  2,  2,  0,  0,  0,  0,  0,  0, -1],
                [ 35,  5, -1, 1, -1, -1, -1, -1, -1,  1, -1,  0,  0,  1,  0],
                [-14, -4, -2, 0,  1, -1,  1, -2,  2,  0, -1,  1,  1,  0,  0],
                [-15, -5,  1, 3, -3,  1,  1,  0, -1,  1, -1,  0,  0,  0, -1],
                [ 14,  6,  2, 2,  2,  0,  2, -1,  0,  0,  0, -1,  1, -1,  0],
                [ -6, -4, -2, 0, -3, -1,  1,  0, -2,  0,  1, -1,  1,  0,  1],
                [  1,  1,  1, 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1]],
                dtype='int')

###############################################################################
#                             CONJUGACY CLASSES
#
# Tables of conjugacy classes taken from CHEVIE. From ranks 2 to 7.
classes2 = [[[1, 2, 1], [], [2]], ["111", "21", "3"]]

classes3 = [[[1, 2, 1, 3, 2, 1], [2], [], [2, 3], [3]],
             ["1111", "211", "22", "31", "4"]]

classes4 = [[[1, 2, 1, 3, 2, 1, 4, 3, 2, 1], [2, 3, 2], [], [2, 3, 2, 4],
                [3], [4], [3, 4]],
            ["11111", "2111", "221", "311", "32", "41", "5"]]

classes5 = [[[1, 2, 1, 3, 2, 1, 4, 3, 2, 1, 5, 4, 3, 2, 1],
                [2, 3, 2, 4, 3, 2], [3], [], [2, 3, 2, 4, 3, 2, 5], [3, 4],
                [2, 3, 2, 4, 5], [3, 5], [5], [3, 4, 5], [4, 5]],
            ["111111", "21111", "2211", "222", "3111", "321", "33", "411",
                "42", "51", "6"]]

classes6 = [[[1, 2, 1, 3, 2, 1, 4, 3, 2, 1, 5, 4, 3, 2, 1, 6, 5, 4, 3, 2, 1],
                [2, 3, 2, 4, 3, 2, 5, 4, 3, 2], [3, 4, 3], [],
                [2, 3, 2, 4, 3, 2, 5, 4, 3, 2, 6], [3, 4, 3, 5], [4],
                [2, 3, 4, 3, 2, 5, 6], [3, 4, 3, 6], [6], [4, 6],
                [3, 4, 3, 5, 6], [4, 5], [5, 6], [4, 5, 6]],
            ["1111111", "211111", "22111", "2221", "31111", "3211", "322",
                "331", "4111", "421", "43", "511", "52", "61", "7"]]

classes7 = [[[1, 2, 1, 3, 2, 1, 4, 3, 2, 1, 5, 4, 3, 2, 1, 6, 5, 4, 3, 2, 1,
                7, 6, 5, 4, 3, 2, 1],
                [2, 3, 2, 4, 3, 2, 5, 4, 3, 2, 6, 5, 4, 3, 2],
                [3, 4, 3, 5, 4, 3], [4], [],
                [2, 3, 2, 4, 3, 2, 5, 4, 3, 2, 6, 5, 4, 3, 2, 7],
                [3, 4, 3, 5, 4, 3, 6], [4, 5], [2, 3, 4, 3, 5, 4, 3, 2, 6, 7],
                [3, 4, 3, 5, 6], [3, 4, 3, 5, 4, 3, 7], [4, 7], [7], [4, 5, 7],
                [5, 7], [3, 4, 3, 5, 4, 3, 6, 7], [4, 5, 6],
                [3, 4, 3, 5, 6, 7], [4, 6, 7], [6, 7], [4, 5, 6, 7], [5, 6, 7]],
            ["11111111", "2111111", "221111", "22211", "2222", "311111",
                "32111", "3221", "3311", "332", "41111", "4211", "422", "431",
                "44", "5111", "521", "53", "611", "62", "71", "8"]] 

class TestAData(unittest.TestCase):

    def test_chartabs(self):
        # Make sure the character tables are being produced correctly.
        res2 = typ2A.chartable(2)
        res3 = typ2A.chartable(3)
        res4 = typ2A.chartable(4)
        res5 = typ2A.chartable(5)
        res6 = typ2A.chartable(6)

        self.assertTrue(np.array_equal(res2, ct2),
                'Rank 2 Character Table not correct!')
        self.assertTrue(np.array_equal(res3, ct3),
                'Rank 3 Character Table not correct!')
        self.assertTrue(np.array_equal(res4, ct4),
                'Rank 4 Character Table not correct!')
        self.assertTrue(np.array_equal(res5, ct5),
                'Rank 5 Character Table not correct!')
        self.assertTrue(np.array_equal(res6, ct6),
                'Rank 6 Character Table not correct!')

    def test_conjclasses(self):
        # Make sure the conjugacy class information is being produced
        # correctly.
        res2 = typ2A.conjclassdata(range(1,3))[0:3:2]
        res3 = typ2A.conjclassdata(range(1,4))[0:3:2]
        res4 = typ2A.conjclassdata(range(1,5))[0:3:2]
        res5 = typ2A.conjclassdata(range(1,6))[0:3:2]
        res6 = typ2A.conjclassdata(range(1,7))[0:3:2]
        res7 = typ2A.conjclassdata(range(1,8))[0:3:2]

        self.assertEqual(res2, classes2,
                'Rank 2 Conjugacy Classes not correct!')
        self.assertEqual(res3, classes3,
                'Rank 3 Conjugacy Classes not correct!')
        self.assertEqual(res4, classes4,
                'Rank 4 Conjugacy Classes not correct!')
        self.assertEqual(res5, classes5,
                'Rank 5 Conjugacy Classes not correct!')
        self.assertEqual(res6, classes6,
                'Rank 6 Conjugacy Classes not correct!')
        self.assertEqual(res7, classes7,
                'Rank 7 Conjugacy Classes not correct!')


def suite():
    tests = ['test_chartabs']

    return unittest.TestLoader().loadTestsFromTestCase(TestAData)

# Run the test suite.
unittest.TextTestRunner(verbosity=2).run(suite())

