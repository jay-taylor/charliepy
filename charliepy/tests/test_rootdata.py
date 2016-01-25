from . import rootdata as rdata
import numpy as np
import unittest

class TestCoreFuncs(unittest.TestCase):

    def test_cartantotype(self):
        # Make sure the Cartan matrices are being classified correctly.
        types = [("A", 0), ("B", 2), ("C", 3), ("D", 4)]
        for typ in types:
            for i in xrange(typ[1], 25):
                cmat = rdata.cartantotype(rdata.cartanmat(typ[0],i))
                self.assertTrue(cmat == [[typ[0], range(i)]],
                    'Cartan matrix ' + typ[0] + str(i) + ' not identified!')
        
        types = [("A~", 1), ("B~", 2), ("C~", 3), ("D~", 4)]
        for typ in types:
            for i in xrange(typ[1], 25):
                cmat = rdata.cartantotype(rdata.cartanmat(typ[0],i))
                self.assertTrue(cmat == [[typ[0], range(i+1)]],
                    'Cartan matrix ' + typ[0] + str(i) + ' not identified!')

        # Exceptional type test.
        cmatE6 = [rdata.cartantotype(rdata.cartanmat("E",6)),
                [["E",range(6)]]]
        cmatE7 = [rdata.cartantotype(rdata.cartanmat("E",7)),
                [["E",range(7)]]]
        cmatE8 = [rdata.cartantotype(rdata.cartanmat("E",8)),
                [["E",range(8)]]]
        cmatF4 = [rdata.cartantotype(rdata.cartanmat("F",4)),
                [["F",range(4)]]]
        cmatG2 = [rdata.cartantotype(rdata.cartanmat("G",2)),
                [["G",range(2)]]]

        cmatEt6 = [rdata.cartantotype(rdata.cartanmat("E~",6)),
                [["E~",range(7)]]]
        cmatEt7 = [rdata.cartantotype(rdata.cartanmat("E~",7)),
                [["E~",range(8)]]]
        cmatEt8 = [rdata.cartantotype(rdata.cartanmat("E~",8)),
                [["E~",range(9)]]]
        cmatFt4 = [rdata.cartantotype(rdata.cartanmat("F~",4)),
                [["F~",range(5)]]]
        cmatGt2 = [rdata.cartantotype(rdata.cartanmat("G~",2)),
                [["G~",range(3)]]]

        types = [cmatE6, cmatE7, cmatE8, cmatF4, cmatG2,
                 cmatEt6, cmatEt7, cmatEt8, cmatFt4, cmatGt2]

        for typ in types:
            T = typ[1][0][0]
            r = typ[1][0][1]
            self.assertTrue(typ[0] == typ[1],
                    'Cartan matrix ' + T + str(r) + ' not identified!')
    

def suite():
    tests = ['test_chartabs']

    return unittest.TestLoader().loadTestsFromTestCase(TestCoreFuncs)

# Run the test suite.
unittest.TextTestRunner(verbosity=2).run(suite())
