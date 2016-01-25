from __future__ import print_function

from . import algebra as alg
from . import data
from . import utils
from . import core

import numpy as np
import itertools
import operator
import functools
import sys

__crystaltypes__ = {"A", "B", "C", "D", "E", "F", "G"}

class CoxeterGroup:
    """Creates an instance of a Coxeter group as a Python class.

    Let (W, S) be a finite Coxeter system with reflection representation
    V. If alpha is a root in V then we define the corresponding
    reflection to be

            s_alpha(x) = x - <x, alpha*>alpha

    where alpha* is the coroot corresponding to alpha. Note that this
    different the definitions used in CHEVIE. In particular, if C is the
    Cartan matrix of (W, S) then we have

            s_i(alpha_j) = alpha_j - C_{ji}*alpha_i

    In particular, we have C_{ij} = <alpha_i, alpha_j*> where alpha_j*
    is the coroot corresponding to alpha_j.

    """
    def __init__(self, *args, **kwargs):
        # We don't allow reinitialisation of the group. This would mean
        # the group could change.
        if hasattr(self, "__name__"):
            return None

        # First, assume we're given the types of the groups.
        if isinstance(args[0], str):
            cmat, ctyp, nonred = core.cartanmat(*args, extras=True)

            self.rank = len(cmat)
            self.cartanmat = cmat
            self.cartantype = ctyp

            for typ in self.cartantype:
                if typ[0] not in __crystaltypes__:
                    raise ValueError("A Coxeter group must be of finite"
                                     "crystallographic type!")

        else:
            self.cartanmat = args[0]
            if not core.is_intcartanmat(self.cartanmat):
                raise ValueError("The input matrix must be a Cartan matrix!")
            self.cartantype = core.cartantotype(self.cartanmat)

            for x in self.cartantype:
                if x[0] not in __crystaltypes__:
                    raise ValueError("The Cartan matrix must be of finite"
                                     "Crystallographic type!")

        # Now get the roots and the number of positive roots.
        self.roots = core.roots_finintcartan(self.cartanmat)
        self.N = len(self.roots)//2

        # Get the reflection degrees.
        self.degrees = []
        for typ in self.cartantype:
            self.degrees += sys.modules[
                'charliepy.data.typ1{}'.format(typ[0])].degrees(len(typ[1]))
        self.degrees.sort()
        self.degrees = tuple(self.degrees)

        # Matrices and permutations of the simple reflections. Note that
        # our matrices are acting on row vectors from the right.
        # Therefore as we have
        #
        #       alpha_j.s_i = alpha_j - c_{ji}*alpha_i
        #
        # we obtain the matrix of s_i acting on the natural module for W
        # by subtracting the ith column of the Cartan matrix from the
        # ith column of the identity matrix. The permutations act on
        # list(range(len(self.roots))) from the left.
        self.matgens = [None]*self.rank
        self.permgens = [None]*self.rank

        for i in range(self.rank):
            mat = np.identity(self.rank, dtype='int')
            mat[:, i] -= self.cartanmat[:, i]
            self.matgens[i] = mat

            # If alpha = sum_i a_i*alpha_i where alpha_i is a simple
            # root then we have alpha.s_j = sum_i a_i*(alpha_i.s_j).
            self.permgens[i] = alg.Perm(
                [self.roots.index(tuple(np.dot(mat.T, alpha).tolist()))
                 for alpha in self.roots]
            )

        # We now construct the groups unique name. This is used for
        # comparisons and hashing. This is a string obtained from the
        # Cartan type, which is first sorted to make it unique. This
        # step signifies succesful initilisation of the group.
        self.__name__ = '.'.join(map(str, utils.flatten(
            sorted(self.cartantype, key=operator.itemgetter(0, 1)))))

    def __repr__(self):
        return "CoxeterGroup(" + ",".join("'" + typ[0] + "'," + str(len(typ[1]))
                                          for typ in self.cartantype) + ")"

    def __len__(self):
        return functools.reduce(operator.mul, self.degrees)

    def __eq__(self, other):
        """
        Two CoxeterGroup instances are considered to be equal if they
        have exactly the same irreducible factors with the same
        labelling of the simple roots. The ordering of these factors
        does not matter.

        Notes
        -----
        This function does not test whether the groups are isomorphic!
        For instance, this function will return False when comparing 


        Example
        -------
        >>> W = clp.CoxeterGroup("A", [0, 1, 2], "B", [3, 4, 5])
        >>> G = clp.CoxeterGroup("B", [3, 4, 5], "A", [0, 1, 2])
        >>> W
        CoxeterGroup('A',3,'B',3)
        >>> G
        CoxeterGroup('B',3,'A',3)
        >>> W == G
        True
        >>> G = clp.CoxeterGroup("A", [3, 4, 5], "B", [0, 1, 2])
        >>> W == G
        False


        """
        return (isinstance(other, CoxeterGroup)
                and self.__name__ == other.__name__)

    def __hash__(self):
        return hash(self.__name__)

    def convert(self, w, s):
        """
        Converts between different representations of elements.

        This converts and element w from any one of the following
        representations to another:

            'c' - coxelem
            'm' - matrix
            'p' - permutation
            'w' - word

        To do this the function excepts the argument form, which takes
        as input a 2 letter string combination. The first character of
        the string denotes the current representation of w and the
        second character denotes the desired representation of w.

        """
        # To obtain a reduced word in any scenario we use Algorithm A
        # from [GP00, pg. 9].

        if s[0] == 'c':
            if s[1] == 'm': # coxelm -> matrix.
                return np.array([self.roots[w[i]] for i in range(self.rank)])

            if s[1] == 'p': # coxelm -> perm.
                # We first convert to a matrix then produce the
                # permutation.
                mat = np.array([self.roots[w[i]] for i in range(self.rank)])

                return alg.Perm(
                    [self.roots.index(tuple(np.dot(mat.T, alpha).tolist()))
                     for alpha in self.roots]
                )

            if s[1] == 'w': # coxelm -> word.
                # Again we first convert to a matrix.
                word = []
                mat = np.array([self.roots[w[i]] for i in range(self.rank)])

                while True:
                    # Indices of simple roots whose image under elm has
                    # a simple root with a negative coefficient.
                    neginds = np.where(np.any(mat < 0, axis=1))[0]
                    if neginds.size:
                        i = neginds[0]
                        word.append(i)
                        mat = np.dot(self.matgens[i], mat)
                    else:
                        break

                return word[::-1]

        if s[0] == 'p':
            if s[1] == 'c': # perm -> coxelm.
                return tuple(w.perm[:self.rank].tolist())

            if s[1] == 'm': # perm -> matrix.
                return np.array([self.roots[w.perm[i]]
                                 for i in range(self.rank)])

            if s[1] == 'w': # perm -> word.
                word = []
                elm = w

                while True:
                    try:
                        i = next(j for j, x in enumerate(elm.perm[:self.rank])
                                 if x >= self.N)
                        word.append(i)
                        elm = self.permgens[i]*elm
                    except StopIteration:
                        break

                return word[::-1]

        if s[0] == 'm':
            if s[1] == 'c': # matrix -> coxelm.
                # The rows of the matrix are the images of the simple
                # roots under the reflection.
                return tuple([self.roots.index(alpha)
                              for alpha in map(tuple, w)])

            if s[1] == 'p': # matrix -> perm.
                # If alpha = sum_i a_i*alpha_i where alpha_i is a simple
                # root then we have alpha.s_j = sum_i a_i*(alpha_i.s_j).
                return alg.Perm(
                    [self.roots.index(tuple(np.dot(w.T, alpha).tolist()))
                     for alpha in self.roots]
                )

            if s[1] == 'w': # matrix -> word.
                # This uses Algorithm A from [GP00, pg. 9].
                word = []
                elm = w

                while True:
                    # Indices of simple roots whose image under elm has
                    # a simple root with a negative coefficient.
                    neginds = np.where(np.any(elm < 0, axis=1))[0]
                    if neginds.size:
                        i = neginds[0]
                        word.append(i)
                        elm = np.dot(self.matgens[i], elm)
                    else:
                        break

                return word[::-1]

        if s[0] == 'w':
            if s[1] == 'c': # word -> coxelm
                elm = np.arange(len(self.roots))
                for i in w:
                    elm = elm[self.permgens[i].perm]
                return tuple(elm[:self.rank].tolist())

            if s[1] == 'p': # word -> perm.
                elm = np.arange(len(self.roots))
                for i in w:
                    elm = elm[self.permgens[i].perm]
                return alg.Perm(elm)

            if s[1] == 'm': # word -> mat.
                elm = np.identity(self.rank, dtype='int')
                for i in w:
                    elm = np.dot(self.matgens[i], elm)
                return elm

        raise ValueError("Types not specified correctly.")

    def length(self, w, s):
        """
        Returns the length of a reduced expression of the element.

        """
        if s == 'c':
            # First convert to matrix then apply the matrix process.
            length = 0
            mat = np.array([self.roots[w[i]] for i in range(self.rank)])

            while True:
                # Find index of lexicographically smallest simple root
                # whose image under elm has a simple root with a
                # negative coefficient.
                neginds = np.where(np.any(elm < 0, axis=1))[0]
                if neginds.size:
                    length += 1
                    elm = np.dot(self.matgens[neginds[0]], elm)
                else:
                    break

            return length

        if s == 'm':
            # This uses Algorithm A from [GP00, pg. 9].
            length = 0
            elm = w

            while True:
                # Find index of lexicographically smallest simple root
                # whose image under elm has a simple root with a
                # negative coefficient.
                neginds = np.where(np.any(elm < 0, axis=1))[0]
                if neginds.size:
                    length += 1
                    elm = np.dot(self.matgens[neginds[0]], elm)
                else:
                    break

            return length

        if s == 'p':
            return np.count_nonzero(w.perm[:self.N] >= self.N)

        if s == 'w':
            return len(w)

        raise ValueError("Type not specified correctly.")

    def reducedword(self, w):
        """
        Converts a word into a reduced word.

        """
        pass

    def leftdescentset(self, w, s):
        """
        This returns a list of indices of the generators s of W
        satisfying l(sw) < l(w).

        """
        # By 1.1.9 of [GP00] we have l(sw) < l(w) iff alpha_s.w <= 0.
        # Thus we need to know which rows of the matrix representing w
        # have all entries <= 0.
        if s == 'c':
            # First convert to matrix then apply the matrix process.
            mat = np.array([self.roots[w[i]] for i in range(self.rank)])
            return np.where(np.all(mat <= 0, axis=1))[0].tolist()

        if s == 'm':
            return np.where(np.all(w <= 0, axis=1))[0].tolist()

        if s == 'p':
            return np.where(w.perm >= self.N)[0].tolist()

        raise ValueError("Type not specified correctly.")

    def rightdescentset(self, w, s):
        """
        This returns a list of indices of the generators s of W
        satisfying l(sw) > l(w).

        """
        if s == 'c':
            # First convert to matrix then apply the matrix process.
            mat = np.array([self.roots[w[i]] for i in range(self.rank)])
            inv = np.linalg.inv(mat)
            return np.where(np.all(inv <= 0, axis=1))[0].tolist()

        if s == 'm':
            # Note that when numpy inverts the matrix it will return a
            # matrix whose dtype is float but this doesn't matter
            # because we only care whether the entry is negative or not.
            inv = np.linalg.inv(w)
            return np.where(np.all(inv <= 0, axis=1))[0].tolist()

        if s == 'p':
            # Numpy's argsort command precisely inverts the permutation.
            inv = np.argsort(w)
            return np.where(inv >= self.N)[0].tolist()

        raise ValueError("Type not specified correctly.")

    def diagram(self):
        """
        Prints the diagram of the Coxeter Group.

        """
        for typ in self.cartantype:
            sys.modules['charliepy.data.typ1{}'.format(typ[0])].diagram(typ[1])











