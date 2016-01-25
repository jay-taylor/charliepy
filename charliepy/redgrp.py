import numpy as np
import itertools
import sys

from . import algebra as alg
from . import core

# Useful global constants
__redgrptypes__ = set(["A","B","C","D","E","F","G"])

# define what is imported under import *
#__all__ = ['RootDatum']


class ReductiveGroup:
    """Either a string giving the type or a pair of matrices X and Y such that
    the rows of X describe the roots in terms of the chosen basis and the rows
    of Y describe the coroots in terms of the corresponding dual basis.
    
    Optional keyword inputs

    trank: an integer giving the size of the central torus
    F: the Frobenius endomorphism. Should be a matrix describing an
        automorphism of the cocharacter group Y on the chosen basis. 
    p: the characteristic of the field defining the reductive group. If not
        specified then this is set to 0.
    cmp: generators for the component group G/G^0. We assume these
        representatives normalise the fixed torus T and the simple roots.
    iso: a string, either 'sc' or 'ad' giving the isogeny type of the group
    """
    def __init__(self, *args, **kwargs):
        global __redgrptypes__

        if isinstance(args[0], str):
            try:
                trank = kwargs['trank']
            except KeyError:
                trank = 0

            try:
                iso = kwargs['iso']
            except KeyError:
                iso = 'ad'

            cmat, ctyp, nonred = core.cartanmat(*args, extras=True)

            self.rankss = len(cmat)
            self.rank = self.rankss + trank
            self.cartanmat = cmat
            self.cartantype = ctyp

            for typ in self.cartantype:
                if typ[0] not in __redgrptypes__:
                    raise ValueError("A root datum must have a "
                                     "crystallographic root system!")

            self.X = alg.Matrix(np.zeros((self.rankss, self.rank),
                dtype='int8'), copy=False)
            self.Y = alg.Matrix(np.zeros((self.rankss, self.rank),
                dtype='int8'), copy=False)
            self.isogenytype = ''

            # Fix the isogeny type.
            if iso == 'sc':
                self.X[:,:self.rankss] = self.cartanmat
                np.fill_diagonal(self.Y, 1)
                self.isogenytype = 'sc'
            else:
                np.fill_diagonal(self.X, 1)
                self.Y.T[:self.rankss,:] = self.cartanmat
                self.isogenytype = 'ad'
        else:
            self.X = alg.Matrix(args[0], copy=False)
            self.Y = alg.Matrix(args[1], copy=False)
            self.rankss = self.X.shape[0]
            self.rank = self.X.shape[1]
            self.isogenytype = ''
            self.cartanmat = np.dot(self.Y, self.X.T)
            if not core.is_intcartanmat(self.cartanmat):
                raise ValueError("The product of X and Y "
                                 "is not a Cartan matrix!")
            self.cartantype = core.cartantotype(self.cartanmat)

            for x in self.cartantype:
                if x[0] not in __redgrptypes__:
                    raise ValueError("The product of X and Y is not "
                                     "a crystallographic Cartan matrix!")

        try:
            self.p = kwargs['p']
        except KeyError:
            self.p = 0

        # Now get the roots and also store them as a dict for fast lookups. This
        # annoyingly seems like a waste of memory but this dramatically improves
        # performance when deciding whether a tuple is a root.
        self.roots = core.roots_finintcartan(self.cartanmat)
        self.N = len(self.roots)/2
        self._dictroots = dict(itertools.izip(self.roots,
                                  xrange(len(self.roots))))

        # Now get the coroots.
        if np.array_equal(self.cartanmat, self.cartanmat.T):
            self.coroots = list(self.roots)
        else:
            self.coroots = core.roots_finintcartan(self.cartanmat.T)

        # Construct the values of the canonical bilinear form for the simple
        # roots. All bilinear form values are determined from this.
        tmp = np.dot(self.coroots, self.cartanmat)
        self._Bform = alg.Matrix(np.dot(tmp.T,tmp), copy=False)
        del tmp

        self.extraspecialsigns = [1 for i in xrange(self.N - self.rankss)]
        core.specialpairs(self)

    def __repr__(self):
        out = 'ReductiveGroup('
        for typ in self.cartantype:
            out += typ[0] + ',' + str(len(typ[1]))
        if self.isogenytype:
            out += ',"' + self.isogenytype + '")'
        else:
            out += ')'

        return out

    def changestructureconstants(self, signs = None):
        # Now get all special pairs.
        if signs:
            self.extraspecialsigns = signs
        else:
            self.extraspecialsigns = [1 for i in xrange(self.N - self.rankss)]
        
        core.specialpairs(self)

# Functions for the structure constant variables.
    def sc_N(self,x1,x2):
        sumisroot = False
        sgn = 1
        facs = [1,1]

        if x2 >= self.N: # beta < 0
            x1 += cmp(self.N, x1)*self.N
            x2 -= self.N
            sgn *= -1
            # N_a,b = -N_-a,-b
        if x1 >= self.N: # alpha < 0
            if x1 - self.N < x2: # alpha < 0 < alpha + beta < beta
                t = tuple([a+b for a,b in
                    itertools.izip(self.roots[x1],self.roots[x2])])
                if t in self._dictroots:
                    # (alpha) + (beta) + (-alpha-beta) = 0
                    # N_a,b = |a+b|^2*N_-a-b,a/|b|^2
                    #       = -|a+b|^2*N_a+b,-a/|b|^2
                    k = self._dictroots[t]
                    facs = [self.B(k,k), self.B(x2,x2)]
                    x1, x2 = k, x1 - self.N
                    sgn *= -1
                    sumisroot = True
                else:
                    return 0
            else: # alpha < alpha + beta < 0 < beta
                t = tuple([-a-b for a,b in
                    itertools.izip(self.roots[x1],self.roots[x2])])
                if t in self._dictroots:
                    # (alpha) + (beta) + (-alpha-beta) = 0
                    # N_a,b = |a+b|^2*N_b,-a-b/|a|^2
                    k = self._dictroots[t]
                    facs = [self.B(k,k), self.B(x1,x1)]
                    x1, x2 = x2, k
                    sumisroot = True
                else:
                    return 0
        if x1 > x2:
            # N_a,b = -N_b,a
            x1, x2 = x2, x1
            sgn *= -1
        if sumisroot or (x1, x2) in self.specialpairs:
            # Note, we're doing integer division here but it doesn't matter as
            # we know that the result must be an integer. However, we cautiously
            # use brackets to avoid premature rounding.
            return sgn*((facs[0]*self.specialpairs[(x1,x2)][3])//facs[1])
        else:
            return 0

    def sc_e(self,x1,x2):
        # This function returns 0 if alpha+beta is not a root.
        sumisroot = False
        sgn = 1

        if x2 >= self.N: # beta < 0
            x1 += cmp(self.N, x1)*self.N
            x2 -= self.N
            sgn *= -1
            # e_a,b = -e_-a,-b
        if x1 >= self.N: # alpha < 0
            if x1 - self.N < x2: # alpha < 0 < alpha + beta < beta
                t = tuple([a+b for a,b in
                    itertools.izip(self.roots[x1],self.roots[x2])])
                if t in self._dictroots:
                    # (alpha) + (beta) + (-alpha-beta) = 0
                    # e_a,b = e_-a-b,a = -e_a+b,-a
                    x1, x2 = self._dictroots[t], x1 - self.N
                    sgn *= -1
                    sumisroot = True
                else:
                    return 0
            else: # alpha < alpha + beta < 0 < beta
                t = tuple([-a-b for a,b in
                    itertools.izip(self.roots[x1],self.roots[x2])])
                if t in self._dictroots:
                    # (alpha) + (beta) + (-alpha-beta) = 0
                    # e_a,b = e_b,-a-b
                    x1, x2 = x2, self._dictroots[t]
                    sumisroot = True
                else:
                    return 0
        if x1 > x2:
            # e_a,b = -e_b,a
            x1, x2 = x2, x1
            sgn *= -1
        if sumisroot or (x1, x2) in self.specialpairs:
            return sgn*cmp(self.specialpairs[(x1,x2)][3],0)
        else:
            return 0

    # The following formulas will be used to compute the p values.
    #
    #   p_a,b = q_-a,b = q_a,-b
    #   p_a,b = p_-a,-b
    #   p_a,b = p_a,b+a - 1 = p_a,b-a + 1
    #   p_a,b = q_-a,b-a + 1 = q_a,-b-a - 1
    #
    def sc_p(self,x1,x2):
        if x2 >= self.N: # b < 0
            # p_a,b = p_-a,-b
            x2 -= self.N
            if x1 >= self.N:
                x1 -= self.N
            else:
                x1 += self.N

        if x1 >= self.N: # a < 0
            x1 -= self.N # Note, we already make alpha positive here.
            if x1 < x2:
                if (x1,x2) in self.specialpairs:
                    # p_a,b = q_-a,b
                    return self.specialpairs[(x1,x2)][2]
            else:
                # Note that x2 < x1 so x2+x1 is a positive root with height
                # strictly greater than x1, so (x1,x2+x1) will be a special
                # pair.
                if (x2,x1) in self.specialpairs:
                    # p_a,b = q_-a,b-a + 1
                    k = self.specialpairs[(x2,x1)][0]
                    if (x1,k) in self.specialpairs:
                        return self.specialpairs[(x1,k)][2] + 1
                    else:
                        return 1
            return 0

        # Can now assume a,b > 0.
        if (x1, x2) in self.specialpairs:
            return self.specialpairs[(x1,x2)][1]
        elif (x2, x1) in self.specialpairs:
            # p_a,b = 1 + p_a,b-a = 1 + q_a,a-b
            # Note that (a-b,a) is a special pair. Also, a-b is a root iff b-a
            # is a root.
            tmp = tuple([a-b for a,b in
                itertools.izip(self.roots[x1],self.roots[x2])])
            if tmp in self._dictroots:
                # We know a+b is a root so p_a,b <= 2.
                k = self._dictroots[tmp]
                if (k,x1) in self.specialpairs:
                    return 2
                else:
                    return 1
        # We now know a+b is not a root but we may still have b-a is a root.
        tmp = tuple([b-a for a,b in
            itertools.izip(self.roots[x1],self.roots[x2])])
        if tmp in self._dictroots:
            # p_a,b = p_a,b-a + 1
            k = self._dictroots[tmp]
            return self.sc_p(x1,k) + 1
        else:
            return 0

    # The following formulas will be used to compute the p values.
    #
    #   q_a,b = p_-a,b = p_a,-b
    #   q_a,b = q_-a,-b
    #   q_a,b = q_a,b+a + 1 = p_a,b-a - 1
    #   q_a,b = p_-a,b-a - 1 = p_a,-b-a + 1
    #
    def sc_q(self,x1,x2):
        if x2 >= self.N: # b < 0
            # q_a,b = q_-a,-b
            x2 -= self.N
            if x1 >= self.N:
                x1 -= self.N
            else:
                x1 += self.N

        if x1 >= self.N: # a < 0
            # q_a,b = p_-a,b
            x1 -= self.N
            if (x1,x2) in self.specialpairs:
                # a+b is a root ---> p_a,b defined
                return self.specialpairs[(x1,x2)][1]
            elif (x2,x1) in self.specialpairs:
                # p_a,b = p_a,b+a - 1
                k = self.specialpairs[(x2,x1)][0]
                if (x1,k) in self.specialpairs:
                    return self.specialpairs[(x1,k)][1] - 1
            # See if we're at the end of a chain.
            tmp = tuple([b-a for a,b in
                itertools.izip(self.roots[x1],self.roots[x2])])
            if tmp in self._dictroots:
                # p_a,b = p_a,b-a + 1
                k = self._dictroots[tmp]
                return self.sc_p(x1,k) + 1
            return 0

        # Can now assume a,b > 0.
        if (x1, x2) in self.specialpairs:
            return self.specialpairs[(x1,x2)][2]
        elif (x2, x1) in self.specialpairs:
            # q_a,b = q_a,b+a + 1
            # Note that (a,b+a) is a special pair.
            k = self.specialpairs[(x2,x1)][0]
            if (x1,k) in self.specialpairs:
                return self.specialpairs[(x1,k)][2] + 1
            else:
                return 1
        else:
            return 0

    def sc_eta(self, x1, x2):
        if x1 == x2 or abs(x1 - x2) == self.N:
            return 0
        p = self.sc_p(x1,x2)
        C = -self.C(x1,x2) # C = q - p
        sgn = 1
        if C > 0:
            for i in xrange(1,C):
                sgn *= self.sc_e(self._dictroots[tuple([b+i*a for a,b in
                    itertools.izip(self.roots[x1],self.roots[x2])])], x1)
        elif C < 0: # q - p < 0
            for i in xrange(C,0):
                sgn *= self.sc_e(self._dictroots[tuple([b+i*a for a,b in
                    itertools.izip(self.roots[x1],self.roots[x2])])], x1)

        return sgn*((-1)**p)



    def B(self, x, y):
        """
        Gives the value of the bilinear form

        (x,y) = sum_{alpha*} <alpha*, x><alpha*, y>
        """
        # alpha = \sum a_ialpha_i and beta = \sum b_jalpha_j
        #
        # (alpha,beta) = \sum_i\sum_j a_i(alpha_i,alpha_j)b_j
        #
        return np.einsum("i,ij,j", self.roots[x], self._Bform, self.roots[y])

    def C(self, x, y):
        """
        Gives the value of the Cartan integer <alpha*,beta> where x is the
        index of the coroot alpha* in the list self.coroots and y is the index
        of the root beta in the list self.roots.
        """
        # alpha = \sum a_ialpha_i* and beta = \sum b_jalpha_j
        #
        # (alpha,beta) = \sum_i\sum_j a_i(alpha_i*, alpha_j)b_j
        #
        return np.einsum("i,ij,j", self.coroots[x], self.cartanmat,
                self.roots[y])

    def applyreflection(self, x, y):
        # Applies the reflection corresponding to the root x to the root y and
        # gives the index of the resulting root in the list self.roots.

        C = np.einsum("i,ij,j", self.coroots[x], self.cartanmat,
                self.roots[y])
        return self._dictroots[tuple([a - C*b for a,b in
            itertools.izip(self.roots[y], self.roots[x])])]



def rootp(G,x,y):
    a = G.roots[x]
    b = G.roots[y]

    p = 0
    tmp = tuple([d-e for d,e in itertools.izip(b,a)])
    while tmp in G._dictroots:
        p += 1
        tmp = tuple([d-e for d,e in itertools.izip(tmp,a)])
    return p

def rootq(G,x,y):
    a = G.roots[x]
    b = G.roots[y]

    q = 0
    tmp = tuple([d+e for d,e in itertools.izip(b,a)])
    while tmp in G._dictroots:
        q += 1
        tmp = tuple([d+e for d,e in itertools.izip(tmp,a)])
    return q










