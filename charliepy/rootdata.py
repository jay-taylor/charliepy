import numpy as np
import itertools
import sys

from . import algebra as alg
from . import core

# global constants
__crystaltypes__ = ["A","B","C","BC","D","E","F","G"]

# define what is imported under import *
__all__ = ['RootDatum']


class RootDatum:
    """Either a string giving the type or a pair of matrices X and Y such that
    the rows of X describe the roots in terms of the chosen basis and the rows
    of Y describe the coroots in terms of the corresponding dual basis.
    
    Optional keyword inputs

    trank: an integer giving the size of the central torus
    iso: a string, either 'sc' or 'ad' giving the isogeny type of the root datum
    """
    def __init__(self, *args, **kwargs):
        global __crystaltypes__
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

            self.rankss = len(cmat) - len([x for x in args if x == "BC"])
            self.rank = self.rankss + trank
            self.cartanmat = cmat
            self.cartantype = ctyp

            for typ in self.cartantype:
                if typ[0] not in __crystaltypes__:
                    raise ValueError("A root datum must be have a "
                                     "crystallographic root system!")

            self.X = alg.Matrix(np.zeros((len(cmat), len(cmat)+trank),
                dtype='int8'), copy=False)
            self.Y = alg.Matrix(np.zeros((len(cmat), len(cmat)+trank),
                dtype='int8'), copy=False)
            self.isogenytype = ''

            # Fix the isogeny type.
            if iso == 'sc':
                self.X[:,:len(cmat)] = self.cartanmat
                np.fill_diagonal(self.Y, 1)
                self.isogenytype = 'sc'
            else:
                np.fill_diagonal(self.X, 1)
                self.Y.T[:len(cmat),:] = self.cartanmat
                self.isogenytype = 'ad'
        else:
            self.X = alg.Matrix(args[0], copy=False)
            self.Y = alg.Matrix(args[1], copy=False)
            self.rankss = self.X.shape[0]
            self.rank = self.X.shape[1]
            self.isogenytype = ''
            self.cartanmat = np.dot(self.X, self.Y.T)
            if not core.is_intcartanmat(self.cartanmat):
                raise ValueError("The product of X and Y "
                                 "is not a Cartan matrix!")
            self.cartantype = core.cartantotype(self.cartanmat)

            nonred = []

            for x in self.cartantype:
                if x[0] not in __crystaltypes__:
                    raise ValueError("The product of X and Y is not "
                                     "a crystallographic Cartan matrix!")
                elif x[0] == "BC":
                    self.rankss -= 1
                    nonred.append([x[1][0],x[1][-1]])

        # Now get the roots.
        self.roots = core.roots_finintcartan(self.cartanmat, nonred=nonred)
        self.N = len(self.roots)/2

        # Now get the coroots.
        if np.array_equal(self.cartanmat, self.cartanmat.T):
            self.coroots = self.roots[:]
        else:
            self.coroots = core.roots_finintcartan(self.cartanmat.T,
                    nonred=[tuple([x[1],x[0]]) for x in nonred])

    def __repr__(self):
        out = 'RootDatum('
        for typ in self.cartantype:
            if typ[0] == "BC" or typ[0][0] == "~":
                out += typ[0] + ',' + str(len(typ[1])-1)
            else:
                out += typ[0] + ',' + str(len(typ[1]))
        if self.isogenytype:
            out += ',"' + self.isogenytype + '")'
        else:
            out += ')'

        return out

    def innerprod(self, x, y):
        """
        Gives the value of the bilinear form

        (x,y) = \sum_{\alpha*} <x,\alpha*><\alpha*,y>
        """
        try:
            r1 = self.roots[x]
            r2 = self.roots[y]
        except TypeError:
            r1 = x
            r2 = y
       
        # The cartan matrix is of the form <alpha_i*,alpha_j> and we have
        #
        #       <x, alpha_k*> = x_i<alpha_j*, alpha_i>a_kj
        #       <alpha_i*, y> = a_ij<alpha_j*, alpha_k>y_k
        #
        # where alpha_k* = a_ij alpha_j* is the kth coroot.
        left = np.einsum("i,ji,kj", r1, self.cartanmat, self.coroots)
        right = np.einsum("ij,jk,k", self.coroots, self.cartanmat, r2)
        return np.dot(left, right)













## Taken from http://code.activestate.com/recipes/576694/
#class OrderedSet(collections.MutableSet):
#
#    def __init__(self, iterable=None):
#        self.end = end = [] 
#        end += [None, end, end]         # sentinel node for doubly linked list
#        self.map = {}                   # key --> [key, prev, next]
#        if iterable is not None:
#            self |= iterable
#
#    def __len__(self):
#        return len(self.map)
#
#    def __contains__(self, key):
#        return key in self.map
#
#    def add(self, key):
#        if key not in self.map:
#            end = self.end
#            curr = end[1]
#            curr[2] = end[1] = self.map[key] = [key, curr, end]
#
#    def discard(self, key):
#        if key in self.map:        
#            key, prev, next = self.map.pop(key)
#            prev[2] = next
#            next[1] = prev
#
#    def __iter__(self):
#        end = self.end
#        curr = end[2]
#        while curr is not end:
#            yield curr[0]
#            curr = curr[2]
#
#    def __reversed__(self):
#        end = self.end
#        curr = end[1]
#        while curr is not end:
#            yield curr[0]
#            curr = curr[1]
#
#    def pop(self, last=True):
#        if not self:
#            raise KeyError('set is empty')
#        key = self.end[1][0] if last else self.end[2][0]
#        self.discard(key)
#        return key
#
#    def __repr__(self):
#        if not self:
#            return '%s()' % (self.__class__.__name__,)
#        return '%s(%r)' % (self.__class__.__name__, list(self))
#
#    def __eq__(self, other):
#        if isinstance(other, OrderedSet):
#            return len(self) == len(other) and list(self) == list(other)
#        return set(self) == set(other)

            










