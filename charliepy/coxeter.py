from . import permutat
from . import utils
from . import braid
from . import cartan

import numpy as np
import itertools
import operator
import weakref
import functools

# define what is imported under import *
__all__ = [
    'CoxeterGroup',
    'reflectionsubgroup',
    'rootinclusion'
]

# Used for testing that Cartan matrices are crystalographic.
_crystaltypes = {"A", "B", "C", "D", "E", "F", "G"}

########################################################################
########################################################################
##                                                                    ##
##                      Main Class definitions                        ##
##                                                                    ##

class CoxeterGroup:
    """
    Creates an instance of a Coxeter group as a Python class.

    Let (W, S) be a finite Coxeter system with reflection representation
    V. If alpha is a root in V then we define the corresponding
    reflection to be

            x.s_alpha = x - <x, alpha*>alpha

    where alpha* is the coroot corresponding to alpha. Note that this is
    different to the definition used in CHEVIE. In particular, if C is the
    Cartan matrix of (W, S) then we have

            alpha_i.s_j = alpha_i - C_{ij}*alpha_j

    In particular, we have C_{ij} = <alpha_i, alpha_j*> where alpha_j*
    is the coroot corresponding to alpha_j.

    """
    def __init__(self, *args, **kwargs):
        # We don't allow reinitialisation of the group. This would mean
        # the group could change.
        if hasattr(self, "_name"):
            return None

        # First, assume we're given the types of the groups.
        if isinstance(args[0], str):
            cmat, ctyp, nonred = cartan.cartanmat(*args, extras=True)

            self.rank = len(cmat)
            self.cartanmat = cmat
            self.cartantype = ctyp

            for typ in self.cartantype:
                if typ[0] not in _crystaltypes:
                    raise ValueError("A Coxeter group must be of finite"
                                     "crystallographic type.")

        else:
            self.cartanmat = np.array(args[0], dtype='int')
            self.rank = self.cartanmat.shape[0]
            if not cartan.is_intcartanmat(self.cartanmat):
                raise ValueError("The input matrix must be a Cartan matrix.")
            self.cartantype = cartan.cartantotype(self.cartanmat)

            for x in self.cartantype:
                if x[0] not in _crystaltypes:
                    raise ValueError("The Cartan matrix must be of finite"
                                     "Crystallographic type.")

        # We give the mapping of the indices of the simple reflections
        # onto the standard indices for their Cartan type. For instance,
        # if
        #
        #       ["A", [6, 10, 3]]
        #
        # is part of the Cartan type of the group then L = _standardinds
        # is a list such that L[6] = 0, L[10] = 1, L[3] = 2. Storing
        # this information is useful when computing class fusions.
        stdinds = [None]*self.rank
        for typ, inds in self.cartantype:
            for i, j in enumerate(inds):
                stdinds[j] = i
        self._standardinds = stdinds

        # Now get the roots and the number of positive roots.
        self.roots = cartan.roots_finintcartan(self.cartanmat)
        self.N = len(self.roots)//2

        # Construct the Coxeter matrix from the Cartan matrix.
        self.coxetermat = np.ones_like(self.cartanmat, dtype='int')
        self.coxetermat *= 2
        np.fill_diagonal(self.coxetermat, 1)
        tmp = self.cartanmat*self.cartanmat.T
        conversion = [None, 3, 4, 6]
        for i in range(1, 4):
            inds = np.where(tmp == i)
            self.coxetermat[inds] = [conversion[i]]*len(inds[0])

        # Get the reflection degrees and the relative root lengths which
        # symmetrize the Cartan matrix.
        self.degrees = []
        self.symform = [None]*self.rank

        for typ, inds in self.cartantype:
            mod = utils.getmodule(typ)
            self.degrees += mod.degrees(len(inds))
            rootlengths = mod.rootlengths(len(inds))
            for ind, elm in zip(inds, rootlengths):
                self.symform[ind] = elm

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
            self.permgens[i] = permutat.Perm(
                [self.roots.index(tuple(np.dot(mat.T, alpha).tolist()))
                 for alpha in self.roots]
            )

        # We now compute for each root the lexicographically smallest
        # simple root which is contained in its W-orbit.
        reps = [None]*len(self.roots)
        self._schreier = dict()

        # Construct the orbit of the ith simple root. This is a basic
        # orbit computation algorithm. We simply keep applying
        # generators until we've seen everything. The dict _schreier is
        # such that for any i in range(self.rank) we have _schreier[i]
        # is a Schreier vector for the simple root alpha_i. This vector
        # allows us to reconstruct for any root alpha in the orbit of
        # alpha_i an element g in W such that alpha.g = alpha_i.
        for i in range(self.rank):
            # If we've already seen the simple root move on.
            if reps[i] is not None:
                continue
            reps[i] = i
            sch = [None]*len(self.roots)
            orb = [i]
            sch[i] = -1
            for item in orb:
                for k, p in enumerate(self.permgens):
                    j = item^p
                    if reps[j] is None:
                        reps[j] = i
                        orb.append(j)
                        sch[j] = k
            self._schreier[i] = tuple(sch)

        self.orbitrepresentatives = reps

        # We now construct the groups unique name. This is used for
        # comparisons and hashing. Here we have to do some sorting and
        # mangling as we want to make the following conversions in the
        # name for uniqueness
        #
        #   D3 -> A3    C2 -> B2   D2 -> A1xA1   B1, C1 -> A1
        #
        # but we don't want to change the naming in the Cartan type.
        # This step signifies succesful initilisation of the group.
        uniquename = [[elm[0][:], elm[1][:]] for elm in self.cartantype]
        for ind, typ in enumerate(uniquename):
            if len(typ[1]) == 3 and typ[0] == "D":
                uniquename[ind] = ["A", [typ[1][1], typ[1][0], typ[1][2]]]
            if len(typ[1]) == 2 and typ[0] == "C":
                uniquename[ind] = ["B", [typ[1][1], typ[1][0]]]
            if len(typ[1]) == 1 and (typ[0] == "B" or typ[0] == "C"):
                uniquename[ind][0] = "A"
            if len(typ[1]) == 2 and typ[0] == "D":
                uniquename[ind] = ["A", [typ[1][0]]]
                uniquename.append(["A", [typ[1][1]]])

        # Now sort the list.
        uniquename.sort(key=operator.itemgetter(0, 1))
        self._name = '*'.join(
                x[0] + "-" + ".".join(map(str, x[1])) for x in uniquename
            )

        # Define the embedding of the group into itself. We use here that the
        # group is hashable so can be a key in a dictionary. Note for this to be
        # the case we must define this after defining self._name.
        self.embeddings = weakref.WeakKeyDictionary(
            {self : tuple(range(self.rank))})

        # This keeps track of the fusion of conjugacy classes into groups with
        # embeddings. This will be computed for the group when the conjugacy
        # classes are computed.
        self.fusion = weakref.WeakKeyDictionary()

    def __repr__(self):
        return "CoxeterGroup(" + ", ".join(
            "'" + typ[0] + "', [" + ", ".join(map(str,typ[1])) + "]"
            for typ in self.cartantype
        ) + ")"

    def __eq__(self, other):
        """
        Two CoxeterGroup instances are considered to be equal if they
        have exactly the same labelled Dynkin diagrams, up to
        automorphisms of the Dynkin diagram. Note that if the connected
        components are isomorphic as Coxeter groups then the type
        labelling of the diagram, e.g., A3 or D3, does not matter; only
        that the vertices of the diagram are labelled in precisely the
        same way.

        Notes
        -----
        This function does not test whether the groups are isomorphic!

        Example
        -------
        >>> W = clp.CoxeterGroup("A", [0, 1, 2], "B", [3, 4, 5])
        >>> G = clp.CoxeterGroup("B", [3, 4, 5], "A", [0, 1, 2])
        >>> W.diagram()
        A3 :  0 - 1 - 2
        B3 :  3 - 4 > 5
        >>> G.diagram()
        B3 :  3 - 4 > 5
        A3 :  0 - 1 - 2
        >>> W == G
        True
        >>> G = clp.CoxeterGroup("A", [3, 4, 5], "B", [0, 1, 2])
        >>> G.diagram()
        A3 :  3 - 4 - 5
        B3 :  0 - 1 > 2
        >>> W == G
        False

        >>> W = clp.CoxeterGroup("A", [0, 1, 2])
        >>> G = clp.CoxeterGroup("D", [0, 1, 2])
        >>> W.diagram()
        A3 :  0 - 1 - 2
        >>> G.diagram()
        D3 :  1 - 0 - 2
        >>> W == G
        False
        >>> G = clp.CoxeterGroup("D", [1, 0, 2])
        D3 :  0 - 1 - 2
        >>> W == G
        True

        """
        return (isinstance(other, CoxeterGroup) and self._name == other._name)

    def __hash__(self):
        return hash(self._name)

    @property
    def size(self):
        """
        The cardinality of the group.

        """
        return functools.reduce(operator.mul, self.degrees, 1)

    def coroot(self, i):
        """
        Gives the coroot of the ith root in the list W.roots. The coroot
        is represented just as the roots are. In other words, this is a
        tuple whose jth entry is the coefficient of the coroot of the
        jth simple root.

        """
        # Assume alpha is the root we are considering and we have
        #
        #       alpha = sum_j a_j.alpha_j
        #
        # then using Lemma 2 of [Bou02, Ch. VI, S1, no. 1] we easily
        # deduce that
        #
        #    alpha* = sum_j a_j.len(alpha_j)//len(alpha).alpha_j*
        #
        # where len denotes the length of the root. If alpha is in the
        # same W-orbit as alpha_k then we have len(alpha) ==
        # len(alpha_k) so we can compute the coroot from the knowledge
        # of the relative root lengths, which are stored in the
        # attribute self.symform.
        k = self.orbitrepresentatives[i]
        return tuple((elm*self.symform[j])//self.symform[k]
                     for j, elm in enumerate(self.roots[i]))

    def __iter__(self):
        # Firstly, we write W as the direct product W_1 x ... x W_k of
        # its irreducible factors. Now for each irreducible factor X of
        # rank n we choose a maximal chain of parabolic subgroups
        #
        #   {1} = X_0 c X_1 c X_2 c ... c X_n = X
        #
        # whose quotients X_{i+1}/X_i are minimal. The choices are
        # obtained from the data module. Note this is better than simply
        # choosing [], [0], [0, 1], ..., etc. because this might yield
        # some parabolics with huge quotients depending on the numbering
        # of the Dynkin diagram.
        #
        # Now for each quotient X_{i+1}/X_i we compute a set of reduced
        # right coset representatives using Algorithm B on pg. 40 of
        # [GP00]. The generator then simply returns the cartesian
        # product of all these coset reps.
        cosetreps = [None]*self.rank
        identity = permutat.Perm(range(2*self.rank))

        c = 0
        for typ, inds in self.cartantype:
            pinds = utils.getmodule(typ).maxparachain(inds)

            for i in range(len(pinds)):
                Y = {identity}
                X = {identity}

                # Get the generators for the larger parabolic subgroup X_{i+1}.
                S = [(j, self.permgens[j]) for j in pinds[:i+1]]

                # Indices for the generators of smaller parabolic X_i.
                T = pinds[:i]

                while Y:
                    Z = set()
                    for x in Y:
                        for k, s in S:
                            # l(xs_k) > l(x) iff alpha_k.x**-1 is positive.
                            if k/x < self.N:
                                z = x*s
                                # l(s_jz) > l(z) iff alpha_j.z is positive.
                                if all(j^z < self.N for j in T):
                                    Z.add(z)
                    X |= Z
                    Y = Z

                # Make sure it returns a consistent result.
                cosetreps[c] = sorted(X)
                c += 1

        return _elms_gen(cosetreps[0], cosetreps[1:])

    def convert(self, w, lab):
        """
        Converts between different representations of elements.

        This converts an element w from any one of the following
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
        if isinstance(w, permutat.Perm):
            if lab == 'w': # perm -> word.
                word = []
                elm = w

                while True:
                    try:
                        i = next(j for j in range(self.rank) if j^elm >= self.N)
                        word.append(i)
                        elm = self.permgens[i]*elm
                    except StopIteration:
                        break

                return word
            
            if lab == 'm': # perm -> matrix.
                return np.array([self.roots[i^w]
                                 for i in range(self.rank)])

            if lab == 'p': # perm -> perm.
                return w

            if lab == 'c': # perm -> coxelm.
                return tuple(i^w for i in range(self.rank))

        if isinstance(w, list):
            if lab == 'p': # word -> perm.
                perm = permutat.Perm(range(2*self.N))
                for i in w:
                    perm *= self.permgens[i]

                return perm

            if lab == 'm': # word -> mat.
                elm = np.identity(self.rank, dtype='int')
                for i in w:
                    elm = np.dot(self.matgens[i], elm)
                return elm

            if lab == 'w': # word -> word.
                return w

            if lab == 'c': # word -> coxelm
                elm = permutat.Perm(range(len(self.roots)))
                for i in w:
                    elm *= self.permgens[i]
                return tuple([i^elm for i in range(self.rank)])

        if isinstance(w, tuple):
            if lab == 'w': # coxelm -> word.
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

            if lab == 'm': # coxelm -> matrix.
                return np.array([self.roots[w[i]] for i in range(self.rank)])

            if lab == 'p': # coxelm -> perm.
                # We first convert to a matrix then produce the
                # permutation.
                mat = np.array([self.roots[w[i]] for i in range(self.rank)])

                return permutat.Perm(
                    [self.roots.index(tuple(np.dot(mat.T, alpha).tolist()))
                     for alpha in self.roots]
                )

            if lab == 'c':
                return w

        if isinstance(w, np.ndarray):
            if lab == 'c': # matrix -> coxelm.
                # The rows of the matrix are the images of the simple
                # roots under the reflection.
                return tuple([self.roots.index(alpha)
                              for alpha in map(tuple, w)])

            if lab == 'p': # matrix -> perm.
                # If alpha = sum_i a_i*alpha_i where alpha_i is a simple
                # root then we have alpha.s_j = sum_i a_i*(alpha_i.s_j).
                return permutat.Perm(
                    [self.roots.index(tuple(np.dot(w.T, alpha).tolist()))
                     for alpha in self.roots]
                )

            if lab == 'w': # matrix -> word.
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

            if lab == 'm': # matrix -> matrix.
                return w

        raise ValueError("Elements cannot be of type {}.".format(type(w)))

    def length(self, w):
        """
        Returns the length of a reduced expression of the element.

        """
        if isinstance(w, permutat.Perm):
            count = 0
            for i in range(self.N):
                if i^w >= self.N:
                    count += 1
            return count

        if isinstance(w, list):
            return len(w)

        if isinstance(w, tuple):
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

        if isinstance(w, np.ndarray):
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

        raise ValueError("Elements cannot be of type {}.".format(type(w)))

    def longestelement(self, J=None, lab='w'):
        """
        This returns the longest element of the Coxeter group as a permutation,
        word, or matrix. The default is to return the longest word.

        """
        if J is None:
            pgens = list(enumerate(self.permgens))
        else:
            pgens = [(i, self.permgens[i]) for i in J]

        N = self.N

        if lab == 'w':
            word = []
            w0 = permutat.Perm(range(2*N))
            flag = True
            while flag:
                flag = False
                for i, s in pgens:
                    if i^w0 < N:
                        flag = True
                        w0 = s*w0
                        word.append(i)
                        break

            return word

        if lab == 'p':
            w0 = permutat.Perm(range(2*N))
            flag = True
            while flag:
                flag = False
                for i, s in pgens:
                    if i^w0 < N:
                        flag = True
                        w0 = s*w0
                        break

            return w0

        raise ValueError("Type not specified correctly.")

    def leftdescentset(self, w):
        """
        This returns a list of indices of the generators s of W
        satisfying l(sw) < l(w).

        """
        # By 1.1.9 of [GP00] we have l(sw) < l(w) iff alpha_s.w <= 0.
        # Thus we need to know which rows of the matrix representing w
        # have all entries <= 0.
        if isinstance(w, tuple):
            # First convert to matrix then apply the matrix process.
            mat = np.array([self.roots[w[i]] for i in range(self.rank)])
            return np.where(np.all(mat <= 0, axis=1))[0].tolist()

        if isinstance(w, np.ndarray):
            return np.where(np.all(w <= 0, axis=1))[0].tolist()

        if isinstance(w, permutat.Perm):
            # l(sw) < l(w) iff alpha_s.w < 0.
            return [i for i in range(self.rank) if i^w >= self.N]

        raise ValueError("Cannot accept elements of type {}.".format(type(w)))

    def rightdescentset(self, w, lab):
        """
        This returns a list of indices of the generators s of W
        satisfying l(ws) < l(w).

        """
        if isinstance(w, tuple):
            # First convert to matrix then apply the matrix process.
            mat = np.array([self.roots[w[i]] for i in range(self.rank)])
            inv = np.linalg.inv(mat)
            return np.where(np.all(inv <= 0, axis=1))[0].tolist()

        if isinstance(w, np.ndarray):
            # Note that when numpy inverts the matrix it will return a
            # matrix whose dtype is float but this doesn't matter
            # because the result will be integral and we only care about the
            # sign.
            inv = np.linalg.inv(w)
            return np.where(np.all(inv <= 0, axis=1))[0].tolist()

        if isinstance(w, permutat.Perm):
            # l(ws) < l(w) iff alpha_s.(w**-1) < 0.
            return [i for i in range(self.rank) if i/w >= self.N]

        raise ValueError("Type not specified correctly.")

    def diagram(self, parent=None):
        """
        Prints the diagram of the Coxeter Group. If the group, say H,
        has an embedding into the group W then setting parent=W will
        print the diagram with the labels provided by the roots in W.

        Examples
        --------
        >>> W = clp.CoxeterGroup("B",7)
        >>> H = clp.reflectionsubgroup(W, [2*W.N-1,0,1,2,3,5,6])
        >>> H.diagram()
                          3
                        /
        D5 :  0 - 1 - 2
                        \
                          6
        B2 :  4 > 5
        >>> H.diagram(parent=W)
                          3
                        /
        D5 :  0 - 1 - 2
                        \
                          36
        B2 :  5 > 6

        """
        for typ, inds in self.cartantype:
            if parent is not None:
                inds = [self.embeddings[parent][i] for i in inds]
            utils.getmodule(typ).diagram(inds)

    def reflection(self, i):
        """
        Returns the reflection of the ith root as a permutation.

        """
        # Assume alpha_j is the simple root representative for the
        # W-orbit containing alpha_i. We first get the Schreier vector
        # for alpha_j and then use a standard algorithm to determine an
        # element g in W such that alpha_i = alpha_j.g. The reflection of
        # alpha_i is then simply
        #
        #           (s_j)^g = g**(-1) * s_j * g,
        #
        # where s_j is the simple reflection of alpha_j.
        if i < self.rank:
            return self.permgens[i]

        j = self.orbitrepresentatives[i]
        vect = self._schreier[j]
        g = permutat.Perm(range(len(self.roots)))
        k = vect[i]

        while k != -1:
            # Assume alpha_i is the current root we're considering and
            # vect[i] = k, then we have k is the index of a simple root
            # and moreover that 
            #
            #       alpha_i = (alpha_i.s_k).s_k 
            #
            # because the generators are involutions. Continuing in this
            # way we obtain g.
            g = self.permgens[k]*g
            i ^= self.permgens[k]
            k = vect[i]

        return self.permgens[j]^g

    def braid(self, w):
        """
        Returns the element of the braid monoid determined by w. We have w is
        either a word in the braid monoid, or a list of permutations.

        """
        return braid.BraidElement(self, w)

def _elms_gen(multd, pool):
    if len(pool):
        return _elms_gen(
            (x*y for x, y in itertools.product(multd, pool[0])), pool[1:])
    else:
        return iter(multd)




########################################################################
########################################################################
##                                                                    ##
##                      Subgroups/Subcosets                           ##
##                                                                    ##
        
def reflectionsubgroup(W, J):
    """
    Constructs an instance of CoxeterGroup as a subgroup.

    """
    # It's easy to get the Cartan matrix by slicing if J defines a
    # standard parabolic subgroup of W.
    gens = set(range(W.rank))
    if gens.issuperset(J):
        P = CoxeterGroup(W.cartanmat[np.ix_(J, J)])
        P.embeddings.update({W : tuple(J)})

        return P

    # Let W' be the subgroup of W generated by the reflections of the
    # roots indexed by J. By 8.2 of [Hum90] we can obtain a set of
    # simple reflections for W' as follows. If T' is the set of
    # reflections of W' then for any w in W' we define a set
    #
    #       N(w) = {t in T' | l(tw) < l(w)},
    #
    # where l denotes the length function of W. A generating set S' of
    # W' is then given by all those t in T' such that N(t) = {t}.

    # First get all the roots and reflections of W'.
    reflsJ = [W.reflection(i) for i in J]
    rootsJ = J[:]
    seen = set(rootsJ)
    for i in rootsJ:
        for s in reflsJ[:len(J)]: # Loop over the generators.
            j = i^s
            if not j in seen:
                rootsJ.append(j)
                seen.add(j)
                reflsJ.append(W.reflection(j))

    simples = [rootsJ[i]
        for i, r in enumerate(reflsJ)
            if len(set([i for i in range(W.N) if i^r >= W.N]) & seen) == 1
               and rootsJ[i] < W.N
    ]
    simples.sort()

    ## We now have to compute the Cartan matrix, which simply
    ## involves computing the Cartan integers. Assume we have
    ##
    ##       alpha = sum_i a_i.alpha_i
    ##       beta* = sum_i b_i.alpha_i*
    ##
    ## then clearly we have
    ##
    ##   <alpha,beta*> = sum_isum_j a_i<alpha_i,alpha_j*>b_j
    ##
    ## so this can easily be computed from the Cartan matrix of W but we
    ## have to get the corresponding coroots first
    roots = np.empty((len(J), W.rank), dtype='int')
    coroots = np.empty((len(J), W.rank), dtype='int')
    for ind, elm in enumerate(simples):
        roots[ind] = W.roots[elm]
        coroots[ind] = W.coroot(elm)
    cmat = np.einsum("ij,jk,kl", roots, W.cartanmat, coroots.T)

    P = CoxeterGroup(cmat)
    P.embeddings.update({W : tuple(simples)})

    return P

def rootinclusion(W, H, i):
    """
    Returns the index of the root H.root[i] in the list W.roots.

    """
    L = list(zip(*[W.roots[j] for j in H.embeddings[W]]))
    r = tuple([sum([k*l for k, l in zip(L[n], H.roots[i])])
               for n in range(W.rank)])
    return W.roots.index(r)


        
