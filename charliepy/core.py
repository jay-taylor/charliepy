from . import permutat

import numpy as np
import itertools
import sys
import operator
import weakref
import functools
import collections

# define what is imported under import *
__all__ = [
    'CoxeterGroup',
    'CoxeterCoset',
    'reflectionsubgroup',
    'cartanmat',
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
    """Creates an instance of a Coxeter group as a Python class.

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
            cmat, ctyp, nonred = cartanmat(*args, extras=True)

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
            if not is_intcartanmat(self.cartanmat):
                raise ValueError("The input matrix must be a Cartan matrix.")
            self.cartantype = cartantotype(self.cartanmat)

            for x in self.cartantype:
                if x[0] not in _crystaltypes:
                    raise ValueError("The Cartan matrix must be of finite"
                                     "Crystallographic type.")

        # Now get the roots and the number of positive roots.
        self.roots = roots_finintcartan(self.cartanmat)
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

        for typ in self.cartantype:
            mod = 'charliepy.data.typ1{}'.format(typ[0])
            self.degrees += sys.modules[mod].degrees(len(typ[1]))
            rootlengths = sys.modules[mod].rootlengths(len(typ[1]))
            for ind, elm in zip(typ[1], rootlengths):
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

    def convert(self, w, lab):
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

        if lab[0] == 'c':
            if lab[1] == 'm': # coxelm -> matrix.
                return np.array([self.roots[w[i]] for i in range(self.rank)])

            if lab[1] == 'p': # coxelm -> perm.
                # We first convert to a matrix then produce the
                # permutation.
                mat = np.array([self.roots[w[i]] for i in range(self.rank)])

                return permutat.Perm(
                    [self.roots.index(tuple(np.dot(mat.T, alpha).tolist()))
                     for alpha in self.roots]
                )

            if lab[1] == 'w': # coxelm -> word.
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

        if lab[0] == 'p':
            if lab[1] == 'c': # perm -> coxelm.
                return tuple([i^w for i in range(self.rank)])

            if lab[1] == 'm': # perm -> matrix.
                return np.array([self.roots[i^w]
                                 for i in range(self.rank)])

            if lab[1] == 'w': # perm -> word.
                word = []
                elm = w

                while True:
                    try:
                        i = next(j for j in range(self.rank) if j^elm >= self.N)
                        word.append(i)
                        elm = self.permgens[i]*elm
                    except StopIteration:
                        break

                return word[::-1]

        if lab[0] == 'm':
            if lab[1] == 'c': # matrix -> coxelm.
                # The rows of the matrix are the images of the simple
                # roots under the reflection.
                return tuple([self.roots.index(alpha)
                              for alpha in map(tuple, w)])

            if lab[1] == 'p': # matrix -> perm.
                # If alpha = sum_i a_i*alpha_i where alpha_i is a simple
                # root then we have alpha.s_j = sum_i a_i*(alpha_i.s_j).
                return permutat.Perm(
                    [self.roots.index(tuple(np.dot(w.T, alpha).tolist()))
                     for alpha in self.roots]
                )

            if lab[1] == 'w': # matrix -> word.
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

        if lab[0] == 'w':
            if lab[1] == 'c': # word -> coxelm
                elm = permutat.Perm(range(len(self.roots)))
                for i in w:
                    elm *= self.permgens[i]
                return tuple([i^elm for i in range(self.rank)])

            if lab[1] == 'p': # word -> perm.
                elm = permutat.Perm(range(len(self.roots)))
                for i in w:
                    elm *= self.permgens[i]
                return elm

            if lab[1] == 'm': # word -> mat.
                elm = np.identity(self.rank, dtype='int')
                for i in w:
                    elm = np.dot(self.matgens[i], elm)
                return elm

        raise ValueError("Types not specified correctly.")

    def length(self, w, lab):
        """
        Returns the length of a reduced expression of the element.

        """
        if lab == 'c':
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

        if lab == 'm':
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

        if lab == 'p':
            count = 0
            for i in range(self.N):
                if i^w >= self.N:
                    count += 1
            return count

        if lab == 'w':
            return len(w)

        raise ValueError("Type not specified correctly.")

    def longestelement(self, J=None, lab='w'):
        """
        This returns the longest element of the Coxeter group as a permutation,
        word, or matrix. The default is to return the longest word.

        """
        if J is None:
            J = list(range(self.rank))
            cartantype = self.cartantype
        else:
            cartantype = cartantotype(self.cartanmat[np.ix_(J, J)])

        if lab == 'w':
            w0 = []
            for typ in cartantype:
                mod = 'charliepy.data.typ1{}'.format(typ[0])
                w0.extend([J[i] for i in sys.modules[mod].longestword(typ[1])])
            return w0

        if lab == 'p':
            w0 = permutat.Perm(range(2*self.N))
            for typ in cartantype:
                mod = 'charliepy.data.typ1{}'.format(typ[0])
                for i in sys.modules[mod].longestword(typ[1]):
                    w0 *= self.permgens[J[i]]
            return w0

        raise ValueError("Type not specified correctly.")

    def leftdescentset(self, w, lab):
        """
        This returns a list of indices of the generators s of W
        satisfying l(sw) < l(w).

        """
        # By 1.1.9 of [GP00] we have l(sw) < l(w) iff alpha_s.w <= 0.
        # Thus we need to know which rows of the matrix representing w
        # have all entries <= 0.
        if lab == 'c':
            # First convert to matrix then apply the matrix process.
            mat = np.array([self.roots[w[i]] for i in range(self.rank)])
            return np.where(np.all(mat <= 0, axis=1))[0].tolist()

        if lab == 'm':
            return np.where(np.all(w <= 0, axis=1))[0].tolist()

        if lab == 'p':
            return [i for i in range(self.rank) if i^w >= self.N]

        raise ValueError("Type not specified correctly.")

    def rightdescentset(self, w, lab):
        """
        This returns a list of indices of the generators s of W
        satisfying l(sw) > l(w).

        """
        if lab == 'c':
            # First convert to matrix then apply the matrix process.
            mat = np.array([self.roots[w[i]] for i in range(self.rank)])
            inv = np.linalg.inv(mat)
            return np.where(np.all(inv <= 0, axis=1))[0].tolist()

        if lab == 'm':
            # Note that when numpy inverts the matrix it will return a
            # matrix whose dtype is float but this doesn't matter
            # because the result will be integral and we only care about the
            # sign.
            inv = np.linalg.inv(w)
            return np.where(np.all(inv <= 0, axis=1))[0].tolist()

        if lab == 'p':
            # Numpy's argsort command precisely inverts the permutation.
            inv = np.argsort(w)
            return np.where(inv[:self.rank] >= self.N)[0].tolist()

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
        for typ in self.cartantype:
            if parent is None:
                inds = typ[1]
            else:
                inds = [self.embeddings[parent][i] for i in typ[1]]
            sys.modules['charliepy.data.typ1{}'.format(typ[0])].diagram(inds)

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


class CoxeterCoset:
    """
    Creates an instance of a Coxeter coset.

    Parameters
    ----------
    W : Coxeter group
        An instance of :class:`.CoxeterGroup`.
    phi : permutation
        Permutation of the Coxeter generators inducing a Coxeter
        automorphism.

    Attributes
    ----------
    phi : permutation
        Permutation of the Coxeter generators defining the automorphism.
    orbits : list
        Orbits of isomorphic factors cyclically permuted by the
        automorphism phi.
    phitype : list
        Gives the type of the automorphism phi.
    order : int
        Cardinality of the coset (same as the atribute group.order).
    
    Notes
    -----
    The permutation ``phi`` may be specified using any list like object
    (for example, an instance of :class:`~.permutat.perm`).  If ``i`` is
    the index of a simple root in ``W.roots`` then ``phi[i]`` may be the
    index for the negative of a simple root.  Using this one may obtain
    natural Coxeter automorphisms using the longest element (see the
    examples).

    Examples
    --------
    The following example constructs the Coxeter coset of the Coxeter
    group of type `\mathrm{A}_4` together with the unique non-trivial
    Coxeter automorphism (induced by the longest element).::

        >>> W = clp.CoxeterGroup("A", 4)
        >>> phi = clp.longestperm(W)
        >>> Wphi = clp.CoxeterCoset(W, phi)
        >>> Wphi.phitype
        [['A', [0, 1, 2, 3], (0, 3)(1, 2)]]

    """
    def __init__(self, W, phi=permutat.Perm([])):
        self.group = W

        # Define the automorphism phi converting any negative roots to
        # positive roots. This means we can use something like the
        # longest element of the Weyl group as a permutation.
        if not isinstance(phi, permutat.Perm):
            raise TypeError("phi must be a permutation.")

        self.phi = phi
        
        # Make sure phi preserves the Cartan matrix. We use here that the
        # Coxeter matrix is symmetric to make things quicker. This should be
        # changed to Coxeter matrix in the future.
        for s in range(self.group.rank):
            for t in range(s):
                if (self.group.cartanmat[s][t]
                        != self.group.cartanmat[s^self.phi][t^self.phi]):
                    print("Permutation does not preserve Cartan matrix!")
                    return None
        
        # Determine the orbits of phi acting cyclically on isomorphic factors.
        typsets = [set(x[1]) for x in self.group.cartantype]
        indperm = self.phi.onsets(typsets)
        self.orbits = indperm.cycles()

        # Determine the twist induced by phi on each orbit and then construct
        # the type of the coxetercoset.
        self.phitype = []

        # If x is an orbit of phi of length n then phi**n stabilises each
        # component of x. We choose the first entry in the orbit as a
        # representative for the automorphism induced by phi**n on each
        # component in that orbit.
        for x in self.orbits:
            self.phitype.append(self.group.cartantype[x[0]][:])
            self.phitype[-1].append(permutat.restrictedperm((self.phi)**len(x),
                    self.group.cartantype[x[0]][1]))

        # Set the order to be the order of the underlying group.
        self.size = self.group.size

    def __repr__(self):
        """
        Lists the type of the Coxeter coset. A collection of factors
        in parentheses indicates that these factors are cyclically
        permuted by the automorphism. A number preceeding the first
        parenthesis denotes that the automorphism induced on that factor
        is of that given order. Orbits of factors are seperated by an x.

        Examples
        --------
        >>> phi = range(8, -1, -1)
        >>> W = clp.CoxeterGroup("A", 3, "A", 3, "A", 3)
        >>> clp.CoxeterCoset(W, phi)
        (A3.A3) x 2A3

        >>> W = clp.CoxeterGroup("A", 2, "A", 2)
        >>> clp.CoxeterCoset(W, clp.cyclestoperm((0, 3, 1, 2)))
        2(A2.A2)

        For Steinberg's triality group we keep track of which version of
        the automorphism occurs.

        >>> W = clp.CoxeterGroup("D", 4)
        >>> clp.CoxeterCoset(W, clp.cyclestoperm([0, 2, 3]))
        3D4
        >>> clp.CoxeterCoset(W, clp.cyclestoperm([0, 3, 2])
        3'D4
        
        """
        orbs = []

        for i, orb in enumerate(self.orbits):
            orbtype = self.group.cartantype[orb[0]]
            ordtwist = self.phitype[i][2].order

            # Get the type of twisting on the orbit.
            if ordtwist == 2:
                orbtxt = '\'2'
            elif ordtwist == 3:
                ind = self.phitype[i][1]
                phi = self.phitype[i][2]
                if phi[ind[0]] == ind[2]:
                    orbtxt = '\'3'
                elif phi[ind[0]] == ind[3]:
                    orbtxt = "\'3'"
            else:
                orbtxt = '\''

            orbtxt += '{}\', ('.format(orbtype[0])
            orbtxt += ', '.join(
                str([i^(self.phi**j) for i in orbtype[1]])
                for j in range(len(orb)))
            orbtxt += ')'

            orbs.append(orbtxt)

        return 'CoxeterCoset(' + ', '.join(orbs) + ')'


class ReductiveGroup:
    """
    Either a string giving the type or a pair of matrices X and Y such
    that the rows of X describe the roots in terms of the chosen basis
    and the rows of Y describe the coroots in terms of the corresponding
    dual basis.
    
    Optional keyword inputs

    trank:  an integer giving the size of the central torus
    F:      the Frobenius endomorphism. Should be a matrix describing an
            automorphism of the cocharacter group Y on the chosen basis. 
    p:      the characteristic of the field defining the reductive
            group. If not specified then this is set to 0.
    cmp:    generators for the component group G/G^0. We assume these
            representatives normalise the fixed torus T and the simple
            roots.
    iso:    a string, either 'sc' or 'ad' giving the isogeny type of the
            group

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

            cmat, ctyp, nonred = cartanmat(*args, extras=True)

            self.rankss = len(cmat)
            self.rank = self.rankss + trank
            self.cartanmat = cmat
            self.cartantype = ctyp

            for typ in self.cartantype:
                if typ[0] not in __redgrptypes__:
                    raise ValueError("A root datum must have a "
                                     "crystallographic root system!")

            self.X = np.zeros((self.rankss, self.rank), dtype='int8')
            self.Y = np.zeros((self.rankss, self.rank), dtype='int8')
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
            self.X = np.array(args[0])
            self.Y = np.array(args[1])
            self.rankss = self.X.shape[0]
            self.rank = self.X.shape[1]
            self.isogenytype = ''
            self.cartanmat = np.dot(self.Y, self.X.T)
            if not is_intcartanmat(self.cartanmat):
                raise ValueError("The product of X and Y "
                                 "is not a Cartan matrix!")
            self.cartantype = cartantotype(self.cartanmat)

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
        self.roots = roots_finintcartan(self.cartanmat)
        self.N = len(self.roots)/2
        self._dictroots = {elm : ind for ind, elm in enumerate(self.roots)}

        # Now get the coroots.
        if np.array_equal(self.cartanmat, self.cartanmat.T):
            self.coroots = list(self.roots)
        else:
            self.coroots = roots_finintcartan(self.cartanmat.T)

        # Construct the values of the canonical bilinear form for the simple
        # roots. All bilinear form values are determined from this.
        tmp = np.dot(self.coroots, self.cartanmat)
        self._Bform = np.dot(tmp.T,tmp)
        del tmp

        self.extraspecialsigns = [1 for i in xrange(self.N - self.rankss)]
        specialpairs(self)

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
        
        specialpairs(self)

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
                    zip(self.roots[x1],self.roots[x2])])
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
                    zip(self.roots[x1],self.roots[x2])])
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
                    zip(self.roots[x1],self.roots[x2])])
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
                    zip(self.roots[x1],self.roots[x2])])
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
                zip(self.roots[x1],self.roots[x2])])
            if tmp in self._dictroots:
                # We know a+b is a root so p_a,b <= 2.
                k = self._dictroots[tmp]
                if (k,x1) in self.specialpairs:
                    return 2
                else:
                    return 1
        # We now know a+b is not a root but we may still have b-a is a root.
        tmp = tuple([b-a for a,b in
            zip(self.roots[x1],self.roots[x2])])
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
                zip(self.roots[x1],self.roots[x2])])
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
                    zip(self.roots[x1],self.roots[x2])])], x1)
        elif C < 0: # q - p < 0
            for i in xrange(C,0):
                sgn *= self.sc_e(self._dictroots[tuple([b+i*a for a,b in
                    zip(self.roots[x1],self.roots[x2])])], x1)

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
        return self._dictroots[tuple([a - C*b for a, b in
            zip(self.roots[y], self.roots[x])])]

        
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


        
########################################################################
########################################################################
##                                                                    ##
##                        Cartan Matrices                             ##
##                                                                    ##

def is_intcartanmat(C):
    """
    Determines whether an integral matrix is a Cartan matrix.

    """
    n = len(C)
    vals = {0, 1, 2, 3, 4}
    for i in range(n):
        if C[i][i] != 2:
            return False
        for j in range(i+1, n):
            if C[i][j] > 0 or C[j][i] > 0:
                return False
            if C[i][j]*C[j][i] not in vals:
                return False
    return True


def cartanmat(*args, **kwargs):
    """
    If extras=True is passed as a keyword argument then the function
    returns a tuple containing: the Cartan matrix, the Cartan type of
    the matrix and the indices of the nonreduced parts of the matrix.

    Notes
    -----
    The Cartan matrices are as in [Bou02]. In particular, the Dynkin
    diagrams of the finite reduced crystallographic Cartan matrices are
    as follows:

         0   1   2        n-1          0   1            n-1
    An   o---o---o-- ... --o      Bn   o---o-- ... --o=>=o
    (n>0)                         (n>0)

                                                       o n-2
         0   1            n-1          0        n-4   /     
    Cn   o---o-- ... --o=<=o      Dn   o-- ... --o---o n-3  
    (n>0)                         (n>1)               \     
                                                       o n-1

         0   1           0   1   2   3          0   2   3   4   5
    G2   o-<-o       F4  o---o=>=o---o      E6  o---o---o---o---o
           6                                            |
                                                        o 1

         0   2   3   4   5   6           0   2   3   4   5   6   7
    E7   o---o---o---o---o---o       E8  o---o---o---o---o---o---o
                 |                               |
                 o 1                             o 1

    The Dynkin diagram of the unique non-reduced irreducible
    crystallographic Cartan matrix is given by:

                        0   1          n-1  n
                   BCn  o=>=o--- ... ---o=>=o
                   (n>2)

    The Dynkin diagrams of the affine Cartan matrices are:

          1   2       n-1  n           1 o
    A~n   o---o-- ...--o---o              \    3        n-1  n
    (n>1)  \      0       /      B~n     2 o---o-- ... --o=>=o
             -----o------        (n>2)    /
                                       0 o
    
                                       1 o                     o n-1
          0   1        n-1  n             \    3        n-3   /
    C~n   o=>=o-- ... --o=<=o    D~n     2 o---o-- ... --o---o n-2
    (n>1)                        (n>3)    /                   \ 
                                       0 o                     o n

          0   2   1           1   2   0           0   1   2   3   4
    B~2   o=>=o=<=o     G~2   o-<-o---o     F~4   o---o---o=>=o---o
                                6

          1   3   4   5   6          0   1   3   4   5   6   7
    E~6   o---o---o---o---o    E~7   o---o---o---o---o---o---o
                  |                              |
                  o 2                            o 2
                  |
                  o 0                1   3   4   5   6   7   8   0
                               E~8   o---o---o---o---o---o---o---o
                                             |
                                             o 2

    """
    rank = 0
    pairs = []
    ctyp = []
    try:
        extras = kwargs['extras']
    except KeyError:
        extras = False

    def sortinds(name, L):
        # We now apply diagram autos so that the labels of the diagram
        # respect the lexicographic ordering given by [0, 1, ...] as
        # much as possible.
        if name == "A":
            if L[0] > L[-1]:
                L.reverse()

        if name == "A~":
            m = L.index(min(L))
            L[:] = L[m:] + L[:m]

        if name == "B~":
            if L[0] > L[1]:
                L[1], L[0] = L[0], L[1]

        if name == "C~":
            if L[0] > L[-1]:
                L.reverse()

        if name == "D":
            if len(L) < 4:
                if L[0] > L[-1]:
                    L.reverse()
            elif len(L) == 4:
                ends = [L[0], L[2], L[3]]
                ends.sort()
                L[:] = [ends[0], L[1], ends[1], ends[2]]
            else:
                if L[-2] > L[-1]:
                    L[-1], L[-2] = L[-2], L[-1]

        if name == "D~":
            if len(L) == 5:
                ends = [L[0], L[1], L[3], L[4]]
                ends.sort()
                L[:] = [ends[0], ends[1], L[2], ends[2], ends[3]]
            else:
                if min(L[0], L[1]) > min(L[-1], L[-2]):
                    L.reverse()
                if L[0] > L[1]:
                    L[1], L[0] = L[0], L[1]
                if L[-2] > L[-1]:
                    L[-1], L[-2] = L[-2], L[-1]

        if name == "E" and len(L) == 6:
            if L[0] > L[-1]:
                L[5], L[0] = L[0], L[5]
                L[4], L[2] = L[2], L[4]

        if name == "E~":
            if len(L) == 6:
                branches = [[L[0], L[2]], [L[1], L[3]], [L[6], L[5]]]
                branches.sort()
                L[:] = [branches[0][0], branches[1][0], branches[0][1],
                        branches[1][1], L[4], branches[2][0], branches[2][1]]

            if len(L) == 7:
                if L[0] > L[7]:
                    L[7], L[0] = L[0], L[7]
                    L[6], L[1] = L[1], L[6]
                    L[5], L[3] = L[3], L[5]

        return None

    for typ in zip(*(iter(args),)*2):
        if typ[0].endswith("~"):
            nam = typ[0][:-1] + "t"
        else:
            nam = typ[0][:]
        try:
            rank += typ[1]
            if typ[0][-1] == "~" or typ[0] == "BC":
                rank += 1
                pairs.append((nam, typ[1]+1, False))
            else:
                pairs.append((nam, typ[1], False))
        except TypeError:
            rank += len(typ[1])
            pairs.append((nam, typ[1], True))

    cmat = np.zeros((rank, rank), dtype='int')

    # Case of the trivial group.
    if not rank:
        if extras:
            return cmat, [["A", []]], []
        else:
            return cmat

    inds = collections.deque(range(rank))
    nonred = []
    for typ in pairs:
        sl = []
        if typ[2]: # Given the indices.
            for i in typ[1]:
                try:
                    inds.remove(i)
                except ValueError:
                    raise ValueError("Indices not aligned!")
                sl.append(i)
        else:
            for i in range(typ[1]):
                sl.append(inds.popleft())

        # Apply diagram automorphisms so that the indices are
        # lexicographically ordered wherever possible.
        sortinds(typ[0], sl)

        ctyp.append([typ[0], sl])

        if sl:
            tmp = np.ix_(sl, sl)
            cmat[tmp] = sys.modules[
                    'charliepy.data.typ1{}'.format(typ[0])].cartanmat(len(sl))
            if typ[0] == 'BC':
                nonred.append([sl[0], sl[-1]])

    if inds:
        raise ValueError("Indices not aligned!")

    if extras:
        return cmat, ctyp, nonred 
    else:
        return cmat


def cartantotype(C):
    """The input is an integral Cartan matrix of finite type and the output is
    the type of the Cartan matrix. It is assumed that the input is genuinely a
    Cartan matrix.

    If the optional flags finite or crystal are set to True then the function
    will return False if the matrix is not finite or crystalographic
    respectively.
    """
    n = len(C)
    if not n:
        return [["A",[]]]

    # nonzero[i] gives all j such that C[i][j] != 0.
    nonzero = [np.nonzero(C[i])[0].tolist() for i in range(n)]
    blocks = []
    tmp = list(range(n))
    while tmp:
        orb = [tmp[0]]
        for i in orb:
            for j in tmp:
                if i in nonzero[j] and j in nonzero[i] and not j in orb:
                    orb.append(j)
        for i in orb: 
            tmp.remove(i)
        orb.sort()
        blocks.append(orb)

    typ = []
    for b in blocks:
        tnonzero = [[b.index(j) for j in nonzero[i]] for i in b]
        ntyp = cartantotype_indec(C[np.ix_(b,b)], tnonzero)
        ntyp[1] = [b[i] for i in ntyp[1]]
        typ.append(ntyp)

    return typ


def cartantotype_indec(C, nonzero):
    """Determines the Cartan type when the matrix is indecomposable. In
    particular the corresponding Dynkin diagram of C is connected.
    """
    n = len(C)
    if n == 1: # A1
        return ["A", [0]]

    if n == 2: # Rank 2.
        #if C[0][1] == -4 and C[1][0] == -1: return ["BC", [0,1]]
        #if C[1][0] == -4 and C[0][1] == -1: return ["BC", [1,0]]
        if C[1][0] == -3 and C[0][1] == -1: return ["G", [0,1]]
        if C[0][1] == -3 and C[1][0] == -1: return ["G", [1,0]]
        if C[0][1] == -2 and C[1][0] == -1: return ["B", [0,1]]
        if C[1][0] == -2 and C[0][1] == -1: return ["B", [1,0]]
        if C[0][1] == -1 and C[1][0] == -1: return ["A", [0,1]]
        if C[0][1] == -2 and C[1][0] == -2: return ["A~", [0,1]]
        return ["U", [0,1]]
    
    # Must now have rank at least 3.
    ends = []
    bpoints = []
    sbpoints = []

    uniq = np.unique(C)
    sbond = np.array_equal(uniq, [-1, 0, 2]) # single bond
    dbond = np.array_equal(uniq, [-2, -1, 0, 2]) # double bond
    tbond = np.array_equal(uniq, [-3, -1, 0, 2]) # triple bond

    for i in range(n):
        if len(nonzero[i]) == 2: # End point
            ends.append(i)
        elif len(nonzero[i]) == 4: # Branch point
            bpoints.append(i)
        elif len(nonzero[i]) == 5: # D4~ branch point
            sbpoints.append(i)
        elif len(nonzero[i]) != 3: # Don't know what it is.
            return ["U", list(range(n))]

    # All vertices now have degree 1, 2, 3 or 4.
    if len(sbpoints) == 1 and sbond: # D4~
        if len(ends) == 4 and n == 5:
            typinds = sorted(ends)
            typinds.insert(2, sbpoints[0])
            return ["D~", typinds]
    if sbpoints:
        return ["U", list(range(n))]

    # All vertices now have degree 1, 2, or 3.
    if len(ends) == 4 and len(bpoints) == 2 and sbond:
        # Could be Dn~. First get the end nodes which are attached to
        # the 2 branch points.
        grpd = []
        for p in bpoints:
            grpd.append([j for j in nonzero[p] if j in ends])
        if (len(grpd[0]), len(grpd[1])) != (2, 2):
            return ["U", list(range(n))]
        # Must now be Dn~
        a, b = bpoints[0], bpoints[1]
        grpd[0].sort()
        grpd[1].sort()
        typinds = grpd[0][:]

        # Now determine the main branch a - ... - b.
        typinds.append(a)
        i = [j for j in nonzero[a] if j != a and j not in grpd[0]][0]
        while i != b:
            typinds.append(i)
            i = [j for j in nonzero[i] if j not in typinds][0]
        typinds.append(b)
        return ["D~", typinds + grpd[1]]
    
    if len(ends) == 3 and len(bpoints) == 1:
        # Note, we must now have precisely 3 straight line branches.
        # Assume a is the unique branch point. If e is the end of a
        # branch then we describe the branch of e as a list [b0, b1,
        # ..., bk] such that the branch is given by
        #
        #       e = b0 - b1 - ... - bk - a
        #
        a = bpoints[0]
        branches = []
        tinds = set([a])
        for e in ends:
            tinds.add(e)
            branch = [e]
            i = e
            while a not in nonzero[i]:
                i = [j for j in nonzero[i] if j not in tinds][0]
                tinds.add(i)
                branch.append(i)
            branches.append(branch)

        # Sort the branches by their lengths so we can classify the
        # diagram easily.
        branches.sort(key=(lambda x: len(x)))
        blens = [len(branch) for branch in branches]

        if sbond: # En, En~, Dn.
            if blens == [1, 2, 2]: # E6
                if branches[1][0] > branches[2][0]:
                    branches[2], branches[1] = branches[1], branches[2]
                branches[2].reverse()
                branches[1].insert(1, branches[0][0])
                typinds = branches[1] + [a] + branches[2]
                return ["E", typinds]
            if blens == [1, 2, 3] or blens == [1, 2, 4]: # E7, E8
                branches[2].reverse()
                branches[1].insert(1, branches[0][0])
                return ["E", branches[1] + [a] + branches[2]]
            if blens == [2, 2, 2]: # E6~
                branches.sort(key=(lambda x: x[0]))
                typinds = [None]*n
                typinds[:3:2] = branches[0]
                typinds[1:4:2] = branches[1]
                typinds[4] = a
                typinds[5:] = branches[2][::-1]
                return ["E~", typinds]
            if blens == [1, 3, 3]: # E7~
                if branches[1][0] > branches[2][0]:
                    branches[2], branches[1] = branches[1], branches[2]
                typinds = [None]*n
                typinds[:2] = branches[1][:2]
                typinds[2] = branches[0][0]
                typinds[3] = branches[1][2]
                typinds[4] = a
                typinds[5:] = branches[2][::-1]
                return ["E~", typinds]
            if blens == [1, 2, 5]: # E8~
                typinds = [None]*n
                typinds[0] = branches[2][0]
                typinds[1] = branches[1][0]
                typinds[2] = branches[0][0]
                typinds[3] = branches[1][1]
                typinds[4] = a
                typinds[5:] = branches[2][-1:0:-1]
                return ["E~", typinds]
            if blens == [1, 1, 1]: # D4
                typinds = branches[0] + branches[1] + branches[2]
                typinds.sort()
                typinds.insert(1, a)
                return ["D", typinds]
            if blens[:2] == [1, 1]: # Dn
                typinds = branches[0] + branches[1]
                typinds.sort()
                return ["D", branches[2] + [a] + typinds]
            return ["U", list(range(n))]

        if dbond and blens == [1, 1, 1]: # B3~
            dbend = []
            for e in ends:
                if -2 in C[:, e]:
                    dbend.append(e)
            if len(dbend) != 1 or -2 in C[:, a]:
                return ["U", list(range(n))]
            e = dbend[0]
            ends.remove(e)
            ends.sort()
            typinds = ends + [a, e]
            return ["B~", typinds]

        if dbond and blens[:2] == [1, 1]: # Bn~
            if np.where(C == -2) == (branches[2][1], branches[2][0]):
                if branches[0][0] > branches[1][0]:
                    typinds = [branches[1][0], branches[0][0], a]
                else:
                    typinds = [branches[0][0], branches[1][0], a]
                return ["B~", typinds + branches[2][::-1]]
            return ["U", list(range(n))]

    if bpoints:
        return ["U", list(range(n))]

    # All vertices now have degree 1 or 2.
    if len(ends) == 2: # An, Bn, Cn, BCn, F4, B2~, Cn~, F4~, G2~.
        # Note, we must now have a straight line.
        ends.sort()
        e = ends[0] # line starts as lex smallest end node.
        line = [e]
        while e != ends[1]:
            e = [j for j in nonzero[e] if j not in line][0]
            line.append(e)

        if sbond: # An.
            return ["A", line]
        if dbond:
            inds = list(zip(*np.where(C == -2)))
            if len(inds) == 1: # Bn, Cn, F4, F4~.
                bond = inds[0]
                if n == 4 and bond == (line[1], line[2]):
                    return ["F", line]
                if n == 4 and bond == (line[2], line[1]):
                    return ["F", line[::-1]]
                if n == 5 and bond == (line[2], line[3]):
                    return ["F~", line]
                if n == 5 and bond == (line[2], line[1]):
                    return ["F~", line[::-1]]
                if bond == (line[n-2], line[n-1]):
                    return ["B", line]
                if bond == (line[1], line[0]):
                    return ["B", line[::-1]]
                if bond == (line[n-1], line[n-2]):
                    return ["C", line]
                if bond == (line[0], line[1]):
                    return ["C", line[::-1]]
            if len(inds) == 2: # B2~, BC2, Cn~.
                S = set(inds)
                if n == 3 and S == set([(line[0],line[1]), (line[2],line[1])]):
                        return ["B~", [line[0], line[2], line[1]]]
                if n > 3 and S == set([(line[0], line[1]),
                                       (line[n-1], line[n-2])]):
                    return ["C~", line]
                if S == set([(line[0], line[1]), (line[n-2], line[n-1])]):
                    return ["BC", line]
                if S == set([(line[1], line[0]), (line[n-1], line[n-2])]):
                    return ["BC", line[::-1]]
        if tbond: # G2~.
            inds = list(zip(*np.where(C == -3)))
            bond = inds[0]
            if n == 3 and len(inds) == 1:
                # Note there is an unusual labelling for G~2.
                if bond == (line[1], line[2]):
                    return ["G~", [line[0], line[2], line[1]]]
                if bond == (line[1], line[0]):
                    return ["G~", [line[2], line[0], line[1]]]

    if ends:
        return ["U", list(range(n))]

    # All vertices now have degree 2. Must be a circuit.
    i = 0
    typinds = [0]
    while len(typinds) < n:
        i = [j for j in nonzero[i] if j not in typinds][0]
        typinds.append(i)
    return ["A~", typinds]


########################################################################
########################################################################
##                                                                    ##
##                              Roots                                 ##
##                                                                    ##

def roots_finintcartan(C, nonred=[]):
    """
    Returns the roots of an integral Cartan matrix of finite type.

    Recall that every root is a Z-linear combination of the simple
    roots. Each root is thus represented by a tuple where the ith entry
    of the tuple gives the coefficient of the simple root alpha_i, with
    the indexing of the simple roots determined by the Cartan matrix.
    """
    # We use the following inductive algorithm to determine the roots.
    # This uses the assertions in Ch. VI, S1, no.  3, Proposition 9 of
    # [Bou02].
    #
    # Consider two roots alpha and beta such that alpha != beta.  The
    # set of integers i such that beta + i.alpha is a root is an
    # interval [-q, p] in Z containing 0. We call [-q, p] the
    # alpha-interval through beta. Clearly we have p-q is contained in
    # [-q, p] and moreover we have
    #
    #     p - q = -<beta, alpha*>.
    #
    # In particular, beta - <beta, alpha*>.alpha is a root. Now if beta
    # + alpha is a root then the alpha-interval through beta + alpha is
    # simply [-q-1, p-1]. As p - q - 1 is contained in this interval we
    # therefore have
    #
    #     beta - (<beta, alpha*> + 1).alpha
    #
    # is a root. Conversely if this is a root then we must have
    #
    #     -q <= p - q - 1  -->  1 <= p
    #
    # so beta + alpha is a root.
    #
    # This observation gives us an inductive approach to determining the
    # positive roots from the Cartan matrix. We proceed by induction on
    # the height of the root. The height 1 and 2 roots are easy to write
    # down. We then use the above characterisation to determine whether
    # adding a simple root gives a new positive root.
    #
    # When dealing with non-reduced root systems we always ensure that
    # each root has a unique expression by never using simple roots of
    # the form 2alpha_i.
    n = C.shape[0]
    droots = list(zip(*nonred))
    if not droots:
        droots = [[], []]
    inds = [i for i in range(n) if i not in droots[1]]
    I = np.identity(n, dtype='int8')

    # We start with the simple roots of height 1.
    roots = [[I[i].tolist() for i in inds]]

    # Now the roots of height 2.
    roots.append(
        [(I[i] + I[j]).tolist() for i in inds for j in range(i, n)
         if (j == i and i in droots[0])
         or (j != i and C[i][j] != 0 and j not in droots[1])]
    )

    # Now proceed by induction to obtain all roots.
    i = 1
    while roots[i]:
        roots.append([])
        prod = np.dot(roots[i], C).tolist()

        for j in range(len(roots[i])):
            r = roots[i][j]
            for k in inds:
                p = prod[j][k]

                # If p is negative then alpha+alpha_k is definitely a
                # root by Ch. VI, S1, no. 3, Corollary to Theorem 1 of
                # [Bou02].
                if p < 0:
                    r[k] += 1
                    if r not in roots[i+1]:
                        roots[i+1].append(list(r))
                    r[k] -= 1

                # If i <=p or r[k] <= p then we cannot have alpha -
                # (p+1)alpha_k is a root. After this all cheap tests are
                # done.
                elif i > p and r[k] > p:
                    r[k] -= p + 1
                    if r in roots[i-p-1]:
                        r[k] += p + 2
                        if r not in roots[i+1]:
                            roots[i+1].append(list(r))
                        r[k] -= 1
                    else:
                        r[k] += p + 1
        i += 1

    # We now just flatten the roots to a single list and add the
    # negative roots.  Note we have to create the tuples here as they
    # are immutable so we wouldn't be able to change their values in the
    # above loop.
    roots = [tuple(x) for x in itertools.chain.from_iterable(roots)]
    return tuple(roots + [tuple(-i for i in r) for r in roots])


def specialpairs(R):
    # For each positive root alpha we find the minimal simple root alpha_j
    # and positive root beta such that alpha = alpha_j + beta. Here R should be
    # a list consisting of the positive roots and n should be the number of
    # simple roots.
    #
    # Note that the ith entry of extspec is the extraspecial pair of the root
    # with position i+n.

    n = len(R.cartanmat)
    hts = [sum(r) for r in R.roots[:R.N]]
    mod = np.einsum("ij,jk,ik->i", R.roots[:R.N], R._Bform,
            R.roots[:R.N]).tolist()
    D = R._dictroots
    extspec = [None]*(R.N-n)
    pairs = dict()

    for i in range(n, R.N):
        first = True
        for j in range(i):
            tmp = tuple([a-b for a, b in zip(R.roots[i],R.roots[j])])
            if tmp in D:
                m = D[tmp]
                if j >= m:
                    if hts[j] > hts[m]:
                        break
                    else:
                        continue
                
                # Get the p and the q values for the root string. Here we use
                # the clasification of such root strings given in Lemma 4.1.1 of
                # [Car92].
                x, y = mod[m], mod[i]
                if x == y:
                    if (tuple([a-b for a,b in
                               zip(R.roots[m],R.roots[j])]) in D):
                        p, q = 1, 2
                    else:
                        p, q = 0, 1
                elif x == 2*y:
                    p, q = 0, 2
                elif 2*x == y:
                    p, q = 1, 1
                elif x == 3*y:
                    p, q = 0, 3
                elif 3*x == y:
                    p, q = 2, 1

                # Get the structure constant for the special pair.
                if first:
                    N = R.extraspecialsigns[i-n]*(p+1)
                    extspec[i-n] = [(j, m), p, q, N]
                    pairs.update({(j, m) : [i, p, q, N]})
                    first = False
                    continue
                else:
                    x, y = extspec[i-n][0][0], extspec[i-n][0][1]
                    tmp = tuple([a-b for a,b in
                        zip(R.roots[m],R.roots[x])])
                    if tmp in D:
                        k = D[tmp]
                        sgn = 1
                        if x > k:
                            p1 = (k,x)
                            sgn *= -1
                        else:
                            p1 = (x,k)
                        if j > k:
                            p2 = (k,j)
                            sgn *= -1
                        else:
                            p2 = (j,k)
                        t1 = sgn*mod[k]*mod[j]*pairs[p1][2]*pairs[p2][2]
                    else:
                        t1 = 0
                    
                    tmp = tuple([a-b for a,b in
                        zip(R.roots[j], R.roots[x])])
                    if tmp in D:
                        k = D[tmp]
                        sgn = 1
                        if x > k:
                            p1 = (k, x)
                            sgn *= -1
                        else:
                            p1 = (x, k)
                        if m > k:
                            p2 = (k, m)
                            sgn *= -1
                        else:
                            p2 = (m, k)
                        t2 = sgn*mod[k]*mod[m]*pairs[p1][2]*pairs[p2][2]
                    else:
                        t2 = 0
                    N = cmp(t1, t2)*(p+1)
                    N = p+1
                    pairs.update({(j, m) : [i, p, q, N]})
    
    R.specialpairs = pairs


def rootp(G, x, y):
    a = G.roots[x]
    b = G.roots[y]

    p = 0
    tmp = tuple([d - e for d, e in zip(b, a)])
    while tmp in G._dictroots:
        p += 1
        tmp = tuple([d - e for d, e in zip(tmp, a)])
    return p


def rootq(G, x, y):
    a = G.roots[x]
    b = G.roots[y]

    q = 0
    tmp = tuple([d + e for d, e in zip(b, a)])
    while tmp in G._dictroots:
        q += 1
        tmp = tuple([d + e for d, e in zip(tmp, a)])
    return q

