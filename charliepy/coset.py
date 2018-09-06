from . import permutat
from . import utils
import weakref

import numpy as np

__all__ = ['CoxeterCoset']


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

        if not isinstance(phi, permutat.Perm):
            raise TypeError("phi must be a permutation.")

        # Get phi as a full permutation on the roots.
        posroots = W.roots[:W.N]
        try:
            pperm = [posroots.index(phi.permute(x)) for x in posroots]
        except IndexError:
            raise ValueError("phi must preserve the positive roots")
        pperm.extend([x + W.N for x in pperm])
        phi = permutat.Perm(pperm)
        self.phi = phi
        
        # Determine the orbits of phi acting cyclically on isomorphic factors.
        typsets = [set(x[1]) for x in W.cartantype]
        indperm = phi.onsets(typsets)
        self.orbits = indperm.cycles()

        # Determine the twist induced by phi on each orbit and then construct
        # the type of the coxetercoset.
        phitype = [None]*len(self.orbits)

        # If x is an orbit of phi of length n then phi**n stabilises each
        # component of x. We choose the first entry in the orbit as a
        # representative for the automorphism induced by phi**n on each
        # component in that orbit.
        for i, orb in enumerate(self.orbits):
            rep = orb[0]
            typ, inds = W.cartantype[rep]
            phi_orb = permutat.restrictedperm(phi**len(orb), inds)
            phitype[i] = [typ, inds, phi_orb]

        self.phitype = phitype

        # Set the order to be the order of the underlying group.
        self.size = self.group.size

        # This keeps track of the fusion of conjugacy classes into groups with
        # embeddings. This will be computed for the group when the conjugacy
        # classes are computed.
        self.fusion = weakref.WeakKeyDictionary()

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
                if ind[0]^phi == ind[2]:
                    orbtxt = '\'3'
                elif ind[0]^phi == ind[3]:
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

    def diagram(self, parent=None):
        """
        Prints the diagram of the Coxeter Coset. If the coset, say H,
        has an embedding into the coset W then setting parent=W will
        print the diagram with the labels provided by the roots in W.

        Examples
        --------
        >>> W = clp.CoxeterGroup("A", 5, "A", 5)
        >>> phi = clp.cyclestoperm((0,9,4,5), (1,8,3,6))
        >>> Wphi = clp.CoxeterCoset(W, phi)
        >>> Wphi.diagram()
        phi acts on the following component as (0, 4)(1, 3)
        A5 :  0 - 1 - 2 - 3 - 4
        A5 :  9 - 8 - 7 - 6 - 5

        """
        for i, orb in enumerate(self.orbits):
            typ, inds = self.group.cartantype[orb[0]]
            phi_orb = self.phitype[i][2]

            print("phi acts on the following component as {}".format(phi_orb))

            for j in range(len(orb)):
                utils.getmodule(typ).diagram([i^(self.phi**j) for i in inds])









