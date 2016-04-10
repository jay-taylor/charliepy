from . import core

import sys
import itertools
import functools
import operator
import collections

__all__ = ['ConjugacyClass',
           'conjugacyclasses',
           'minlengthrep',
           'cyclicshiftclass']


########################################################################
########################################################################
##                                                                    ##
##                        Conjugacy Classes                           ##
##                                                                    ##

class ConjugacyClass(collections.namedtuple("ConjugacyClass",
            ["parent", "name", "centsize", "rep"])):
    # Setting slots to an empty tuple prevents a dict for the class
    # being created.
    __slots__ = () 
    def __repr__(self):
        return "ConjugacyClass( {}, {} )".format(self.parent,
                                                 self.name)

    @property
    def length(self):
        return self.parent.size//self.centsize


def conjugacyclasses(W):
    """
    Computes the conjugacy classes of a Coxeter group.

    Parameters
    ----------
    W : :class:`core.CoxeterGroup` or :class:`core.CoxeterCoset`

    Returns
    -------
    out : list
        Each entry of the list is an instance of :class:`ConjugacyClass'
        which contains relevant information concerning the class.

    Notes
    -----
        * For the conventions concerning the raw data for each
          irreducible type see the databook.
        * For a coset this computes the phi-conjugacy classes of the
          coset.
        * For a group the representatives are *very good* in the
          sense of [GM97]_.
        * Once computed the classes are stored as the atribute
          ``conjugacyclasses`` of W.

    Examples
    --------
    >>> conjugacyclasses(coxeter("G",2))
    {'classnames'   : [('A_0',), ('~A_1',), ('A_1',), ('G_2',),
                       ('A_2',), ('A_1+~A_1',)],
     'classlengths' : [1, 3, 3, 2, 2, 1],
     'reps'         : [[], [1], [0], [0, 1], [0, 1, 0, 1],
                       [0, 1, 0, 1, 0, 1]]}

    >>> W = coxeter([[2,0,-1,0],[0,2,0,-2],[-1,0,2,0],[0,-1,0,2]])
    >>> W.cartantype
    [['A', [0, 2]], ['C', [3, 1]]]
    >>> conjugacyclasses(W)['reps']
    [[], [3], [3, 1, 3, 1], [1], [3, 1], [0], [0, 3],
     [0, 3, 1, 3, 1], [0, 1], [0, 3, 1], [0, 2], [0, 2, 3],
     [0, 2, 3, 1, 3, 1], [0, 2, 1], [0, 2, 3, 1]]

    """
    try:
        return W.conjugacyclasses
    except AttributeError:
        pass

    # First get the data for each irreducible component. This dynamically
    # calls the correct module for the type. Note that, here, the function
    # returns centraliser orders not class lengths.
    if isinstance(W, core.CoxeterGroup):
        if len(W.cartantype) == 1:
            gen = sys.modules[
                'charliepy.data.typ1{}'.format(W.cartantype[0][0])
                ].conjclasses(W.cartantype[0][1])

            W.conjugacyclasses = [ConjugacyClass(W, elm[0], elm[1], elm[2])
                                  for elm in gen]

        else:
            W.conjugacyclasses = []

            gen = itertools.product(*(sys.modules[
                'charliepy.data.typ1{}'.format(typ[0])].conjclasses(typ[1])
                 for typ in W.cartantype))

            for elm in gen:
                name, clen, rep = zip(*elm)
                W.conjugacyclasses.append(
                    ConjugacyClass(
                        W, ",".join(name),
                        functools.reduce(operator.mul, clen, 1),
                        list(itertools.chain.from_iterable(rep))
                    )
                )

    elif isinstance(W, core.CoxeterCoset):
        if len(W.phitype) == 1:
            gen = sys.modules[
                'charliepy.data.typ{}{}'.format(
                    W.phitype[0][2].order, W.phitype[0][0])
            ].conjclasses(W.phitype[0][1], phi=W.phitype[0][2])

            W.conjugacyclasses = [ConjugacyClass(W, *elm) for elm in gen]

        else:
            W.conjugacyclasses = []

            gen = itertools.product(*(
                sys.modules[
                    'charliepy.data.typ{}{}'.format(typ[2].order, typ[0])
                ].conjclasses(typ[1], phi=typ[2]) for typ in W.phitype)
            )

            for elm in gen:
                name, clen, rep = zip(*elm)
                W.conjugacyclasses.append(
                    ConjugacyClass(
                        W, ",".join(name),
                        functools.reduce(operator.mul, clen, 1),
                        list(itertools.chain.from_iterable(rep))
                    )
                )

    else:
        raise NotImplemented("W must be a CoxeterGroup or CoxeterCoset.")

    return W.conjugacyclasses

def cyclicshiftclass(W, w):
    """
    Tests if an element of W, given as a permutation on the roots, has
    minimal length in its conjugacy class. If this is the case () is
    returned and if not then a tuple (x, i) is returned where x is an
    element of W and i is the index of a simple reflection such that
    l(x) < l(s_ixs_i) and s_ixs_i lies in the same cyclic shift class as
    w.

    """
    # Here we follow Algorithm G from pg. 80 of [GP00] to construct the
    # cyclic shift class of the element w. Note we keep the result as
    # both a set and a list, this means we have the advantage of fast
    # lookups but also predictable return results.
    Y = {w}
    X = set()
    while Y:
        y = next(iter(Y))
        l = W.length(y, 'p')
        X.add(y)
        Y.remove(y)
        for s in W.permgens:
            z = y^s
            if W.length(z, 'p') == l and not z in X:
                Y.add(z)
    return X

def minlengthrep(W, w):
    """
    This returns the elements of minimal length in the conjugacy class of w as a
    permutation on the roots.

    """
    # Here we follow Algorithm H from pg. 83 of [GP00] to construct the
    # cyclic shift class of the element w.
    Y = {w}
    X = set()
    while Y:
        y = next(iter(Y))
        l = W.length(y, 'p')
        X.add(y)
        Y.remove(y)
        for s in W.permgens:
            z = y^s
            lz = W.length(z, 'p')
            if lz < l:
                Y = {z}
                X = set()
                break
            elif lz == l and not z in X:
                Y.add(z)
    return next(iter(X))#min(X, key=(lambda x: W.convert(x, 'pw')))









