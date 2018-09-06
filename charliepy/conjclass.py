from . import utils
from . import coxeter
from . import coset

import sys
import itertools
import functools
import operator
import collections

__all__ = ['ConjugacyClass',
           'conjugacyclasses',
           'minlengthelms',
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
    W : :class:`coxeter.CoxeterGroup` or :class:`coset.CoxeterCoset`

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
    if isinstance(W, coxeter.CoxeterGroup):
        if len(W.cartantype) == 1:
            typ = W.cartantype[0]
            gen = utils.getmodule(typ[0]).conjclasses(typ[1])

            W.conjugacyclasses = [
                ConjugacyClass(W, *elm) for elm in gen]

        else:
            W.conjugacyclasses = []

            gen = itertools.product(*(utils.getmodule(typ).conjclasses(inds)
                    for (typ, inds) in W.cartantype))

            for elm in gen:
                name, clen, rep = zip(*elm)
                W.conjugacyclasses.append(
                    ConjugacyClass(
                        W, ",".join(name),
                        functools.reduce(operator.mul, clen, 1),
                        list(itertools.chain.from_iterable(rep))
                    )
                )

    elif isinstance(W, coset.CoxeterCoset):
        if len(W.phitype) == 1:
            typ, inds, phi = W.phitype[0]
            gen = utils.getmodule(typ, phi.order).conjclasses(inds, phi=phi)

            W.conjugacyclasses = [ConjugacyClass(W, *elm) for elm in gen]

        else:
            W.conjugacyclasses = []

            gen = itertools.product(*(
                utils.getmodule(typ, phi.order).conjclasses(inds, phi=phi)
                for (typ, inds, phi) in W.phitype)
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

    # Describe the fusion of conjugacy classes into itself.
    W.fusion[W] = [i for i in range(len(W.conjugacyclasses))]

    return W.conjugacyclasses

def identifyclass(W, w):
    """
    Given an element of W this returns the name of the class containing
    w as a string.

    """
    if isinstance(W, coxeter.CoxeterGroup):
        # First make sure the element is a word.
        w = W.convert(w, 'w')

        # If the group has more than one irreducible component then we
        # have to seperate out the word into its irreducible
        # constituents.
        if len(W.cartantype) == 1:
            typ, inds = W.cartantype[0]

            return utils.getmodule(typ).wordtoclass(
                    len(inds), [W._standardinds[i] for i in w])
        else:
            inds = [I for _, I in W.cartantype]
            words = [[] for _ in W.cartantype]
            labs = []

            for a in w:
                ind = next(i for i, I in enumerate(inds) if a in I)
                words[ind].append(a)

            for word, (typ, I) in zip(words, W.cartantype):
                labs.append(utils.getmodule(typ).wordtoclass(
                    len(I), [W._standardinds[i] for i in word]))

            return ",".join(labs)

def fusionclasses(H, W):
    """
    This computes the fusion of the conjugacy classes of H into W. For this to
    work we must have H.embeddings[W] is defined.

    """
    # See if we already computed this before.
    try:
        return H.fusion[W]
    except KeyError:
        pass

    # Make sure H embeds into W.
    try:
        embed = H.embeddings[W]
    except KeyError:
        raise KeyError("No embedding of H into W.")

    # Make sure we have the classes of H and W.
    Wclasses = conjugacyclasses(W)
    Hclasses = conjugacyclasses(H)

    # We construct a dict matching class names of W to the position of the
    # class in W.conjugacyclasses. This will make the lookups below much
    # more efficient.
    pos_Wclass = dict(zip((c.name for c in W.conjugacyclasses),
        range(len(W.conjugacyclasses))))

    # This will store the matching of the classes.
    res = [None]*len(Hclasses)

    # For each conjugacy class of H we obtain a word in the standard
    # generators of W representing that class. We then use the corresponding
    # algorithm from the data package to determine the class.
    if len(W.cartantype) == 1:
        typ = W.cartantype[0][0]
        n = len(W.cartantype[0][1])

        # The generators of H are reflections in W. Here we obtain these
        # reflections as words in the standard generators of W. This allows us
        # to easily build a word in the standard generators from the
        # representatives of the H-classes.
        refword = dict(zip(embed,
               ([W._standardinds[j] for j in W.convert(W.reflection(x), 'w')]
                   for x in embed)))

        for i, cls in enumerate(Hclasses):
            word = []
            for x in cls.rep:
                word.extend(refword[embed[x]])
            res[i] = pos_Wclass[utils.getmodule(typ).wordtoclass(n, word)]

    else:
        inds = [set(I) for _, I in W.cartantype]

        for i, cls in enumerate(Hclasses):
            words = [[] for _ in W.cartantype]
            labs = [None]*len(W.cartantype)

            # The generators of H are reflections in W. Here we obtain these
            # reflections as words in the generators of W. Note, we cannot
            # write them as words in the standard generators because we have to
            # test which irreducible component each letter in the word
            # corresponds to.
            refword = dict(zip(embed,
                   (W.convert(W.reflection(x), 'w') for x in embed)))

            # We break the word up as a product of words in the irreducible
            # components of W. Note, as in the irreducible case, we write the
            # word in terms of the standard generators of W.
            for x in cls.rep:
                rword = refword[embed[x]]
                for a in rword:
                    ind = next(j for j, I in enumerate(inds) if a in I)
                    words[ind].append(W._standardinds[a])

            for j, (word, (typ, rank)) in enumerate(zip(words, W.cartantype)):
                labs[j] = utils.getmodule(typ).wordtoclass(len(rank), word)
            
            res[i] = pos_Wclass[",".join(labs)]

    # Store the fusion of the classes.
    H.fusion[W] = res

    return res

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

def minlengthelms(W, w):
    """
    This returns the elements of minimal length in the conjugacy
    class of w as a permutation on the roots.

    """
    # Here we follow Algorithm H from pg. 83 of [GP00].
    if isinstance(W, coset.CoxeterCoset):
        phi = W.phi
        gens = W.group.permgens
        phigens = [(i^phi, gens[i^phi]) for i in range(W.group.rank)]
        N = W.group.N
        w = W.group.convert(w, 'p')
    else:
        gens = W.permgens
        phigens = list(enumerate(gens))
        N = W.N
        w = W.convert(w, 'p')

    Y = {w}
    X = set()
    while Y:
        y = Y.pop()
        X.add(y)
        for (i, s), (j, t) in zip(enumerate(gens), phigens):
            z = y*t
            # l(syt) < l(y) iff alpha_s.z < 0 and alpha_t.y**-1 < 0.
            if i^z >= N:
                v = s*z
                if j/y >= N:
                    Y = {v}
                    X = set()
                    break
                elif not v in X:
                    Y.add(v)
            # l(syt) == l(y) if alpha_s.z > 0 and alpha_t.y**-1 < 0.
            elif j/y >= N:
                v = s*z
                if not v in X:
                    Y.add(v)
    return sorted(X)









