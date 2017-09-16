#############################################################################
#  This module is designed to be used with CharLiePy. It implements permutations
#  using numpy arrays. This allows us fast access for permutation arithmetic
import numpy as np
import operator
import itertools

__all__ = [
    'Perm',
    'cycle',
    'cycledecomp',
    'cycles',
    'cyclestoperm',
    'cycletype',
    'orderperm',
    'permonsets',
    'restrictedperm'
]

bitsize = "16"
intsize = "uint" + bitsize

class Perm:
    """
    A permutation on the integers `\{0,\dots,n-1\}`.

    Parameters
    ----------
    perm : list
        The entry of ``perm[i]`` should give the image of ``i`` under
        the permutation.

    Attributes
    ----------
    perm : numpy array
        Stores the input perm as a numpy array.
    deg : int
        If the permutation is on `\{0,\dots,n-1\}` then this is `n`.

    Notes
    -----
    To create more efficient memory management, it is assumed that the
    degree of a permutation does not exceed 2\ :sup:`16`. This should be
    more than sufficient for the purposes we treat in CharLiePy. If for some
    reason one needs to consider permutations of larger degree then one
    can change the global attribute ``bitsize``. For example, setting::
    
        >>> clp.permutat.bitsize = "32"

    will allow permutations of degree up to 2\ :sup:`32`.

    Examples
    --------
    >>> p = clp.permutat.Perm([6, 5, 4, 3, 2, 1, 0])
    >>> p
    (0, 6)(1, 5)(2, 4)
    >>> p.deg
    7
    >>> p.perm
    array([6, 5, 4, 3, 2, 1, 0], dtype=uint16)

    """
    def __init__(self, perm):
        # Note that this always makes a new permutation by default, this
        # matches behaviour of Python's builtins.
        self.perm = np.array(perm, dtype = intsize, copy=True)
        self.deg = len(perm)

    def __repr__(self):
        """
        Returns the cycle decomposition of the permutation as a string.
        
        """
        return 'Perm(' + ''.join(str(t) for t in cycledecomp(self)) + ')'

    def __eq__(self, other):
        """
        Returns true if the images of the permutations are the same and
        false otherwise.

        """
        # Checks first that the perms agree on the parts which are defined
        # and then makes sure that the remaining parts of the longer perm
        # are just the identity.
        if self.deg <= other.deg:
            if not np.array_equal(self.perm, other.perm[:self.deg]):
                return False
            if not np.array_equal(other.perm[self.deg],
                                  range(self.deg, other.deg)):
                return False
        else:
            if not np.array_equal(self.perm[:other.perm], other.perm):
                return False
            if not np.array_equal(self.perm[other.perm:],
                                  range(other.deg, self.deg)):
                return False

        return True

    def __mul__(p, q):
        """
        Returns the product of two permutations.

        The multiplication of permutations is defined so that (p*q)(i) =
        q(p(i)) for all i. In other words, the composition of the
        functions is read from left to right. This is compatible with our
        convention that the symmetric group acts on the right of the set
        `\{1,\dots,n\}`.
        
        Examples
        --------
        >>> p = clp.permutat.cyclestoperm([0, 1])
        >>> q = clp.permutat.cyclestoperm([1, 2])
        >>> p*q
        Perm((0, 2, 1))
        >>> q*p
        Perm((0, 1, 2))

        """
        # If the inner permutation has smaller degree than the outer permutation
        # then this is easy. Simply do i --> q(p(i)) if i < p.deg, otherwise
        # set i --> q(i). If the degree of the inner permutation is bigger then
        # it is essentially the same but we must first extend the list by the
        # identity.
        try:
            if p.deg == q.deg:
                new = q.perm[p.perm]
            elif p.deg < q.deg:
                new = np.append(q.perm[p.perm], q.perm[p.deg:])
            else:
                new = np.append(q.perm, range(q.deg, p.deg))
                new = new[p.perm]

            return Perm(new)

        except:
            return NotImplemented("Multiplication is only defined for"
                                  "permutations.")

    def __pow__(self, other):
        """
        The power of a permutation is only defined for integers. In this
        case it returns the i-fold composition of the permutation with
        itself.
        
        Examples
        --------
        >>> w = clp.permutat.cyclestoperm([4, 2, 3, 1])
        >>> w
        (1, 4, 2, 3)
        >>> w**2
        (1, 2)(3, 4)
        >>> w**3
        (1, 3, 2, 4)
        >>> w**4
        ()
        >>> w**-1
        (1, 3, 2, 4)
        >>> w**-2
        (1, 2)(3, 4)

        """
        try:
            if other == 0:
                return Perm(range(self.deg))

            elif other > 0:
                q = self.perm
                for i in range(other-1):
                    q = self.perm[q]
                return Perm(q)

            elif other < 0:
                # First have to construct the inverse of self. This makes
                # inverting slower than regular multiplication but still quick
                # enough.
                inv = np.argsort(self.perm)

                # Now take appropriate powers.
                q = inv
                for i in range(-other-1):
                    q = inv[q]
                return Perm(q)

        except:
            return NotImplemented("Powers are only defined for integers.")

    def __xor__(left, right):
        """
        For permutations x and y we have x^y is defined to be the
        conjugation of x by y, i.e. y^-1 * x * y. In all other cases it
        is undefined.

        Examples
        --------
        >>> p = clp.permutat.cyclestoperm([0, 1, 2, 3])
        >>> q = clp.permutat.cyclestoperm([0, 1], [2, 3])
        >>> p^q
        Perm((0, 3, 2, 1))
        >>> q^p
        Perm((0, 3)(1, 2))

        """
        try:
            # Now carry out the product.
            return (right**(-1) * left) * right

        except:
            return NotImplemented("Xor only defined for permutations.")

    def __getitem__(self, key):
        """
        If key is an integer then this returns the image of key under
        the permutation.

        .. warning::
            It is **always** faster to do self.perm[key] than self[key]
            and this should be the prefered method in CharLiePy
            implementation.

        If key is equivalent to slice(i, j) then this is equivalent to::

            >>> clp.permutat.restrictedperm(self, range(i, j))
        
        If self does not preserve the slice then this will throw an
        error. Note also that key must be a slice object. To use NumPy's
        advanced indexing then use self.perm[key].

        Examples
        --------
        >>> w = clp.permutat.cyclestoperm([0, 2, 3], [6, 4, 5])
        >>> w[0]
        2
        >>> w[0:4]
        (0, 2, 3)

        """
        if isinstance(key, int):
            return self.perm[key]
        elif isinstance(key, slice):
            return Perm(np.append(range(key.start), self.perm[key]))
        else:
            return NotImplemented

    # len returns the degree of the permutation. 
    def __len__(self):
        """Returns the degree of the permutation."""
        return self.deg

    # The remaining simply pass the request to the underlying numpy array.
    def __setitem__(self, key, value):
        """Changes the permutation so that the image of key is value."""
        self.perm[key] = value

    def __iter__(self):
        """Makes an iterable out of the underlying numpy array."""
        return iter(self.perm)

    def __reversed__(self):
        """Reverses the underlying numpy array."""
        return reversed(self.perm)

    def __contains__(self, item):
        """Returns True if the permutation is defined on item."""
        return item in self.perm

def permonsets(S, w):
    """
    Determines the permutation induced on a collection of sets.

    Parameters
    ----------
    S : list
        Each set in the list should be preserved by the permutation w.
    w : :class:`perm`
        A permutation inducing a permutation of the sets in S.

    Returns
    -------
    out : :class:`perm`
        The permutation on the entries in S induced by w.

    Examples
    --------
    >>> p = clp.permutat.perm(list(range(9))[::-1])
    >>> p
    (0, 8)(1, 7)(2, 6)(3, 5)
    >>> S = [set([0, 1, 2]), set([3, 4, 5]), set([6, 7, 8])]
    >>> clp.permutat.permonsets(S, p)
    (0, 2)

    """
    numsets = len(S)
    out = list(range(numsets))

    for i in range(numsets):
        image = set([w.perm[item] for item in S[i]])
        out[i] = S.index(image)

    return Perm(out)

def cyclestoperm(*cycles):
    """
    Constructs a permutation from a product of *disjoint* cycles.

    Parameters
    ----------
    cycles : list of lists
        Each entry of cycles should be a list like object describing a
        cycle.

    Returns
    -------
    out : :class:`perm`
        A permutation whose cycle decomposition is given by cycles.

    Notes
    -----
        * This function cannot accept cycles which are not disjoint.
        * One can force the degree of the resulting permutation to be as
          large as one likes by adding trivial cycles.

    Examples
    --------
    >>> p = clp.permutat.cyclestoperm([0, 2, 1], [8, 5, 3], [7, 6])
    >>> p
    (0, 2, 1)(3, 8, 5)(6, 7)
    >>> p.deg
    9

    By adding a trivial cycle we may now change the degree of the
    resulting permutation.

    >>> q = clp.permutat.cyclestoperm([0, 2, 1], [8, 5, 3], [7, 6], [20])
    >>> q
    (0, 2, 1)(3, 8, 5)(6, 7)
    >>> q.deg
    21

    """
    m = [max(x) for x in cycles]
    # If cycles is an empty list then return the identity.
    if not m:
        return Perm([])
    else:
        # Find the largest number occuring in cycles, this will be the degree
        # of the permutation. Then construct the identity of that degree.
        m = max(m)
        out = list(range(m+1))
        # We read the cycles from right to left.
        for item in cycles:
            # Makes i go to the next thing in the cycle.
            for i, j in pairwise(item):
                out[i] = j
            # Close the cycle.
            out[item[-1]] = item[0]
    
    return Perm(out)

def cycle(w, i):
    """
    Returns the cycle of w containing the point i.

    Notes
    -----
        * The cycle will always begin with the smallest entry.
        * The function modifies the variable ``i``.

    Examples
    --------
    >>> w = clp.permutat.cyclestoperm([1, 3, 2], [4, 5], [7, 8, 9])
    >>> clp.permutat.cycle(w, 3)
    (1, 3, 2)

    """
    # This function changes i.
    cycle = list()
    while not i in cycle:
        cycle.append(i)
        i = w.perm[i]

    # Return the cycle starting at the smallest entry.
    ind = cycle.index(min(cycle))
    return tuple(cycle[ind:] + cycle[:ind])

def cycles(w, L=None):
    """
    Gives the cycles of ``w`` acting on the entries in ``L``.

    Parameters
    ----------
    w : :class:`perm`
    L : iterable, optional
        The entries of ``L`` must form a subset of
        ``list(range(w.deg))`` such that ``w`` restricted to ``L`` is a
        permutation on ``L``.

    Returns
    -------
    out : tuple of tuples
        Each entry of ``out`` is a tuple describing a cycle of ``w``
        acting on ``L``.

    Notes
    -----
        * If ``L`` is not specified then ``range(w.deg)`` is used and
          the function returns almost the same as :func:`cycledecomp`,
          however the ouput of ``cycles`` includes trivial cycles of
          length 1.
        * The ordering of the cycles depends upon the ordering of ``L``.
          However, each cycle will start with the smallest entry.

    Examples
    --------
    >>> w = clp.permutat.cyclestoperm([1, 3, 2], [4, 5], [7, 8, 9])
    >>> clp.permutat.cycles(w, range(4, 10))
    ((4, 5), (6,), (7, 8, 9))
    >>> clp.permutat.cycles(w, reversed(range(4, 10)))
    ((7, 8, 9), (6,), (4, 5))
    
    """
    # If L is not specified then set it to range(w.deg)
    if L is None:
        L = range(w.deg)

    cycles = list()
    # Keeps track of which cycles we have already computed.
    cyclesum = list()
    for iItem in L:
        # If we have already computed the cycle containing iItem then move on.
        if iItem in cyclesum:
            continue
        icyc = cycle(w, iItem)
        cycles.append(icyc)
        cyclesum += icyc
     
    return tuple(cycles)

def cycledecomp(w):
    """
    Determines the decomposition of a permutation as a product of
    disjoint cycles.

    Returns
    -------
    out : tuple of tuples
        Each entry of ``out`` is a tuple describing a cycle in the cycle
        decomposition of the permutation.

    Notes
    -----
    This is similar to ``cycles(w)`` except cycles of length 1 are
    ignored.

    Examples
    --------
    >>> w = clp.permutat.cyclestoperm([3, 1, 2], [4, 5], [9, 8, 7])
    >>> clp.permutat.cycledecomp(w)
    ((1, 2, 3), (4, 5), (7, 9, 8))
    >>> clp.permutat.cycles(w)
    ((0,), (1, 2, 3), (4, 5), (6,), (7, 9, 8))

    """
    cycles = list()

    # Keeps track of which cycles we have already computed.
    cyclesum = list()
    for iItem in w.perm:
        # If we have already computed the cycle containing iItem then move on.
        if iItem in cyclesum:
            continue
        # Compute the cycle containg iItem and only include it if it is non
        # trivial
        icyc = cycle(w, iItem)
        if len(icyc) != 1:
            cycles.append(icyc)
        cyclesum += icyc

    # If w is the identity then cycles will be empty, so fix this.
    if cycles == []:
        cycles.append(tuple())
        
    return tuple(cycles) 

def cycletype(w, part=False):
    """
    Returns the cycletype of a permutation.

    Parameters
    ----------
    w : :func:`perm`
    part : bool, optional
        If ``True`` then the function returns a list whose entries
        consist of the lengths of the cycles in the cycle decomposition
        of ``w``. This is naturally a partition of ``w.deg``.

    Returns
    -------
    out : list of ints
        For any ``i in range(1, w.deg+1)`` we have ``out[i]`` is the
        number of cycles of length ``i`` occurring in the cycle
        decomposition of ``w``. As a convention we set ``out[0] = 0``.

    Notes
    -----
    If ``part`` is set to ``True`` then the resulting list is sorted and
    weakly decreasing.

    Examples
    --------
    >>> w = clp.permutat.cyclestoperm([1, 2, 3], [4, 5], [7, 9, 8])
    >>> clp.permutat.cycletype(w)
    [0, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0]
    >>> clp.permutat.cycletype(w, part=True)
    [3, 3, 2, 1, 1]

    >>> w = clp.permutat.cyclestoperm([0, 1, 2, 3, 4, 5, 6])
    >>> clp.permutat.cycletype(w)
    [0, 0, 0, 0, 0, 0, 0, 1]
    >>> clp.permutat.cycletype(w, part=True)
    [7]

    """
    # This is essentially the same as cycledecomp but only keeps track of
    # the length of the cycles and not the cycles themselves.
    if part:
        cyctype = []
        cyclesum = set() 
        for iItem in w.perm:
            # Check if we've already seen it.
            if iItem in cyclesum:
                continue

            # Compute the cycle containing iItem.
            cycle = set()
            while not i in cycle:
                cycle.append(i)
                i = w.perm[i]
            cyctype.append(len(icyc))
            cyclesum |= icyc
        cyctype.sort(reverse=True)
    else:
        cyctype = [0]*(w.deg+1)
        cyclesum = set() 
        for iItem in w.perm:
            # Check if we've already seen it.
            if iItem in cyclesum:
                continue

            # Compute the cycle containing iItem.
            cycle = set()
            while not i in cycle:
                cycle.append(i)
                i = w.perm[i]
            cyctype.append(len(icyc))
            cyclesum |= icyc
            icyc = cycle(w, iItem)
            cyctype[len(icyc)] += 1
            cyclesum += icyc

    return cyctype

def restrictedperm(w, L):
    """
    Gives a permutation obtained by restriction to a subset.

    Parameters
    ----------
    w : :func:`perm`
    L : list of ints
        The entries of ``L`` must form a subset of
        ``list(range(w.deg))`` such that the restriction of ``w`` to
        ``L`` is a permutation on ``L``.

    Returns
    -------
    out : :func:`perm`
        A permutation such that ``out[i] == w[i]`` for all ``i in L``.

    Notes
    -----
    The degree of the resulting permutation is determined by the largest
    entry in ``L``.
    
    Examples
    --------
    >>> w = clp.permutat.cyclestoperm([1, 2, 3], [4, 5], [7, 9, 8])
    >>> clp.permutat.restrictedperm(w, [1, 2, 3, 4, 5, 6])
    (1, 2, 3)(4, 5)
    >>> clp.permutat.restrictedperm(w, [4, 5, 7, 8, 9])
    (4, 5)(7, 9, 8)
    >>> (clp.permutat.restrictedperm(w, [1, 2, 3, 4, 5, 6]).deg ==
    ...  clp.permutat.restrictedperm(w, [1, 2, 3, 4, 5]).deg)
    False

    """
    out = list(range(max(L)+1))
    for i in L:
        out[i] = w.perm[i]

    return Perm(out)

def orderperm(w):
    """
    Returns the order of a permutation.

    Examples
    --------
    >>> w = cyclestoperm([0, 1, 2], [3, 4], [5, 6, 7])
    >>> w
    (0, 1, 2)(3, 4)(5, 6, 7)
    >>> orderperm(w)
    6

    """
    i = 1
    idperm = np.arange(w.deg)
    p = w.perm
    while not np.array_equal(p, idperm):
        p = w.perm[p]
        i += 1
    return i


##### Utility Functions ######

# The following is from the itertools documentation.
def pairwise(iterable):
    "s = (s0, s1, ...) -> (s0, s1), (s1, s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)
