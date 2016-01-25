#from pycox import core
from __future__ import print_function

from . import algebra as alg
from .algebra import permutat
from . import data
from . import utils
from . import core

import numpy as np
import itertools
import sys

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

    >>> W = pcx.coxeter("A", 4)
    >>> phi = pcx.longestperm(W)
    >>> Wphi = pcx.coset(W, phi)
    >>> Wphi.phitype
    [['A', [0, 1, 2, 3], (0, 3)(1, 2)]]

    """
    def __init__(self, W, phi=tuple()):
        self.group = W

        # Convert negative to positive roots.
        def flipint(N, n):
            if n > N-1:
                return n - N
            else:
                return n

        # define the automorphism phi
        self.phi = permutat.perm([flipint(self.group.N, i) for i in phi])
        
        # Make sure phi preserves the Coxeter matrix. We use here that the
        # Coxeter matrix is symmetric to make things quicker.
        for s in self.group.rank:
            for t in xrange(s):
                if (self.group.coxetermat[s][t]
                        != self.group.coxetermat[self.phi[s]][self.phi[t]]):
                    print("Permutation does not preserve Coxeter matrix!")
                    return None
        
        # Determine the orbits of phi acting cyclically on isomorphic factors.
        typsets = [set(x[1]) for x in self.group.cartantype]
        indperm = permutat.permonsets(typsets, self.phi)
        self.orbits = permutat.cycles(indperm)

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
        self.order = self.group.order

    def __repr__(self):
        """Lists the type of the Coxeter coset. A collection of factors
        in parentheses indicates that these factors are cyclically
        permuted by the automorphism. A number preceeding the first
        parenthesis denotes that the automorphism induced on that factor
        is of that given order. Orbits of factors are seperated by an x.

        Examples
        --------
        >>> phi = range(9)[::-1]
        >>> coset(coxeter("A",3,"A",3,"A",3), phi)
        (A3.A3) x 2A3

        >>> W = coxeter("A",2,"A",2)
        >>> coset(W, permutat.cyclestoperm([0, 3, 1, 2]))
        2(A2.A2)

        For Steinberg's triality group we keep track of which version of
        the automorphism occurs.

        >>> phi = permutat.cyclestoperm([0,1,3])
        >>> W = coxeter("D",4)
        >>> coset(W, permutat.cyclestoperm([0, 1, 3])
        3D4
        >>> coset(W, permutat.cyclestoperm([0, 3, 1])
        3'D4
        
        """
        pout = ''

        for i, orb in enumerate(self.orbits):
            orbtype = self.group.cartantype[orb[0]]
            ordtwist = permutat.orderperm(self.phitype[i][2])

            # Add the order of the twist if it's non-trivial.
            if ordtwist == 2:
                pout += '2'
            elif ordtwist == 3:
                ind = self.phitype[i][1]
                phi = self.phitype[i][2]
                if phi[ind[0]] == ind[1]:
                    pout += '3'
                elif phi[ind[0]] == ind[3]:
                    pout += "3'"

            if len(orb) != 1:
                pout += '('
            pout += (len(orb)-1)*(orbtype[0] + str(len(orbtype[1])) + '.')
            pout += orbtype[0] + str(len(orbtype[1]))
            if len(orb) != 1:
                pout += ')'
            pout += ' x '

        # Take off the final space and x.
        return pout[:-2]


# This overwrites the conjugacyclasses function in core.
def conjugacyclasses(W):
    """Computes the conjugacy classes of a Coxeter group/coset.

    Parameters
    ----------
    W : :class:`~pycox.coxeter` or :class:`coset`

    Returns
    -------
    out : dictionary
        Contains three entries ``reps``, ``classlengths`` and
        ``classnames`` which are all lists. The first contains
        representatives of minimal length, the second the size of the
        conjugacy classes and finally the names of the conjugacy
        classes.

    See Also
    --------
    core.conjugacyclass
    core.conjtomin

    Notes
    -----
        * The conventions are the same as in `CHEVIE`_; in particular
          the ordering of the classes matches that in CHEVIE.
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

    >>> W = coxeter("A",3)
    >>> phi = range(3)[::-1]
    >>> Wphi = coset(W,phi)
    >>> Wphi.phitype
    [['A', [0, 1, 2], (0, 2)]]
    >>> conjugacyclasses(Wphi)
    {'classnames'   : [('1111',), ('211',), ('22',), ('31',), ('4',)],
     'classlengths' : [3, 6, 1, 8, 6],
     'reps': [[0, 1, 0, 2, 1, 0], [1], [], [1, 2], [2]]}

    """
    if 'conjugacyclasses' in dir(W):
        return W.conjugacyclasses

    reps = list()
    clens = list()
    names = list()

    if isinstance(W, core.coxeter):
        # First get the data for each irreducible component. This dynamically
        # calls the correct module for the type. Note that, here, the function
        # returns centraliser orders not class lengths.
        for typ in W.cartantype:
            treps, tclens, tnames = sys.modules['pycox.data.typ1'
                + typ[0]].conjclassdata(typ[1])
            reps.append(treps)
            clens.append(tclens)
            names.append(tnames)

        # Only need to take cartesian products of class reps and class lengths
        # if there is more than one irreducible component.
        if len(W.cartantype) > 1:
            reps = [list(itertools.chain.from_iterable(rep))
                    for rep in itertools.product(*reps)]
            clens = [W.order//sum(rep) for rep in itertools.product(*clens)]
        else:
            reps = reps[0]
            clens = [W.order//cl for cl in clens[0]]
 
        names = list(itertools.product(*names))

    elif isinstance(W, coset):
        for typ in W.phitype:
            m = 0
            if typ[0][0] == 'I':
                m = int(typ[0][1])
            tord = permutat.orderperm(typ[2])
            treps, tclens, tnames = sys.modules['pycox.data.typ'
                + str(tord) + typ[0]].conjclassdata(typ[1], m=m, phi=typ[2])
            reps.append(treps)
            clens.append(tclens)
            names.append(tnames)

        # Only need to take cartesian products of class reps and class lengths
        # if there is more than one irreducible component.
        if len(W.phitype) > 1:
            reps = [list(itertools.chain.from_iterable(rep))
                    for rep in itertools.product(*reps)]
            clens = [W.order//sum(rep) for rep in itertools.product(*clens)]
        else:
            reps = reps[0]
            clens = [W.order//cl for cl in clens[0]]
 
        names = list(itertools.product(*names))

    W.conjugacyclasses = {'reps':reps, 'classnames':names,
                          'classlengths':clens}

    return W.conjugacyclasses

def chartable(W, chars=True):
    """returns the  ordinary  character table  of W,  together with related 
    additional information.
   
    The result is a dictionary with at least the following entries:
  
      classlengths    sizes of the conjugacy classes
      classnames      see help to 'conjugacyclasses'
      classreps       representatives of minimal length
      charnames       tuples of names for the characters put together 
                                    from the irreducible components of W
      irreducibles    square matrix of character values
      position_id     position of trivial character
      position_sgn    position of sign character
      position_refl   position of reflection character (W irreducible)
      permsgn         permutation induced by tensoring with sign
      b               b-invariants (see also 'fakedegree')
      a               a-invariants (with respect to equal parameters; 
                                                see also 'ainvariants').
  
    The most expensive part of this function  is the  computation of the 
    character table.  If the optional argument  'chars' is set to False, 
    then the entry  'irreducibles'  will be omitted,  in which case  the 
    result  of the  function  is roughly  equivalent  to the  gap-chevie 
    function 'ChevieCharInfo'.
  
    The  raw data  for the  various types of  irreducible finite Coxeter 
    groups  are   explicitly  stored  in  this   module  and  called via 
    'irrchardata(typ,n)'.  For  a  general  W  the  data  are then built  
    together  from the irreducible components using  'W.cartantype'. See 
    also 'chartablesymmetric', 'chartableB', 'chartableD'.
  
    >>> chartable(coxeter("G",2))
    {'irreducibles'  : [[1, 1, 1, 1, 1, 1], 
                        [1,-1,-1, 1, 1, 1], 
                        [1, 1,-1,-1, 1,-1], 
                        [1,-1, 1,-1, 1,-1], 
                        [2, 0, 0, 1,-1,-2], 
                        [2, 0, 0,-1,-1, 2]], 
      'classnames'   : [(' ',),('~A_1',),('A_1',),('G_2',),('A_2',),
                        ('A_1+~A_1',)], 
      'classreps'    : [[],[1],[0],[0,1],[0,1,0,1],[0,1,0,1,0,1]],
      'classlengths' : [1,3,3,2,2,1],
      'b'            : [0,6,3,3,1,2], 
      'a'            : [0,6,1,1,1,1], 
      'charnames'    : [('phi_{1,0}',),('phi_{1,6}',),("phi_{1,3}'",),
                        ("phi_{1,3}'',)",('phi_{2,1}',),('phi_{2,2}',)]], 
      'position_id'  : 0,              
      'position_sgn' : 1,              
      'position_refl': 4,              
      'permsgn'      : [1,0,3,2,4,5]}  
  
    See also 'displaychartable'.
    """
    ctabtest = 'chartable' in dir(W)
    
    if ctabtest and (chars == False or 'irreducibles' in W.chartable.keys()):
        return W.chartable

    names, avals, bvals = [], [], []

    if isinstance(W, core.coxeter):
        # First get the data for each irreducible component. This
        # dynamically calls the correct module for the type.
        for typ in W.cartantype:
            m = 0
            if typ[0][0] == 'I':
                m = int(typ[0][1])

            tnames, taval, tbval = sys.modules['pycox.data.typ1'
                + typ[0]].irrchardata(len(typ[1]), m=m)
            names.append(tnames)
            avals.append(taval)
            bvals.append(tbval)

        # Only need to take cartesian products of a values and b values if
        # there is more than one irreducible component.
        if len(W.cartantype) > 1:
            avals = [sum(aval) for aval in itertools.product(*avals)]
            bvals = [sum(bval) for bval in itertools.product(*bvals)]
        else:
            avals, bvals = avals[0], bvals[0]
 
        names = list(itertools.product(*names))

        # Get conjugacy class information.
        if 'conjugacyclasses' in dir(W):
            conjclasses = W.conjugacyclasses
        else:
            conjclasses = conjugacyclasses(W)

        # If it's not there then construct the chartable entry for W.
        if not ctabtest:
            W.chartable = {'charnames':names,
                           'classnames':conjclasses['classnames'],
                           'classlengths':conjclasses['classlengths'],
                           'classreps':conjclasses['reps'],
                           'b':bvals, 'a':avals}

        try:
            W.irreducibles
        except:
            if chars:
                # Get the character table of the first irreducible component.
                typ = W.cartantype[0]
                m = 0
                if typ[0][0] == 'I':
                    m = int(typ[0][1])
                mat = sys.modules['pycox.data.typ1'
                        + typ[0]].chartable(len(typ[1]), m=m)

                # Take the successive kronecker product with all remaining
                # ones.
                for typ in W.cartantype[1:]:
                    m = 0
                    if typ[0][0] == 'I':
                        m = int(typ[0][1])
                    mat = np.kron(mat, sys.modules['pycox.data.typ1'
                            + typ[0]].chartable(len(typ[1]), m=m))
                
                # Determine position of identity and sign characters.
                idchar = np.array([1]*len(conjclasses['reps']))
                sgn = np.array([(-1)**len(w)
                        for w in conjclasses['reps']], dtype='int')

                # Position of reflection character.
                if len(W.cartantype) > 1:
                    prf = False
                else:
                    Wrank = len(W.rank)
                    if Wrank == 0:
                        prf = 0
                    else:
                        L = zip(mat[:,0], bvals)
                        prf = L.index((Wrank, 1))

                # Now add the data we've just collected.
                W.chartable['irreducibles'] = mat
                W.chartable['position_id'] = utils.npindex(mat, idchar)
                W.chartable['position_sgn'] = utils.npindex(mat, sgn)
                W.chartable['position_refl'] = prf
                #W.chartable['permsgn'] = psgn
    elif isinstance(W, coset):
        # First get the data for each irreducible component. This
        # dynamically calls the correct module for the type.
        for typ in W.phitype:
            m = 0
            if typ[0][0] == 'I':
                m = int(typ[0][1])

            # If phi is the identity on the factor then irrchardata will return
            # the avals and bvals as well, so just take the first entry.
            tord = permutat.orderperm(typ[2])
            tnames = sys.modules['pycox.data.typ' + str(tord)
                + typ[0]].irrchardata(len(typ[1]), m=m, phi=typ[2])[0]
            names.append(tnames)

        names = list(itertools.product(*names))

        # Get conjugacy class information.
        if 'conjugacyclasses' in dir(W):
            conjclasses = W.conjugacyclasses
        else:
            conjclasses = conjugacyclasses(W)

        # If it's not already there then construct the chartable entry for W.
        if not ctabtest:
            W.chartable = {'charnames':names,
                           'classnames':conjclasses['classnames'],
                           'classlengths':conjclasses['classlengths'],
                           'classreps':conjclasses['reps']}

        try:
            W.irreducibles
        except:
            if chars:
                # Get the character table of the first irreducible component.
                typ = W.phitype[0]
                m = 0
                if typ[0][0] == 'I':
                    m = int(typ[0][1])
                tord = permutat.orderperm(typ[2])
                mat = sys.modules['pycox.data.typ' + str(tord)
                        + typ[0]].chartable(len(typ[1]), m=m, phi=typ[2])

                # Take the successive kronecker product with all remaining
                # ones.
                for typ in W.phitype[1:]:
                    m = 0
                    if typ[0][0] == 'I':
                        m = int(typ[0][1])
                    tord = permutat.orderperm(typ[2])
                    mat = np.kron(mat, sys.modules['pycox.data.typ' + str(tord)
                            + typ[0]].chartable(len(typ[1]), m=m, phi=typ[2]))
                
                ## Determine position of identity and sign characters.
                #idchar = np.array([1]*len(conjclasses['reps']))
                #sgn = np.array([(-1)**len(w)
                #        for w in conjclasses['reps']], dtype='int')

                ## Position of reflection character.
                #if len(W.cartantype) > 1:
                #    prf = False
                #else:
                #    Wrank = len(W.rank)
                #    if Wrank == 0:
                #        prf = 1
                #    else:
                #        L = zip(mat[:,0], bvals)
                #        prf = L.index((Wrank, 1))

                # Now add the data we've just collected.
                W.chartable['irreducibles'] = mat
                #W.chartable['position_id'] = utils.npindex(mat, idchar)
                #W.chartable['position_sgn'] = utils.npindex(mat, sgn)
                #W.chartable['position_refl'] = prf

    return W.chartable

def linearperm(W, lchar=-1):
    """Computes the permutation on the character table of a Coxeter
    group W induced by tensoring all irreducible characters with a
    linear character.

    If called with no optional arguments the function uses the sign
    character of W. Additionally one may choose which linear character
    to use by setting lchar to be the index of the linear character in
    the character table of W.
    
    Note that this function will need to compute the character table of
    W if it has not already been computed.

    >>> chartable(coxeter("G",2))['irreducibles']
    array([[ 1,  1,  1,  1,  1,  1],
           [ 1, -1, -1,  1,  1,  1],
           [ 1,  1, -1, -1,  1, -1],
           [ 1, -1,  1, -1,  1, -1],
           [ 2,  0,  0,  1, -1, -2],
           [ 2,  0,  0, -1, -1,  2]])
    >>> linearperm(coxeter("G",2))
    (0, 1)(2, 3)
    >>> linearperm(coxeter("G",2), lchar=2)
    (0, 2)(1, 3)(4, 5)
    """
    charinfo = chartable(W)
    chartab = charinfo['irreducibles']

    if lchar == -1:
        lchar = charinfo['position_sgn']

    return permutat.perm(utils.npperm(chartab, chartab*chartab[lchar])) 









