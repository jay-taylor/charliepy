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


########################################################################
########################################################################
##                                                                    ##
##              Character Tables and Conjugacy Classes                ##
##                                                                    ##

#def conjugacyclasses(W):
#    """
#    Computes the conjugacy classes of a Coxeter coset.
#
#    Parameters
#    ----------
#    W : :class:`CoxeterCoset`
#
#    Returns
#    -------
#    out : dictionary
#        Contains three entries ``reps``, ``classlengths`` and
#        ``classnames`` which are all lists. The first contains
#        representatives of minimal length, the second the size of the
#        conjugacy classes and finally the names of the conjugacy
#        classes.
#
#    See Also
#    --------
#    core.conjugacyclass
#    core.conjtomin
#
#    Notes
#    -----
#        * The conventions are the same as in `CHEVIE`_; in particular
#          the ordering of the classes matches that in CHEVIE.
#        * For the conventions concerning the raw data for each
#          irreducible type see the databook.
#        * For a coset this computes the phi-conjugacy classes of the
#          coset.
#        * For a group the representatives are *very good* in the
#          sense of [GM97]_.
#        * Once computed the classes are stored as the atribute
#          ``conjugacyclasses`` of W.
#
#    Examples
#    --------
#    >>> W = clp.CoxeterGroup("A", 3)
#    >>> phi = range(2, -1, -1)
#    >>> Wphi = clp.CoxeterCoset(W, phi)
#    >>> Wphi.phitype
#    [['A', [0, 1, 2], (0, 2)]]
#    >>> conjugacyclasses(Wphi)
#    {'classnames'   : [('1111',), ('211',), ('22',), ('31',), ('4',)],
#     'classlengths' : [3, 6, 1, 8, 6],
#     'reps': [[0, 1, 0, 2, 1, 0], [1], [], [1, 2], [2]]}
#
#    """
#    try:
#        return W.conjugacyclasses
#    except AttributeError:
#        pass
#
#    # First get the data for each irreducible component. This dynamically
#    # calls the correct module for the type. Note that, here, the function
#    # returns centraliser orders not class lengths.
#
#    return W.conjugacyclasses
#
#def chartable(W):
#    """returns the  ordinary  character table  of W,  together with related 
#    additional information.
#   
#    The result is a dictionary with at least the following entries:
#  
#      classlengths    sizes of the conjugacy classes
#      classnames      see help to 'conjugacyclasses'
#      classreps       representatives of minimal length
#      charnames       tuples of names for the characters put together 
#                                    from the irreducible components of W
#      irreducibles    square matrix of character values
#      position_id     position of trivial character
#      position_sgn    position of sign character
#      position_refl   position of reflection character (W irreducible)
#      permsgn         permutation induced by tensoring with sign
#      b               b-invariants (see also 'fakedegree')
#      a               a-invariants (with respect to equal parameters; 
#                                                see also 'ainvariants').
#  
#    The most expensive part of this function  is the  computation of the 
#    character table.  If the optional argument  'chars' is set to False, 
#    then the entry  'irreducibles'  will be omitted,  in which case  the 
#    result  of the  function  is roughly  equivalent  to the  gap-chevie 
#    function 'ChevieCharInfo'.
#  
#    The  raw data  for the  various types of  irreducible finite Coxeter 
#    groups  are   explicitly  stored  in  this   module  and  called via 
#    'irrchardata(typ,n)'.  For  a  general  W  the  data  are then built  
#    together  from the irreducible components using  'W.cartantype'. See 
#    also 'chartablesymmetric', 'chartableB', 'chartableD'.
#  
#    >>> chartable(coxeter("G",2))
#    {'irreducibles'  : [[1, 1, 1, 1, 1, 1], 
#                        [1,-1,-1, 1, 1, 1], 
#                        [1, 1,-1,-1, 1,-1], 
#                        [1,-1, 1,-1, 1,-1], 
#                        [2, 0, 0, 1,-1,-2], 
#                        [2, 0, 0,-1,-1, 2]], 
#      'classnames'   : [(' ',),('~A_1',),('A_1',),('G_2',),('A_2',),
#                        ('A_1+~A_1',)], 
#      'classreps'    : [[],[1],[0],[0,1],[0,1,0,1],[0,1,0,1,0,1]],
#      'classlengths' : [1,3,3,2,2,1],
#      'b'            : [0,6,3,3,1,2], 
#      'a'            : [0,6,1,1,1,1], 
#      'charnames'    : [('phi_{1,0}',),('phi_{1,6}',),("phi_{1,3}'",),
#                        ("phi_{1,3}'',)",('phi_{2,1}',),('phi_{2,2}',)]], 
#      'position_id'  : 0,              
#      'position_sgn' : 1,              
#      'position_refl': 4,              
#      'permsgn'      : [1,0,3,2,4,5]}  
#  
#    See also 'displaychartable'.
#    """
#    ctabtest = 'chartable' in dir(W)
#    
#    if ctabtest and (chars == False or 'irreducibles' in W.chartable.keys()):
#        return W.chartable
#
#    names, avals, bvals = [], [], []
#
#    # First get the data for each irreducible component. This
#    # dynamically calls the correct module for the type.
#    for typ in W.phitype:
#        m = 0
#        if typ[0][0] == 'I':
#            m = int(typ[0][1])
#
#        # If phi is the identity on the factor then irrchardata will return
#        # the avals and bvals as well, so just take the first entry.
#        tord = typ[2].order
#        tnames = sys.modules['pycox.data.typ' + str(tord)
#            + typ[0]].irrchardata(len(typ[1]), m=m, phi=typ[2])[0]
#        names.append(tnames)
#
#    names = list(itertools.product(*names))
#
#    # Get conjugacy class information.
#    if 'conjugacyclasses' in dir(W):
#        conjclasses = W.conjugacyclasses
#    else:
#        conjclasses = conjugacyclasses(W)
#
#    # If it's not already there then construct the chartable entry for W.
#    if not ctabtest:
#        W.chartable = {'charnames':names,
#                       'classnames':conjclasses['classnames'],
#                       'classlengths':conjclasses['classlengths'],
#                       'classreps':conjclasses['reps']}
#
#    try:
#        W.irreducibles
#    except:
#        if chars:
#            # Get the character table of the first irreducible component.
#            typ = W.phitype[0]
#            m = 0
#            if typ[0][0] == 'I':
#                m = int(typ[0][1])
#            tord = typ[2].order
#            mat = sys.modules['pycox.data.typ' + str(tord)
#                    + typ[0]].chartable(len(typ[1]), m=m, phi=typ[2])
#
#            # Take the successive kronecker product with all remaining
#            # ones.
#            for typ in W.phitype[1:]:
#                m = 0
#                if typ[0][0] == 'I':
#                    m = int(typ[0][1])
#                tord = typ[2].order
#                mat = np.kron(mat, sys.modules['pycox.data.typ' + str(tord)
#                        + typ[0]].chartable(len(typ[1]), m=m, phi=typ[2]))
#            
#            ## Determine position of identity and sign characters.
#            #idchar = np.array([1]*len(conjclasses['reps']))
#            #sgn = np.array([(-1)**len(w)
#            #        for w in conjclasses['reps']], dtype='int')
#
#            ## Position of reflection character.
#            #if len(W.cartantype) > 1:
#            #    prf = False
#            #else:
#            #    Wrank = len(W.rank)
#            #    if Wrank == 0:
#            #        prf = 1
#            #    else:
#            #        L = zip(mat[:,0], bvals)
#            #        prf = L.index((Wrank, 1))
#
#            # Now add the data we've just collected.
#            W.chartable['irreducibles'] = mat
#            #W.chartable['position_id'] = utils.npindex(mat, idchar)
#            #W.chartable['position_sgn'] = utils.npindex(mat, sgn)
#            #W.chartable['position_refl'] = prf
#
#    return W.chartable
#
#def linearperm(W, lchar=-1):
#    """Computes the permutation on the character table of a Coxeter
#    group W induced by tensoring all irreducible characters with a
#    linear character.
#
#    If called with no optional arguments the function uses the sign
#    character of W. Additionally one may choose which linear character
#    to use by setting lchar to be the index of the linear character in
#    the character table of W.
#    
#    Note that this function will need to compute the character table of
#    W if it has not already been computed.
#
#    >>> chartable(coxeter("G",2))['irreducibles']
#    array([[ 1,  1,  1,  1,  1,  1],
#           [ 1, -1, -1,  1,  1,  1],
#           [ 1,  1, -1, -1,  1, -1],
#           [ 1, -1,  1, -1,  1, -1],
#           [ 2,  0,  0,  1, -1, -2],
#           [ 2,  0,  0, -1, -1,  2]])
#    >>> linearperm(coxeter("G",2))
#    (0, 1)(2, 3)
#    >>> linearperm(coxeter("G",2), lchar=2)
#    (0, 2)(1, 3)(4, 5)
#    """
#    charinfo = chartable(W)
#    chartab = charinfo['irreducibles']
#
#    if lchar == -1:
#        lchar = charinfo['position_sgn']
#
#    return permutat.Perm(utils.npperm(chartab, chartab*chartab[lchar])) 









