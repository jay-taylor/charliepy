from . import permutat
from . import utils
from . import core

import numpy as np
import itertools
import sys


########################################################################
########################################################################
##                                                                    ##
##              Character Tables and Conjugacy Classes                ##
##                                                                    ##

def conjugacyclasses(W):
    """
    Computes the conjugacy classes of a Coxeter coset.

    Parameters
    ----------
    W : :class:`CoxeterCoset`

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
    >>> W = clp.CoxeterGroup("A", 3)
    >>> phi = range(2, -1, -1)
    >>> Wphi = clp.CoxeterCoset(W, phi)
    >>> Wphi.phitype
    [['A', [0, 1, 2], (0, 2)]]
    >>> conjugacyclasses(Wphi)
    {'classnames'   : [('1111',), ('211',), ('22',), ('31',), ('4',)],
     'classlengths' : [3, 6, 1, 8, 6],
     'reps': [[0, 1, 0, 2, 1, 0], [1], [], [1, 2], [2]]}

    """
    try:
        return W.conjugacyclasses
    except AttributeError:
        pass

    # First get the data for each irreducible component. This dynamically
    # calls the correct module for the type. Note that, here, the function
    # returns centraliser orders not class lengths.

    return W.conjugacyclasses

def chartable(W):
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

    # First get the data for each irreducible component. This
    # dynamically calls the correct module for the type.
    for typ in W.phitype:
        m = 0
        if typ[0][0] == 'I':
            m = int(typ[0][1])

        # If phi is the identity on the factor then irrchardata will return
        # the avals and bvals as well, so just take the first entry.
        tord = typ[2].order
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
            tord = typ[2].order
            mat = sys.modules['pycox.data.typ' + str(tord)
                    + typ[0]].chartable(len(typ[1]), m=m, phi=typ[2])

            # Take the successive kronecker product with all remaining
            # ones.
            for typ in W.phitype[1:]:
                m = 0
                if typ[0][0] == 'I':
                    m = int(typ[0][1])
                tord = typ[2].order
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

    return permutat.Perm(utils.npperm(chartab, chartab*chartab[lchar])) 









