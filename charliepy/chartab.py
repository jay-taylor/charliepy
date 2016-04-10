from . import core
from . import permutat
from . import utils
from . import conjclass

import sys
import itertools
import functools
import operator
import collections
import numpy as np

__all__ = ['CharacterTable']


########################################################################
########################################################################
##                                                                    ##
##                           Characters                               ##
##                                                                    ##

class Character(collections.namedtuple("Character",
            ["parent", "name", "param", "degree", "a", "b"])):
    # Setting slots to an empty tuple prevents a dict for the class
    # being created.
    __slots__ = () 
    def __repr__(self):
        return "Character( {}, {} )".format(self.parent,
                                            self.name)

    @property
    def length(self):
        return self.parent.size//self.centsize

class CharacterTable:
    """
    returns the  ordinary  character table  of W,  together with related 
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
    def __init__(self, *args, **kwargs):
        W = args[0]
        self.parent = W

        # Get the conjugacy classes first.
        self.classes = conjclass.conjugacyclasses(W)

        if isinstance(W, core.CoxeterGroup):
            # We treat the irreducible case differently because it's simpler.
            if len(W.cartantype) == 1:
                typ, inds = W.cartantype[0]
                n = len(inds)

                # Get the correct module from the data package for the
                # irreducible factor.
                mod = utils.getmodule(typ)

                # Irreducible character data.
                self.charnames, self.a, self.b = zip(*mod.irrchars(n))

                # Finally the irreducible characters.
                self.irreducibles = mod.chartable(n)

            else:
                names, avals, bvals = [], [], []

                # Assume W is a direct product of n irreducible Coxeter
                # groups. The following generators work as follows. The call
                # to itertools.product constructs tuples of the form
                #
                #       (x_0, ..., x_n)
                #
                # where x_i denotes the information for irreducible
                # character/conjugacy class of the ith factor of W. The call
                # to zip then puts all of this information together into a
                # single tuple.

                # Irreducible character data.
                gen = (zip(*elm)
                   for elm in itertools.product(
                       *(utils.getmodule(typ).irrchars(len(inds))
                           for (typ, inds) in W.cartantype)))

                self.charnames, self.a, self.b = [], [], []

                for name, a, b in gen:
                    self.charnames.append(",".join(name))
                    self.a.append(sum(a))
                    self.b.append(sum(b))

                # Finally the values of the irreducible characters.
                typ, inds = self.parent.cartantype[0]
                mat = utils.getmodule(typ).chartable(len(inds))

                # Take the successive tensor product with all remaining
                # character tables.
                for (typ, inds) in W.cartantype[1:]:
                    new =  utils.getmodule(typ).chartable(len(inds))
                    mat = [[a*b for a, b in itertools.product(mrow, nrow)]
                           for mrow, nrow in itertools.product(mat, new)]

                self.irreducibles = mat

        elif isinstance(self.parent, core.CoxeterCoset):
            # We treat the irreducible case differently because it's simpler.
            if len(self.parent.phitype) == 1:
                # Get the correct module from the data package for the
                # irreducible factor.
                module = sys.modules[
                    'charliepy.data.typ{}{}'.format(
                        self.parent.phitype[0][2],
                        self.parent.phitype[0][0])]
                inds = self.parent.phitype[0][1]

                # First the conjugacy class data.
                self.classnames, self.centralisers = list(
                    zip(*module.conjclasses_min(inds)))

                # Now the irreducible character data. We don't need the
                # a-values and b-values for cosets.
                self.charnames, a, b = list(
                    zip(*module.irrchars(len(inds))))
                del a, b

                # Finally the irreducible characters.
                self.irreducibles = module.chartable(len(inds))

            else:
                names, avals, bvals = [], [], []

                # Assume W is a direct product of n irreducible Coxeter
                # groups. The following generators work as follows. The call
                # to itertools.product constructs tuples of the form
                #
                #       (x_0, ..., x_n)
                #
                # where x_i denotes the information for irreducible
                # character/conjugacy class of the ith factor of W. The call
                # to zip then puts all of this information together into a
                # single tuple.

                # First the conjugacy class data.
                gen = (zip(*elm) for elm in 
                    itertools.product(*(
                        sys.modules[
                            'charliepy.data.typ{}{}'.format(
                                typ[2].order, typ[0])
                        ].conjclasses_min(typ[1])
                        for typ in self.parent.phitype)))

                self.classnames, self.centralisers = [], []

                for name, cent in gen:
                    self.classnames.append(','.join(name))
                    self.centralisers.append(cent)

                # Now the irreducible character data.
                gen = (zip(*elm) for elm in 
                    itertools.product(*(
                        sys.modules[
                            'charliepy.data.typ{}{}'.format(
                                typ[2].order, typ[0])
                        ].irrchars(len(typ[1]))
                        for typ in self.parent.phitype)))

                self.charnames = []

                for name, a, b in gen:
                    self.charnames.append(",".join(name))

                # Finally the values of the irreducible characters.
                typ = self.parent.phitype[0]
                mat = sys.modules[
                    'charliepy.data.typ{}{}'.format(
                        typ[2].order, typ[0])
                ].chartable(len(typ[1]))

                # Take the successive kronecker product with all
                # remaining ones. Note that the definition of the
                # Kronecker product in numpy is compatible with the
                # Cartesian product defined in itertools.
                for typ in self.parent.phitype[1:]:
                    new =  sys.modules[
                        'charliepy.data.typ{}{}'.format(typ[2].order, typ[0])
                    ].chartable(len(typ[1]))
                    mat = [[a*b for a, b in itertools.product(mrow, nrow)]
                           for mrow, nrow in itertools.product(mat, new)]

                self.irreducibles = mat

        else:
            raise NotImplemented("Only implemented for Coxeter Groups and"
                "Coxeter Cosets.")

    def __repr__(self):
        return "CharacterTable({})".format(self.parent)

    def __str__(self):
        return utils.arraytostring(
            self.irreducibles,
            rowlabels = self.charnames,
            collabels = self.classnames,
            maxwidth = 78
        )

    @property
    def classnames(self):
        return [cls.name for cls in self.classes]



    def isconsistent(self):
        """
        Checks whether the row and column orthogonality relations are
        satisfied for the character table. This doesn't rule out the
        possibility of an incorrect table but gives some indication that
        things are working correctly. This is mostly here for testing
        purposes and can be time consuming for large tables, say with
        more than 1000 rows and columns.

        """
        size = self.parent.size
        cents = self.centralisers
        X = self.irreducibles

        # Row orthogonality relations.
        for i in range(len(X)):
            if sum((size*a*a)//c for a, c in zip(X[i], cents)) != size:
                return False
            for j in range(i+1, len(X)):
                if sum((size*a*b)//c for a, b, c in zip(X[i], X[j], cents)):
                    return False

        # Column orthogonality relations.
        for i in range(len(X)):
            if sum(chi[i]**2 for chi in X) != self.centralisers[i]:
                return False
            for j in range(i+1, len(X)):
                if sum(chi[i]*chi[j] for chi in X):
                    return False

        return True



class InductionTable:
    """
    Takes as input the table of a 

    """
    pass


