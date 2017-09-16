from . import coxeter
from . import coset
from . import permutat
from . import utils
from . import conjclass

import numpy as np
import itertools
import collections

__all__ = ['CharacterTable',
           'Character']


########################################################################
########################################################################
##                                                                    ##
##                           Characters                               ##
##                                                                    ##

class Character(collections.namedtuple("Character",
            ["parent", "name", "a", "b", "values"])):
    # Setting slots to an empty tuple prevents a dict for the class
    # being created.
    __slots__ = () 
    def __repr__(self):
        return "Character( {}, {} )".format(self.parent,
                                            self.name)

    #@property
    #def length(self):
    #    return self.parent.size//self.centsize

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

        if isinstance(W, coxeter.CoxeterGroup):
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
                typ, inds = W.cartantype[0]
                mat = utils.getmodule(typ).chartable(len(inds))

                # Take the successive tensor product with all remaining
                # character tables.
                for (typ, inds) in W.cartantype[1:]:
                    new = utils.getmodule(typ).chartable(len(inds))
                    mat = [[a*b for a, b in itertools.product(mrow, nrow)]
                           for mrow, nrow in itertools.product(mat, new)]

                self.irreducibles = mat

        elif isinstance(W, coset.CoxeterCoset):
            # We treat the irreducible case differently because it's simpler.
            if len(W.phitype) == 1:
                typ, inds, phi = W.phitype[0]
                n = len(inds)

                # Get the correct module from the data package for the
                # irreducible factor.
                mod = utils.getmodule(typ, phi.order)

                # Irreducible character data. We don't need the
                # a-values and b-values for cosets.
                self.charnames, a, b = zip(*mod.irrchars(n))

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
                gen = (zip(*elm) for elm in itertools.product(
                       *(utils.getmodule(typ, phi.order).irrchars(len(inds))
                           for (typ, inds, phi) in W.phitype)))

                self.charnames = list(','.join(name) for name, a, b in gen)

                # Finally the values of the irreducible characters.
                (typ, inds, phi) = W.phitype[0]
                mat = utils.getmodule(typ).chartable(len(inds))

                # Take the successive tensor product with all remaining
                # character tables.
                for (typ, inds, phi) in W.phitype[1:]:
                    new = utils.getmodule(typ, phi.order).chartable(len(inds))
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

def inducecharacter(chi, H, W):
    """
    Induce the character chi of the subgroup H to W.

    """
    Hclasses = conjclass.conjugacyclasses(H)
    Wclasses = conjclass.conjugacyclasses(W)

    # Get the fusion of the conjugacy classes of H into W. This will raise an
    # error if H.embeddings[W] is undefined.
    fusion = conjclass.fusionclasses(H, W)

    # Construct the induced character I = Ind_H^W(chi).
    indchar = [0]*len(Wclasses)
    for i, j in enumerate(fusion):
        indchar[j] += (Wclasses[j].centsize*chi[i])//Hclasses[i].centsize

    return indchar

def restrictcharacter(chi, H, W):
    """
    Restrict the character chi to the subgroup H of W.

    """
    Hclasses = conjclass.conjugacyclasses(H)
    Wclasses = conjclass.conjugacyclasses(W)

    # Get the fusion of the conjugacy classes of H into W. This will raise an
    # error if H.embeddings[W] is undefined.
    fusion = conjclass.fusionclasses(H, W)

    # Construct the restricted character R = Res_H^W(chi).
    reschar = [None]*len(Hclasses)
    for i, j in enumerate(fusion):
        reschar[i] = chi[j]

    return reschar

def inducedecompose(chi, S, T):
    """
    Induce the character chi from a subgroup H to W and decompose the induced
    character into irreducible constituents.

    """
    H = S.parent
    W = T.parent

    Hclasses = S.classes
    Wclasses = T.classes

    # Get the fusion of the conjugacy classes of H into W. This will raise an
    # error if H.embeddings[W] is undefined.
    fusion = conjclass.fusionclasses(H, W)

    # Construct the induced character I = Ind_H^W(chi).
    Ichar = [0]*len(Wclasses)
    for i, j in enumerate(fusion):
        Ichar[j] += (Wclasses[j].centsize*chi[i])//Hclasses[i].centsize

    res = [0]*len(Wclasses)

    # Degree of Ind_H^W(chi).
    d = (W.size*chi[0])//H.size
    for i, psi in enumerate(T.irreducibles):
        if psi[0] <= d:
            mult = 0
            for j, (a, b) in enumerate(zip(Ichar, psi)):
                mult += (W.size*a*b)//Wclasses[j].centsize
            mult = mult//W.size
            if mult:
                res[i] = mult
                d -= mult*psi[0]
                if not d:
                    break

    return res


def inductiontable(S, T):
    """
    Takes as input the character table S of a group H and the character table
    T of a group W for which H.embeddings[W] is defined. What is produced is a
    matrix of multiplicites decomposing the induced characters of H into
    irreducible constituents.

    """
    H = S.parent
    W = T.parent
    Hclasses = S.classes
    Wclasses = T.classes

    # Get the fusion of the conjugacy classes of H into W. This will raise an
    # error if H.embeddings[W] is undefined.
    fusion = conjclass.fusionclasses(H, W)
    splitting = []

    # We need to know how the classes of W split upon restriction to H. The only
    # classes of W we consider are those whose intersection with H is non-empty.
    # If the intersection is empty then the corresponding value of the induced
    # character will be 0.
    for i in range(len(Wclasses)):
        tmp = [j for j, y in enumerate(fusion) if y == i]
        if tmp:
            splitting.append([i, tmp])


    index = W.size//H.size

    res = [[0]*len(Wclasses) for _ in range(len(Hclasses))]

    # Now get the multiplicities. Assume chi is an irreducible character of H
    # and psi is an irreducible character of W. Set I = Ind_H^W(chi) then the
    # multiplicity of psi in I is computed by the formula
    #
    # <I, psi> = \sum_{[w]} I(w)psi(w)/|C_W(w)| 
    #
    # where the sum is over all the W-conjugacy classes; note the character
    # values are integers here. The values of I are given by the formula 
    #
    # I(w) = |C_W(w)|\sum_{[x]} chi(x)/|C_H(x)|
    #
    # Here the sum is over all the H-conjugacy classes contained in the
    # W-conjugacy class containing w. This means we can write this as
    #
    # <I, psi> = 1/|W|\sum_{[w]} psi(w)(\sum_{[x]} |W|chi(x)/|C_H(x)|).
    #
    # By arranging things in this way we avoid integer division errors because
    # |C_H(x)| divides |W|.
    for i, chi in enumerate(S.irreducibles):
        # Degree of Ind_H^W(chi).
        d = index*chi[0]
        for j, psi in enumerate(T.irreducibles):
            if psi[0] > d:
                continue
            else:
                mult = sum(sum((W.size*chi[a])//Hclasses[a].centsize
                                for a in inds)*psi[k]
                                    for k, inds in splitting)//W.size
                if mult:
                    res[i][j] = mult
                    d -= mult*psi[0]
                    if not d:
                        break

    return res









