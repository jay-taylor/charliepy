# We describe here the basic multiplication rules for the braid monoid B
# corresponding to a finite Coxeter group W. In what follows we identify
# an element of W with an element of B by taking some (any) reduced
# expression for the element. For the relevant background material see
# Chapter 4 of [GP00]. It is likely that using Garside theory one could
# implement more efficient algorithms but we won't need this here.
#
# Deligne Normal Form
# -------------------
# For each element g of the braid monoid there exists a unique element
# alpha(g) of W of maximal possible length such that
#
#       g = alpha(g)|omega(g)
#
# with omega(g) a uniquely determined element of B. Here we use | to
# denote multiplication in the braid monoid, i.e., concatenation. The
# uniqueness of omega(g) follows from the fact that B is left
# simplifiable.
#
# Now every element g of B may be written as a unique product
#
#       g = g_1|...|g_n
#
# such that g_i = alpha(g_i|...|g_n). This is called the Deligne
# normal form of the element. Another way of writing this is that we
# have
#
#   g_1 = alpha(g),
#   g_2 = alpha(omega(g)),
#   g_3 = alpha(omega(omega(g))),
#       ...
#
# and so on. Each element of the Braid monoid is represented in
# CharLiePy by its Deligne normal form. In other words, each element of
# B is expressed as a list of elements of W which we take to be
# permutations on the roots.
#
# Multiplication
# --------------
# To multiply elements of the Braid monoid we have to explain how to
# multiply with the Deligne normal form. We start with the case of two
# elements x, w of W. There exist unique elements y, z in W such that
#
#       xw = xyz   and   l(xyz) = l(xy) + l(z) = l(x) + l(y) + l(z)
#
# with the multiplication performed inside W. With this we have
#
#       alpha(x|w) = alpha2(x, w) = xy,
#       omega(x|w) = omega2(x, w) = z.
#
# Now consider the case where g = a|b is an element of B in its Deligne
# normal form, which has two factors, and w is in W. Note that a and b
# are elements of W. Then we have
#
#   alpha(g|w) = alpha2(a, alpha(b|w)) = alpha2(a, alpha2(b, w))
#   omega(g|w) = omega2(a, alpha(b|w))omega(b|w)
#              = omega2(a, alpha2(b, w))omega2(b, w)
#
# Repeated applications of these rules allow us to consider the case of
# g|w for any g in B and w in W. Furthermore, looping again through the
# factors in the Deligne-Normal form of an element h allows us to
# compute the normal form of the product g|h for any elements g, h in B.

from . import permutat
from . import coxeter
from . import coset

import collections
import itertools

__all__ = ['BraidElement',
           'Brieskornform']

class BraidElement:
    """
    Each instance of this class represents an element in the Braid monoid
    associated to the parent group.

    """
    # Don't create a dict to save memory.
    __slots__ = ("parent", "nform")

    def __init__(self, W, elms, allperms=False):
        self.parent = W

        if allperms:
            self.nform = Deligneform(W, elms)
        else:
            g = [None]*len(elms)

            for i, x in enumerate(elms):
                if isinstance(x, int):
                    g[i] = W.permgens[x]
                else:
                    g[i] = x

            self.nform = Deligneform(W, g)

    def __repr__(self):
        return ".".join(
            "".join(str(i) for i in self.parent.convert(x, 'w'))
            for x in self.nform)

    def __mul__(self, r): # self*r
        return BraidElement(self.parent, self.nform + r.nform, True)

    def __xor__(self, i): # self^i
        if isinstance(i, int) and i > 0:
            return BraidElement(self.parent, self.nform*i, True)
        else:
            return NotImplemented


def alpha2(W, x, w):
    """
    Given two elements of a weyl group W, as permutations on the roots, this
    function returns the element alpha_2(x, w) in the Braid monoid associated to
    W.

    """
    # Recall from above that there exist unique elements y, z in W such
    # that
    #
    #   xw = xyz   and   l(xyz) = l(xy) + l(z) = l(x) + l(y) + l(z).
    #
    # We then have alpha2(x, w) = xy and omega2(x, w) = z.
    #
    # To obtain these elements we simply loop through the simple
    # reflections of W and do a replacement
    #
    #   x, w -> xs, sw
    #
    # whenever we have l(xs) > l(x) and l(sw) < l(w). The fact that what
    # we wind up with are the desired elements follows easily from the
    # cancellation law, see Theorem 1.2.5 of [GP00].
    N = W.N

    flag = True
    while flag:
        flag = False
        for i, s in enumerate(W.permgens):
            if i^w >= N and i/x < N:
                flag = True
                x = x*s
                w = s*w
                break

    return x, w


def alpha(W, g):
    """
    Given a list g of elements of W, as permutations on the roots, this function
    returns alpha(g) and omega(g), i.e., we have g[0]*...*g[k] =
    alpha(g)omega(g) in the Braid monoid and alpha(g) is maximal with this
    property.

    """
    if g:
        alph = g.pop()
        omeg = collections.deque()

        for x in reversed(g):
            alph, o = alpha2(W, x, alph)
            # We only add o if it's not the identity.
            if o:
                omeg.appendleft(o)

        return (alph, omeg)
    else:
        return (permutat.Perm(range(2*W.N)), [])


def Deligneform(W, g):
    """
    Returns the Deligne normal form of the sequence of elements g.

    """
    nform = []

    alph, omeg = alpha(W, g)
    nform.append(alph)

    while omeg:
        a, omeg = alpha(W, omeg)
        nform.append(a)

    return nform

def Brieskornform(b):
    """
    Returns the Brieskorn--Saito normal form of the element in the braid monoid.

    """
    W = b.parent
    N = W.N
    bform = []

    # See pg. 109 of [GP00] for a description of the Brieskorn--Saito
    # normal form.
    for w in b.nform:
        bform_w = []

        x = w
        while x:
            J = W.leftdescentset(x)
            bform_w.append(J)

            # We now replace x with wJ*x. The method for constructing
            # wJ*x is based on Lemma 1.1.8 of [GP00], see also Lemma
            # 1.5.2 of [GP00].
            flag = True
            while flag:
                flag = False
                for j in J:
                    if j^x >= N:
                        flag = True
                        x = W.permgens[j]*x

        bform.extend(bform_w)

    return bform

def isverygood(W, rep):
    """
    Returns True if the element w is a very good element and False otherwise. If
    the element has odd order then this simply tests if it is a good element.

    """
    # Treat the coset case first.
    if isinstance(W, coset.CoxeterCoset):
        w = W.group.convert(rep, 'p')
        phi = W.phi

        # Special case for the identity.
        if not w:
            return True

        wphi = w*W.phi
        d = wphi.order

        # Get the braid element b = w|phi(w)|phi^2(w)|...|phi^{d-1}(w).
        b = [w]
        for i in range(1, d):
            b.append(b[-1]^phi)

        b = BraidElement(W.group, b)
        bform = Brieskornform(b)

        # Group consecutive like terms together in the Brieskorn form.
        bform_mult = [[set(x), len(list(g))]
                      for x, g in itertools.groupby(bform)]

        # Check it's a good element.
        J = bform_mult[0][0]
        if bform_mult[0][1] % 2:
            return False

        for x in bform_mult[1:]:
            K, m = x
            if not (K.issubset(J) and K != J):
                return False
            if m % 2:
                return False
            J = K

        return True

    # Now the group case.
    w = W.convert(rep, 'p')
    d = w.order

    # Special case for the identity.
    if not w:
        return True

    # Get the Brieskorn normal form for w^d.
    b = BraidElement(W, [w])
    bform = Brieskornform(b^d)

    # Group consecutive like terms together in the Brieskorn form.
    bform_mult = [[set(x), len(list(g))] for x, g in itertools.groupby(bform)]

    # Check it's a good element.
    J = bform_mult[0][0]
    for x in bform_mult[1:]:
        K = x[0]
        if not (K.issubset(J) and K != J):
            return False
        J = K

    # Now check if it's a very good element.
    if d % 2:
        return True

    # Get the Brieskorn normal form for w^(d/2).
    bform2 = Brieskornform(b^(d//2))

    # Get the Brieskorn form with multiplicities.
    bform2_mult = [[set(x), len(list(g))] for x, g in itertools.groupby(bform2)]

    # Check it's a very good element.
    if len(bform_mult) != len(bform2_mult):
        return False
    for x, y in zip(bform_mult, bform2_mult):
        if x[0] != y[0] or x[1] != 2*y[1]:
            return False

    return True
