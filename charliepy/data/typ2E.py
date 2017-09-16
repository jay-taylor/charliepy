########################################################
# Data for type 2E_6 Coxeter groups and reductive groups
#
from . import rawdata as _rd

def conjclasses(inds, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type E_6 and phi is the unique graph automorphism
    of order 2. The data is adapted from the corresponding files in GAP-Chevie.

    """
    return (
        (nam, cent, [inds[i] for i in rep])
        for (nam, cent, rep) in _rd.E6_2conjclasses
    )

def conjclasses_min(inds, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type E_6 and phi is the unique graph automorphism
    of order 2. The data is adapted from the corresponding files in GAP-Chevie.

    """
    return (cls[:2] for cls in _rd.E6_2conjclasses)

def irrchars(n, labels=(), **kwargs):
    if "E6frame" in labels or "lusztig" in labels:
        return (
            ('{d}{p}'.format(d=deg, p=frame), None, None)
            for (deg, a, b, frame) in _rd.E6irrchardata
        )
    else:
        return (
            ('phi{{{d},{b}}}'.format(d=deg, b=b), None, None)
            for (deg, a, b, frame) in _rd.E6irrchardata
        )


def chartable(n, **kwargs):
    # Here we use Lusztig's preferred extension defined in [L85, IV,
    # 17.2]. In particular w.phi acts on the representation E as
    # w(-1)^a_E. This involves multiplying all rows of the character
    # table by appropriate signs.
    signs = [(-1)**a for (deg, a, b, frame) in _rd.E6irrchardata]
    return [[s*x for x in row] for s, row in zip(signs, _rd.E6chartable)]
