########################################################
# Data for type 2F_4 Coxeter groups and reductive groups
#

def conjclassdata(ind, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type D_n and phi is the unique graph automorphism
    of order 2. The data is adapted from the corresponding files in GAP-Chevie.
  
    """
    classes = (
        ([], 16, '2a'), ([1, 2, 1], 8, '8a'), ([0], 4, '4a'),
        ([0, 1, 0, 2, 1, 0, 3, 2, 1, 0], 12, '24a'), ([0, 1], 12, '24b'),
        ([1], 8, '8a'), ([0, 1, 0, 2, 1, 0, 2, 1], 16, '8b'),
        ([0, 1, 2, 1], 6, '12a'),
        ([0, 1, 0, 2, 1, 0, 2, 1, 3, 2, 1, 0], 48, '4b'),
        ([0, 1, 0, 2, 1, 0, 2, 1, 3, 2, 1, 0, 2, 1, 3, 2, 1, 0], 96, '8c'),
        ([0, 1, 0, 2, 1, 0], 96, '8d')
    )

    return (([ind[j] for j in cls[0]], cls[1], cls[2]) for cls in classes)

