########################################################
# Data for type 2I_2(m) Coxeter groups
#

def conjclassdata(ind, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type I_2(m) and phi is the unique graph
    automorphism of order 2. The data is adapted from the corresponding files
    in GAP-Chevie.
    """
    m = kwargs['m']
    mpar = m%2

    # Stores the data: representatives, centraliser orders, names.
    repcentnam = [[[]], [2], [' ']]

    # Create words based on indices of the irreducible component.
    for i in range((m+1)//2):
        rep = ind[0:1] + ind[::-1]*i
        repcentnam[0] += [rep]
        repcentnam[1] += [m]
        repcentnam[2] += [''.join(str(x) for x in rep)]

    # If m is odd then the last representative is in the centre of W.
    if mpar == 1:
        repcentnam[1][-1] = 2*m

    return repcentnam
