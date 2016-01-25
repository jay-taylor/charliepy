########################################################
# Data for type 2F_4 Coxeter groups and reductive groups
#

def conjclassdata(ind, **kwargs):
    """
    Returns the conjugacy class data for the coset W.phi where W is an
    irreducible Weyl group of type D_n and phi is the unique graph automorphism
    of order 2. The data is adapted from the corresponding files in GAP-Chevie.
  
    """
    # stores the data: representatives, centraliser orders, names
    repcentnam = list() 

    repcentnam.append([[], [1, 2, 1], [0], [0, 1, 0, 2, 1, 0, 3, 2, 1, 0],
                       [0, 1], [1], [0, 1, 0, 2, 1, 0, 2, 1], [0, 1, 2, 1],
                       [0, 1, 0, 2, 1, 0, 2, 1, 3, 2, 1, 0],
                       [0, 1, 0, 2, 1, 0, 2, 1, 3, 2, 1, 0, 2, 1, 3, 2, 1, 0],
                       [0, 1, 0, 2, 1, 0]])

    # redefine the reps based on the indices of the irreducible component
    repcentnam[0] = [[ind[j] for j in rep] for rep in repcentnam[0]]

    repcentnam.append([16, 8, 4, 12, 12, 8, 16, 6, 48, 96, 96])

    repcentnam.append(["2a", "8a", "4a", "24a", "24b", "8a", "8b", "12a",
                       "4b", "8c", "8d"])

    return repcentnam

