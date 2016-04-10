import charliepy as clp
import numpy as np
import itertools
import subprocess
import collections

def coxeterclasses(W):
    """
    Returns representatives of the Coxeter classes as defined in [GP00]. Each
    such representative is a list of generators which is the lexicographically
    smallest list in its class.

    """
    # This uses Algorithm F on pg. 55 of [GP00]. We first get the powerset of
    # the generating set. Note that the elements are sorted lexicographiclly,
    # see the itertools documentation.
    P = [list(x)
         for x in itertools.chain.from_iterable(
              itertools.combinations(range(W.rank), r)
                  for r in range(W.rank + 1)
         )
    ]

    coxclasses = []

    while P:
        V = [P[0]]
        E = []
        i, k = 0, 0

        while i <= k:
            K = V[i]
            w1 = W.longestelement(K, 'p')
            for s in (i for i in range(W.rank) if i not in K):
                d = w1 * W.longestelement(K + [s], 'p')
                L = sorted(W.convert(W.permgens[i]^d, 'pw')[0] for i in K)
                if L not in V:
                    V.append(L)
                    k += 1
            i += 1

        V.sort()
        coxclasses.append(V)
        P = [x for x in P if x not in V]

    return coxclasses

def normaliser(W, J):
    """
    This returns the index [N_W(W_J) : W_J] of the parabolic subgroup in its
    normaliser.

    """
    # We start by adding to J those simple reflections which centralise J. The
    # subgroup W_L is contained in the noramliser of W_J and to determine the
    # size of [N_W(W_J) : W_J] we need only determine
    #
    #     [N_W(W_J) : W_J] = [N_W(W_J) : W_L] . [W_L : W_J]
    #
    L = J + [i for i in range(W.rank) if all(W.cartanmat[i,j] == 0 for j in J)]

    M = clp.reflectionsubgroup(W, L).size
    M = M//clp.reflectionsubgroup(W, J).size
    N = 1

    dL = W.longestelement(L, 'p') * W.longestelement(lab='p')
    k, Y = W.length(dL, 'p'), {dL}

    S = {W.permgens[i] for i in J}

    # We use Algorithm C on pg. 46 of [GP00] to construct the right cosets of
    # W_L in W and check which representatives normalise J. We treat the trivial
    # subgroup as a special case. Note this algorithm misses the identity but
    # we've already included it above.
    while k:
        # Count the number of elements of Y which normalise J.
        N += len([x for x in Y if all(s^x in S for s in S)])

        # Replace Y by those of the next lowest length.
        Y = {x*s for x in Y for s in W.permgens if W.length(x*s, 'p') < k}
        k -= 1

    return M*N





###############################################################
##                                                           ##
##                          Type G2                          ##
##                                                           ##
###############################################################
W = clp.CoxeterGroup("G", 2)

G2coxclasses = []

# We now generate the labels for the Coxeter classes. We sort based on the
# conditions: rank, type, long/short roots. Note that because we're actually
# sorting in reverse we prefix the ~ notation with a 0/1 so that long roots are
# preferred over short roots.
for orb in coxeterclasses(W):
    lab = []
    J = orb[0]
    cartantype = clp.core.cartantotype(W.cartanmat[np.ix_(J, J)])
    for typ in cartantype:
        if typ[1] and typ[0] == "A" and W.symform[J[typ[1][0]]] == 1:
            lab.append([len(typ[1]), typ[0], "0~"])
        else:
            lab.append([len(typ[1]), typ[0], "1"])
    lab.sort(reverse=True)

    # We now get rid of the 0/1 label because we'll sort forwards below.
    for typ in lab:
        typ[2] = typ[2][1:]

    lab.insert(0, len(J))
    lab.append(J)
    G2coxclasses.append([orb, lab])

G2coxclasses.sort(key=(lambda pair: pair[1]))

G2coxclasses = [
    ('+'.join([typ[2] + typ[1] + str(typ[0]) for typ in x[1:-1]]), len(o), o[0])
        for o, x in G2coxclasses
]

# Check each Coxeter class is uniquely labelled.
assert len(set(x[0] for x in G2coxclasses)) == len(G2coxclasses)





###############################################################
##                                                           ##
##                          Type F4                          ##
##                                                           ##
###############################################################
W = clp.CoxeterGroup("F", 4)

F4coxclasses = []

# We now generate the labels for the Coxeter classes. We sort based on the
# conditions: rank, type, long/short roots. Note that because we're actually
# sorting in reverse we prefix the ~ notation with a 0/1 so that long roots are
# preferred over short roots.
for orb in coxeterclasses(W):
    lab = []
    J = orb[0]
    cartantype = clp.core.cartantotype(W.cartanmat[np.ix_(J, J)])
    for typ in cartantype:
        if typ[1] and typ[0] == "A" and W.symform[J[typ[1][0]]] == 1:
            lab.append([len(typ[1]), typ[0], "0~"])
        else:
            lab.append([len(typ[1]), typ[0], "1"])
    lab.sort(reverse=True)

    # We now get rid of the 0/1 label because we'll sort forwards below.
    for typ in lab:
        typ[2] = typ[2][1:]

    lab.insert(0, len(J))
    lab.append(J)
    F4coxclasses.append([orb, lab])

F4coxclasses.sort(key=(lambda pair: pair[1]))

F4coxclasses = [
    ('+'.join([typ[2] + typ[1] + str(typ[0]) for typ in x[1:-1]]), len(o), o[0])
        for o, x in F4coxclasses
]

# Check each Coxeter class is uniquely labelled.
assert len(set(x[0] for x in F4coxclasses)) == len(F4coxclasses)






###############################################################
##                                                           ##
##                          Type E6                          ##
##                                                           ##
###############################################################
W = clp.CoxeterGroup("E", 6)

E6coxclasses = []

for orb in coxeterclasses(W):
    lab = []
    J = list(orb[0])
    cartantype = clp.core.cartantotype(W.cartanmat[np.ix_(J, J)])
    for typ in cartantype:
        lab.append([len(typ[1]), typ[0]])
    lab.sort(reverse=True)
    lab.insert(0, len(J))
    E6coxclasses.append([orb, lab])

E6coxclasses.sort(key=(lambda pair: pair[1]))

# We now make the labels of the Coxeter classes as we wish.
tmp = [None]*len(E6coxclasses)

for i, pair in enumerate(E6coxclasses):
    o, x = pair
    lab_list = [typ[1] + str(typ[0]) for typ in x[1:]]

    # We now get the unique elements in this list, note it maintains the order.
    seen = set()
    unique = [z for z in lab_list if z not in seen and not seen.add(z)]
    count = collections.Counter(lab_list)
    lab = "+".join(str(count[z]) + z if count[z] != 1 else z for z in unique)

    tmp[i] = (lab, len(o), o[0])

E6coxclasses = tmp

# Check each Coxeter class is uniquely labelled.
assert len(set(x[0] for x in E6coxclasses)) == len(E6coxclasses)




###############################################################
##                                                           ##
##                          Type E7                          ##
##                                                           ##
###############################################################
W = clp.CoxeterGroup("E", 7)

E7coxclasses = []

for orb in coxeterclasses(W):
    lab = []
    J = list(orb[0])
    cartantype = clp.core.cartantotype(W.cartanmat[np.ix_(J, J)])
    for typ in cartantype:
        lab.append([len(typ[1]), typ[0]])
    lab.sort(reverse=True)
    lab.insert(0, len(J))
    E7coxclasses.append([orb, lab])

E7coxclasses.sort(key=(lambda pair: pair[1]))

# We now make the labels of the Coxeter classes as we wish.
tmp = [None]*len(E7coxclasses)

for i, pair in enumerate(E7coxclasses):
    o, x = pair
    lab_list = [typ[1] + str(typ[0]) for typ in x[1:]]

    # We now get the unique elements in this list, note it maintains the order.
    seen = set()
    unique = [z for z in lab_list if z not in seen and not seen.add(z)]
    count = collections.Counter(lab_list)
    lab = "+".join(str(count[z]) + z if count[z] != 1 else z for z in unique)

    tmp[i] = [lab, len(o), o[0]]

E7coxclasses = tmp

# Now apply the prime labellings to the classes of the same type.
for lab in ["3A1", "A3+A1", "A5"]:
    a, b = [i for i, x in enumerate(E7coxclasses) if x[0] == lab]

    # Make sure the smallest class comes first.
    if E7coxclasses[a][2] > E7coxclasses[b][2]:
        E7coxclasses[a], E7coxclasses[b] = E7coxclasses[b], E7coxclasses[a]

    E7coxclasses[a][0] += "'"
    E7coxclasses[b][0] += "''"

# We have to create the tuples here because otherwise we can't modify the labels
# in the previous loop.
E7coxclasses = [tuple(x) for x in E7coxclasses]

# Check each Coxeter class is uniquely labelled.
assert len(set(x[0] for x in E7coxclasses)) == len(E7coxclasses)





###############################################################
##                                                           ##
##                          Type E8                          ##
##                                                           ##
###############################################################
W = clp.CoxeterGroup("E", 8)

E8coxclasses = []

for orb in coxeterclasses(W):
    lab = []
    J = list(orb[0])
    cartantype = clp.core.cartantotype(W.cartanmat[np.ix_(J, J)])
    for typ in cartantype:
        lab.append([len(typ[1]), typ[0]])
    lab.sort(reverse=True)
    lab.insert(0, len(J))
    E8coxclasses.append([orb, lab])

E8coxclasses.sort(key=(lambda pair: pair[1]))

# We now make the labels of the Coxeter classes as we wish.
tmp = [None]*len(E8coxclasses)

for i, pair in enumerate(E8coxclasses):
    o, x = pair
    lab_list = [typ[1] + str(typ[0]) for typ in x[1:]]

    # We now get the unique elements in this list, note it maintains the order.
    seen = set()
    unique = [z for z in lab_list if z not in seen and not seen.add(z)]
    count = collections.Counter(lab_list)
    lab = "+".join(str(count[z]) + z if count[z] != 1 else z for z in unique)

    tmp[i] = (lab, len(o), o[0])

E8coxclasses = tmp

# Check each Coxeter class is uniquely labelled.
assert len(set(x[0] for x in E8coxclasses)) == len(E8coxclasses)






# Write the irreducible character data to file.
def writedata(typ, classes):
    """
    Writes a function to the corresponding data file returning the Coxeter
    classes.

    """
    ## Construct the TextWrapper class.
    #wrapper = textwrap.TextWrapper(initial_indent=" "*4,
    #                               subsequent_indent=" "*4,
    #                               width=78)

    with open("../rawdata.py", 'a') as datafile:
        datafile.write(
            "{}coxeterclasses = (\n".format(typ)
        )
        
        datafile.write(" "*4)
        linelength = 4

        # Write out the list.
        for tup in classes[:-1]:
            s = str(tup) + ', '
            if linelength + len(s) <= 78:
                datafile.write(s)
                linelength += len(s)
            else:
                datafile.write("\n    " + s)
                linelength = 4 + len(s)

        # Treat the last item differently.
        tup = classes[-1]
        s = str(classes[-1]) + '\n'
        if linelength + len(s) <= 78:
            datafile.write(s)
        else:
            datafile.write("\n    " + s)

        datafile.write(")\n\n")

head_str = "{0:#<74}\n##{0:^70}##\n##{1:^70}##\n##{0:^70}##\n{0:#<74}\n"

# Write a header to the section.
with open("../rawdata.py", 'a') as datafile:
    datafile.write(head_str.format("", "Coxeter Classes"))

# Write the data to the data files.
writedata('G2', G2coxclasses)
writedata('F4', F4coxclasses)
writedata('E6', E6coxclasses)
writedata('E7', E7coxclasses)
writedata('E8', E8coxclasses)





