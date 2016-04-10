# In the labelling of [Car72] each element of the Weyl group w is
# written as a product of reflections. Below we give a list containing:
# Carter's label, a list of reflections whose product corresponds to
# that label, and the characteristic polynomial of the class.
#
# Note that most of the time the element produced by Carter is simply a
# Coxeter element in a reflection subgroup of the specified type. In
# this instance the order of the product of the reflections clearly does
# not matter. However when we are not in this case the order is more
# important. In these cases, which we refer to as exceptional, we give a
# realisation of Carter's admissible diagram inside the root system.
#
# The characteristic polynomial is given as an idiocy check to make sure
# we have matched the classes up correctly. For Carter's labelling the
# characteristic polynomial can be computed from the information given
# in §6 of the paper. Here we describe the polynomial by its coefficient
# list. This makes it easy to compare with GAP output.
#
# When comparing the classes given by the Carter labelling with those
# produced by GAP we use the fact that a conjugacy class can be uniquely
# identified by using cycle types of permutations associated to the
# class. Usually this is a representative of the class as a permutation
# but for instance in G2 and F4 one needs to look at the restricted
# partitions obtained from the action on the two orbits of the roots.

import charliepy as clp
import numpy as np
import collections
import textwrap
import G2ctab, F4ctab, E6ctab, E7ctab, E8ctab



###########################################################################
##                                                                       ##
##                         Auxiliiary Functions                          ##
##                                                                       ##
def identifycoxclass(H, C, J):
    """
    Returns representatives of the Coxeter classes as defined in [GP00]. Each
    such representative is a list of generators which is the lexicographically
    smallest list in its class.

    """
    # This uses Algorithm F on pg. 55 of [GP00].
    try:
        return C.index(J)
    except ValueError:
        pass

    V = [J]
    E = []
    i, k = 0, 0

    while i <= k:
        K = V[i]
        w1 = H.longestelement(K, 'p')
        for s in (j for j in range(H.rank) if j not in K):
            d = w1 * H.longestelement(K + [s], 'p')
            L = sorted(H.convert(H.permgens[l]^d, 'pw')[0] for l in K)
            if L not in V:
                try:
                    return C.index(L)
                except ValueError:
                    V.append(L)
                    k += 1
        i += 1

def minlengthrep(W, w):
    """
    This returns an element of minimal length in the conjugacy class of w as a
    permutation on the roots.

    """
    # Here we follow Algorithm H from pg. 83 of [GP00].
    Y = {w}
    X = set()
    while Y:
        y = next(iter(Y))
        l = W.length(y, 'p')
        X.add(y)
        Y.remove(y)
        for s in W.permgens:
            z = y^s
            lz = W.length(z, 'p')
            if lz < l:
                Y = {z}
                X = set()
                break
            elif lz == l and not z in X:
                Y.add(z)
    return next(iter(X))




###############################################################
##                                                           ##
##                          Type G2                          ##
##                                                           ##
###############################################################
W = clp.CoxeterGroup("G", 2)

G2reps = [W.convert(np.array(x), 'mp') for x in G2ctab.classreps]

# We reconstruct the orbits of the Weyl group acting on the roots.
o1 = [i for i in range(2*W.N) if W.orbitrepresentatives[i] == 0]
o2 = [i for i in range(2*W.N) if W.orbitrepresentatives[i] == 1]

# Each class is uniquely determined by the pair of cycletypes of the
# restricted permutation on each orbit.
G2labs = [(clp.restrictedperm(w, o1).cycletype(True),
           clp.restrictedperm(w, o2).cycletype(True)) for w in G2reps]

# Check that each of these labels is unique.
assert len(set(G2labs)) == len(G2labs)

G2carter = [
['A0',     [[], []],     [1, -2, 1]],
['~A1',    [[0], []],    [-1, 0, 1]],
['A1',     [[1], []],    [-1, 0, 1]],
['G2',     [[0], [1]],   [1, -1, 1]],
['A2',     [[1], [5]],   [1, 1, 1] ],
['A1+~A1', [[0], [5]],   [1, 2, 1] ]
]


# We check that the labels given above in the Carter labelling are correct.
# Namely, we check that the listed elements generate a reflection subgroup of
# the specified type.
for nam, J, dum in G2carter:
    orbreps = W.orbitrepresentatives

    K = J[0] + J[1]

    H = clp.reflectionsubgroup(W, K)

    # We get the list consisting of the rank, type, and "~" label
    # indicating if a type A subsystem is made up of long or short roots
    # and then sort the list in reverse. Note, we use 0s and 1s in front of the
    # "~" labelling to ensure that long roots are preferred over short roots
    # (we're sorting in reverse).
    lab_list = []
    for typ in H.cartantype:
        if (typ[1] and typ[0] == "A"
            and W.symform[orbreps[H.embeddings[W][typ[1][0]]]] == 1):
            lab_list.append([len(typ[1]), typ[0], "0~"])
        else:
            lab_list.append([len(typ[1]), typ[0], "1"])
    lab_list.sort(reverse=True)

    lab = "+".join(typ[2][1:] + typ[1] + str(typ[0]) for typ in lab_list)

    assert lab == nam, (lab, nam)

# We now match the Carter classes above with those produced by GAP.
G2classmatch = []

for k, x in enumerate(G2carter):
    w1 = clp.Perm(range(2*W.N))
    w2 = clp.Perm(range(2*W.N))

    # Check that the two elements are involutions.
    for i in x[1][0]:
        w1 *= W.reflection(i)
    assert w1.order <= 2, x[0]

    for i in x[1][1]:
        w2 *= W.reflection(i)
    assert w2.order <= 2, x[0]

    w = w1*w2

    param = (clp.restrictedperm(w, o1).cycletype(True),
             clp.restrictedperm(w, o2).cycletype(True))
    j = G2labs.index(param)

    # Check the characteristic polynomials match.
    assert x[2] == G2ctab.classcharpols[j], x[0]
    G2classmatch.append(j)

# Check we have a bijection between the labellings.
assert len(set(G2classmatch)) == len(G2classmatch)

# These are very good representatives from CHEVIE.
G2vgwords = [
[],
[0],
[1],
[0, 1],
[0, 1, 0, 1],
[0, 1, 0, 1, 0, 1]
]

# Now sort the conjugacy classes based on their cuspidal label.
G2coxclasses = list(clp.data.typ1G.coxeterclasses(list(range(W.rank))))
tmp = [x[0] for x in G2coxclasses]
G2sortedclasses = []

for i, j in enumerate(G2classmatch):
    # Get a permutation rep of the GAP class.
    w = W.convert(np.array(G2ctab.classreps[j]), 'mp')

    # Get a minimal length element in the class and determine J(w_min).
    # We then determine the representative of the Coxeter class
    # containing J(w_min).
    w_min = W.convert(minlengthrep(W, w), 'pw')
    J = list(set(w_min))
    a = identifycoxclass(W, tmp, J)
    G2sortedclasses.append(((i, j), (a, len(w_min), G2vgwords[i])))

    # Check that the very good representative is of minimal length.
    assert len(w_min) == len(G2vgwords[i])

G2sortedclasses.sort(key = (lambda pair: pair[1]))
G2sortedclasses = [x[0] for x in G2sortedclasses]


# Get the pairs consisting of the degree and b-value.
G2degb = [(G2ctab.irrchars[i][G2classmatch[0]], G2ctab.bvals[i])
          for i in range(len(G2reps))]

G2irrdata = [
[1, 0, 0, ""  ],
[1, 6, 6, ""  ],
[1, 1, 3, "'" ],
[1, 1, 3, "''"],
[2, 1, 1, ""  ],
[2, 1, 2, ""  ]
]

G2irrmatch = []

for x in G2irrdata:
    if x[3] == "":
        G2irrmatch.append(G2degb.index((x[0], x[2])))
    elif x[3] == "'":
        db = (x[0], x[2])
        for i in range(len(G2reps)):
            if G2degb[i] == db and G2ctab.irrchars[i][G2classmatch[1]] == 1:
                G2irrmatch.append(i)
    else:
        db = (x[0], x[2])
        for i in range(len(G2reps)):
            if G2degb[i] == db and G2ctab.irrchars[i][G2classmatch[1]] == -1:
                G2irrmatch.append(i)

assert len(set(G2irrmatch)) == len(G2irrmatch)

# Sort the irreducible characters by degree and b-value.
G2sortedirrs = sorted(enumerate(G2irrmatch),
                      key=(lambda x: [G2irrdata[x[0]][y] for y in [0, 2, 3]]))










###############################################################
##                                                           ##
##                          Type F4                          ##
##                                                           ##
###############################################################
W = clp.CoxeterGroup("F", 4)

F4reps = [W.convert(np.array(x), 'mp') for x in F4ctab.classreps]

# We reconstruct the orbits of the Weyl group acting on the roots.
o1 = [i for i in range(2*W.N) if W.orbitrepresentatives[i] == 0]
o2 = [i for i in range(2*W.N) if W.orbitrepresentatives[i] == 2]

# Each class is uniquely determined by the pair of cycletypes of the
# restricted permutation on each orbit.
F4labs = [(clp.restrictedperm(w, o1).cycletype(True),
           clp.restrictedperm(w, o2).cycletype(True)) for w in F4reps]

# Check that each of these labels is unique.
assert len(set(F4labs)) == len(F4labs)

# The exceptional diagrams are realised by the following configurations
# of roots
#
#     D4(a1):  0  - 1       F4(a1):  1 = 2
#              |    |                |   |
#              15 - 13               0 = 11
#
F4carter = [
['A0',      [[], []],             [1, -4, 6, -4, 1]],
['4A1',     [[1, 8], [15, 23]],   [1, 4, 6, 4, 1]  ],
['2A1',     [[1], [23]],          [1, 0, -2, 0, 1] ],
['A2',      [[0], [1]],           [1, -1, 0, -1, 1]],
['D4',      [[0], [1, 8, 15]],    [1, 1, 0, 1, 1]  ],
['D4(a1)',  [[0, 13], [1, 15]],   [1, 0, 2, 0, 1]  ],
['~A2',     [[2], [3]],           [1, -1, 0, -1, 1]],
['C3+A1',   [[1, 3], [2, 23]],    [1, 1, 0, 1, 1]  ],
['A2+~A2',  [[0, 2], [3, 22]],    [1, 2, 3, 2, 1]  ],
['F4(a1)',  [[0, 2], [1, 11]],    [1, -2, 3, -2, 1]],
['F4',      [[0, 2], [1, 3]],     [1, 0, -1, 0, 1] ],
['A1',      [[0], []],            [-1, 2, 0, -2, 1]],
['3A1',     [[1, 8], [15]],       [-1, -2, 0, 2, 1]],
['~A2+A1',  [[0, 3], [2]],        [-1, -1, 0, 1, 1]],
['C3',      [[1, 3], [2]],        [-1, 1, 0, -1, 1]],
['A3',      [[0, 21], [1]],       [-1, 0, 0, 0, 1] ],
['~A1',     [[2], []],            [-1, 2, 0, -2, 1]],
['2A1+~A1', [[1, 3], [23]],       [-1, -2, 0, 2, 1]],
['A2+~A1',  [[0, 3], [1]],        [-1, -1, 0, 1, 1]],
['B3',      [[0, 2], [1]],        [-1, 1, 0, -1, 1]],
['B2+A1',   [[1, 23], [2]],       [-1, 0, 0, 0, 1] ],
['A1+~A1',  [[0], [2]],           [1, 0, -2, 0, 1] ],
['B2',      [[1], [2]],           [1, -2, 2, -2, 1]],
['A3+~A1',  [[0, 21], [1, 3]],    [1, 2, 2, 2, 1]  ],
['B4',      [[0, 2], [1, 15]],    [1, 0, 0, 0, 1]  ]
]

# We check that the labels given above in the Carter labelling are correct
# whenever the element is a Coxeter element in a reflection subgroup.
for nam, J, dum in F4carter:
    orbreps = W.orbitrepresentatives

    if "(" in nam:
        continue

    K = J[0] + J[1]

    H = clp.reflectionsubgroup(W, K)

    # We get the list consisting of the rank, type, and "~" label
    # indicating if a type A subsystem is made up of long or short roots
    # and then sort the list in reverse. Note, we use 0s and 1s in front of the
    # "~" labelling to ensure that long roots are preferred over short roots
    # (we're sorting in reverse).
    lab_list = []
    for typ in H.cartantype:
        if (typ[1] and typ[0] == "A"
            and W.symform[orbreps[H.embeddings[W][typ[1][0]]]] == 1):
            lab_list.append([len(typ[1]), typ[0], "0~"])
        else:
            lab_list.append([len(typ[1]), typ[0], "1"])
    lab_list.sort(reverse=True)

    lab_list = [typ[2][1:] + typ[1] + str(typ[0]) for typ in lab_list]

    # We now get the unique elements in this list, note it maintains the order.
    seen = set()
    unique = [z for z in lab_list if z not in seen and not seen.add(z)]
    count = collections.Counter(lab_list)
    lab = "+".join(str(count[z]) + z if count[z] != 1 else z for z in unique)

    assert lab == nam, (lab, nam)

# We now match the above classes to those produced by GAP.
F4classmatch = []

for k, x in enumerate(F4carter):
    w1 = clp.Perm(range(2*W.N))
    w2 = clp.Perm(range(2*W.N))

    for i in x[1][0]:
        w1 *= W.reflection(i)
    assert w1.order <= 2, x[0]

    for i in x[1][1]:
        w2 *= W.reflection(i)
    assert w2.order <= 2, x[0]

    w = w1*w2

    param = (clp.restrictedperm(w, o1).cycletype(True),
             clp.restrictedperm(w, o2).cycletype(True))
    j = F4labs.index(param)
    assert x[2] == F4ctab.classcharpols[j]
    F4classmatch.append(j)

# Check we have a bijection between the labellings.
assert len(set(F4classmatch)) == len(F4classmatch)

# These are very good representatives from CHEVIE.
F4vgwords = [
[],
[0, 1, 0, 2, 1, 0, 2, 1, 2, 3, 2, 1, 0, 2, 1, 2, 3, 2, 1, 0, 2, 1, 2, 3],
[1, 2, 1, 2],
[1, 0],
[0, 1, 2, 1, 2, 3, 2, 1, 2, 3],
[0, 1, 0, 2, 1, 0, 2, 3, 2, 1, 2, 3],
[3, 2],
[0, 1, 0, 2, 1, 0, 2, 1, 2, 3],
[0, 1, 0, 2, 1, 0, 2, 3, 2, 1, 0, 2, 1, 2, 3, 2],
[0, 1, 0, 2, 1, 2, 3, 2],
[0, 1, 2, 3],
[0],
[1, 2, 1, 2, 3, 2, 1, 2, 3],
[0, 3, 2],
[3, 2, 1],
[1, 2, 1, 0, 2],
[2],
[0, 1, 0, 2, 1, 0, 2, 1, 2],
[1, 0, 3],
[2, 1, 0],
[1, 3, 2, 1, 2],
[0, 2],
[2, 1],
[0, 1, 0, 2, 1, 0, 2, 1, 2, 3, 2, 1, 2, 3],
[0, 1, 2, 1, 3, 2]
]

# We now check that these representatives genuinely belong to the matched class.
for i, word in zip(F4classmatch, F4vgwords):
    w = W.convert(word, 'wp')
    param = (clp.restrictedperm(w, o1).cycletype(True),
             clp.restrictedperm(w, o2).cycletype(True))
    assert param == F4labs[i]

# Now sort the conjugacy classes based on their cuspidal label.
F4coxclasses = list(clp.data.typ1F.coxeterclasses(list(range(W.rank))))
tmp = [x[0] for x in F4coxclasses]
F4sortedclasses = []

for i, j in enumerate(F4classmatch):
    # Get a permutation rep of the GAP class.
    w = W.convert(np.array(F4ctab.classreps[j]), 'mp')

    # Get a minimal length element in the class and determined J(w_min).
    # We then determine the representative of the Coxeter class
    # containing J(w_min).
    w_min = W.convert(minlengthrep(W, w), 'pw')
    J = list(set(w_min))
    a = identifycoxclass(W, tmp, J)
    F4sortedclasses.append(((i, j), (a, len(w_min), F4vgwords[i])))

    # Check that the very good representative is of minimal length.
    assert len(w_min) == len(F4vgwords[i])

F4sortedclasses.sort(key = (lambda pair: pair[1]))
F4sortedclasses = [x[0] for x in F4sortedclasses]


# Get the index of the trivial class.
ind_A0 = F4classmatch[[i for i, x in enumerate(F4carter) if x[0] == 'A0'][0]]

# Get the pairs consisting of the degree and b-value from the characters
# computed by GAP.
F4degb = [(F4ctab.irrchars[i][ind_A0], F4ctab.bvals[i])
          for i in range(len(F4reps))]

F4irrdata = [
[1,  0,  0 , "",   '_1'],
[1,  4,  12, "''", '_2'],
[1,  4,  12, "'",  '_3'],
[1,  24, 24, "",   '_4'],
[2,  1,  4,  "''", '_1'],
[2,  13, 16, "'",  '_2'],
[2,  1,  4,  "'",  '_3'],
[2,  13, 16, "''", '_4'],
[4,  4,  8,  "",   '_1'],
[9,  2,  2,  "",   '_1'],
[9,  4,  6,  "''", '_2'],
[9,  4,  6,  "'",  '_3'],
[9,  10, 10, "",   '_4'],
[6,  4,  6,  "'" , '_1'],
[6,  4,  6,  "''", '_2'],
[12, 4,  4,  "",   ''  ],
[4,  1,  1,  "",   '_2'],
[4,  4,  7,  "''", '_3'],
[4,  4,  7,  "'",  '_4'],
[4,  13, 13, "",   '_5'],
[8,  3,  3,  "''", '_1'],
[8,  9,  9,  "'",  '_2'],
[8,  3,  3,  "'",  '_3'],
[8,  9,  9,  "''", '_4'],
[16, 4,  5,  "",   ''  ]
]

# The following table gives the values of the ' and '' characters on the
# conjugacy classes A1 and B2. The first entry gives the value of the '
# character and the second entry gives the value of the '' character.
# The prime labelling is such that phi'(x) <= phi''(x) whenever x is
# contained in one of these two classes. This condition uniquely
# determines the prime labelling. Note the first column gives the pairs
# (d, b) where d is the degree of the character and b is its b-value.
#
#              A1      B2
# (1, 12)  | -1  1 | -1  -1
# (2, 4)   |  0  2 |  0   0
# (2, 16)  | -2  0 |  0   0
# (4, 7)   | -2  2 | -2  -2
# (6, 6)   |  0  0 | -2   2
# (8, 3)   |  0  4 |  0   0
# (8, 9)   | -4  0 |  0   0
# (9, 6)   | -3  3 | -1  -1

F4irrmatch = []

F4dbdict = {x : list() for x in F4degb}
for i, x in enumerate(F4degb):
    F4dbdict[x].append(i)

ind_A1 = F4classmatch[[i for i, x in enumerate(F4carter) if x[0] == 'A1'][0]]
ind_B2 = F4classmatch[[i for i, x in enumerate(F4carter) if x[0] == 'B2'][0]]

# Make sure the lists are ordered by ' then ''.
for x in [(1, 12), (2, 4), (2, 16), (4, 7), (8, 3), (8, 9), (9, 6)]:
    assert len(F4dbdict[x]) == 2
    a, b = F4dbdict[x]
    if F4ctab.irrchars[a][ind_A1] > F4ctab.irrchars[b][ind_A1]:
        F4dbdict[x] = [b, a]

assert len(F4dbdict[(6, 6)]) == 2
a, b = F4dbdict[(6, 6)]
if F4ctab.irrchars[a][ind_B2] > F4ctab.irrchars[b][ind_B2]:
    F4dbdict[(6, 6)] = [b, a]

# F4irrmatch[i] is the index of the irreducible character in the character table
# computed by GAP which corresopnds to the character at the ith index of
# F4irrdata.
for x in F4irrdata:
    if x[3] == "":
        F4irrmatch.append(F4dbdict[(x[0], x[2])][0])
    elif x[3] == "'":
        F4irrmatch.append(F4dbdict[(x[0], x[2])][0])
    else:
        F4irrmatch.append(F4dbdict[(x[0], x[2])][1])

assert len(set(F4irrmatch)) == len(F4irrmatch)

# Sort the irreducible characters by degree and b-value.
F4sortedirrs = sorted(enumerate(F4irrmatch),
                      key=(lambda x: [F4irrdata[x[0]][y] for y in [0, 2, 4]]))




###############################################################
##                                                           ##
##                          Type E6                          ##
##                                                           ##
###############################################################
W = clp.CoxeterGroup("E", 6)

E6reps = [W.convert(np.array(x), 'mp') for x in E6ctab.classreps]

E6labs = [w.cycletype(True) for w in E6reps]

# Check that each of these labels is unique.
assert len(set(E6labs)) == len(E6labs)

# The exceptional diagrams are realised by the following configurations
# of roots
#
#     E6(a1):  3 - 4 - 5      E6(a2):  3 - 4  - 5
#              |   |                   |   |    |
#              1 - 12 - 0              1 - 12 - 28 
#
# The diagrams D4(a1) and D5(a1) are easily obtained from these diagrams
# by deleting nodes.
E6carter = [
['A0',     [[], []],                 [1, -6, 15, -20, 15, -6, 1]],
['4A1',    [[0, 3], [5, 35]],        [1, 2, -1, -4, -1, 2, 1]   ],
['2A1',    [[0], [3]],               [1, -2, -1, 4, -1, -2, 1]  ],
['3A2',    [[0, 1, 4], [2, 5, 34]],  [1, 3, 6, 7, 6, 3, 1]      ],
['A2',     [[0], [2]],               [1, -3, 3, -2, 3, -3, 1]   ],
['2A2',    [[0, 4], [2, 5]],         [1, 0, 0, -2, 0, 0, 1]     ],
['D4(a1)', [[1, 4], [3, 12]],        [1, -2, 3, -4, 3, -2, 1]   ],
['A3+A1',  [[0, 3], [2, 5]],         [1, 0, -1, 0, -1, 0, 1]    ],
['A4',     [[0, 3], [2, 4]],         [1, -1, 0, 0, 0, -1, 1]    ],
['E6(a2)', [[1, 4, 28], [3, 12, 5]], [1, -1, 2, -1, 2, -1, 1]   ],
['D4',     [[1, 2, 4], [3]],         [1, -1, -1, 2, -1, -1, 1]  ],
['A5+A1',  [[0, 3, 5], [2, 4, 35]],  [1, 2, 2, 2, 2, 2, 1]      ],
['A2+2A1', [[0, 1], [2, 4]],         [1, 1, -1, -2, -1, 1, 1]   ],
['E6(a1)', [[0, 1, 4], [3, 5, 12]],  [1, 0, 0, 1, 0, 0, 1]      ],
['E6',     [[0, 3, 5], [1, 2, 4]],   [1, 1, 0, -1, 0, 1, 1]     ],
['A1',     [[0], []],                [-1, 4, -5, 0, 5, -4, 1]   ],
['3A1',    [[0, 3], [5]],            [-1, 0, 3, 0, -3, 0, 1]    ],
['A3+2A1', [[0, 3, 5], [2, 35]],     [-1, -2, -1, 0, 1, 2, 1]   ],
['A3',     [[0, 3], [2]],            [-1, 2, -1, 0, 1, -2, 1]   ],
['A2+A1',  [[0, 4], [2]],            [-1, 1, 1, 0, -1, -1, 1]   ],
['2A2+A1', [[0, 1, 4], [2, 5]],      [-1, -2, -2, 0, 2, 2, 1]   ],
['A5',     [[0, 3, 5], [2, 4]],      [-1, 0, 0, 0, 0, 0, 1]     ],
['D5',     [[1, 2, 4], [3, 5]],      [-1, 0, 1, 0, -1, 0, 1]    ],
['A4+A1',  [[0, 3, 35], [2, 4]],     [-1, -1, 0, 0, 0, 1, 1]    ],
['D5(a1)', [[1, 4], [3, 5, 12]],     [-1, 1, -1, 0, 1, -1, 1]   ]
]

# We check that the labels given above in the Carter labelling are correct
# whenever the element is a Coxeter element in a reflection subgroup.
for nam, J, dum in E6carter:
    if "(" in nam:
        continue

    K = J[0] + J[1]

    H = clp.reflectionsubgroup(W, K)

    # We get the list consisting of the rank and type.
    lab_list = []
    for typ in H.cartantype:
        lab_list.append([len(typ[1]), typ[0]])
    lab_list.sort(reverse=True)

    lab_list = [typ[1] + str(typ[0]) for typ in lab_list]

    # We now get the unique elements in this list, note it maintains the order.
    seen = set()
    unique = [z for z in lab_list if z not in seen and not seen.add(z)]
    count = collections.Counter(lab_list)
    lab = "+".join(str(count[z]) + z if count[z] != 1 else z for z in unique)

    assert lab == nam, (lab, nam)

E6classmatch = []

for k, x in enumerate(E6carter):
    w1 = clp.Perm(range(2*W.N))
    w2 = clp.Perm(range(2*W.N))

    for i in x[1][0]:
        w1 *= W.reflection(i)
    assert w1.order <= 2, x[0]

    for i in x[1][1]:
        w2 *= W.reflection(i)
    assert w2.order <= 2, x[0]

    w = w1*w2
    
    param = w.cycletype(True)
    j = E6labs.index(param)
    assert x[2] == E6ctab.classcharpols[j], x[0]
    E6classmatch.append(j)

# Check we have a bijection between the labellings.
assert len(set(E6classmatch)) == len(E6classmatch)

# These are very good representatives from CHEVIE.
E6vgwords = [
[],
[2, 3, 2, 1, 3, 2, 4, 3, 2, 1, 3, 4],
[0, 3],
[0, 2, 0, 3, 2, 0, 1, 3, 4, 3, 2, 0, 1, 3, 2, 4, 5, 4, 3, 2, 1, 3, 4, 5],
[0, 2],
[0, 2, 4, 5],
[2, 3, 2, 1, 3, 4],
[0, 3, 2, 5],
[0, 3, 2, 1],
[0, 1, 2, 0, 4, 3, 5, 4, 3, 1, 2, 3],
[2, 3, 1, 4],
[0, 1, 2, 3, 1, 2, 3, 5, 4, 3, 1, 2, 3, 4],
[0, 2, 1, 4],
[0, 2, 3, 2, 1, 3, 4, 5],
[0, 3, 5, 1, 2, 4],
[0],
[0, 3, 5],
[0, 2, 3, 2, 1, 3, 2, 4, 3, 2, 1, 3, 4],
[0, 3, 2],
[0, 2, 1],
[0, 2, 1, 4, 5],
[0, 3, 5, 2, 4],
[0, 2, 3, 1, 4],
[0, 3, 2, 1, 5],
[0, 3, 1, 4, 3, 1, 2],
]

# We now check that these representatives genuinely belong to the matched class.
for i, word in zip(E6classmatch, E6vgwords):
    w = W.convert(word, 'wp')
    assert w.cycletype(True) == E6labs[i]

# Now sort the conjugacy classes based on their cuspidal label.
E6coxclasses = list(clp.data.typ1E.coxeterclasses(list(range(W.rank))))
tmp = [x[0] for x in E6coxclasses]
E6sortedclasses = []

for i, j in enumerate(E6classmatch):
    # Get a permutation rep of the GAP class.
    w = W.convert(np.array(E6ctab.classreps[j]), 'mp')

    # Get a minimal length element in the class and determined J(w_min).
    # We then determine the representative of the Coxeter class
    # containing J(w_min).
    w_min = W.convert(minlengthrep(W, w), 'pw')
    J = list(set(w_min))
    a = identifycoxclass(W, tmp, J)
    E6sortedclasses.append(((i, j), (a, len(w_min), E6vgwords[i])))

    # Check that the very good representative is of minimal length.
    assert len(w_min) == len(E6vgwords[i])

E6sortedclasses.sort(key = (lambda pair: pair[1]))
E6sortedclasses = [x[0] for x in E6sortedclasses]

# Get the pairs consisting of the degree and b-value.
E6degb = [(E6ctab.irrchars[i][ind_A0], E6ctab.bvals[i])
          for i in range(len(E6reps))]

E6irrdata = [
[1,  0,  0,  '_p' ],
[1,  36, 36, "_p'"],
[10, 7,  9,  '_s' ],
[6,  1,  1,  '_p' ],
[6,  25, 25, "_p'"],
[20, 7,  10, '_s' ],
[15, 3,  5,  '_p' ],
[15, 15, 17, "_p'"],
[15, 3,  4,  '_q' ],
[15, 15, 16, "_q'"],
[20, 2,  2,  '_p' ],
[20, 20, 20, "_p'"],
[24, 6,  6,  '_p' ],
[24, 12, 12, "_p'"],
[30, 3,  3,  '_p' ],
[30, 15, 15, "_p'"],
[60, 7,  8,  '_s' ],
[80, 7,  7,  '_s' ],
[90, 7,  8,  '_s' ],
[60, 5,  5,  '_p' ],
[60, 11, 11, "_p'"],
[64, 4,  4,  '_p' ],
[64, 13, 13, "_p'"],
[81, 6,  6,  '_p' ],
[81, 10, 10, "_p'"]
]

ind_A0 = E6classmatch[[i for i, x in enumerate(E6carter) if x[0] == 'A0'][0]]

E6irrmatch = [E6degb.index((x[0], x[2])) for x in E6irrdata]

assert len(set(E6irrmatch)) == len(E6irrmatch)

# Sort the irreducible characters by degree and b-value.
E6sortedirrs = sorted(enumerate(E6irrmatch),
                      key=(lambda x: [E6irrdata[x[0]][y] for y in [0, 2]]))




###############################################################
##                                                           ##
##                          Type E7                          ##
##                                                           ##
###############################################################
W = clp.CoxeterGroup("E", 7)

E7reps = [W.convert(np.array(x), 'mp') for x in E7ctab.classreps]

# Each class is uniquely detetmined by the cycle type of w and w*w0.
w0 = W.longestelement(lab='p')
E7labs = [(w.cycletype(True),
           (w*w0).cycletype(True)) for w in E7reps]

# Check that each of these labels is unique.
assert len(set(E7labs)) == len(E7labs)

# The exceptional diagrams are realised by the following configurations
# of roots
#
#   E7(a1):  3 - 4  - 5 - 6    E7(a2):  60 - 3 - 4  - 5
#            |   |                           |   | 
#            1 - 14 - 0                      1 - 14 - 0
#
#   E7(a3):  3 - 4  - 5  - 6   E7(a4):      1 - 3
#            |   |    |                   /   /   \
#            1 - 14 - 36                14 - 4    54
#                                         \   \  /
#                                          36 - 5
#
# The diagrams D4(a1), D5(a1), D6(a1), D6(a2), E6(a1) and E6(a2) can
# easily be obtained from the above diagrams by deleting nodes.
E7carter = [
['A0',         [[], []],                     [-1, 7, -21, 35, -35, 21, -7, 1]],
['6A1',        [[1, 2, 4], [6, 27, 48]],     [-1, -5, -9, -5, 5, 9, 5, 1]    ],
["4A1''",      [[2, 4], [6, 62]],            [-1, -1, 3, 3, -3, -3, 1, 1]    ],
['2A1',        [[0], [1]],                   [-1, 3, -1, -5, 5, 1, -3, 1]    ],
["4A1'",       [[1, 4], [6, 62]],            [-1, -1, 3, 3, -3, -3, 1, 1]    ],
['A2',         [[0], [2]],                   [-1, 4, -6, 5, -5, 6, -4, 1]    ],
['3A2',        [[0, 1, 5], [3, 6, 61]],      [-1, -2, -3, -1, 1, 3, 2, 1]    ],
['2A2',        [[0, 4], [2, 5]],             [-1, 1, 0, 2, -2, 0, -1, 1]     ],
['D4(a1)',     [[1, 4], [3, 14]],            [-1, 3, -5, 7, -7, 5, -3, 1]    ],
["A3+A1'",     [[1, 4], [3, 6]],             [-1, 1, 1, -1, 1, -1, -1, 1]    ],
['A3+3A1',     [[1, 2, 5], [4, 6, 62]],      [-1, -3, -3, -1, 1, 3, 3, 1]    ],
['D4(a1)+2A1', [[1, 4, 6], [3, 14, 48]],     [-1, -1, -1, -1, 1, 1, 1, 1]    ],
["A3+A1''",    [[2, 4], [3, 6]],             [-1, 1, 1, -1, 1, -1, -1, 1]    ],
['A4',         [[0, 3], [2, 4]],             [-1, 2, -1, 0, 0, 1, -2, 1]     ],
['D4+2A1',     [[1, 2, 4], [3, 6, 62]],      [-1, -2, 0, 1, -1, 0, 2, 1]     ],
['D4',         [[1, 2, 4], [3]],             [-1, 2, 0, -3, 3, 0, -2, 1]     ],
['E6(a2)',     [[1, 4, 36], [3, 5, 14]],     [-1, 2, -3, 3, -3, 3, -2, 1]    ],
['A2+2A1',     [[0, 6], [2, 4]],             [-1, 0, 2, 1, -1, -2, 0, 1]     ],
['D6(a2)',     [[0, 1, 4, 60], [3, 14]],     [-1, 1, 0, -2, 2, 0, -1, 1]     ],
["A5+A1''",    [[2, 4, 6], [3, 5, 62]],      [-1, -1, 0, 0, 0, 0, 1, 1]      ],
["A5+A1'",     [[1, 4, 6], [3, 5, 62]],      [-1, -1, 0, 0, 0, 0, 1, 1]      ],
['A6',         [[0, 3, 5], [2, 4, 6]],       [-1, 0, 0, 0, 0, 0, 0, 1]       ],
['D5+A1',      [[1, 2, 4], [3, 5, 62]],      [-1, -1, 1, 1, -1, -1, 1, 1]    ],
['D6(a1)',     [[1, 4, 6], [3, 14, 5]],      [-1, 1, -1, 1, -1, 1, -1, 1]    ],
['E6(a1)',     [[0, 1, 4], [3, 5, 14]],      [-1, 1, 0, -1, 1, 0, -1, 1]     ],
['D6',         [[1, 2, 4, 6], [3, 5]],       [-1, 0, 1, 0, 0, -1, 0, 1]      ],
['A3+A2+A1',   [[1, 2, 5], [3, 6, 62]],      [-1, -2, -2, -1, 1, 2, 2, 1]    ],
['D5(a1)+A1',  [[1, 4, 62], [3, 5, 14]],     [-1, 0, 0, -1, 1, 0, 0, 1]      ],
['E6',         [[0, 3, 5], [1, 2, 4]],       [-1, 0, 1, 1, -1, -1, 0, 1]     ],
['A4+A2',      [[0, 3, 5], [4, 6, 61]],      [-1, -1, -1, 0, 0, 1, 1, 1]     ],
['7A1',        [[1, 2, 4, 6], [27, 48, 62]], [1, 7, 21, 35, 35, 21, 7, 1]    ],
['A1',         [[0], []],                    [1, -5, 9, -5, -5, 9, -5, 1]    ],
["3A1'",       [[1, 4], [6]],                [1, -1, -3, 3, 3, -3, -1, 1]    ],
['5A1',        [[1, 2, 4], [6, 27]],         [1, 3, 1, -5, -5, 1, 3, 1]      ],
["3A1''",      [[2, 4], [6]],                [1, -1, -3, 3, 3, -3, -1, 1]    ],
['D4+3A1',     [[1, 2, 4, 6], [3, 48, 62]],  [1, 4, 6, 5, 5, 6, 4, 1]        ],
['E7(a4)',     [[1, 4, 36, 54], [3, 5, 14]], [1, -2, 3, -1, -1, 3, -2, 1]    ],
['D6(a2)+A1',  [[0, 1, 4, 60], [3, 6, 14]],  [1, 1, 0, 2, 2, 0, 1, 1]        ],
['2A3+A1',     [[0, 4, 6, 60], [1, 2, 5]],   [1, 3, 5, 7, 7, 5, 3, 1]        ],
["A3+2A1''",   [[2, 4, 6], [5, 62]],         [1, 1, -1, -1, -1, -1, 1, 1]    ],
['A3',         [[0, 3], [2]],                [1, -3, 3, -1, -1, 3, -3, 1]    ],
['D4(a1)+A1',  [[1, 4, 6], [3, 14]],         [1, -1, 1, -1, -1, 1, -1, 1]    ],
["A3+2A1'",    [[1, 4, 6], [5, 62]],         [1, 1, -1, -1, -1, -1, 1, 1]    ],
['D6+A1',      [[1, 2, 4, 6], [3, 5, 62]],   [1, 2, 1, 0, 0, 1, 2, 1]        ],
['A2+A1',      [[0, 4], [2]],                [1, -2, 0, 1, 1, 0, -2, 1]      ],
['A2+3A1',     [[0, 1, 4], [2, 6]],          [1, 2, 0, -3, -3, 0, 2, 1]      ],
['A5+A2',      [[0, 1, 4, 6], [3, 5, 61]],   [1, 2, 3, 3, 3, 3, 2, 1]        ],
['D4+A1',      [[1, 2, 4], [3, 6]],          [1, 0, -2, 1, 1, -2, 0, 1]      ],
['2A2+A1',     [[0, 1, 5], [3, 6]],          [1, 1, 0, -2, -2, 0, 1, 1]      ],
["A5'",        [[1, 4, 6], [3, 5]],          [1, -1, 0, 0, 0, 0, -1, 1]      ],
["A5''",       [[2, 4, 6], [3, 5]],          [1, -1, 0, 0, 0, 0, -1, 1]      ],
['E7(a1)',     [[0, 1, 4, 6], [3, 5, 14]],   [1, 0, 0, 0, 0, 0, 0, 1]        ],
['D5',         [[1, 2, 4], [3, 5]],          [1, -1, -1, 1, 1, -1, -1, 1]    ],
['A7',         [[0, 3, 5, 52], [2, 4, 6]],   [1, 1, 1, 1, 1, 1, 1, 1]        ],
['E7',         [[1, 2, 4, 6], [0, 3, 5]],    [1, 1, 0, -1, -1, 0, 1, 1]      ],
['A4+A1',      [[0, 3, 6], [2, 4]],          [1, 0, -1, 0, 0, -1, 0, 1]      ],
['D5(a1)',     [[3, 14, 5], [1, 4]],         [1, -2, 2, -1, -1, 2, -2, 1]    ],
['A3+A2',      [[0, 3, 5], [2, 6]],          [1, 0, 0, -1, -1, 0, 0, 1]      ],
['E7(a2)',     [[0, 1, 4, 60], [3, 5, 14]],  [1, 0, -1, 1, 1, -1, 0, 1]      ],
['E7(a3)',     [[1, 4, 6, 36], [3, 5, 14]],  [1, -1, 1, 0, 0, 1, -1, 1]      ]
]

# We check that the labels given above in the Carter labelling are correct
# whenever the element is a Coxeter element in a reflection subgroup.
for nam, J, dum in E7carter:
    if "(" in nam:
        continue

    K = J[0] + J[1]

    H = clp.reflectionsubgroup(W, K)

    # We get the list consisting of the rank and type and sort it in reverse
    # order.
    lab_list = []
    for typ in H.cartantype:
        lab_list.append([len(typ[1]), typ[0]])
    lab_list.sort(reverse=True)

    lab_list = [typ[1] + str(typ[0]) for typ in lab_list]

    # We now get the unique elements in this list, note it maintains the order.
    # From this we can reconstruct the Carter label.
    seen = set()
    unique = [z for z in lab_list if z not in seen and not seen.add(z)]
    count = collections.Counter(lab_list)
    lab = "+".join(str(count[z]) + z if count[z] != 1 else z for z in unique)

    assert lab == nam.rstrip("'"), (lab, nam)

# We now match the classes to those produced by GAP.
E7classmatch = []

for k, x in enumerate(E7carter):
    w1 = clp.Perm(range(2*W.N))
    w2 = clp.Perm(range(2*W.N))

    for i in x[1][0]:
        w1 *= W.reflection(i)
    assert w1.order <= 2, x[0]

    for i in x[1][1]:
        w2 *= W.reflection(i)
    assert w2.order <= 2, x[0]

    w = w1*w2

    param = (w.cycletype(True), (w*w0).cycletype(True))
    j = E7labs.index(param)
    assert x[2] == E7ctab.classcharpols[j], x[0]
    E7classmatch.append(j)

# Check we have a bijection between the labellings.
assert len(set(E7classmatch)) == len(E7classmatch)

# We need to make sure we got the prime labellings right. The two primed
# classes can be distinguished by their class sizes. We record this
# information here, it is taken directly from the tables in [Car72].
E7primetest = [
[4,  "4A1'",     2**2 * 3**3 * 5 * 7],
[2,  "4A1''",    2**0 * 3**2 * 5 * 7],
[9,  "A3+A1'",   2**3 * 3**3 * 5 * 7],
[12, "A3+A1''",  2**4 * 3**4 * 5 * 7],
[20, "A5+A1'",   2**7 * 3**3 * 5 * 7],
[19, "A5+A1''",  2**7 * 3**2 * 5 * 7],
[32, "3A1'",     2**0 * 3**2 * 5 * 7],
[34, "3A1''",    2**2 * 3**3 * 5 * 7],
[42, "A3+2A1'",  2**4 * 3**4 * 5 * 7],
[39, "A3+2A1''", 2**3 * 3**3 * 5 * 7],
[49, "A5'",      2**7 * 3**2 * 5 * 7],
[50, "A5''",     2**7 * 3**3 * 5 * 7]
]

for x in E7primetest:
    assert E7carter[x[0]][0] == x[1], (E7carter[x[0]][0], x[1])
    assert x[2] == W.size//E7ctab.cents[E7classmatch[x[0]]], x[1]

# These are very good representatives taken from CHEVIE.
E7vgwords = [
[],
[6, 5, 6, 4, 5, 6, 3, 4, 5, 6, 1, 3, 4, 5, 6, 2, 3, 4, 5, 6, 1, 3, 4, 5, 2, 3,
    4, 1, 3, 2],
[4, 3, 4, 1, 3, 4, 2, 3, 4, 1, 3, 2],
[6, 4],
[6, 4, 1, 2],
[6, 5],
[5, 4, 5, 3, 4, 5, 1, 3, 2, 3, 4, 5, 1, 3, 4, 2, 0, 2, 3, 4, 1, 3, 2, 0],
[6, 5, 3, 1],
[4, 3, 4, 1, 3, 2],
[6, 4, 5, 1],
[6, 4, 3, 4, 1, 3, 4, 2, 3, 4, 1, 3, 2, 0],
[6, 5, 6, 4, 5, 3, 4, 1, 3, 4, 2, 3, 4, 1, 3, 2],
[6, 4, 5, 2],
[6, 5, 4, 3],
[6, 5, 4, 3, 4, 1, 3, 4, 2, 3, 4, 1, 3, 2],
[4, 3, 1, 2],
[0, 1, 2, 0, 4, 3, 5, 4, 3, 1, 2, 3],
[6, 5, 3, 0],
[6, 5, 6, 4, 5, 3, 4, 1, 3, 2],
[0, 1, 2, 3, 1, 2, 3, 5, 4, 3, 1, 2, 3, 4],
[6, 4, 1, 5, 3, 0],
[6, 5, 4, 3, 2, 0],
[6, 4, 1, 2, 3, 0],
[2, 3, 1, 2, 3, 6, 5, 4],
[5, 4, 3, 4, 1, 3, 2, 0],
[6, 5, 4, 3, 1, 2],
[6, 4, 5, 1, 2, 0],
[6, 4, 3, 4, 1, 3, 2, 0],
[5, 3, 0, 4, 2, 1],
[6, 5, 3, 1, 2, 0],
[6, 5, 6, 4, 5, 6, 3, 4, 5, 6, 1, 3, 4, 5, 6, 2, 3, 4, 5, 6, 1, 3, 4, 5, 2, 3,
    4, 1, 3, 2, 0, 2, 3, 4, 5, 6, 1, 3, 4, 5, 2, 3, 4, 1, 3, 2, 0, 2, 3, 4, 5,
    6, 1, 3, 4, 5, 2, 3, 4, 1, 3, 2, 0],
[6],
[6, 4, 1],
[6, 4, 3, 4, 1, 3, 4, 2, 3, 4, 1, 3, 2],
[6, 4, 2],
[6, 5, 6, 4, 5, 6, 3, 4, 5, 6, 1, 3, 4, 5, 6, 2, 3, 4, 5, 6, 1, 3, 4, 5, 2, 3,
    4, 1, 3, 2, 0],
[6, 5, 6, 4, 5, 6, 3, 4, 5, 1, 3, 2, 3, 4, 1, 3, 2, 0, 2, 3, 1],
[6, 5, 6, 4, 5, 3, 4, 5, 1, 3, 4, 5, 2, 3, 4, 5, 0, 2, 3, 4, 1, 3, 2],
[6, 5, 6, 4, 5, 6, 3, 4, 5, 6, 1, 3, 4, 5, 6, 2, 3, 4, 5, 6, 1, 3, 4, 5, 2, 0,
    2, 3, 4, 1, 3, 2, 0],
[5, 4, 3, 4, 1, 3, 4, 2, 3, 4, 1, 3, 2],
[6, 4, 5],
[6, 4, 3, 4, 1, 3, 2],
[6, 4, 5, 1, 2],
[6, 5, 4, 3, 4, 1, 3, 4, 2, 3, 4, 1, 3, 2, 0],
[6, 5, 3],
[6, 4, 1, 2, 0],
[6, 5, 4, 5, 3, 4, 5, 1, 3, 2, 3, 4, 5, 1, 3, 4, 2, 0, 2, 3, 4, 1, 3, 2, 0],
[6, 4, 3, 1, 2],
[6, 5, 3, 1, 0],
[6, 4, 1, 5, 3],
[6, 4, 2, 5, 3],
[6, 5, 4, 3, 4, 1, 3, 2, 0],
[5, 4, 3, 1, 2],
[6, 5, 6, 4, 5, 3, 4, 1, 3, 4, 2, 3, 4, 1, 3, 2, 0],
[6, 5, 4, 3, 1, 2, 0],
[6, 4, 5, 3, 0],
[2, 3, 1, 2, 3, 5, 4],
[6, 4, 5, 2, 0],
[6, 5, 6, 4, 5, 3, 4, 1, 3, 2, 0],
[1, 3, 1, 2, 4, 3, 1, 6, 5, 4, 3, 2, 0],
]

# We now check that these representatives genuinely belong to the matched class.
for i, word in zip(E7classmatch, E7vgwords):
    w = W.convert(word, 'wp')
    param = (w.cycletype(True), (w*w0).cycletype(True))
    assert param == E7labs[i]

# Now sort the conjugacy classes based on their cuspidal label.
E7coxclasses = list(clp.data.typ1E.coxeterclasses(list(range(W.rank))))
tmp = [x[0] for x in E7coxclasses]
E7sortedclasses = []

for i, j in enumerate(E7classmatch):
    # Get a permutation rep of the GAP class.
    w = W.convert(np.array(E7ctab.classreps[j]), 'mp')

    # Get a minimal length element in the class and determined J(w_min).
    # We then determine the representative of the Coxeter class
    # containing J(w_min).
    w_min = W.convert(minlengthrep(W, w), 'pw')
    J = list(set(w_min))
    a = identifycoxclass(W, tmp, J)
    E7sortedclasses.append(((i, j), (a, len(w_min), E7vgwords[i])))

    # Check that the very good representative is of minimal length.
    assert len(w_min) == len(E7vgwords[i])

E7sortedclasses.sort(key = (lambda pair: pair[1]))
E7sortedclasses = [x[0] for x in E7sortedclasses]



# Get the pairs consisting of the degree and b-value.
E7degb = [(E7ctab.irrchars[i][ind_A0], E7ctab.bvals[i])
          for i in range(len(E7reps))]

E7irrdata = [
[1,   0,  0,  '_a' ],
[1,   63, 63, "_a'"],
[7,   46, 46, '_a' ],
[7,   1,  1,  "_a'"],
[15,  25, 28, '_a' ],
[15,  4,  7,  "_a'"],
[21,  3,  6,  '_a' ],
[21,  30, 33, "_a'"],
[21,  36, 36, '_b' ],
[21,  3,  3,  "_b'"],
[27,  2,  2,  '_a' ],
[27,  37, 37, "_a'"],
[35,  16, 22, '_a' ],
[35,  7,  13, "_a'"],
[35,  3,  4,  '_b' ],
[35,  30, 31, "_b'"],
[56,  30, 30, '_a' ],
[56,  3,  3,  "_a'"],
[70,  16, 18, '_a' ],
[70,  7,  9,  "_a'"],
[84,  10, 12, '_a' ],
[84,  13, 15, "_a'"],
[105, 25, 26, '_a' ],
[105, 4,  5,  "_a'"],
[105, 6,  6,  '_b' ],
[105, 21, 21, "_b'"],
[105, 12, 12, '_c' ],
[105, 15, 15, "_c'"],
[120, 4,  4,  '_a' ],
[120, 25, 25, "_a'"],
[168, 6,  6,  '_a' ],
[168, 21, 21, "_a'"],
[189, 8,  10, '_a' ],
[189, 15, 17, "_a'"],
[189, 22, 22, '_b' ],
[189, 5,  5,  "_b'"],
[189, 20, 20, '_c' ],
[189, 7,  7,  "_c'"],
[210, 6,  6,  '_a' ],
[210, 21, 21, "_a'"],
[210, 10, 10, '_b' ],
[210, 13, 13, "_b'"],
[216, 15, 16, '_a' ],
[216, 8,  9,  "_a'"],
[280, 16, 18, '_a' ],
[280, 7,  9,  "_a'"],
[280, 7,  8,  '_b' ],
[280, 16, 17, "_b'"],
[315, 16, 16, '_a' ],
[315, 7,  7,  "_a'"],
[336, 13, 14, '_a' ],
[336, 10, 11, "_a'"],
[378, 14, 14, '_a' ],
[378, 9,  9,  "_a'"],
[405, 8,  8,  '_a' ],
[405, 15, 15, "_a'"],
[420, 10, 10, '_a' ],
[420, 13, 13, "_a'"],
[512, 11, 12, '_a' ],
[512, 11, 11, "_a'"]
]

ind_A0 = E7classmatch[[i for i, x in enumerate(E7carter) if x[0] == 'A0'][0]]

E7irrmatch = [E7degb.index((x[0], x[2])) for x in E7irrdata]

assert len(set(E7irrmatch)) == len(E7irrmatch)

tmp = [[E7ctab.irrchars[E7irrmatch[i]][E7classmatch[j]]
           for j in range(len(E7reps))] for i in range(len(E7reps))]

assert tmp == clp.data.typ1E.chartable(7)

# Sort the irreducible characters by degree and b-value.
E7sortedirrs = sorted(enumerate(E7irrmatch),
                      key=(lambda x: [E7irrdata[x[0]][y] for y in [0, 2]]))





###############################################################
##                                                           ##
##                          Type E8                          ##
##                                                           ##
###############################################################
W = clp.CoxeterGroup("E", 8)

E8reps = [W.convert(np.array(x), 'mp') for x in E8ctab.classreps]

# Construct the longest element of the Weyl group.
#w0 = W.convert(-np.identity(8, dtype='int'), 'mp')
E8labs = [w.cycletype(True) for w in E8reps]

# Check that each of these labels is unique.
assert len(set(E8labs)) == len(E8labs)

# The exceptional diagrams are realised by the following configurations
# of roots
#
#   2D4(a1):  3 - 4       7 - 96
#             |   |       |   |
#             1 - 16      6 - 73
#
#   D8(a1):  3 - 4  - 5 - 6 - 7 - 96
#            |   |
#            1 - 16
#
#   D8(a2):      4  - 5  - 6 - 7 - 96
#                |    |
#            1 - 3  - 31
#
#   D8(a3):          5  - 6 - 7 - 96
#                    |    |
#            1 - 3 - 4  - 47
#
#   E8(a1):  3 - 4  - 5 - 6 - 7      E8(a2):       3 - 4  - 5 - 6
#            |   |                                 |   | 
#            1 - 16 - 0                       41 - 1 - 16 - 0
#
#   E8(a3):  83 - 3 - 4  - 5         E8(a4):  3  - 4  - 5  - 6 - 7
#                 |   |                       |    |    |
#            33 - 1 - 16 - 0                  1  - 16 - 43
#
#   E8(a5):       3 - 4  - 5  - 6    E8(a6):  83 - 3 - 4  - 5
#                 |   |    |                  |    |   |    |
#            42 - 1 - 16 - 43                 34 - 1 - 16 - 43
#
#   E8(a7):      1 - 3               E8(a8):      37 ---  54
#              /   /   \                        / |     / |
#            16 - 4    70 - 7                 3 --|-- 4   |
#              \   \  /                       |   25 -|-- 46
#               43 - 5                        | /     | /
#                                             2 ----  16
#   
# The diagrams D4(a1), D5(a1), D6(a1), D6(a2), D7(a1), D7(a2), E6(a1),
# E6(a2), E7(a1), E7(a2), E7(a3), and E7(a4) can easily be obtained from
# the above diagrams by deleting nodes.
E8carter = [
['A0',        [[], []],                    [1,-8,28,-56,70,-56,28,-8,1]],
['8A1',       [[1,2,4,6], [31,60,96,119]], [1,8,28,56,70,56,28,8,1]    ],
["4A1'",      [[2,4], [6,96]],             [1,0,-4,0,6,0,-4,0,1]       ],
['2A1',       [[0], [1]],                  [1,-4,4,4,-10,4,4,-4,1]     ],
['6A1',       [[1,2,4], [6,31,60]],        [1,4,4,-4,-10,-4,4,4,1]     ],
['2D4(a1)',   [[1,4,6,96], [3,7,16,73]],   [1,0,4,0,6,0,4,0,1]         ],
["4A1''",     [[2,4], [6,119]],            [1,0,-4,0,6,0,-4,0,1]       ],
['A2',        [[0], [2]],                  [1,-5,10,-11,10,-11,10,-5,1]],
['D4+4A1',    [[1,2,4,6], [3,60,96,119]],  [1,5,10,11,10,11,10,5,1]    ],
['4A2',       [[0,1,4,7], [2,5,62,118]],   [1,4,10,16,19,16,10,4,1]    ],
['E8(a8)',    [[2,4,37,46], [3,16,25,54]], [1,-4,10,-16,19,-16,10,-4,1]],
['3A2',       [[0,1,4], [2,5,62]],         [1,1,1,-2,-2,-2,1,1,1]      ],
['E7(a4)+A1', [[1,4,43,70], [3,5,16,119]], [1,-1,1,2,-2,2,1,-1,1]      ],
['2A2',       [[0,4], [2,5]],              [1,-2,1,-2,4,-2,1,-2,1]     ],
['2D4',       [[1,2,4,7], [3,6,60,96]],    [1,2,1,2,4,2,1,2,1]         ],
['D4(a1)',    [[1,4], [3,16]],             [1,-4,8,-12,14,-12,8,-4,1]  ],
['2A3+2A1',   [[2,4,6,96], [3,7,95,111]],  [1,4,8,12,14,12,8,4,1]      ],
["2A3'",      [[2,4,7], [3,6,96]],         [1,0,0,0,-2,0,0,0,1]        ],
['A3+A1',     [[0,3], [2,5]],              [1,-2,0,2,-2,2,0,-2,1]      ],
['A3+3A1',    [[1,2,6], [4,7,117]],        [1,2,0,-2,-2,-2,0,2,1]      ],
['D8(a3)',    [[1,4,6,96], [3,5,7,47]],    [1,0,0,0,2,0,0,0,1]         ],
["2A3''",     [[2,4,6], [3,7,117]],        [1,0,0,0,-2,0,0,0,1]        ],
['A4',        [[0,3], [2,4]],              [1,-3,3,-1,0,-1,3,-3,1]     ],
['D6+2A1',    [[1,2,4,6], [3,5,96,119]],   [1,3,3,1,0,1,3,3,1]         ],
['2A4',       [[0,3,5,7], [1,2,6,116]],    [1,2,3,4,5,4,3,2,1]         ],
['E8(a6)',    [[1,4,43,83], [3,5,16,34]],  [1,-2,3,-4,5,-4,3,-2,1]     ],
['A2+4A1',    [[0,1,4], [2,6,119]],        [1,3,2,-3,-6,-3,2,3,1]      ],
['D4',        [[1,2,4], [3]],              [1,-3,2,3,-6,3,2,-3,1]      ],
['E6(a2)+A2', [[1,4,7,43], [3,5,16,119]],  [1,0,2,0,3,0,2,0,1]         ],
['A2+2A1',    [[0,1], [2,4]],              [1,-1,-2,1,2,1,-2,-1,1]     ],
['D4+2A1',    [[1,2,4], [3,6,119]],        [1,1,-2,-1,2,-1,-2,1,1]     ],
['E8(a3)',    [[0,1,4,83], [3,5,16,33]],   [1,0,-2,0,3,0,-2,0,1]       ],
['E6(a2)',    [[1,4,43], [3,5,16]],        [1,-3,5,-6,6,-6,5,-3,1]     ],
['A5+A2+A1',  [[0,5,7,115], [1,2,4,6]],    [1,3,5,6,6,6,5,3,1]         ],
["A5+A1'",    [[2,4,6], [3,5,96]],         [1,0,-1,0,0,0,-1,0,1]       ],
['D4+A2',     [[1,2,4,6], [3,7]],          [1,0,-1,0,0,0,-1,0,1]       ],
['2A2+2A1',   [[0,1,4], [2,5,7]],          [1,2,1,-2,-4,-2,1,2,1]      ],
['D6(a2)',    [[1,4,6,31], [3,5]],         [1,-2,1,2,-4,2,1,-2,1]      ],
['D8(a1)',    [[1,4,6,96], [3,5,7,16]],    [1,0,1,0,0,0,1,0,1]         ],
["A5+A1''",   [[2,4,6], [3,5,119]],        [1,0,-1,0,0,0,-1,0,1]       ],
['A6',        [[0,3,5], [2,4,6]],          [1,-1,0,0,0,0,0,-1,1]       ],
['D8',        [[1,2,4,6,119], [3,5,7]],    [1,1,0,0,0,0,0,1,1]         ],
['D5+A1',     [[1,2,4], [3,5,7]],          [1,0,-2,0,2,0,-2,0,1]       ],
['D6(a1)',    [[1,4,6], [3,5,16]],         [1,-2,2,-2,2,-2,2,-2,1]     ],
['A7+A1',     [[2,4,6,96], [3,5,7,111]],   [1,2,2,2,2,2,2,2,1]         ],
['E6(a1)',    [[0,1,4], [3,5,16]],         [1,-2,1,1,-2,1,1,-2,1]      ],
['E7+A1',     [[1,2,4,6], [0,3,5,119]],    [1,2,1,-1,-2,-1,1,2,1]      ],
['A8',        [[0,3,5,7], [2,4,6,119]],    [1,1,1,1,1,1,1,1,1]         ],
['E8(a4)',    [[1,4,6,43], [3,5,7,16]],    [1,-1,1,-1,1,-1,1,-1,1]     ],
['A4+2A1',    [[2,4,6], [1,5,7]],          [1,1,-1,-1,0,-1,-1,1,1]     ],
['D6',        [[1,2,4,6], [3,5]],          [1,-1,-1,1,0,1,-1,-1,1]     ],
['E8(a2)',    [[0,1,4,6], [3,5,16,41]],    [1,0,-1,0,1,0,-1,0,1]       ],
['D4(a1)+A2', [[1,4,6], [3,16,7]],         [1,-1,2,-3,2,-3,2,-1,1]     ],
['D5(a1)+A3', [[1,4,7,73], [3,5,16,96]],   [1,1,2,3,2,3,2,1,1]         ],
['E6+A2',     [[0,3,5,7], [1,2,4,119]],    [1,2,2,0,-1,0,2,2,1]        ],
['E8(a7)',    [[1,4,43,70], [3,5,7,16]],   [1,-2,2,0,-1,0,2,-2,1]      ],
['E6',        [[0,3,5], [1,2,4]],          [1,-1,-1,0,2,0,-1,-1,1]     ],
['E7(a2)+A1', [[0,1,4,83], [3,5,7,16]],    [1,1,-1,0,2,0,-1,1,1]       ],
['A3+A2+A1',  [[0,3,6], [1,5,7]],          [1,1,0,-1,-2,-1,0,1,1]      ],
['D5(a1)+A1', [[1,4,7], [3,5,16]],         [1,-1,0,1,-2,1,0,-1,1]      ],
['E8(a1)',    [[0,1,4,6], [3,5,7,16]],     [1,0,0,0,-1,0,0,0,1]        ],
['A4+A2',     [[0,4,6], [2,5,7]],          [1,0,0,-1,0,-1,0,0,1]       ],
['D8(a2)',    [[1,4,6,31,96], [3,5,7]],    [1,0,0,1,0,1,0,0,1]         ],
['E8(a5)',    [[1,4,6,43], [3,5,16,42]],   [1,-1,0,1,-1,1,0,-1,1]      ],
['E8',        [[0,3,5,7], [1,2,4,6]],      [1,1,0,-1,-1,-1,0,1,1]      ],
['A1',        [[0], []],                   [-1,6,-14,14,0,-14,14,-6,1] ],
['7A1',       [[1,2,4,6], [31,60,119]],    [-1,-6,-14,-14,0,14,14,6,1] ],
['3A1',       [[1,2], [4]],                [-1,2,2,-6,0,6,-2,-2,1]     ],
['5A1',       [[1,2,4], [6,31]],           [-1,-2,2,6,0,-6,-2,2,1]     ],
['A3',        [[0,3], [2]],                [-1,4,-6,4,0,-4,6,-4,1]     ],
['A3+4A1',    [[2,4,6,60], [3,96,119]],    [-1,-4,-6,-4,0,4,6,4,1]     ],
["A3+2A1'",   [[2,4,6], [3,96]],           [-1,0,2,0,0,0,-2,0,1]       ],
['D4(a1)+A1', [[1,4,6], [3,16]],           [-1,2,-2,2,0,-2,2,-2,1]     ],
['2A3+A1',    [[0,3,6,119], [1,4,7]],      [-1,-2,-2,-2,0,2,2,2,1]     ],
['D4(a1)+A3', [[1,4,6,96], [3,16,7]],      [-1,0,-2,0,0,0,2,0,1]       ],
["A3+2A1''",  [[2,4,6], [3,119]],          [-1,0,2,0,0,0,-2,0,1]       ],
['A2+A1',     [[0,4], [2]],                [-1,3,-2,-1,0,1,2,-3,1]     ],
['D4+3A1',    [[1,2,4,6], [3,60,96]],      [-1,-3,-2,1,0,-1,2,3,1]     ],
['3A2+A1',    [[0,1,4,7], [2,5,118]],      [-1,-3,-5,-4,0,4,5,3,1]     ],
['E7(a4)',    [[1,4,43,70], [3,5,16]],     [-1,3,-5,4,0,-4,5,-3,1]     ],
['A2+3A1',    [[1,2,6], [4,7]],            [-1,-1,2,3,0,-3,-2,1,1]     ],
['D4+A1',     [[1,2,4], [3,6]],            [-1,1,2,-3,0,3,-2,-1,1]     ],
['2A2+A1',    [[0,4,7], [2,5]],            [-1,0,1,2,0,-2,-1,0,1]      ],
['D6(a2)+A1', [[1,4,6,31], [3,5,96]],      [-1,0,1,-2,0,2,-1,0,1]      ],
['A5+A2',     [[0,4,6,119], [2,5,7]],      [-1,-1,-1,0,0,0,1,1,1]      ],
['E6(a2)+A1', [[1,4,43,7], [3,5,16]],      [-1,1,-1,0,0,0,1,-1,1]      ],
['A5',        [[0,3,5], [2,4]],            [-1,2,-1,0,0,0,1,-2,1]      ],
['A5+2A1',    [[1,4,6,119], [2,5,7]],      [-1,-2,-1,0,0,0,1,2,1]      ],
['D5',        [[1,2,4], [3,5]],            [-1,2,0,-2,0,2,0,-2,1]      ],
['D5+2A1',    [[1,2,4,7], [3,5,73]],       [-1,-2,0,2,0,-2,0,2,1]      ],
["A7'",       [[2,4,6,96], [3,5,7]],       [-1,0,0,0,0,0,0,0,1]        ],
["A7''",      [[1,4,6,96], [3,5,7]],       [-1,0,0,0,0,0,0,0,1]        ],
['A4+A1',     [[0,3,6], [2,4]],            [-1,1,1,-1,0,1,-1,-1,1]     ],
['D6+A1',     [[1,2,4,6], [3,5,119]],      [-1,-1,1,1,0,-1,-1,1,1]     ],
['A3+A2+2A1', [[0,3,5,81], [2,6,119]],     [-1,-3,-4,-3,0,3,4,3,1]     ],
['D5(a1)',    [[1,4], [3,5,16]],           [-1,3,-4,3,0,-3,4,-3,1]     ],
['A3+A2',     [[0,3,5], [2,6]],            [-1,1,0,1,0,-1,0,-1,1]      ],
['D4+A3',     [[1,2,4,7], [3,6,119]],      [-1,-1,0,-1,0,1,0,1,1]      ],
['D5(a1)+A2', [[1,4,7], [3,5,16,96]],      [-1,0,-1,0,0,0,1,0,1]       ],
['D7',        [[1,2,4,6], [3,5,7]],        [-1,0,1,0,0,0,-1,0,1]       ],
['E6+A1',     [[0,3,5,7], [1,2,4]],        [-1,-1,1,2,0,-2,-1,1,1]     ],
['E7(a2)',    [[0,1,4,83], [3,5,16]],      [-1,1,1,-2,0,2,-1,-1,1]     ],
['A6+A1',     [[2,4,6,231], [3,5,7]],      [-1,-1,0,0,0,0,0,1,1]       ],
['E7(a1)',    [[0,1,4,6], [3,5,16]],       [-1,1,0,0,0,0,0,-1,1]       ],
['E6(a1)+A1', [[0,1,4,7], [3,5,16]],       [-1,0,1,-1,0,1,-1,0,1]      ],
['E7',        [[0,3,5], [1,2,4,6]],        [-1,0,1,1,0,-1,-1,0,1]      ],
['A4+A3',     [[0,3,6,119], [2,4,7]],      [-1,-1,-1,-1,0,1,1,1,1]     ],
['D7(a1)',    [[1,4,6], [3,5,7,16]],       [-1,1,-1,1,0,-1,1,-1,1]     ],
['D5+A2',     [[1,2,4,7], [3,5,119]],      [-1,-1,0,1,0,-1,0,1,1]      ],
['D7(a2)',    [[1,4,6,31], [3,5,7]],       [-1,1,0,-1,0,1,0,-1,1]      ],
['A4+A2+A1',  [[0,3,5,119], [1,2,6]],      [-1,-2,-2,-1,0,1,2,2,1]     ],
['E7(a3)',    [[1,4,6,43], [3,5,16]],      [-1,2,-2,1,0,-1,2,-2,1]     ]
]

# We check that the labels given above in the Carter labelling are correct
# whenever the element is a Coxeter element in a reflection subgroup.
for nam, J, dum in E8carter:
    if "(" in nam:
        continue

    K = J[0] + J[1]

    H = clp.reflectionsubgroup(W, K)

    # We get the list consisting of the rank and type.
    lab_list = []
    for typ in H.cartantype:
        lab_list.append([len(typ[1]), typ[0]])
    lab_list.sort(reverse=True)

    lab_list = [typ[1] + str(typ[0]) for typ in lab_list]

    # We now get the unique elements in this list, note it maintains the order.
    seen = set()
    unique = [z for z in lab_list if z not in seen and not seen.add(z)]
    count = collections.Counter(lab_list)
    lab = "+".join(str(count[z]) + z if count[z] != 1 else z for z in unique)

    assert lab == nam.rstrip("'"), (lab, nam)

# We now match the classes with those produced by GAP.
E8classmatch = []

for k, x in enumerate(E8carter):
    w1 = clp.Perm(range(2*W.N))
    w2 = clp.Perm(range(2*W.N))

    for i in x[1][0]:
        w1 *= W.reflection(i)
    assert w1.order <= 2, x[0]

    for i in x[1][1]:
        w2 *= W.reflection(i)
    assert w2.order <= 2, x[0]

    w = w1*w2

    param = w.cycletype(True)
    j = E8labs.index(param)
    assert x[2] == E8ctab.classcharpols[j], x[0]
    E8classmatch.append(j)

## Check we have a bijection between the labellings.
assert len(set(E8classmatch)) == len(E8classmatch)

# We need to make sure we got the prime labellings right. The two primed
# classes can be distinguished by their class sizes. We record this
# information here, it is taken directly from the tables in [Car72].
E8primetest = [
[2,  "4A1'",     2**1  * 3**2 * 5**2 * 7],
[6,  "4A1''",    2**3  * 3**4 * 5**2 * 7],
[17, "2A3'",     2**4  * 3**5 * 5**2 * 7],
[21, "2A3''",    2**7  * 3**5 * 5**2 * 7],
[34, "A5+A1'",   2**10 * 3**2 * 5**2 * 7],
[39, "A5+A1''",  2**10 * 3**4 * 5**2 * 7],
[71, "A3+2A1'",  2**5  * 3**3 * 5**2 * 7],
[75, "A3+2A1''", 2**6  * 3**5 * 5**2 * 7],
[90, "A7'",      2**8  * 3**5 * 5**2 * 7],
[91, "A7''",     2**10 * 3**5 * 5**2 * 7],
]

for x in E8primetest:
    assert E8carter[x[0]][0] == x[1], (x[0], E8carter[x[0]][0], x[1])
    assert x[2] == W.size//E8ctab.cents[E8classmatch[x[0]]], x[1]

# These are very good representatives taken from CHEVIE.
E8vgwords = [
[],
[7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 1, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6,
    7, 1, 3, 4, 5, 6, 2, 3, 4, 5, 1, 3, 4, 2, 3, 1, 0, 2, 3, 4, 5, 6, 7, 1, 3,
    4, 5, 6, 2, 3, 4, 5, 1, 3, 4, 2, 3, 1, 0, 2, 3, 4, 5, 6, 7, 1, 3, 4, 5, 6,
    2, 3, 4, 5, 1, 3, 4, 2, 3, 1, 0, 2, 3, 4, 5, 6, 7, 1, 3, 4, 5, 2, 3, 4, 1,
    3, 2, 0, 2, 3, 4, 5, 6, 1, 3, 4, 5, 2, 3, 4, 1, 3, 2, 0],
[4, 3, 4, 1, 3, 4, 2, 3, 4, 1, 3, 2],
[5, 0],
[6, 5, 6, 4, 5, 6, 3, 4, 5, 6, 1, 3, 4, 5, 6, 2, 3, 4, 5, 6, 1, 3, 4, 5, 2, 3,
    4, 1, 3, 2],
[7, 6, 5, 6, 4, 5, 3, 4, 5, 6, 1, 3, 4, 5, 2, 3, 4, 5, 6, 7, 1, 3, 4, 5, 2, 3,
    4, 0, 2, 3, 4, 5, 6, 7, 1, 3, 4, 5, 6, 2, 3, 4, 5, 1, 3, 4, 2, 3, 1, 0, 2,
    3, 4, 5, 6, 1, 3, 4, 2, 0],
[6, 4, 1, 2],
[5, 6],
[0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 1, 2, 0,
    3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3,
    2, 0, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7],
[7, 6, 7, 5, 6, 4, 5, 6, 3, 1, 3, 4, 2, 3, 4, 5, 6, 7, 1, 3, 4, 5, 2, 3, 4, 1,
    0, 2, 3, 4, 5, 6, 7, 1, 3, 4, 5, 6, 2, 3, 4, 5, 1, 3, 4, 2, 3, 0, 2, 3, 4,
    5, 6, 7, 1, 3, 4, 5, 6, 2, 3, 4, 5, 1, 3, 4, 2, 3, 1, 0, 2, 3, 4, 5, 1, 3,
    4, 2, 3, 0],
[4, 5, 3, 2, 3, 4, 5, 6, 7, 1, 3, 4, 5, 2, 3, 1, 0, 2, 3, 4, 5, 6, 7, 1, 3, 4,
    5, 6, 2, 3, 4, 1, 0, 2, 3, 4, 5, 6, 7, 1],
[5, 4, 1, 3, 4, 5, 2, 3, 4, 1, 3, 2, 0, 2, 3, 4, 5, 1, 3, 4, 2, 3, 0, 2],
[1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 0, 3, 4, 5, 7, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3,
    1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6],
[6, 7, 1, 3],
[1, 3, 1, 4, 3, 1, 5, 4, 3, 1, 2, 3, 4, 5, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7, 6,
    5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 7],
[3, 1, 3, 2, 3, 4],
[0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 1, 2, 0,
    3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3,
    2, 0, 6, 7, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7],
[1, 2, 3, 1, 2, 3, 5, 4, 6, 5, 4, 3, 1, 2, 3, 4],
[2, 6, 5, 7],
[0, 1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 3, 4, 6],
[3, 4, 5, 6, 2, 3, 4, 1, 3, 2, 0, 2, 3, 4, 5, 6, 7, 1, 3, 4, 2, 3, 1, 0, 2, 3,
    4, 5, 6, 7],
[6, 5, 7, 0, 3, 2],
[4, 5, 2, 3],
[2, 3, 2, 4, 3, 2, 5, 4, 3, 2, 7, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3,
    2, 6, 5, 4, 3, 1],
[6, 4, 5, 1, 3, 4, 2, 3, 4, 5, 6, 7, 1, 3, 4, 5, 6, 0, 2, 3, 4, 5, 6, 7, 1, 3,
    4, 5, 6, 2, 3, 4, 1, 3, 2, 0, 2, 3, 4, 5, 6, 7, 1, 3, 4, 5, 2, 3],
[7, 6, 5, 6, 4, 3, 4, 5, 1, 3, 2, 3, 4, 5, 6, 0, 2, 3, 4, 5, 1, 3, 4, 2],
[7, 6, 4, 3, 4, 1, 3, 4, 2, 3, 4, 1, 3, 2],
[4, 3, 1, 2],
[0, 1, 2, 0, 3, 1, 2, 0, 3, 4, 3, 1, 2, 0, 3, 2, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3,
    1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7, 6],
[6, 3, 4, 0],
[4, 3, 4, 1, 3, 4, 2, 3, 4, 5, 6, 1, 3, 2],
[7, 6, 5, 6, 4, 1, 3, 2, 3, 4, 1, 3, 0, 2, 3, 4, 5, 1, 3, 4],
[2, 3, 1, 2, 4, 3, 1, 2, 0, 3, 4, 5],
[3, 2, 4, 3, 1, 2, 3, 4, 5, 4, 3, 1, 2, 3, 4, 5, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6,
    7, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 7],
[1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 0, 3, 4, 5],
[7, 6, 4, 1, 3, 2],
[6, 7, 4, 1, 2, 0],
[5, 4, 5, 6, 3, 1, 3, 4, 2, 3],
[4, 3, 5, 4, 3, 1, 2, 6, 5, 4, 3, 1, 2, 7, 6, 5, 4, 3, 1, 2, 0, 3],
[5, 3, 4, 1, 6, 0],
[7, 4, 5, 6, 1, 3],
[1, 6, 5, 4, 3, 1, 2, 3, 7, 6, 5, 4, 3, 1, 2, 0, 3, 2],
[7, 4, 5, 1, 2, 3],
[1, 3, 1, 2, 3, 4, 5, 6],
[1, 3, 1, 4, 3, 1, 2, 5, 4, 3, 1, 2, 3, 4, 5, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3,
    1, 5, 4, 3, 2, 0, 6, 7],
[4, 2, 3, 4, 5, 1, 3, 0],
[2, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 7, 6, 5],
[0, 2, 3, 2, 0, 4, 3, 1, 2, 3, 4, 7, 6, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4,
    6, 5],
[1, 2, 3, 2, 0, 4, 3, 1, 2, 3, 4, 5, 6, 7],
[5, 6, 7, 4, 1, 0],
[6, 4, 5, 1, 2, 3],
[5, 6, 7, 4, 5, 6, 3, 4, 2, 3, 1, 0],
[6, 7, 3, 1, 2, 3, 4, 1],
[0, 2, 0, 3, 4, 3, 2, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 6, 5, 4, 3, 1, 2, 0, 3,
    2, 4, 3, 1, 5, 4, 3, 2, 0, 6, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7],
[0, 1, 2, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 4, 5, 4, 3, 1, 2, 3, 4, 5, 7, 6],
[0, 1, 2, 0, 3, 1, 4, 3, 1, 2, 0, 3, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7],
[5, 3, 0, 4, 2, 1],
[0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 2, 0, 5, 4, 3, 1, 2, 3, 4, 5, 6, 7],
[4, 6, 5, 1, 0, 2],
[7, 0, 1, 2, 3, 1, 4, 3],
[7, 6, 5, 1, 3, 4, 2, 3, 0, 2],
[7, 6, 5, 4, 2, 0],
[2, 3, 1, 2, 4, 3, 1, 2, 3, 5, 4, 3, 1, 2, 6, 5, 4, 3, 1, 2, 0, 3, 4, 5, 6, 7],
[7, 5, 4, 3, 1, 2, 3, 4, 5, 6, 0, 2, 3, 4, 1, 3],
[6, 7, 5, 1, 0, 2, 3, 4],
[2],
[6, 5, 6, 4, 5, 6, 3, 4, 5, 6, 1, 3, 4, 5, 6, 2, 3, 4, 5, 6, 1, 3, 4, 5, 2, 3,
    4, 1, 3, 2, 0, 2, 3, 4, 5, 6, 1, 3, 4, 5, 2, 3, 4, 1, 3, 2, 0, 2, 3, 4, 5,
    6, 1, 3, 4, 5, 2, 3, 4, 1, 3, 2, 0],
[7, 5, 0],
[6, 4, 3, 4, 1, 3, 4, 2, 3, 4, 1, 3, 2],
[3, 0, 2],
[1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 3, 4, 5, 4, 3, 1, 2, 3, 4, 5, 6, 5, 4, 3, 1, 2,
    3, 4, 5, 6, 7],
[1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 3, 4, 5],
[7, 4, 3, 1, 2, 3, 1],
[1, 3, 1, 2, 3, 4, 3, 1, 2, 0, 3, 2, 4, 5, 4, 3, 1, 2, 3, 4, 5, 6, 5, 4, 3, 1,
    2, 0, 3, 2, 4, 5, 6],
[1, 2, 3, 4, 3, 1, 2, 3, 4, 6, 5, 7, 6, 5, 4, 3, 1, 2, 3, 4, 5],
[7, 4, 1, 3, 0],
[3, 1, 0],
[1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 3, 4, 5, 4, 3, 1, 2, 3, 4, 5, 6, 5, 4, 3, 1, 2,
    0, 3, 4, 5, 6],
[7, 4, 5, 3, 4, 5, 2, 3, 4, 5, 1, 3, 2, 0, 2, 3, 4, 5, 1, 3, 4, 2, 3, 1, 0],
[4, 1, 3, 4, 0, 2, 3, 4, 5, 1, 3, 4, 2, 3, 1, 0, 2, 3, 4, 5, 6],
[7, 4, 1, 2, 0],
[6, 4, 3, 1, 2],
[7, 5, 4, 2, 0],
[1, 2, 3, 1, 2, 0, 3, 4, 3, 1, 5, 4, 3, 1, 2, 0, 3, 2, 4, 3, 5, 4, 6],
[0, 1, 2, 0, 3, 1, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 5, 4, 3, 1, 2, 3, 4, 5, 6],
[0, 1, 2, 3, 1, 5, 4, 3, 1, 2, 3, 4, 7],
[6, 4, 1, 5, 3],
[7, 1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 0, 3, 4, 5],
[1, 2, 4, 0, 3],
[1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 3, 4, 6, 5, 7],
[1, 2, 4, 5, 4, 3, 1, 2, 3, 6, 5, 4, 3, 1, 2, 0, 3],
[3, 5, 7, 0, 2, 4, 6],
[4, 5, 1, 3, 0],
[5, 6, 4, 3, 4, 1, 3, 4, 0, 2, 3, 4, 1, 3, 2],
[1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 0, 3, 4, 7, 6],
[3, 1, 4, 3, 1, 2, 0],
[6, 5, 3, 1, 2],
[1, 2, 3, 1, 2, 3, 4, 6, 5, 4, 3, 1, 2, 3, 4, 5, 7],
[0, 3, 1, 4, 3, 1, 2, 6, 7],
[5, 6, 7, 4, 3, 1, 2],
[7, 5, 3, 0, 1, 2, 4],
[3, 1, 4, 3, 2, 0, 5, 4, 6, 5, 4],
[7, 5, 6, 1, 3, 4, 0],
[4, 2, 3, 4, 5, 6, 1, 3, 0],
[7, 4, 5, 2, 3, 1, 0, 2, 3],
[5, 6, 4, 1, 3, 2, 0],
[7, 5, 6, 3, 0, 1, 2],
[6, 7, 4, 5, 3, 1, 3, 2, 3],
[6, 7, 3, 4, 1, 2, 0],
[3, 2, 4, 3, 1, 2, 3, 4, 5, 7, 6],
[7, 5, 6, 4, 1, 2, 0],
[1, 3, 1, 2, 0, 5, 4, 3, 1, 2, 3, 4, 6]
]

# We now check that these representatives genuinely belong to the matched class.
for i, word in zip(E8classmatch, E8vgwords):
    w = W.convert(word, 'wp')
    assert w.cycletype(True) == E8labs[i]

# Now sort the conjugacy classes based on their cuspidal label.
E8coxclasses = list(clp.data.typ1E.coxeterclasses(list(range(W.rank))))
tmp = [x[0] for x in E8coxclasses]
E8sortedclasses = []

for i, j in enumerate(E8classmatch):
    # Get a permutation rep of the GAP class.
    w = W.convert(np.array(E8ctab.classreps[j]), 'mp')

    # Get a minimal length element in the class and determined J(w_min).
    # We then determine the representative of the Coxeter class
    # containing J(w_min).
    w_min = W.convert(minlengthrep(W, w), 'pw')
    J = list(set(w_min))
    a = identifycoxclass(W, tmp, J)
    E8sortedclasses.append(((i, j), (a, len(w_min), E8vgwords[i])))

    # Check that the very good representative is of minimal length.
    assert len(w_min) == len(E8vgwords[i])

E8sortedclasses.sort(key = (lambda pair: pair[1]))
E8sortedclasses = [x[0] for x in E8sortedclasses]



# Get the pairs consisting of the degree and b-value.
E8degb = [(E8ctab.irrchars[i][ind_A0], E8ctab.bvals[i])
          for i in range(len(E8reps))]

E8irrdata = [
[1,    0,   0,   '_x'  ],
[1,    120, 120, "_x'" ],
[28,   3,   8,   '_x'  ],
[28,   63,  68,  "_x'" ],
[35,   2,   2,   '_x'  ],
[35,   74,  74,  "_x'" ],
[70,   16,  32,  '_y'  ],
[50,   4,   8,   '_x'  ],
[50,   52,  56,  "_x'" ],
[84,   3,   4,   '_x'  ],
[84,   63,  64,  "_x'" ],
[168,  16,  24,  '_y'  ],
[175,  8,   12,  '_x'  ],
[175,  32,  36,  "_x'" ],
[210,  4,   4,   '_x'  ],
[210,  52,  52,  "_x'" ],
[420,  16,  20,  '_y'  ],
[300,  6,   8,   '_x'  ],
[300,  42,  44,  "_x'" ],
[350,  8,   14,  '_x'  ],
[350,  32,  38,  "_x'" ],
[525,  12,  12,  '_x'  ],
[525,  36,  36,  "_x'" ],
[567,  6,   6,   '_x'  ],
[567,  46,  46,  "_x'" ],
[1134, 16,  20,  '_y'  ],
[700,  13,  16,  '_xx' ],
[700,  25,  28,  "_xx'"],
[700,  6,   6,   '_x'  ],
[700,  42,  42,  "_x'" ],
[1400, 16,  20,  '_y'  ],
[840,  12,  14,  '_x'  ],
[840,  24,  26,  "_x'" ],
[1680, 16,  22,  '_y'  ],
[972,  10,  12,  '_x'  ],
[972,  30,  32,  "_x'" ],
[1050, 8,   10,  '_x'  ],
[1050, 32,  34,  "_x'" ],
[2100, 20,  20,  '_y'  ],
[1344, 7,   8,   '_x'  ],
[1344, 37,  38,  "_x'" ],
[2688, 16,  20,  '_y'  ],
[1400, 8,   8,   '_x'  ],
[1400, 32,  32,  "_x'" ],
[1575, 8,   10,  '_x'  ],
[1575, 32,  34,  "_x'" ],
[3150, 16,  18,  '_y'  ],
[2100, 13,  16,  '_x'  ],
[2100, 25,  28,  "_x'" ],
[4200, 16,  18,  '_y'  ],
[2240, 10,  10,  '_x'  ],
[2240, 28,  28,  "_x'" ],
[4480, 16,  16,  '_y'  ],
[2268, 10,  10,  '_x'  ],
[2268, 30,  30,  "_x'" ],
[4536, 16,  18,  '_y'  ],
[2835, 14,  14,  '_x'  ],
[2835, 22,  22,  "_x'" ],
[5670, 16,  18,  '_y'  ],
[3200, 15,  16,  '_x'  ],
[3200, 21,  22,  "_x'" ],
[4096, 11,  12,  '_x'  ],
[4096, 26,  26,  "_x'" ],
[4200, 12,  12,  '_x'  ],
[4200, 24,  24,  "_x'" ],
[6075, 14,  14,  '_x'  ],
[6075, 22,  22,  "_x'" ],
[8,    1,   1,   '_z'  ],
[8,    91,  91,  "_z'" ],
[56,   7,   19,  '_z'  ],
[56,   37,  49,  "_z'" ],
[112,  3,   3,   '_z'  ],
[112,  63,  63,  "_z'" ],
[160,  4,   7,   '_z'  ],
[160,  52,  55,  "_z'" ],
[448,  16,  25,  '_w'  ],
[400,  6,   7,   '_z'  ],
[400,  42,  43,  "_z'" ],
[448,  7,   9,   '_z'  ],
[448,  37,  39,  "_z'" ],
[560,  5,   5,   '_z'  ],
[560,  47,  47,  "_z'" ],
[1344, 16,  19,  '_w'  ],
[840,  10,  13,  '_z'  ],
[840,  28,  31,  "_z'" ],
[1008, 7,   9,   '_z'  ],
[1008, 37,  39,  "_z'" ],
[2016, 16,  19,  '_w'  ],
[1296, 10,  13,  '_z'  ],
[1296, 30,  33,  "_z'" ],
[1400, 10,  11,  '_zz' ],
[1400, 28,  29,  "_zz'"],
[1400, 7,   7,   '_z'  ],
[1400, 37,  37,  "_z'" ],
[2400, 15,  17,  '_z'  ],
[2400, 21,  23,  "_z'" ],
[2800, 13,  13,  '_z'  ],
[2800, 25,  25,  "_z'" ],
[5600, 16,  19,  '_w'  ],
[3240, 9,   9,   '_z'  ],
[3240, 31,  31,  "_z'" ],
[3360, 12,  13,  '_z'  ],
[3360, 24,  25,  "_z'" ],
[7168, 16,  17,  '_w'  ],
[4096, 11,  11,  '_z'  ],
[4096, 26,  27,  "_z'" ],
[4200, 15,  15,  '_z'  ],
[4200, 21,  21,  "_z'" ],
[4536, 13,  13,  '_z'  ],
[4536, 23,  23,  "_z'" ],
[5600, 15,  15,  '_z'  ],
[5600, 21,  21,  "_z'" ]
]

ind_A0 = E8classmatch[[i for i, x in enumerate(E8carter) if x[0] == 'A0'][0]]

E8irrmatch = [E8degb.index((x[0], x[2])) for x in E8irrdata]

assert len(set(E8irrmatch)) == len(E8irrmatch)

# Sort the irreducible characters by degree and b-value.
E8sortedirrs = sorted(enumerate(E8irrmatch),
                      key=(lambda x: [E8irrdata[x[0]][y] for y in [0, 2]]))


###########################################################################
##                                                                       ##
##                          Write Data to File                           ##
##                                                                       ##
###########################################################################
head_str = "{0:#<74}\n##{0:^70}##\n##{1:^70}##\n##{0:^70}##\n{0:#<74}\n"

# Write the conjugacy class data to file.
def writeclassdata(typ, classes):
    """
    Writes a function to the corresponding data file returning the Coxeter
    classes.

    """
    # Construct the TextWrapper class.
    wrapper = textwrap.TextWrapper(initial_indent=" "*4,
                                   subsequent_indent=" "*8,
                                   width=78)

    with open("../rawdata.py", 'a') as datafile:
        datafile.write(
            "{}conjclasses = (\n".format(typ)
        )

        # Write out the list.
        for nam, cent, rep in classes[:-1]:
            s = '("{}", {}, {}),'.format(nam, cent, rep)
            datafile.write("\n".join(wrapper.wrap(s) + [""]))

        # Treat the last item differently.
        nam, cent, rep = classes[-1]
        s = '("{}", {}, {})\n'.format(nam, cent, rep)
        datafile.write("\n".join(wrapper.wrap(s) + [""]))

        datafile.write(")\n\n")


G2clsdata = [(G2carter[i][0], G2ctab.cents[j], G2vgwords[i])
             for i, j in G2sortedclasses]
F4clsdata = [(F4carter[i][0], F4ctab.cents[j], F4vgwords[i])
             for i, j in F4sortedclasses]
E6clsdata = [(E6carter[i][0], E6ctab.cents[j], E6vgwords[i])
             for i, j in E6sortedclasses]
E7clsdata = [(E7carter[i][0], E7ctab.cents[j], E7vgwords[i])
             for i, j in E7sortedclasses]
E8clsdata = [(E8carter[i][0], E8ctab.cents[j], E8vgwords[i])
             for i, j in E8sortedclasses]

# Write a header to the section.
with open("../rawdata.py", 'a') as datafile:
    datafile.write(head_str.format("", "Conjugacy Classes"))

writeclassdata("G2", G2clsdata)
writeclassdata("F4", F4clsdata)
writeclassdata("E6", E6clsdata)
writeclassdata("E7", E7clsdata)
writeclassdata("E8", E8clsdata)

# Write the irreducible character data to file.
def writeirrdata(typ, irrs):
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
            "{}irrchardata = (\n".format(typ)
        )
        
        datafile.write(" "*4)
        linelength = 4

        # Write out the list.
        for tup in irrs[:-1]:
            s = str(tup) + ', '
            if linelength + len(s) <= 78:
                datafile.write(s)
                linelength += len(s)
            else:
                datafile.write("\n    " + s)
                linelength = 4 + len(s)

        # Treat the last item differently.
        tup = irrs[-1]
        s = str(irrs[-1]) + '\n'
        if linelength + len(s) <= 78:
            datafile.write(s)
        else:
            datafile.write("\n    " + s)

        datafile.write(")\n\n")


G2irrdatafinal = [tuple(G2irrdata[i]) for i, j in G2sortedirrs]
F4irrdatafinal = [tuple(F4irrdata[i]) for i, j in F4sortedirrs]
E6irrdatafinal = [tuple(E6irrdata[i]) for i, j in E6sortedirrs]
E7irrdatafinal = [tuple(E7irrdata[i]) for i, j in E7sortedirrs]
E8irrdatafinal = [tuple(E8irrdata[i]) for i, j in E8sortedirrs]

# Write a header to the section.
with open("../rawdata.py", 'a') as datafile:
    datafile.write(head_str.format("", "Irreducible Characters"))

writeirrdata("G2", G2irrdatafinal)
writeirrdata("F4", F4irrdatafinal)
writeirrdata("E6", E6irrdatafinal)
writeirrdata("E7", E7irrdatafinal)
writeirrdata("E8", E8irrdatafinal)


# Write the character table to file.
def writechartable(typ, table):
    """
    Writes the character table computed by GAP to the rawdata file.

    """
    # Construct the TextWrapper class.
    wrapper = textwrap.TextWrapper(initial_indent=" "*4,
                                   subsequent_indent=" "*8,
                                   width=78)

    with open("../rawdata.py", 'a') as datafile:
        datafile.write(
            "{}chartable = [\n".format(typ)
        )

        # Write out the list.
        for row in table[:-1]:
            datafile.write("\n".join(wrapper.wrap(str(row))))
            datafile.write(",\n")

        # Treat the last item differently.
        datafile.write("\n".join(wrapper.wrap(str(table[-1]))))
        datafile.write("\n]\n\n")


G2chartab = [[G2ctab.irrchars[row][col] for i, col in G2sortedclasses]
             for k, row in G2sortedirrs]
F4chartab = [[F4ctab.irrchars[row][col] for i, col in F4sortedclasses]
             for k, row in F4sortedirrs]
E6chartab = [[E6ctab.irrchars[row][col] for i, col in E6sortedclasses]
             for k, row in E6sortedirrs]
E7chartab = [[E7ctab.irrchars[row][col] for i, col in E7sortedclasses]
             for k, row in E7sortedirrs]
E8chartab = [[E8ctab.irrchars[row][col] for i, col in E8sortedclasses]
             for k, row in E8sortedirrs]

# Write a header to the section.
with open("../rawdata.py", 'a') as datafile:
    datafile.write(head_str.format("", "Character Tables"))

writechartable("G2", G2chartab)
writechartable("F4", F4chartab)
writechartable("E6", E6chartab)
writechartable("E7", E7chartab)
writechartable("E8", E8chartab)


