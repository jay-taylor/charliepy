# Here we determine the phi-conjugacy classes for 2E6 and 3D4.

import charliepy as clp

import functools
import itertools
import operator
import textwrap


###############################################################
##                                                           ##
##                         Type 2E6                          ##
##                                                           ##
###############################################################

# The non-trivial automorphism phi : W -> W can be identified with
# conjugation by the longest element w0. Hence we obtain a W-equivariant
# bijection W.phi -> W defined by w.phi -> w.w0. This gives an
# identification of phi-conjugacy classes of W with the conjugacy
# classes of W. The phi-conjugacy class containing w is mapped to the
# conjugacy class containing w.w0.
W = clp.CoxeterGroup("E", 6)
E6conjclasses = clp.conjugacyclasses(W)

# Construct the coset.
phi = clp.cyclestoperm((0, 5), (2, 4))
Wphi = clp.CoxeterCoset(W, phi)

# Get phi as a full permutation of the roots.
phi = Wphi.phi

# Longest element as a permutation.
w0 = W.longestelement(lab='p')

E6cycstructs = [W.convert(c.rep, 'p').cyclestructure for c in E6conjclasses]

# Check the cycle structure determines the class uniquely.
assert len(set(E6cycstructs)) == len(E6cycstructs)

# These are very good class representatives taken from CHEVIE.
E6greps = (
[0, 1, 2, 0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 2, 4, 3, 1, 5, 4, 3, 1, 2, 0,
    3, 2, 4, 3, 1, 5, 4, 3, 2, 0],
[],
[2, 3, 2, 4, 3, 2],
[0, 1, 3, 2, 0, 4, 3, 2, 5, 4, 3, 2],
[0, 1, 2, 0, 3, 2, 0, 4, 3, 2, 0, 5, 4, 3, 2, 0],
[1, 2, 3, 1, 2, 3, 5, 4, 3, 1, 2, 3, 4, 5],
[0, 3, 1, 2, 0, 3, 2, 4, 3, 1, 2, 0, 3, 5, 4, 3, 2, 0],
[0, 1],
[3, 4, 3, 1, 2, 0, 3, 4],
[3, 1, 4, 3, 1, 2, 3, 4, 5, 4, 3, 1, 2, 3, 4, 5],
[1, 3],
[0, 4],
[4, 3],
[0, 1, 4, 3],
[0, 1, 2, 0, 3, 2],
[0, 2, 0, 3, 2, 0, 4, 3, 2, 0, 5, 4, 3, 2, 0],
[1],
[0],
[1, 2, 3, 2, 4, 3, 2],
[0, 2, 3, 2, 4, 3, 2],
[0, 2, 0, 3, 2],
[0, 1, 4],
[1, 4, 3],
[0, 4, 3],
[0, 1, 3]
)

E6phiclasses = [None]*len(E6conjclasses)

for rep in E6greps:
    w = W.convert(rep, 'p')
    wphi = w*w0

    # Get the class containing w*w0.
    i = E6cycstructs.index(wphi.cyclestructure)
    c = E6conjclasses[i]

    # Check w is a good representative.
    assert clp.braid.isverygood(Wphi, w)

    # Check w is a minimal length representative.
    wmin = clp.minlengthelms(Wphi, w)[0]
    assert len(rep) == W.length(wmin)

    E6phiclasses[i] = (i, (c.name, c.centsize, rep))

# Sort so they're in the same order as the conjugacy classes of E6.
E6phiclasses.sort(key=operator.itemgetter(0))
E6phiclasses = [tup for _, tup in E6phiclasses]



###############################################################
##                                                           ##
##                         Type 3D4                          ##
##                                                           ##
###############################################################

# The Dynkin diagram is labelled as follows:
#
#                               2
#                             /
#                       0 - 1
#                             \
#                               3
#
# We assume that phi is the automorphism of order 3 which is
#
#           phi:  0 -> 2 -> 3 -> 0
#

W = clp.CoxeterGroup("D", 4)
Wphi = clp.CoxeterCoset(W, clp.cyclestoperm((0, 2, 3)))
phi = Wphi.phi

# Get all the elements of W as permutations on the roots.
Welms = {x for x in W}

# We now use a basic orbit algorithm to work out the phi-conjugacy classes.
genpairs = [(W.permgens[i], W.permgens[i^phi]) for i in range(4)]
D4orbs = []

while Welms:
    orb = [Welms.pop()]

    for w in orb:
        for (s, t) in genpairs:
            z = s*w*t
            if z not in orb:
                orb.append(z)
                Welms.remove(z)

    D4orbs.append(orb)

# Get the minimal length of the elements in an orbit.
D4minlens = [min(map(W.length, orb)) for orb in D4orbs]

# Get all minimal length elements in each orbit.
Cmin = [[W.convert(x, 'w') for x in orb if W.length(x) == mlen]
        for mlen, orb in zip(D4minlens, D4orbs)]

# We define an embedding of W into F4.
H = clp.CoxeterGroup("F", 4)
#embed = (2, 3, 5, 7) # short roots
embed = (1, 0, 8, 15) # long roots
W.embeddings[H] = embed

s0, s1, s2, s3 = embed

# There is a unique element of H which induces the triality automorphism
# on the underlying four roots (not just on the reflections). This gives us an
# embedding of the coset W.phi into H.
autos = []
for x in H:
    if (s0^x, s1^x, s2^x, s3^x) == (s2, s1, s3, s0):
        autos.append(x)

assert len(autos) == 1

Hphi = autos[0]

# Just for paranoias sake, let's also check that conjugation with phi
# permutes the generators as expected.
WgensH = [H.reflection(i) for i in embed]
s0, s1, s2, s3 = WgensH
assert (s0^Hphi, s1^Hphi, s2^Hphi, s3^Hphi) == (s2, s1, s3, s0)

# We now consider how to label the phi-conjugacy classes in 3D4. Below
# we check that if O is a phi-conjugacy class of D4 then it is uniquely
# determined by the conjugacy classes of F4 which meet the set Omin.phi
# of minimal length elements in the coset. This essentially allows us to
# distinguish the classes.
F4conjclasses = clp.conjugacyclasses(H)

# Get the orbits of short roots and long roots.
o1 = [i for i, x in enumerate(H.orbitrepresentatives) if x == 0]
o2 = [i for i, x in enumerate(H.orbitrepresentatives) if x == 2]

# Determine the cyclestructures of the permutations on each orbit.
F4cycstructs = [None]*len(F4conjclasses)
for i, c in enumerate(F4conjclasses):
    w = H.convert(c.rep, 'p')
    param = (clp.restrictedperm(w, o1).cyclestructure,
             clp.restrictedperm(w, o2).cyclestructure)
    F4cycstructs[i] = param

# Check the cyclestructures completely distinguish the classes.
assert len(set(F4cycstructs)) == len(F4cycstructs)

# Get the parameters of the conjugacy classes of F4 meeting Omin.
D4labels = []
allparams = set()

identity = clp.Perm(range(2*H.N))

for orb in Cmin:
    S = set()
    for g in orb:
        w = functools.reduce(operator.mul, (WgensH[i] for i in g), identity)
        w = w*Hphi
        param = (clp.restrictedperm(w, o1).cyclestructure,
                 clp.restrictedperm(w, o2).cyclestructure)
        S.add(param)
    D4labels.append(S)
    allparams.update(S)

# Remove duplicates.
for i, lab in enumerate(D4labels):
    if len(lab) > 1:
        others = set.union(*(D4labels[:i] + D4labels[i+1:]))
        lab -= others
        assert len(lab) == 1

# Check that each class is uniquely labelled this way.
assert all(D4labels.count(x) == 1 for x in D4labels)


D4greps = (
[0],
[],
[0, 2, 1, 0, 2, 1],
[1],
[0, 1],
[0, 2, 1, 0, 2, 3, 1, 2],
[0, 2, 1, 2],
)


D4phiclasses = []

for rep in D4greps:
    # Get the element in D4 and F4.
    x = W.convert(rep, 'p')
    w = functools.reduce(operator.mul, (WgensH[i] for i in rep), identity)

    # Determine the class parameter of w.phi in F4.
    wphi = w*Hphi
    param = (clp.restrictedperm(wphi, o1).cyclestructure,
             clp.restrictedperm(wphi, o2).cyclestructure)
    
    # Check it's a good representative.
    assert clp.braid.isverygood(Wphi, x)

    # Get the orbit it's contained in.
    i = D4labels.index({param})
    cent = W.size//len(D4orbs[i])

    # Check it's minimal length.
    assert len(rep) == D4minlens[i]

    # Get its label as a class of F4.
    i = F4cycstructs.index(param)
    D4phiclasses.append((i, (F4conjclasses[i].name, cent, rep)))

# Sort the classes by the order of the classes in F4.
#D4phiclasses.sort(key=operator.itemgetter(0))
#D4phiclasses = [tup for _, tup in D4phiclasses]

# Now the irreducible characters of the coset. According to Lusztig all
# of these can be obtained by restricting representations of F4. It
# suffices to show that for each irreducible character chi of D4
# invariant under phi there exists an irreducible character psi of F4
# such that Res(psi) = chi.
D4ct = clp.CharacterTable(W)

# Now determine the permutation of the conjugacy classes induced by phi.
D4classes = clp.conjugacyclasses(W)
reps = [W.convert(C.rep, 'p') for C in D4classes]
repsperm = [clp.conjclass.identifyclass(W, x^phi) for x in reps]
p = clp.Perm([D4ct.classnames.index(nam) for nam in repsperm])

# Check which characters are fixed by phi. Note it doesn't matter whether we
# consider phi or phi' = phi^-1 here.
D4irrs = D4ct.irreducibles
irrsperm = [[chi[i^p] for i in range(len(D4ct.classes))] for chi in D4irrs]
irrsfix = [i for i in range(len(D4irrs)) if D4irrs[i] == irrsperm[i]]

# Now get the characters of the coset.
F4ct = clp.CharacterTable(H)
F4irrs = F4ct.irreducibles
indtb = clp.chartab.inductiontable(D4ct, F4ct)
restb = [list(x) for x in zip(*indtb)]

# Find a character restricting irreducibly.
charpairs = [None]*len(irrsfix)

for i, j in enumerate(irrsfix):
    chi = [0]*len(D4irrs)
    chi[j] = 1
    psi = restb.index(chi)
    charpairs[i] = (j, psi)

## Now restrict these characters to the coset.
D4_3ct = [[F4irrs[j][cls[0]] for cls in D4phiclasses] for (i, j) in charpairs]

D4_3irrdata = [(D4ct.charnames[i], None, None) for (i, j) in charpairs]


###########################################################################
##                                                                       ##
##                          Write Data to File                           ##
##                                                                       ##
###########################################################################
head_str = "{0:#<74}\n##{0:^70}##\n##{1:^70}##\n##{0:^70}##\n{0:#<74}\n"

# Write the conjugacy class data to file.
def writeclassdata(typ, twist, classes):
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
            "{}_{}conjclasses = (\n".format(typ, twist)
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


# Write a header to the section.
with open("../rawdata.py", 'a') as datafile:
    datafile.write(head_str.format("", "Phi-Conjugacy Classes"))

writeclassdata("E6", "2", E6phiclasses)
writeclassdata("D4", "3", [cls[1] for cls in D4phiclasses])

# Write the character table to file.
def writechartable(typ, twist, table):
    """
    Writes the character table computed by GAP to the rawdata file.

    """
    # Construct the TextWrapper class.
    wrapper = textwrap.TextWrapper(initial_indent=" "*4,
                                   subsequent_indent=" "*8,
                                   width=78)

    with open("../rawdata.py", 'a') as datafile:
        datafile.write(
            "{}_{}chartable = [\n".format(typ, twist)
        )

        # Write out the list.
        for row in table[:-1]:
            datafile.write("\n".join(wrapper.wrap(str(row))))
            datafile.write(",\n")

        # Treat the last item differently.
        datafile.write("\n".join(wrapper.wrap(str(table[-1]))))
        datafile.write("\n]\n\n")

# Write the irreducible character data to file.
def writeirrdata(typ, twist, irrs):
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
            "{}_{}irrchardata = (\n".format(typ, twist)
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

# Write a header to the section.
with open("../rawdata.py", 'a') as datafile:
    datafile.write(head_str.format("", "Phi-Irreducible Characters"))

writeirrdata("D4", 3, D4_3irrdata)

# Write a header to the section.
with open("../rawdata.py", 'a') as datafile:
    datafile.write(head_str.format("", "Phi-Character Tables"))

writechartable("D4", "3", D4_3ct)

