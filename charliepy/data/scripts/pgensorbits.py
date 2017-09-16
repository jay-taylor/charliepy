import charliepy as clp
import textwrap

# Header string to signify start of a new section.
head_str = "{0:#<74}\n##{0:^70}##\n##{1:^70}##\n##{0:^70}##\n{0:#<74}\n"

# Construct the TextWrapper class.
wrapper = textwrap.TextWrapper(initial_indent=" "*4,
                               subsequent_indent=" "*8,
                               width=78)

# We first compute the standard permutation generators.
with open("../rawdata.py", 'a') as datafile:
    # Write the header file for the section.
    datafile.write(head_str.format("", "Permutation Generators"))

    for s, i in [["G", 2], ["F", 4], ["E", 6], ["E", 7], ["E", 8]]:
        permgens = clp.CoxeterGroup(s, i).permgens
        datafile.write("{}{}permgens = [\n".format(s, i))
        for gen in permgens[:-1]:
            datafile.write("\n".join(wrapper.wrap(str(gen.cycledecomp()))))
            datafile.write(",\n")

        # Treat the last item differently.
        datafile.write("\n".join(wrapper.wrap(str(permgens[-1].cycledecomp()))))
        datafile.write("\n]\n\n")

# Now the orbits of the permutations on the roots.
with open("../rawdata.py", 'a') as datafile:
    # Write the header file for the section.
    datafile.write(head_str.format("", "Root Orbits"))

    for s, i in [["G", 2], ["F", 4]]:
        W = clp.CoxeterGroup(s, i)

        if s == "G":
            a, b = 0, 1
        elif s == "F":
            a, b = 0, 2

        o1 = [i for i, x in enumerate(W.orbitrepresentatives) if x == a]
        o2 = [i for i, x in enumerate(W.orbitrepresentatives) if x == b]
        datafile.write("{}{}rootorbits = [\n".format(s, i))
        datafile.write("\n".join(wrapper.wrap(str(o1))))
        datafile.write(",\n")

        # Treat the last item differently.
        datafile.write("\n".join(wrapper.wrap(str(o2))))
        datafile.write("\n]\n\n")
