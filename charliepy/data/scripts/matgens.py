import charliepy as clp

with open("matgens.gap", 'a') as datafile:
    for s, i in [["G", 2], ["F", 4], ["E", 6], ["E", 7], ["E", 8]]:
        matgens = [gen.tolist() for gen in clp.CoxeterGroup(s, i).matgens]
        datafile.write("{}{}gens := [\n".format(s, i))
        for gen in matgens:
            datafile.write(" [{},\n".format(gen[0]))
            for row in gen[1:-1]:
                datafile.write("  {},\n".format(row))
            datafile.write("  {}],\n".format(gen[-1]))
        datafile.write("];;\n\n")




