# First we give matrix generators for all the groups. These were
# produced from CharLiePy by first constructing the group and then
# looking up the attribute .matgens.

loop := [
    [G2gens, "./G2ctab.py"],
    [F4gens, "./F4ctab.py"],
    [E6gens, "./E6ctab.py"],
    [E7gens, "./E7ctab.py"],
    [E8gens, "./E8ctab.py"]
];;

for pair in loop do
    gens := pair[1];;
    file := pair[2];;

    # Construct the group given by the generators and construct the
    # character table.
    W := Group(gens);;
    T := CharacterTable(W);;

    # Get a representative of each conjugacy class.
    classreps := List(ConjugacyClasses(T), x -> Representative(x));;

    # Compute the characteristic polynomial of each such representative.
    classcharpols := List(classreps,
         x -> CoefficientsOfUnivariatePolynomial(CharacteristicPolynomial(x)));;

    # Get the values of the character table.
    irrchars := List(Irr(T), x -> ValuesOfClassFunction(x));;

    # We first compute the character of the natural module for W. Then for
    # each irreducible character we compute the Molien series with the
    # character of the natural module. From this we can obtain the b-value
    # of the irreducible character. It is the degree of the least term
    # occurring in the Molien series.
    natchar := Character(T,
                   List(ConjugacyClasses(T), x -> Trace(Representative(x))));;
    bvals := [];;
    for i in [1..Length(Irr(T))] do
        mol := MolienSeries(natchar, Irr(T)[i]);;
        b := 0;;
        while true do
            if ValueMolienSeries(mol, b) <> 0 then 
                break;
            fi;
            b := b + 1;;
        od;
        bvals[i] := b;;
    od;

    # We now print the computed data to a file. First the class reps.
    PrintTo(file, "classreps = [\n");

    for i in [1..Length(classreps)-1] do
      AppendTo(file, classreps[i], ",\n");
    od;

    AppendTo(file, classreps[Length(classreps)], "\n");
    AppendTo(file, "]\n\n");

    # Next the characteristic polynomials of the class reps.
    AppendTo(file, "classcharpols = [\n");

    for i in [1..Length(classreps)-1] do
      AppendTo(file, classcharpols[i], ",\n");
    od;

    AppendTo(file, classcharpols[Length(classreps)], "\n");
    AppendTo(file, "]\n\n");

    # Next the centraliser orders.
    AppendTo(file, "cents = ", SizesCentralizers(T), "\n\n");

    # Next the b-values.
    AppendTo(file, "bvals = ", bvals, "\n\n");

    # Finally the values of the irreducible characters.
    AppendTo(file, "irrchars = [\n");

    for i in [1..Length(irrchars)-1] do
      AppendTo(file, irrchars[i], ",\n");
    od;

    AppendTo(file, irrchars[Length(irrchars)], "\n");
    AppendTo(file, "]");
od;

quit;
