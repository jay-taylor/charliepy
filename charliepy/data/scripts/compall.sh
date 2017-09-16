#!/usr/bin/env bash

# Use these to indicate the commands for GAP and Python in the chosen
# shell.
GAP="/usr/local/bin/gap4.sh"
PYTHON=python3

# This will delete the contents of the rawdata file!
> ../rawdata.py

# First get the matrix generators of the exceptional Weyl groups.
$PYTHON matgens.py

# Use GAP4 to compute the character tables and conjugacy classes.
cat matgens.gap ctab.gap | $GAP -q

# Compute the Coxeter classes of the generators of parabolic subgroups.
$PYTHON coxclasses.py

# Compute the information concerning characters and classes.
$PYTHON chartab.py

# Write out the standard generators as permutations and the orbits of
# the roots.
$PYTHON pgensorbits.py

# Clean up the generated files.
rm matgens.gap
rm *ctab.py
