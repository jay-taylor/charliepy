#!/usr/bin/env bash

# Use these to indicate the commands for GAP and Python in the chosen shell.
GAP="/usr/local/bin/gap4.sh"
PYTHON=python3

# First get the matrix generators of the exceptional Weyl groups.
$PYTHON matgens.py

# Use GAP4 to compute the character tables.
cat matgens.gap ctab.gap | $GAP -q

# Compute the Coxeter classes of the generators of parabolic subgroups.
$PYTHON coxclasses.py

# Compute the information concerning characters and classes.
$PYTHON chartab.py

# Clean up the generated files.
rm matgens.gap
rm *ctab.py
