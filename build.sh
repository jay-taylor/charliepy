#!/usr/bin/env bash

# Use these to indicate the command for Python in the chosen shell.
PYTHON=python3

# Update the timestamps on extension files.
touch ./charliepy/src/permutat/permutat.c
touch ./charliepy/data/src/_cdata.c

# Compile the extensions.
$PYTHON setup.py build_ext --inplace
