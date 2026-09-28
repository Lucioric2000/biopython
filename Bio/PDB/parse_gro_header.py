# Copyright 2004 Kristian Rother.
# Revisions copyright 2004 Thomas Hamelryck.
# Revisions copyright 2009 Lucio Montero.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Parse the header of a GRO file."""

import sys


def _nice_case(line):
    """Make a lowercase string with capitals."""
    line_lower = line.lower()
    s = ""
    i = 0
    next_cap = 1
    while i < len(line_lower):
        c = line_lower[i]
        if "a" <= c <= "z" and next_cap:
            c = c.upper()
            next_cap = 0
        elif c in (" ", ".", ",", ";", ":", "\t", "-", "_"):
            next_cap = 1
        s += c
        i += 1
    return s


def parse_gro_header(file):
    """Return the header lines of a gro file as a dictionary.

    Dictionary keys are: name, atom_count, box.
    """
    header = []
    if isinstance(file, str):
        f = open(file)
        close_file = True
    else:
        f = file
        close_file = False
    header.append(f.readline())
    line_count = 0
    atoms_count = (
        -10
    )  # Negative value to initialize without participating in a wrong true evaluation
    for line in f:
        line_count += 1
        if line_count == 1:
            atoms_count = int(line)
            header.append(atoms_count)
        if line_count == (atoms_count + 2):
            header.append(line)
    if close_file:
        f.close()
    return _parse_gro_header_list(header)


def _parse_gro_header_list(header):
    return {"name": header[0], "atom_count": header[1], "box": header[2]}


if __name__ == "__main__":
    filename = sys.argv[1]
    with open(filename) as f:
        dct = parse_gro_header(f)

    # print the dictionary
    for k, v in dct.items():
        print("-" * 40)
        print(k)
        print(v)
