# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parser for GRO files."""

import re
import warnings

import numpy

from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.parse_gro_header import _parse_gro_header_list

water_names = ["HOH", "WAT", "SOL", "T3P", "T4P"]


def ElementGuesser(atom_name, residue_name):
    """Guess the element from atom name and residue name."""
    el_start = 0
    if atom_name[0] == "A":  # Handle leading 'A' seen in ATP
        el_start = 1
    while atom_name[el_start].isdigit() or atom_name[el_start] in "'*":
        el_start += 1
    el_end = el_start + 1
    while (el_end < len(atom_name)) and (not atom_name[el_end].isdigit()):
        el_end += 1
    resname_upper = residue_name.upper()
    atom_root_name_upper = atom_name[el_start:el_end].upper()
    atom_root_name_upper_0to2 = atom_root_name_upper[:2]
    if resname_upper in ("NA", "NA+") and atom_root_name_upper in ("NA", "NA+"):
        element = "NA"
    elif resname_upper in ("CA", "CA2+", "CAL", "CA++") and atom_root_name_upper in (
        "CA",
        "CA2+",
        "CAL",
        "CA++",
    ):
        element = "CA"
    elif atom_root_name_upper_0to2 in ("SE", "SR", "MG", "CU", "CL", "BR"):
        element = atom_root_name_upper_0to2
    elif atom_root_name_upper[0] in "CHONPSK":
        element = atom_root_name_upper[0]
    else:
        raise PDBConstructionException(
            f"Could not assign element to atom name {atom_name}"
        )
    return element


class GROParser:
    """Parse a GRO file and return a Structure object."""

    def __init__(self, PERMISSIVE=1, structure_builder=None):
        """Create a GROParser object.

        The GRO parser calls a number of standard methods in an aggregated
        StructureBuilder object. Normally this object is instantiated by the
        GROParser object itself, but if the user provides his own StructureBuilder
        object, the latter is used instead.

        Arguments:
         - PERMISSIVE - int, if this is 0 exceptions in constructing the
           SMCRA data structure are fatal. If 1 (DEFAULT), the exceptions are
           caught, but some residues or atoms will be missing. THESE EXCEPTIONS
           ARE DUE TO PROBLEMS IN THE GRO FILE!.
         - structure_builder - an optional user implemented StructureBuilder class.

        """
        if structure_builder is not None:
            self.structure_builder = structure_builder
        else:
            self.structure_builder = StructureBuilder()
        self.header = None
        self.trailer = None
        self.line_counter = 0
        self.PERMISSIVE = PERMISSIVE

    # Public methods

    def get_structure(self, id, file):
        """Return the structure.

        Arguments:
         - id - string, the id that will be used for the structure
         - file - name of the GRO file OR an open filehandle

        """
        self.header = None
        self.trailer = None
        self.filename_in_read = file
        # Make a StructureBuilder instance (pass id of structure as parameter)
        self.structure_builder.init_structure(id)
        if isinstance(file, str):
            with open(file) as f:
                lines = f.readlines()
        else:
            lines = file.readlines()
        if not lines:
            raise ValueError("Empty file.")
        self._parse(lines)
        self.structure_builder.set_header(self.header)
        # Return the Structure instance
        getstruct = self.structure_builder.get_structure()
        self.filename_in_read = None
        return getstruct

    def get_header(self):
        """Return the header."""
        return self.header

    def get_trailer(self):
        """Return the trailer."""
        return self.trailer

    # Private methods

    def _parse(self, header_coords_trailer):
        """Parse the GRO file (PRIVATE)."""
        # In the case of GRO files we don't extract the header; because it is a specific header/trailer line set for each
        # time frame
        # self.header, coords_trailer=self._get_header(header_coords_trailer)
        # Parse the atomic data; return the PDB file trailer
        self.trailer = self._parse_coordinates(header_coords_trailer)

    def _get_header(self, header_coords_trailer):
        """Get the header of the GRO file, return the rest (PRIVATE)."""
        structure_builder = self.structure_builder
        atoms_count = (
            -10
        )  # Negative value to initialize the variable without participating in a wrong true evaluation
        header = []
        for i in range(0, len(header_coords_trailer)):
            structure_builder.set_line_counter(i + 1)
            if (i + 1) == 1:
                header.append(header_coords_trailer[0])
            if (i + 1) == 2:
                atoms_count = int(header_coords_trailer[1])
                header.append(atoms_count)
            if (i + 1) == (atoms_count + 2):
                header.append(header_coords_trailer[i])
        # Return the rest of the coords+trailer for further processing
        self.line_counter = 2
        coords_trailer = header_coords_trailer[2 : atoms_count + 2]
        header_dict = _parse_gro_header_list(header)
        return header_dict, coords_trailer

    def _parse_coordinates(self, coords_trailer):
        """Parse the atomic data in the GRO file (PRIVATE)."""
        local_line_counter = 0
        structure_builder = self.structure_builder
        # Flag we have an open model
        model_open = 0  # 0: closed, 1: description line, 2: atom count, 3: coords, 4: box, 5: after box
        current_chain_id = None
        current_segid = None
        current_residue_id = None
        current_resname = None
        thismodel_atoms_line_counter = 0
        atoms_count_for_this_model = 0
        for i in range(0, len(coords_trailer)):
            line = coords_trailer[i].rstrip()
            global_line_counter = self.line_counter + local_line_counter + 1
            structure_builder.set_line_counter(global_line_counter)
            if model_open in (0, 5):
                # Initialize the Model - there was no explicit MODEL record
                current_model_id = line
                structure_builder.init_model(current_model_id)
                model_open = 1
                thismodel_atoms_line_counter = 0
            elif model_open == 1:
                atoms_count_for_this_model = int(line)
                model_open = 2
            elif thismodel_atoms_line_counter == atoms_count_for_this_model:
                model_open = 4
                splut = line.split()
                try:
                    box = [float(val) for val in splut]
                except ValueError:
                    raise PDBConstructionException(
                        "Invalid box line: " + str(line)
                    ) from None
                structure_builder.model.box = box
                model_open = 5
            elif model_open in (2, 3):
                if model_open == 2:
                    model_open = 3
                resseq = int(line[0:5].split()[0])  # sequence identifier
                fullresname = line[5:10]
                # get rid of whitespace in residue names
                split_list_resname = fullresname.split()
                if len(split_list_resname) != 1:
                    # residue name has internal spaces, e.g. " N B ", so
                    # we do not strip spaces
                    resname = fullresname
                else:
                    # atom name is like " CA ", so we can strip spaces
                    resname = split_list_resname[0]
                fullname = line[10:15]
                # get rid of whitespace in atom names
                split_list = fullname.split()
                if len(split_list) != 1:
                    # atom name has internal spaces, e.g. " N B ", so
                    # we do not strip spaces
                    name = fullname
                else:
                    # atom name is like " CA ", so we can strip spaces
                    name = split_list[0]
                try:
                    serial_number = int(line[15:20])
                except ValueError:
                    serial_number = 0
                if resname in water_names:
                    hetero_flag = "W"
                else:
                    hetero_flag = " "

                altloc = " "
                chainid = " "
                icode = " "
                residue_id = (hetero_flag, resseq, icode)
                # atomic coordinates
                try:
                    x = float(line[20:28]) * 10
                    y = float(line[28:36]) * 10
                    z = float(line[36:44]) * 10
                except ValueError:
                    raise PDBConstructionException(
                        "Invalid or missing coordinate(s) at line %i."
                        % global_line_counter
                    ) from None
                coord = numpy.array((x, y, z), "f")
                # occupancy & B factor (absent in GRO files!)
                occupancy = 1.0
                bfactor = 0.0

                element = ElementGuesser(name, resname)
                segid = ""

                if current_segid != segid:
                    current_segid = segid
                    structure_builder.init_seg(current_segid)
                if current_chain_id != chainid:
                    current_chain_id = chainid
                    structure_builder.init_chain(current_chain_id)
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(
                            resname, hetero_flag, resseq, icode
                        )
                    except PDBConstructionException as message:
                        self._handle_GRO_exception(message, global_line_counter)
                elif current_residue_id != residue_id or current_resname != resname:
                    current_residue_id = residue_id
                    current_resname = resname
                    try:
                        structure_builder.init_residue(
                            resname, hetero_flag, resseq, icode
                        )
                    except PDBConstructionException as message:
                        self._handle_GRO_exception(message, global_line_counter)
                # init atom
                try:
                    structure_builder.init_atom(
                        name,
                        coord,
                        bfactor,
                        occupancy,
                        altloc,
                        fullname,
                        serial_number,
                        element,
                    )
                except PDBConstructionException as message:
                    self._handle_GRO_exception(message, global_line_counter)
                local_line_counter += 1
                thismodel_atoms_line_counter += 1
        # EOF (does not end in END or CONECT)
        self.line_counter = self.line_counter + local_line_counter
        return []

    def _handle_GRO_exception(self, message, line_counter):
        """Handle exception in StructureBuilder (PRIVATE).

        This method catches an exception that occurs in the StructureBuilder
        object (if PERMISSIVE==1), or raises it again, this time adding the
        PDB line number to the error message.
        """
        message = "%s at line %i of file %s." % (
            message,
            line_counter,
            self.filename_in_read,
        )
        if self.PERMISSIVE:
            # just print a warning - some residues/atoms may be missing
            warnings.warn(
                "PDBConstructionException: %s\n"
                "Exception ignored.\n"
                "Some atoms or residues may be missing in the data structure."
                % message,
                PDBConstructionWarning,
                stacklevel=2,
            )
        else:
            # exceptions are fatal - raise again with new message (including line nr)
            raise PDBConstructionException(message)


if __name__ == "__main__":
    import sys

    p = GROParser(PERMISSIVE=1)

    filename = sys.argv[1]
    s = p.get_structure("scr", filename)

    for m in s:
        p_obj = m.get_parent()
        assert p_obj is s
        for c in m:
            p_obj = c.get_parent()
            assert p_obj is m
            for r in c:
                print(r)
                p_obj = r.get_parent()
                assert p_obj is c
                for a in r:
                    p_obj = a.get_parent()
                    if p_obj is not r:
                        print(p_obj, r)
