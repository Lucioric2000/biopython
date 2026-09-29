# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Output of GRO files."""

from .PDBIO import Select

_ATOM_FORMAT_STRING = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"

_select = Select()


class GROIO:
    """Write a Structure object (or a subset of a Structure object) as a GRO file.

    Examples
    --------
        >>> from Bio.PDB.GROParser import GROParser
        >>> from Bio.PDB.GROIO import GROIO
        >>> p = GROParser()
        >>> s = p.get_structure("1fat", "1fat.gro")
        >>> io = GROIO()
        >>> io.set_structure(s)
        >>> io.save("out.gro")
        >>> import os
        >>> os.remove("out.gro")  # tidy up

    """

    # private methods
    def _get_atom_line(
        self,
        atom,
        hetfield,
        segid,
        atom_number,
        resname,
        resseq,
        icode,
        chain_id,
        element="  ",
        charge="  ",
    ):
        """Return an atom GRO string."""
        name = atom.get_name()
        x, y, z = atom.get_coord()
        gro_x = x / 10
        gro_y = y / 10
        gro_z = z / 10
        args = (resseq, resname, name, atom_number, gro_x, gro_y, gro_z)
        return _ATOM_FORMAT_STRING % args

    # Public methods

    def set_structure(self, structure):
        """Set the structure to be written out."""
        self.structure = structure

    def renumber(self, sort_inside_residues=False):
        """Renumber atoms and residues."""
        atoms_written = 0
        for model in self.structure.get_list():
            model_residues_written = 0
            new_res_list = []
            new_res_dict = {}
            for chain in model.get_list():
                for residue in chain:
                    if sort_inside_residues:
                        residue.sort()
                    hetfield, resseq, icode = residue.get_id()
                    resname = residue.get_resname()
                    segid = residue.get_segid()
                    residuerepresented = False
                    for atom in residue.get_unpacked_list():
                        if not residuerepresented:
                            model_residues_written += 1
                            residuerepresented = True
                            new_res_id = (
                                residue.get_id()[0],
                                model_residues_written,
                                residue.get_id()[2],
                            )
                            residue.id = new_res_id
                            new_res_list.append(residue)
                            new_res_dict[residue.get_id()] = residue
                        atoms_written += 1
                        atom.set_serial_number(atoms_written)
            chains_to_remove = model.child_list[1:]
            for chain in chains_to_remove:
                model.detach_child(chain.id)
            chain0 = model.child_list[0]
            chain0.id = " "
            # update parent child_dict with chain identifier change
            chain0.parent.child_dict.clear()
            chain0.parent.child_dict[chain0.id] = chain0
            chain0.child_list[:] = new_res_list
            chain0.child_dict.clear()
            chain0.child_dict.update(new_res_dict)

    def save(
        self,
        file,
        select=_select,
        sort_inside_residues=False,
        update_structobj_numbering=False,
    ):
        """Save the structure to a file.

        :param file: output file
        :type file: string or filehandle

        :param select: selects which entities will be written.
        :type select: object

        Typically select is a subclass of L{Select}, it should
        have the following methods:

         - accept_model(model)
         - accept_chain(chain)
         - accept_residue(residue)
         - accept_atom(atom)

        These methods should return 1 if the entity is to be
        written out, 0 otherwise.

        :param sort_inside_residues: whether atoms are sorted inside residues.
        :type sort_inside_residues: boolean

        :param update_structobj_numbering: whether to update structure object numbering.
        :type update_structobj_numbering: boolean
        """
        get_atom_line = self._get_atom_line
        if isinstance(file, str):
            fp = open(file, "w")
            close_file = 1
        else:
            # filehandle, I hope :-)
            fp = file
            close_file = 0
        atoms_count = {}
        for model in self.structure.get_list():
            if not select.accept_model(model):
                continue
            atoms_count[model.id] = 0
            for chain in model.get_list():
                if not select.accept_chain(chain):
                    continue
                for residue in chain.get_unpacked_list():
                    if not select.accept_residue(residue):
                        continue
                    for atom in residue.get_unpacked_list():
                        if select.accept_atom(atom):
                            atoms_count[model.id] += 1
        atoms_written = 0
        struct_get_list = self.structure.get_list()
        for i_model, model in enumerate(struct_get_list):
            if not select.accept_model(model):
                continue
            fp.write(f"{model.id}\n")
            fp.write(f"    {atoms_count[model.id]}\n")
            model_residues_written = 0
            for chain in model.get_list():
                if not select.accept_chain(chain):
                    continue
                chain_id = chain.get_id()
                for residue in chain.get_unpacked_list():
                    if sort_inside_residues:
                        residue.sort()
                    if not select.accept_residue(residue):
                        continue
                    hetfield, resseq, icode = residue.get_id()
                    resname = residue.get_resname()
                    segid = residue.get_segid()
                    residuerepresented = False
                    for atom in residue.get_unpacked_list():
                        if select.accept_atom(atom):
                            if not residuerepresented:
                                model_residues_written += 1
                                residuerepresented = True
                                if update_structobj_numbering:
                                    old_residue_id = residue.id
                                    new_residue_id = (
                                        residue.id[0],
                                        model_residues_written,
                                        residue.id[2],
                                    )
                                    if new_residue_id != old_residue_id:
                                        residue.id = new_residue_id
                                        del residue.parent.child_dict[old_residue_id]
                                        residue.parent.child_dict[new_residue_id] = (
                                            residue
                                        )
                            atoms_written += 1
                            s = get_atom_line(
                                atom,
                                hetfield,
                                segid,
                                atoms_written,
                                resname,
                                model_residues_written,
                                icode,
                                chain_id,
                            )
                            fp.write(s)
                            if update_structobj_numbering:
                                atom.set_serial_number(atoms_written)

            if hasattr(model, "box"):
                box = model.box
            else:
                box = (0.0, 0.0, 0.0)
            fp.write(" ".join(["%9.5f"] * len(box)) % tuple(box) + "\n")
        if close_file:
            fp.close()


if __name__ == "__main__":
    import sys
    from Bio.PDB.GROParser import GROParser

    g = GROParser(PERMISSIVE=1)

    s = g.get_structure("test", sys.argv[1])

    io = GROIO()
    io.set_structure(s)
    io.save("out1.gro")

    fg = open("out2.gro", "w")
    s1 = g.get_structure("test1", sys.argv[1])
    s2 = g.get_structure("test2", sys.argv[2])
    io = GROIO()
    io.set_structure(s1)
    io.save(fg)
    io.set_structure(s2)
    io.save(fg)
    fg.close()
