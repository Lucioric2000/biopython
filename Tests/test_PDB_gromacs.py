# Copyright 2026 by Biopython contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Unit tests for the Bio.PDB GROParser and GROIO modules."""

from io import StringIO
import os
import tempfile
import unittest

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None

from Bio.PDB.GROIO import GROIO
from Bio.PDB.GROParser import GROParser


class ParseSmallGRO_tests(unittest.TestCase):
    """Testing with real GRO files."""

    @classmethod
    def setUpClass(cls):
        cls.permissive = GROParser()
        cls.strict = GROParser(PERMISSIVE=0)

    def test_SMCRA(self):
        """Walk down the structure hierarchy and test parser reliability."""
        s = self.permissive.get_structure("scr", "PDB/waters.gro")
        for m in s:
            p = m.get_parent()
            self.assertEqual(s, p)
            for c in m:
                p = c.get_parent()
                self.assertEqual(m, p)
                for r in c:
                    p = r.get_parent()
                    self.assertEqual(c, p)
                    for a in r:
                        p = a.get_parent()
                        self.assertEqual(r.get_resname(), p.get_resname())

    def test_waters_strict(self):
        """Parse waters.gro file in strict mode."""
        structure = self.strict.get_structure("example", "PDB/waters.gro")
        self.assertEqual(len(structure), 1)
        model = structure["MD of 2 waters, t= 0.0"]
        self.assertEqual(model.id, "MD of 2 waters, t= 0.0")
        self.assertEqual(model.level, "M")
        self.assertEqual(len(model), 1)
        chain = next(model.get_chains())
        self.assertEqual(chain.level, "C")
        self.assertEqual(len(chain), 2)
        self.assertEqual(
            " ".join(residue.resname for residue in chain),
            "WATER WATER",
        )
        self.assertEqual(
            " ".join(atom.name for atom in chain.get_atoms()),
            "OW1 HW2 HW3 OW1 HW2 HW3",
        )
        self.assertEqual(
            " ".join(atom.element for atom in chain.get_atoms()),
            "O H H O H H",
        )

        # Asserts that all the atoms are read with the correct occupancy (1.0) and bfactor (0.0)
        for atom in chain.get_atoms():
            assert np.isclose(atom.bfactor, 0.0)
            assert np.isclose(atom.occupancy, 1.0)

        self.assertEqual(len(model.box), 3)
        self.assertAlmostEqual(model.box[0], 1.82060)
        self.assertAlmostEqual(model.box[1], 1.82060)
        self.assertAlmostEqual(model.box[2], 1.82060)

        # Tests the chain residue 2, atom 2
        residue2 = list(chain.get_residues())[1]
        atom2res2 = list(residue2.get_atoms())[1]
        assert atom2res2.name == "HW2"
        assert atom2res2.parent.resname == "WATER"
        residue_id = atom2res2.parent.get_id()
        assert residue_id[1] == 2
        assert np.isclose(atom2res2.coord, [13.37, 0.02, 6.80]).all()

    def test_gro_io_save(self):
        """Save and reload structure using GROIO."""
        structure = self.strict.get_structure("example", "PDB/waters.gro")
        with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".gro") as tmp:
            tmp_path = tmp.name
        try:
            io = GROIO()
            io.set_structure(structure)
            io.save(tmp_path)
            with open(tmp_path) as f2, open("PDB/waters.gro") as f1:
                struct1 = list(f1.readlines())
                struct2 = list(f2.readlines())
                assert struct1[0].strip() == struct2[0].strip()  # Title string
                assert (
                    struct1[1].strip() == struct2[1].strip()
                )  # Number of atoms string
                # assert the box numbers
                box1_numbers = [float(x) for x in struct1[-1].strip().split()]
                box2_numbers = [float(x) for x in struct2[-1].strip().split()]
                assert np.isclose(box1_numbers, box2_numbers).all()
                structure1_atoms = list(structure.get_atoms())
                assert len(structure1_atoms) == len(struct1) - 3
                for i_atom_line in range(2, len(struct1) - 1):
                    structure1_atom = structure1_atoms[i_atom_line - 2]
                    items_atom_struct_1 = struct1[i_atom_line].strip().split()
                    items_atom_struct_2 = struct2[i_atom_line].strip().split()
                    # we convert the coordinates to angstroms for testing
                    residue_id = structure1_atom.parent.get_id()
                    assert (
                        f"{residue_id[1]}{structure1_atom.parent.resname}"
                        == items_atom_struct_2[0]
                    )
                    assert structure1_atom.name == items_atom_struct_2[1]
                    coords = [float(c) * 10 for c in items_atom_struct_1[3:6]]
                    assert np.isclose(structure1_atom.coord, coords).all(), (
                        structure1_atom.coord,
                        coords,
                    )
                    # Checks that the coordinates are saved in the re-created structure, but not
                    # the velocities (list elements from 6 onwards)
                    assert items_atom_struct_1[0:3] == items_atom_struct_2[0:3]
                    coords_struct_1 = [float(c) for c in items_atom_struct_1[3:6]]
                    coords_struct_2 = [float(c) for c in items_atom_struct_2[3:6]]
                    assert np.isclose(coords_struct_1, coords_struct_2).all()
                    # Check that the velocities were not saved in the output GROMACS file
                    assert len(items_atom_struct_2) == 6
            reloaded = self.strict.get_structure("reloaded", tmp_path)
            self.assertEqual(len(reloaded), 1)
            reloaded_model = reloaded["MD of 2 waters, t= 0.0"]
            self.assertEqual(len(reloaded_model), 1)
            reloaded_chain = next(reloaded_model.get_chains())
            self.assertEqual(len(reloaded_chain), 2)
            self.assertEqual(
                " ".join(atom.name for atom in reloaded_chain.get_atoms()),
                "OW1 HW2 HW3 OW1 HW2 HW3",
            )
            # Tests the chain residue 2, atom 2
            residue2 = list(reloaded_chain.get_residues())[1]
            atom2res2 = list(residue2.get_atoms())[1]
            assert atom2res2.name == "HW2"
            assert atom2res2.parent.resname == "WATER"
            residue_id = atom2res2.parent.get_id()
            assert residue_id[1] == 2
            assert np.isclose(atom2res2.coord, [13.37, 0.02, 6.80]).all()
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
