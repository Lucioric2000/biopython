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
from Bio.PDB.parse_gro_header import parse_gro_header


class ParseSmallGRO_tests(unittest.TestCase):
    """Testing with real GRO files."""

    @classmethod
    def setUpClass(cls):
        cls.permissive = GROParser()
        cls.strict = GROParser(PERMISSIVE=0)

    def test_empty(self):
        """Parse an empty file."""
        handle = StringIO()
        with self.assertRaises(ValueError) as context_manager:
            _ = self.permissive.get_structure("MT", handle)
        self.assertEqual(str(context_manager.exception), "Empty file.")

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
        model = structure[0.0]
        self.assertEqual(model.id, 0.0)
        self.assertEqual(model.level, "M")
        self.assertEqual(len(model), 1)
        chain = model[" "]
        self.assertEqual(chain.id, " ")
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

    def test_parse_gro_header(self):
        """Parse GRO header."""
        header = parse_gro_header("PDB/waters.gro")
        self.assertEqual(header["name"].strip(), "MD of 2 waters, t= 0.0")
        self.assertEqual(header["atom_count"], 6)
        self.assertEqual(header["box"].strip(), "1.82060   1.82060   1.82060")

    def test_gro_io_save(self):
        """Save and reload structure using GROIO."""
        structure = self.strict.get_structure("example", "PDB/waters.gro")
        with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".gro") as tmp:
            tmp_path = tmp.name
        try:
            io = GROIO()
            io.set_structure(structure)
            io.save(tmp_path)

            reloaded = self.strict.get_structure("reloaded", tmp_path)
            self.assertEqual(len(reloaded), 1)
            reloaded_model = reloaded[0.0]
            self.assertEqual(len(reloaded_model), 1)
            reloaded_chain = reloaded_model[" "]
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
