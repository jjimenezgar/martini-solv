import io
import unittest
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory

from martini_solv.short_md import archive


class TestArchive(unittest.TestCase):
    def test_simulation_files_layout_and_topology_rewrite(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "system.gro").write_text("test\n1\n    1PROT   BB    1   0.0     0.0     0.0\n   1.0 1.0 1.0\n")
            (root / "system.top").write_text(
                '#include "martini_v3.0.0.itp"\n'
                '#include "Protein.itp"\n\n'
                '[ system ]\nTest\n\n'
                '[ molecules ]\nProtein 1\n'
            )
            (root / "martini_v3.0.0.itp").write_text("[ defaults ]\n")
            (root / "Protein.itp").write_text("[ moleculetype ]\nProtein 1\n")
            (root / "minimization.mdp").write_text("integrator = steep\n")
            (root / "manifest.json").write_text("{}\n")
            (root / "build.log").write_text("ok\n")

            payload = archive(root)
            with zipfile.ZipFile(io.BytesIO(payload)) as zipped:
                names = set(zipped.namelist())
                self.assertIn("Simulation_Files/0_topology/system.top", names)
                self.assertIn("Simulation_Files/0_topology/system_itp/Protein.itp", names)
                self.assertIn("Simulation_Files/1_mdp/minimization.mdp", names)
                self.assertIn("Simulation_Files/1_mdp/nvt.mdp", names)
                self.assertIn("Simulation_Files/1_mdp/npt.mdp", names)
                self.assertIn("Simulation_Files/1_mdp/production.mdp", names)
                self.assertIn("Simulation_Files/2_system/system.gro", names)
                self.assertIn("Simulation_Files/metadata/manifest.json", names)
                self.assertIn("Simulation_Files/README.txt", names)

                top = zipped.read("Simulation_Files/0_topology/system.top").decode()
                self.assertIn('#include "system_itp/Protein.itp"', top)
                self.assertIn('#include "system_itp/martini_v3.0.0.itp"', top)
