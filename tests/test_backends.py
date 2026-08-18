from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from mutflow.backends import discover_pymol, discover_schrodinger, inspect_pymol


class BackendDiscoveryTests(unittest.TestCase):
    def test_explicit_schrodinger_root_resolves_run_exe(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "run.exe").touch()
            found = discover_schrodinger(str(root))
            self.assertIsNotNone(found)
            assert found is not None
            self.assertEqual(found.source, "workflow")
            self.assertEqual(found.launcher, (root / "run.exe").resolve())

    def test_explicit_conda_pymol_root_resolves_python(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "conda-meta").mkdir()
            (root / "conda-meta" / "pymol-3.0.0-test.json").write_text("{}", encoding="utf-8")
            (root / "python.exe").touch()
            found = discover_pymol(str(root))
            self.assertIsNotNone(found)
            assert found is not None
            self.assertEqual(found.launcher, (root / "python.exe").resolve())

    def test_direct_pymolwin_is_refused_without_launching(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            launcher = Path(temp) / "PyMOLWin.exe"
            launcher.touch()
            probe = inspect_pymol(str(launcher))
            self.assertEqual(probe.status, "FAILED")
            self.assertIn("unsafe", probe.message.lower())


if __name__ == "__main__":
    unittest.main()
