from __future__ import annotations

import unittest
from types import SimpleNamespace

from mutflow.pymol_worker import _clash_metrics, _minimum_distance


def atom(index: int, coord: tuple[float, float, float], symbol: str, name: str):
    return SimpleNamespace(
        index=index,
        coord=coord,
        symbol=symbol,
        name=name,
        chain="A",
        resn="TST",
        resi=str(index),
    )


class FakeCmd:
    def __init__(self) -> None:
        self.selections = {
            "left": [atom(1, (0.0, 0.0, 0.0), "C", "C1")],
            "right": [atom(2, (1.5, 0.0, 0.0), "O", "O1")],
        }

    def get_model(self, selection: str):
        return SimpleNamespace(atom=self.selections[selection])


class WorkerGeometryTests(unittest.TestCase):
    def test_minimum_distance_uses_distinct_heavy_atoms(self) -> None:
        value = _minimum_distance(
            FakeCmd(),
            {"name": "test_distance", "selection_a": "left", "selection_b": "right"},
        )
        self.assertEqual(value, 1.5)

    def test_clash_uses_configured_radii_and_thresholds(self) -> None:
        result = _clash_metrics(
            FakeCmd(),
            {
                "name": "test_clash",
                "selection_a": "left",
                "selection_b": "right",
                "vdw_radii": {"C": 1.70, "O": 1.52},
                "clash_overlap_angstrom": 0.4,
                "severe_overlap_angstrom": 0.8,
            },
        )
        self.assertEqual(result["clash_count"], 1)
        self.assertEqual(result["severe_clash_count"], 1)
        self.assertAlmostEqual(result["max_overlap_angstrom"], 1.72)
        self.assertEqual(result["min_contact_distance_angstrom"], 1.5)
        self.assertIn("C1", result["worst_pair"])


if __name__ == "__main__":
    unittest.main()
