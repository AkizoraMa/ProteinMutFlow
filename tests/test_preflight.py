from __future__ import annotations

import unittest
from pathlib import Path
from unittest import mock

from mutflow.preflight import PlannedVariant, run_preflight


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


class PreflightTests(unittest.TestCase):
    def test_existing_structure_example_is_ready_and_writes_nothing(self) -> None:
        output = EXAMPLES / "demo_output"
        self.assertFalse(output.exists())
        report = run_preflight(EXAMPLES / "workflow_existing_structures.yaml")
        self.assertEqual(report.status, "READY")
        self.assertEqual(len(report.variants), 2)
        self.assertEqual(report.max_variants, 500)
        self.assertEqual(report.variants[0].mutation_name, "WT")
        self.assertEqual(report.variants[1].mutation_name, "A:S10A")
        self.assertFalse(output.exists())

    def test_saturation_expands_wt_plus_nineteen_without_test_backend_calls(self) -> None:
        report = run_preflight(
            EXAMPLES / "workflow_saturation.yaml",
            run_backend_probes=False,
            auto_discover_backends=False,
        )
        self.assertEqual(report.status, "CHECK")
        self.assertEqual(len(report.variants), 20)
        self.assertEqual(sum(item.mutation_count == 0 for item in report.variants), 1)
        self.assertEqual(sum(item.mutation_count == 1 for item in report.variants), 19)
        self.assertIn("SCHRODINGER_NOT_PROBED", {issue.code for issue in report.issues})

    def test_default_variant_limit_refuses_an_explicitly_oversized_batch(self) -> None:
        expanded = [
            PlannedVariant(f"V{index:04d}", f"synthetic-{index}", ())
            for index in range(501)
        ]
        with mock.patch(
            "mutflow.preflight._read_explicit_variants",
            return_value=expanded,
        ):
            report = run_preflight(EXAMPLES / "workflow_existing_structures.yaml")

        self.assertEqual(report.status, "REFUSED")
        self.assertIn("VARIANT_LIMIT_EXCEEDED", {issue.code for issue in report.issues})


if __name__ == "__main__":
    unittest.main()
