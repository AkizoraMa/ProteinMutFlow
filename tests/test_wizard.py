from __future__ import annotations

import tempfile
import unittest
import shutil
from pathlib import Path
from unittest import mock

import yaml

from mutflow.cli import build_parser
from mutflow.configuration import load_workflow, validate_workflow_data
from mutflow.preflight import run_preflight
from mutflow.wizard import InitError, run_init


ROOT = Path(__file__).resolve().parents[1]


class InitWizardTests(unittest.TestCase):
    @staticmethod
    def _input_from(responses: list[str]):
        iterator = iter(responses)

        def answer(prompt: str) -> str:
            try:
                return next(iterator)
            except StopIteration as exc:
                raise AssertionError(f"unexpected prompt: {prompt}") from exc

        return answer

    def test_minimal_explicit_metrics_only_workflow_is_created_once(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            destination = Path(temp) / "workflow.yaml"
            messages: list[str] = []
            responses = [
                "./wt.pdb",
                "",
                "",
                "./mutations.csv",
                "no",
                "",
                "",
                "",
            ]
            created = run_init(
                destination,
                input_fn=self._input_from(responses),
                output_fn=messages.append,
            )

            self.assertEqual(created, destination.resolve())
            config = load_workflow(destination)
            self.assertFalse(config["modeling"]["enabled"])
            self.assertEqual(config["variants"]["mode"], "explicit")
            self.assertEqual(config["metrics"]["sasa"], [])
            self.assertEqual(config["execution"]["max_variants"], 500)
            self.assertNotIn("backends", yaml.safe_load(destination.read_text(encoding="utf-8")))
            self.assertTrue(any(message.startswith("next: mutflow preflight") for message in messages))

    def test_saturation_workflow_collects_metric_qc_and_optional_discovery(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            destination = Path(temp) / "workflow.yaml"
            responses = [
                "./prepared.pdb",
                "",
                "saturation",
                "A:10:G, B:20A:D",
                "",
                "",
                "sasa",
                "site_sasa",
                "chain A and resi 10",
                "yes",
                "none",
                "yes",
                "expected",
                "site_backbone",
                "chain A and resi 10 and name N+CA+C+O",
                "4",
                "",
                "no",
                "",
                "",
                "",
            ]
            run_init(
                destination,
                input_fn=self._input_from(responses),
                output_fn=lambda message: None,
            )

            raw = yaml.safe_load(destination.read_text(encoding="utf-8"))
            config = load_workflow(destination)
            self.assertTrue(config["modeling"]["enabled"])
            self.assertEqual(config["variants"]["sites"][1]["insertion_code"], "A")
            self.assertTrue(config["metrics"]["sasa"][0]["calculate_wt_delta"])
            self.assertEqual(
                config["quality_control"]["required_selections"][0]["expected_count"],
                4,
            )
            self.assertNotIn("backends", raw)

    def test_existing_destination_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            destination = Path(temp) / "workflow.yaml"
            destination.write_text("owned by user\n", encoding="utf-8")
            with self.assertRaisesRegex(InitError, "refusing to overwrite"):
                run_init(
                    destination,
                    input_fn=lambda prompt: "",
                    output_fn=lambda message: None,
                )
            self.assertEqual(destination.read_text(encoding="utf-8"), "owned by user\n")

    def test_cli_exposes_init_with_default_destination(self) -> None:
        args = build_parser().parse_args(["init"])
        self.assertEqual(args.command, "init")
        self.assertEqual(args.destination, Path("workflow.yaml"))


class PublicExampleContractTests(unittest.TestCase):
    def test_public_1ubq_workflow_and_mutations_match_selected_fixture(self) -> None:
        example = ROOT / "examples" / "public_1ubq"
        raw = yaml.safe_load((example / "workflow.yaml").read_text(encoding="utf-8"))
        config = validate_workflow_data(raw)
        mutations = (example / "mutations.csv").read_text(encoding="utf-8").splitlines()

        self.assertEqual(config["input"]["structure"], "./1UBQ_prepared.pdb")
        self.assertEqual(config["variants"]["file"], "./mutations.csv")
        self.assertEqual(mutations[1:], ["M001,A:G10A", "M002,A:G10A_A:K11A"])
        self.assertEqual(len(config["metrics"]["sasa"]), 2)
        self.assertEqual(len(config["metrics"]["minimum_distance"]), 1)
        self.assertEqual(len(config["metrics"]["clashes"]), 1)
        self.assertEqual(
            config["quality_control"]["required_selections"][0]["expected_count"],
            8,
        )
        self.assertEqual(config["output"]["structure_formats"], ["pdb"])

    def test_public_1ubq_configuration_reaches_ready_with_non_scientific_probe_mock(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            case = Path(temp) / "public_1ubq"
            shutil.copytree(ROOT / "examples" / "public_1ubq", case)
            (case / "1UBQ_prepared.pdb").write_text(
                """ATOM      1  N   GLY A  10       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  GLY A  10       1.000   0.000   0.000  1.00 20.00           C
ATOM      3  C   GLY A  10       2.000   0.000   0.000  1.00 20.00           C
ATOM      4  O   GLY A  10       3.000   0.000   0.000  1.00 20.00           O
ATOM      5  N   LYS A  11       2.000   1.000   0.000  1.00 20.00           N
ATOM      6  CA  LYS A  11       3.000   1.000   0.000  1.00 20.00           C
ATOM      7  C   LYS A  11       4.000   1.000   0.000  1.00 20.00           C
ATOM      8  O   LYS A  11       5.000   1.000   0.000  1.00 20.00           O
ATOM      9  CB  LYS A  11       3.000   2.000   0.000  1.00 20.00           C
END
""",
                encoding="utf-8",
            )
            with mock.patch(
                "mutflow.preflight._backend_inspection",
                return_value=([], ["non-scientific mock"]),
            ):
                report = run_preflight(
                    case / "workflow.yaml",
                    run_backend_probes=False,
                    auto_discover_backends=False,
                )

            self.assertEqual(report.status, "READY")
            self.assertEqual(len(report.variants), 3)
            self.assertEqual(len(report.metric_names), 4)


if __name__ == "__main__":
    unittest.main()
