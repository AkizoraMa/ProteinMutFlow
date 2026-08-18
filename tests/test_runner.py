from __future__ import annotations

import csv
import json
import shutil
import tempfile
import unittest
from importlib import resources
from pathlib import Path
from unittest import mock

import yaml
from jsonschema import Draft202012Validator

from mutflow.backends import BackendLocation, BackendProbe
from mutflow.runner import RunExecutionError, run_workflow
from mutflow import runner as runner_module


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


class RunnerStateTests(unittest.TestCase):
    def _make_case(self, root: Path) -> Path:
        for name in ("synthetic_wt.pdb", "synthetic_S10A.pdb", "mutations_existing.csv"):
            shutil.copyfile(EXAMPLES / name, root / name)
        config = yaml.safe_load(
            (EXAMPLES / "workflow_existing_structures.yaml").read_text(encoding="utf-8")
        )
        config["output"]["directory"] = "./output"
        workflow = root / "workflow.yaml"
        workflow.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        return workflow

    def _make_modeling_case(self, root: Path) -> Path:
        shutil.copyfile(EXAMPLES / "synthetic_wt.pdb", root / "synthetic_wt.pdb")
        (root / "mutations_modeling.csv").write_text(
            "variant_id,mutations\nV001,A:S10A\nV002,A:S10A_A:D11E\n",
            encoding="utf-8",
        )
        config = yaml.safe_load(
            (EXAMPLES / "workflow_existing_structures.yaml").read_text(encoding="utf-8")
        )
        config["variants"] = {
            "mode": "explicit",
            "include_wt": True,
            "file": "./mutations_modeling.csv",
        }
        config["modeling"] = {
            "enabled": True,
            "backend": "schrodinger",
            "strategy": "sequential_from_wt",
            "local_minimization": {
                "radius_angstrom": 5.0,
                "include_nearby_nonprotein": True,
            },
        }
        config["metrics"] = {}
        config["quality_control"] = {"required_selections": []}
        config["output"]["directory"] = "./output"
        workflow = root / "workflow.yaml"
        workflow.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        return workflow

    @staticmethod
    def _successful_model_result(variant, input_structure: Path, staging: Path):
        output = staging / f"{variant.variant_id}.pdb"
        shutil.copyfile(input_structure, output)
        return {
            "variant_id": variant.variant_id,
            "status": "OK",
            "mutation_readback_status": (
                "NOT_REQUESTED" if variant.mutation_count == 0 else "OK"
            ),
            "issue_codes": [],
            "outputs": {"pdb": str(output)},
            "details": {},
            "duration_seconds": 0.1,
            "notes": "mock model completed",
        }

    def test_existing_structure_run_produces_only_compact_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            workflow = self._make_case(root)
            store = run_workflow(workflow)

            self.assertEqual(
                {item.name for item in store.output.iterdir()},
                {"structures", "results.csv", "run.json", "run.log"},
            )
            self.assertEqual(
                {item.name for item in store.structures.iterdir()},
                {"WT.pdb", "V001.pdb"},
            )
            with store.results_path.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertTrue(all(row["overall_status"] == "OK" for row in rows))
            self.assertEqual(store.run_document["state"], "COMPLETED")
            self.assertEqual(store.run_document["counts"]["ok"], 2)
            self.assertTrue(
                Path(
                    store.run_document["input"]["workflow"]["normalized_config"]["input"]["structure"]
                ).is_absolute()
            )
            self.assertFalse(any(root.glob(".*.tmp")))

            schema_root = resources.files("mutflow.schemas")
            run_schema = json.loads(
                schema_root.joinpath("run.schema.json").read_text(encoding="utf-8")
            )
            results_schema = json.loads(
                schema_root.joinpath("results.schema.json").read_text(encoding="utf-8")
            )
            Draft202012Validator(run_schema).validate(store.run_document)
            Draft202012Validator(results_schema).validate(store.rows)

    def test_second_run_refuses_existing_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            workflow = self._make_case(Path(temp))
            run_workflow(workflow)
            with self.assertRaisesRegex(RunExecutionError, "OUTPUT_NOT_EMPTY"):
                run_workflow(workflow)

    def test_one_copy_failure_is_isolated_and_reconciled(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            workflow = self._make_case(Path(temp))
            original_copy = runner_module.atomic_copy

            def fail_wt_only(source: Path, target: Path) -> None:
                if target.name == "WT.pdb":
                    raise OSError("synthetic copy failure")
                original_copy(source, target)

            with mock.patch("mutflow.runner.atomic_copy", side_effect=fail_wt_only):
                store = run_workflow(workflow)

            by_id = {row["variant_id"]: row for row in store.rows}
            self.assertEqual(by_id["WT"]["overall_status"], "FAILED")
            self.assertEqual(by_id["V001"]["overall_status"], "OK")
            self.assertEqual(store.run_document["state"], "FAILED")
            self.assertEqual(store.run_document["counts"]["failed"], 1)
            self.assertEqual(store.run_document["counts"]["ok"], 1)
            self.assertEqual(
                sum(store.run_document["counts"][key] for key in ("ok", "check", "failed", "not_run")),
                store.run_document["counts"]["requested"],
            )

    def test_mocked_metric_batch_populates_values_deltas_clashes_and_qc(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            workflow = self._make_case(root)
            config = yaml.safe_load(workflow.read_text(encoding="utf-8"))
            config["metrics"] = {
                "sasa": [
                    {
                        "name": "ligand_sasa",
                        "selection": "resn LIG",
                        "calculate_wt_delta": True,
                        "dot_solvent": True,
                        "solvent_radius_angstrom": 1.4,
                        "dot_density": 3,
                    }
                ],
                "minimum_distance": [
                    {
                        "name": "ligand_distance",
                        "selection_a": "resn LIG",
                        "selection_b": "chain A and resi 10",
                        "calculate_wt_delta": True,
                    }
                ],
                "clashes": [
                    {
                        "name": "ligand_clash",
                        "selection_a": "resn LIG",
                        "selection_b": "polymer.protein",
                    }
                ],
            }
            config["quality_control"] = {
                "required_selections": [
                    {
                        "name": "ligand",
                        "selection": "resn LIG",
                        "expected_count": 1,
                        "severity": "fail",
                    }
                ]
            }
            workflow.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

            location = BackendLocation("pymol", root / "python.exe", root, "mock")
            probe = BackendProbe(
                "pymol", "OK", location, {"pymol": "test", "python": "test"}
            )
            metric_results = {
                "WT": {
                    "variant_id": "WT",
                    "status": "OK",
                    "metrics": {
                        "ligand_sasa": 10.0,
                        "ligand_distance": 2.0,
                        "ligand_clash": {
                            "clash_count": 1,
                            "severe_clash_count": 0,
                            "max_overlap_angstrom": 0.5,
                            "min_contact_distance_angstrom": 2.2,
                            "worst_pair": "A--B",
                        },
                    },
                    "selection_counts": {"ligand": 1},
                    "error": "",
                },
                "V001": {
                    "variant_id": "V001",
                    "status": "OK",
                    "metrics": {
                        "ligand_sasa": 12.5,
                        "ligand_distance": 1.5,
                        "ligand_clash": {
                            "clash_count": 2,
                            "severe_clash_count": 1,
                            "max_overlap_angstrom": 1.1,
                            "min_contact_distance_angstrom": 1.3,
                            "worst_pair": "C--D",
                        },
                    },
                    "selection_counts": {"ligand": 1},
                    "error": "",
                },
            }
            with (
                mock.patch("mutflow.preflight._backend_inspection", return_value=([], ["mock"])),
                mock.patch("mutflow.runner.inspect_pymol", return_value=probe),
                mock.patch("mutflow.runner.run_pymol_metrics", return_value=metric_results),
            ):
                store = run_workflow(workflow)

            by_id = {row["variant_id"]: row for row in store.rows}
            self.assertEqual(by_id["WT"]["metric__ligand_sasa__delta_angstrom2"], 0.0)
            self.assertEqual(by_id["V001"]["metric__ligand_sasa__angstrom2"], 12.5)
            self.assertEqual(by_id["V001"]["metric__ligand_sasa__delta_angstrom2"], 2.5)
            self.assertEqual(by_id["V001"]["metric__ligand_distance__angstrom"], 1.5)
            self.assertEqual(by_id["V001"]["metric__ligand_distance__delta_angstrom"], -0.5)
            self.assertEqual(by_id["V001"]["metric__ligand_clash__clash_count"], 2)
            self.assertEqual(by_id["V001"]["metric__ligand_clash__worst_pair"], "C--D")
            self.assertEqual(by_id["V001"]["metrics_status"], "OK")
            self.assertEqual(by_id["V001"]["qc_status"], "OK")
            self.assertEqual(store.run_document["metrics"][0]["type"], "sasa")
            self.assertEqual(
                [metric["type"] for metric in store.run_document["metrics"]],
                ["sasa", "minimum_distance", "clash"],
            )
            self.assertEqual(store.run_document["environment"]["pymol"]["pymol"], "test")

    def test_mocked_schrodinger_modeling_promotes_only_final_structures(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            workflow = self._make_modeling_case(root)
            location = BackendLocation("schrodinger", root / "run.exe", root, "mock")
            probe = BackendProbe(
                "schrodinger",
                "OK",
                location,
                {"suite": "test", "python": "test"},
            )

            def fake_models(location, config, variants, input_structure, staging):
                for variant in variants:
                    output = staging / f"{variant.variant_id}.pdb"
                    shutil.copyfile(input_structure, output)
                    is_check = variant.variant_id == "V002"
                    yield {
                        "variant_id": variant.variant_id,
                        "status": "CHECK" if is_check else "OK",
                        "mutation_readback_status": (
                            "NOT_REQUESTED"
                            if variant.mutation_count == 0
                            else "CHECK"
                            if is_check
                            else "OK"
                        ),
                        "issue_codes": ["MUTATION_READBACK_HIS_ALIAS"] if is_check else [],
                        "outputs": {"pdb": str(output)},
                        "details": {
                            "selected_residue_count": 4,
                            "selected_atom_count": 20,
                            "selected_fraction": 0.25,
                            "fixed_heavy_atom_max_displacement_angstrom": 0.0,
                        }
                        if variant.mutation_count
                        else {},
                        "duration_seconds": 0.1,
                        "notes": "mock model completed",
                    }

            with (
                mock.patch("mutflow.preflight._backend_inspection", return_value=([], ["mock"])),
                mock.patch("mutflow.runner.inspect_schrodinger", return_value=probe),
                mock.patch("mutflow.runner.iter_schrodinger_models", side_effect=fake_models),
            ):
                store = run_workflow(workflow)

            self.assertEqual(
                {item.name for item in store.output.iterdir()},
                {"structures", "results.csv", "run.json", "run.log"},
            )
            self.assertEqual(
                {item.name for item in store.structures.iterdir()},
                {"WT.pdb", "V001.pdb", "V002.pdb"},
            )
            by_id = {row["variant_id"]: row for row in store.rows}
            self.assertEqual(by_id["V001"]["model_status"], "OK")
            self.assertEqual(by_id["V001"]["mutation_readback_status"], "OK")
            self.assertEqual(by_id["V002"]["overall_status"], "CHECK")
            self.assertEqual(store.run_document["counts"]["ok"], 2)
            self.assertEqual(store.run_document["counts"]["check"], 1)
            self.assertEqual(store.run_document["state"], "COMPLETED_WITH_CHECKS")
            self.assertEqual(
                store.run_document["execution"]["modeling_strategy"],
                "sequential_from_wt",
            )
            self.assertEqual(store.run_document["environment"]["schrodinger"]["suite"], "test")

    def test_one_mocked_model_failure_does_not_hide_other_variants(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            workflow = self._make_modeling_case(root)
            location = BackendLocation("schrodinger", root / "run.exe", root, "mock")
            probe = BackendProbe("schrodinger", "OK", location, {"suite": "test"})

            def fake_models(location, config, variants, input_structure, staging):
                for variant in variants:
                    if variant.variant_id == "V001":
                        yield {
                            "variant_id": variant.variant_id,
                            "status": "FAILED",
                            "mutation_readback_status": "FAILED",
                            "issue_codes": ["SCHRODINGER_MODEL_FAILED"],
                            "outputs": {},
                            "details": {},
                            "duration_seconds": 0.1,
                            "notes": "synthetic modeling failure",
                        }
                        continue
                    output = staging / f"{variant.variant_id}.pdb"
                    shutil.copyfile(input_structure, output)
                    yield {
                        "variant_id": variant.variant_id,
                        "status": "OK",
                        "mutation_readback_status": (
                            "NOT_REQUESTED" if variant.mutation_count == 0 else "OK"
                        ),
                        "issue_codes": [],
                        "outputs": {"pdb": str(output)},
                        "details": {},
                        "duration_seconds": 0.1,
                        "notes": "mock model completed",
                    }

            with (
                mock.patch("mutflow.preflight._backend_inspection", return_value=([], ["mock"])),
                mock.patch("mutflow.runner.inspect_schrodinger", return_value=probe),
                mock.patch("mutflow.runner.iter_schrodinger_models", side_effect=fake_models),
            ):
                store = run_workflow(workflow)

            by_id = {row["variant_id"]: row for row in store.rows}
            self.assertEqual(by_id["V001"]["overall_status"], "FAILED")
            self.assertEqual(by_id["V002"]["overall_status"], "OK")
            self.assertEqual(store.run_document["counts"]["failed"], 1)
            self.assertEqual(store.run_document["counts"]["ok"], 2)

    def test_interrupted_modeling_resumes_only_not_run_variants(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            workflow = self._make_modeling_case(root)
            location = BackendLocation("schrodinger", root / "run.exe", root, "mock")
            probe = BackendProbe("schrodinger", "OK", location, {"suite": "test"})
            calls: list[list[str]] = []

            def fake_models(location, config, variants, input_structure, staging):
                calls.append([variant.variant_id for variant in variants])
                attempt = len(calls)
                for variant in variants:
                    if attempt == 1 and variant.variant_id != "WT":
                        raise KeyboardInterrupt
                    yield self._successful_model_result(
                        variant, input_structure, staging
                    )

            with (
                mock.patch("mutflow.preflight._backend_inspection", return_value=([], ["mock"])),
                mock.patch("mutflow.runner.inspect_schrodinger", return_value=probe),
                mock.patch("mutflow.runner.iter_schrodinger_models", side_effect=fake_models),
            ):
                with self.assertRaises(KeyboardInterrupt):
                    run_workflow(workflow)

                output = root / "output"
                interrupted = json.loads(
                    (output / "run.json").read_text(encoding="utf-8")
                )
                self.assertEqual(interrupted["state"], "INTERRUPTED")
                with (output / "results.csv").open(
                    "r", encoding="utf-8", newline=""
                ) as handle:
                    before = {row["variant_id"]: row for row in csv.DictReader(handle)}
                self.assertEqual(before["WT"]["overall_status"], "OK")
                self.assertEqual(before["V001"]["overall_status"], "NOT_RUN")

                store = run_workflow(workflow, resume=True)

            after = {row["variant_id"]: row for row in store.rows}
            self.assertEqual(calls, [["WT", "V001", "V002"], ["V001", "V002"]])
            self.assertEqual(after["WT"], {
                **before["WT"],
                "mutation_count": 0,
                "duration_seconds": 0.1,
            })
            self.assertTrue(all(row["overall_status"] == "OK" for row in store.rows))
            self.assertEqual(store.run_document["state"], "COMPLETED")
            self.assertEqual(store.run_document["execution"]["new_or_resumed"], "resumed")
            self.assertEqual(store.run_document["execution"]["resume_count"], 1)
            self.assertEqual(
                {item.name for item in store.output.iterdir()},
                {"structures", "results.csv", "run.json", "run.log"},
            )

    def test_resume_refuses_changed_workflow_or_input_without_touching_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            workflow = self._make_modeling_case(root)
            location = BackendLocation("schrodinger", root / "run.exe", root, "mock")
            probe = BackendProbe("schrodinger", "OK", location, {"suite": "test"})

            def interrupt_after_wt(location, config, variants, input_structure, staging):
                for variant in variants:
                    if variant.variant_id != "WT":
                        raise KeyboardInterrupt
                    yield self._successful_model_result(
                        variant, input_structure, staging
                    )

            with (
                mock.patch("mutflow.preflight._backend_inspection", return_value=([], ["mock"])),
                mock.patch("mutflow.runner.inspect_schrodinger", return_value=probe),
                mock.patch(
                    "mutflow.runner.iter_schrodinger_models",
                    side_effect=interrupt_after_wt,
                ),
            ):
                with self.assertRaises(KeyboardInterrupt):
                    run_workflow(workflow)

            output = root / "output"
            original_workflow = workflow.read_text(encoding="utf-8")
            before = {
                name: (output / name).read_bytes()
                for name in ("results.csv", "run.json", "run.log")
            }
            config = yaml.safe_load(workflow.read_text(encoding="utf-8"))
            config["modeling"]["local_minimization"]["radius_angstrom"] = 6.0
            workflow.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
            with (
                mock.patch("mutflow.preflight._backend_inspection", return_value=([], ["mock"])),
                mock.patch("mutflow.runner.inspect_schrodinger", return_value=probe),
            ):
                with self.assertRaisesRegex(
                    RunExecutionError, "resume compatibility check failed"
                ):
                    run_workflow(workflow, resume=True)
            after = {
                name: (output / name).read_bytes()
                for name in ("results.csv", "run.json", "run.log")
            }
            self.assertEqual(after, before)

            workflow.write_text(original_workflow, encoding="utf-8")
            input_structure = root / "synthetic_wt.pdb"
            input_structure.write_text(
                input_structure.read_text(encoding="utf-8") + "REMARK HASH CHANGE\n",
                encoding="utf-8",
            )
            with (
                mock.patch("mutflow.preflight._backend_inspection", return_value=([], ["mock"])),
                mock.patch("mutflow.runner.inspect_schrodinger", return_value=probe),
            ):
                with self.assertRaisesRegex(
                    RunExecutionError, "resume compatibility check failed"
                ):
                    run_workflow(workflow, resume=True)
            final = {
                name: (output / name).read_bytes()
                for name in ("results.csv", "run.json", "run.log")
            }
            self.assertEqual(final, before)

    def test_resume_does_not_retry_terminal_failed_variants(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            workflow = self._make_modeling_case(root)
            location = BackendLocation("schrodinger", root / "run.exe", root, "mock")
            probe = BackendProbe("schrodinger", "OK", location, {"suite": "test"})
            calls: list[list[str]] = []

            def fake_models(location, config, variants, input_structure, staging):
                calls.append([variant.variant_id for variant in variants])
                if len(calls) == 1:
                    yield self._successful_model_result(
                        variants[0], input_structure, staging
                    )
                    yield {
                        "variant_id": "V001",
                        "status": "FAILED",
                        "mutation_readback_status": "FAILED",
                        "issue_codes": ["SCHRODINGER_MODEL_FAILED"],
                        "outputs": {},
                        "details": {},
                        "duration_seconds": 0.1,
                        "notes": "synthetic terminal failure",
                    }
                    raise KeyboardInterrupt
                for variant in variants:
                    yield self._successful_model_result(
                        variant, input_structure, staging
                    )

            with (
                mock.patch("mutflow.preflight._backend_inspection", return_value=([], ["mock"])),
                mock.patch("mutflow.runner.inspect_schrodinger", return_value=probe),
                mock.patch("mutflow.runner.iter_schrodinger_models", side_effect=fake_models),
            ):
                with self.assertRaises(KeyboardInterrupt):
                    run_workflow(workflow)
                store = run_workflow(workflow, resume=True)

            by_id = {row["variant_id"]: row for row in store.rows}
            self.assertEqual(calls[1], ["V002"])
            self.assertEqual(by_id["V001"]["overall_status"], "FAILED")
            self.assertEqual(by_id["V001"]["notes"], "synthetic terminal failure")
            self.assertEqual(by_id["V002"]["overall_status"], "OK")
            self.assertEqual(store.run_document["state"], "FAILED")

    def test_interrupted_metric_batch_resumes_without_recopying_structures(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            workflow = self._make_case(root)
            config = yaml.safe_load(workflow.read_text(encoding="utf-8"))
            config["metrics"] = {
                "sasa": [
                    {
                        "name": "protein_sasa",
                        "selection": "polymer.protein",
                        "calculate_wt_delta": True,
                    }
                ]
            }
            workflow.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
            location = BackendLocation("pymol", root / "python.exe", root, "mock")
            probe = BackendProbe("pymol", "OK", location, {"pymol": "test"})
            metric_results = {
                variant_id: {
                    "variant_id": variant_id,
                    "status": "OK",
                    "metrics": {"protein_sasa": value},
                    "selection_counts": {},
                    "error": "",
                }
                for variant_id, value in (("WT", 100.0), ("V001", 104.0))
            }
            with (
                mock.patch("mutflow.preflight._backend_inspection", return_value=([], ["mock"])),
                mock.patch("mutflow.runner.inspect_pymol", return_value=probe),
                mock.patch(
                    "mutflow.runner.run_pymol_metrics",
                    side_effect=[KeyboardInterrupt(), metric_results],
                ),
            ):
                with self.assertRaises(KeyboardInterrupt):
                    run_workflow(workflow)
                output = root / "output"
                structures_before = {
                    item.name: item.read_bytes()
                    for item in (output / "structures").iterdir()
                }
                store = run_workflow(workflow, resume=True)

            structures_after = {
                item.name: item.read_bytes() for item in store.structures.iterdir()
            }
            by_id = {row["variant_id"]: row for row in store.rows}
            self.assertEqual(structures_after, structures_before)
            self.assertEqual(
                by_id["V001"]["metric__protein_sasa__delta_angstrom2"], 4.0
            )
            self.assertEqual(store.run_document["state"], "COMPLETED")


if __name__ == "__main__":
    unittest.main()
