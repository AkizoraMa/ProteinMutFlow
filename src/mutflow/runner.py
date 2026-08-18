from __future__ import annotations

import hashlib
import sys
import tempfile
import time
import uuid
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from mutflow import __version__
from mutflow.backends import BackendProbe, inspect_pymol, inspect_schrodinger
from mutflow.configuration import load_workflow, resolve_workflow_paths
from mutflow.metrics import (
    MetricExecutionError,
    metric_definitions,
    metric_result_columns,
    normalized_metric_name,
    run_pymol_metrics,
)
from mutflow.modeling import ModelingExecutionError, iter_schrodinger_models
from mutflow.preflight import PlannedVariant, PreflightError, run_preflight
from mutflow.state import RESULT_FIELDS, RunStore, atomic_copy, utc_now


class RunExecutionError(RuntimeError):
    """Raised when the requested workflow is outside the implemented run slice."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve(base: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _run_id() -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}-{uuid.uuid4().hex[:8]}"


def _initial_row(variant: PlannedVariant, metric_columns: list[str]) -> dict[str, Any]:
    row = {field: "" for field in [*RESULT_FIELDS, *metric_columns]}
    row.update(
        {
            "schema_version": "1.0",
            "variant_id": variant.variant_id,
            "mutation_name": variant.mutation_name,
            "mutation_count": variant.mutation_count,
            "model_status": "NOT_RUN",
            "mutation_readback_status": "NOT_RUN",
            "metrics_status": "NOT_RUN",
            "qc_status": "NOT_RUN",
            "overall_status": "NOT_RUN",
            "duration_seconds": None,
        }
    )
    return row


def _source_for_variant(variant: PlannedVariant, wt_input: Path) -> Path:
    if variant.mutation_count == 0:
        return wt_input
    if variant.source_structure is None:
        raise RunExecutionError(f"variant has no existing source structure: {variant.variant_id}")
    return variant.source_structure


def _input_record(
    workflow: Path,
    config: dict[str, Any],
    variants: list[PlannedVariant],
    wt_input: Path,
) -> dict[str, Any]:
    base = workflow.parent
    variants_file: dict[str, str] | None = None
    if config["variants"]["mode"] == "explicit":
        path = _resolve(base, config["variants"]["file"])
        variants_file = {"path": str(path), "sha256": _sha256(path)}
    variant_structures = [
        {
            "variant_id": variant.variant_id,
            "path": str(variant.source_structure),
            "sha256": _sha256(variant.source_structure),
        }
        for variant in variants
        if variant.source_structure is not None
    ]
    return {
        "structure": {
            "path": str(wt_input),
            "format": "pdb",
            "sha256": _sha256(wt_input),
        },
        "workflow": {
            "path": str(workflow),
            "sha256": _sha256(workflow),
            "normalized_config": config,
        },
        "variants_file": variants_file,
        "variant_structures": variant_structures,
    }


def _run_document(
    workflow: Path,
    config: dict[str, Any],
    variants: list[PlannedVariant],
    wt_input: Path,
    command: str,
    schrodinger_probe: BackendProbe | None,
    pymol_probe: BackendProbe | None,
) -> dict[str, Any]:
    by_count = Counter(variant.mutation_count for variant in variants)
    return {
        "schema_version": "1.0",
        "run_id": _run_id(),
        "mutflow_version": __version__,
        "state": "RUNNING",
        "started_at": utc_now(),
        "finished_at": None,
        "command": command,
        "input": _input_record(workflow, config, variants, wt_input),
        "variants": {
            "mode": config["variants"]["mode"],
            "include_wt": config["variants"]["include_wt"],
            "requested": len(variants),
            "by_mutation_count": {str(key): by_count[key] for key in sorted(by_count)},
            "items": [
                {
                    "variant_id": variant.variant_id,
                    "mutation_name": variant.mutation_name,
                    "mutation_count": variant.mutation_count,
                }
                for variant in variants
            ],
        },
        "environment": {
            "operating_system": sys.platform,
            "core_python": sys.version.split()[0],
            "schrodinger": (
                {
                    "available": True,
                    "requested": True,
                    "launcher": str(schrodinger_probe.location.launcher),
                    "source": schrodinger_probe.location.source,
                    "headless_probe": "OK",
                    **schrodinger_probe.versions,
                }
                if schrodinger_probe is not None
                and schrodinger_probe.location is not None
                else {"available": False, "requested": False}
            ),
            "pymol": (
                {
                    "available": True,
                    "requested": True,
                    "launcher": str(pymol_probe.location.launcher),
                    "source": pymol_probe.location.source,
                    "headless_probe": "OK",
                    **pymol_probe.versions,
                }
                if pymol_probe is not None and pymol_probe.location is not None
                else {"available": False, "requested": False}
            ),
        },
        "execution": {
            "modeling_backend": "schrodinger" if config["modeling"]["enabled"] else None,
            "modeling_strategy": (
                config["modeling"]["strategy"] if config["modeling"]["enabled"] else None
            ),
            "local_minimization": (
                dict(config["modeling"]["local_minimization"])
                if config["modeling"]["enabled"]
                else None
            ),
            "new_or_resumed": "new",
        },
        "metrics": metric_definitions(config),
        "counts": {
            "requested": len(variants),
            "ok": 0,
            "check": 0,
            "failed": 0,
            "not_run": len(variants),
        },
        "outputs": {
            "structures_directory": "structures",
            "results_csv": "results.csv",
            "log": "run.log",
        },
    }


def _schema_major(value: Any) -> str:
    return str(value).split(".", maxsplit=1)[0]


def _validate_resume(
    store: RunStore,
    expected_document: dict[str, Any],
    expected_fieldnames: list[str],
    variants: list[PlannedVariant],
    structure_formats: list[str],
) -> list[PlannedVariant]:
    errors: list[str] = []
    existing = store.run_document
    if _schema_major(existing.get("schema_version")) != _schema_major(
        expected_document["schema_version"]
    ):
        errors.append("run schema major version differs")
    if existing.get("mutflow_version") != expected_document["mutflow_version"]:
        errors.append("MutFlow version differs")
    if existing.get("input") != expected_document["input"]:
        errors.append("input, workflow, or normalized configuration differs")
    if existing.get("variants") != expected_document["variants"]:
        errors.append("expanded variant mapping differs")
    if existing.get("metrics") != expected_document["metrics"]:
        errors.append("metric definitions differ")
    if store.fieldnames != expected_fieldnames:
        errors.append("results.csv header differs")
    if existing.get("state") not in {"INTERRUPTED", "RUNNING"}:
        errors.append(f"run state {existing.get('state')!r} is not resumable")

    expected_by_id = {variant.variant_id: variant for variant in variants}
    existing_ids = [str(row.get("variant_id", "")) for row in store.rows]
    if len(existing_ids) != len(set(existing_ids)):
        errors.append("results.csv contains duplicate variant IDs")
    if set(existing_ids) != set(expected_by_id):
        errors.append("results.csv variant IDs differ")

    allowed_structure_names = {
        f"{variant.variant_id}.{output_format}"
        for variant in variants
        for output_format in structure_formats
    }
    actual_structure_names = {item.name for item in store.structures.iterdir()}
    unexpected_structures = actual_structure_names - allowed_structure_names
    if unexpected_structures:
        errors.append(
            f"structures directory contains unexpected files: {sorted(unexpected_structures)}"
        )
    nonfiles = [item.name for item in store.structures.iterdir() if not item.is_file()]
    if nonfiles:
        errors.append(f"structures directory contains non-files: {sorted(nonfiles)}")

    pending: list[PlannedVariant] = []
    structures_root = store.structures.resolve()
    for row in store.rows:
        variant_id = str(row.get("variant_id", ""))
        variant = expected_by_id.get(variant_id)
        if variant is None:
            continue
        if _schema_major(row.get("schema_version")) != "1":
            errors.append(f"{variant_id}: results schema major version differs")
        try:
            mutation_count = int(row.get("mutation_count"))
        except (TypeError, ValueError):
            mutation_count = -1
        if (
            row.get("mutation_name") != variant.mutation_name
            or mutation_count != variant.mutation_count
        ):
            errors.append(f"{variant_id}: mutation mapping differs")

        overall = row.get("overall_status")
        if overall == "NOT_RUN":
            model_status = row.get("model_status")
            metrics_status = row.get("metrics_status")
            model_pending = model_status == "NOT_RUN"
            metrics_pending = (
                model_status in {"OK", "NOT_REQUESTED"}
                and metrics_status == "NOT_RUN"
                and bool(row.get("structure_path"))
            )
            if not (model_pending or metrics_pending):
                errors.append(f"{variant_id}: NOT_RUN stage state is inconsistent")
            pending.append(variant)
        elif overall not in {"OK", "CHECK", "FAILED"}:
            errors.append(f"{variant_id}: unknown overall status {overall!r}")

        structure_path = str(row.get("structure_path") or "")
        if not structure_path:
            if row.get("model_status") in {"OK", "NOT_REQUESTED"}:
                errors.append(f"{variant_id}: completed modeling has no structure path")
            continue
        candidate = (store.output / structure_path).resolve()
        try:
            candidate.relative_to(structures_root)
        except ValueError:
            errors.append(f"{variant_id}: structure path escapes structures directory")
            continue
        if not candidate.is_file() or candidate.stat().st_size == 0:
            errors.append(f"{variant_id}: referenced structure is missing or empty")
        if row.get("model_status") in {"OK", "NOT_REQUESTED"}:
            for output_format in structure_formats:
                expected_structure = store.structures / f"{variant_id}.{output_format}"
                if not expected_structure.is_file() or expected_structure.stat().st_size == 0:
                    errors.append(
                        f"{variant_id}: expected {output_format} structure is missing or empty"
                    )

    if errors:
        raise RunExecutionError("resume compatibility check failed: " + "; ".join(errors))
    if not pending:
        raise RunExecutionError("RESUME_NOT_NEEDED: no NOT_RUN variants remain")
    return pending


def _selection_qc(
    counts: dict[str, Any], checks: list[dict[str, Any]]
) -> tuple[str, list[str], list[str]]:
    if not checks:
        return "NOT_REQUESTED", [], []
    status = "OK"
    codes: list[str] = []
    notes: list[str] = []
    for check in checks:
        name = check["name"]
        count = int(counts.get(name, -1))
        passed = True
        if "expected_count" in check:
            passed = count == check["expected_count"]
        if "minimum_count" in check:
            passed = passed and count >= check["minimum_count"]
        if "maximum_count" in check:
            passed = passed and count <= check["maximum_count"]
        if passed:
            continue
        slug = normalized_metric_name(name).upper()
        codes.append(f"QC_SELECTION_{slug}_COUNT")
        notes.append(f"selection {name} count={count}")
        if check["severity"] == "fail":
            status = "FAILED"
        elif status != "FAILED":
            status = "CHECK"
    return status, codes, notes


def run_workflow(
    workflow_path: Path,
    command: str | None = None,
    *,
    resume: bool = False,
) -> RunStore:
    workflow = workflow_path.expanduser().resolve()
    try:
        report = run_preflight(workflow, resume=resume)
    except PreflightError as exc:
        raise RunExecutionError(str(exc)) from exc
    if report.status != "READY":
        issue_text = "; ".join(f"{issue.code}: {issue.message}" for issue in report.issues)
        raise RunExecutionError(f"preflight status is {report.status}: {issue_text}")

    config = load_workflow(workflow)
    resolved_config = resolve_workflow_paths(config, workflow)
    modeling_enabled = config["modeling"]["enabled"]
    pymol_requested = bool(
        report.metric_names
        or config["quality_control"]["required_selections"]
    )
    schrodinger_probe: BackendProbe | None = None
    if modeling_enabled:
        explicit = resolved_config["backends"].get("schrodinger", {}).get("home")
        schrodinger_probe = inspect_schrodinger(explicit)
        if schrodinger_probe.status != "OK" or schrodinger_probe.location is None:
            raise RunExecutionError(
                f"Schrödinger backend is not ready: {schrodinger_probe.summary}"
            )
    pymol_probe: BackendProbe | None = None
    if pymol_requested:
        explicit = resolved_config["backends"].get("pymol", {}).get("launcher")
        pymol_probe = inspect_pymol(explicit)
        if pymol_probe.status != "OK" or pymol_probe.location is None:
            raise RunExecutionError(f"PyMOL backend is not ready: {pymol_probe.summary}")

    metric_columns = metric_result_columns(config)
    fieldnames = [*RESULT_FIELDS, *metric_columns]
    rows = [_initial_row(variant, metric_columns) for variant in report.variants]
    run_document = _run_document(
        workflow,
        resolved_config,
        report.variants,
        report.input_path,
        command or f"mutflow run {workflow}",
        schrodinger_probe,
        pymol_probe,
    )
    command_text = command or (
        f"mutflow run --resume {workflow}" if resume else f"mutflow run {workflow}"
    )
    if resume:
        try:
            store = RunStore.open_existing(report.output_path)
            pending_variants = _validate_resume(
                store,
                run_document,
                fieldnames,
                report.variants,
                config["output"]["structure_formats"],
            )
        except RunExecutionError:
            raise
        except Exception as exc:
            raise RunExecutionError(f"cannot open resume output: {exc}") from exc
        store.begin_resume(command_text)
    else:
        try:
            store = RunStore.create(
                report.output_path,
                rows,
                run_document,
                fieldnames=fieldnames,
            )
        except Exception as exc:
            raise RunExecutionError(f"cannot initialize output: {exc}") from exc
        pending_variants = list(report.variants)

    row_by_id = {str(row["variant_id"]): row for row in store.rows}
    model_variants = [
        variant
        for variant in pending_variants
        if row_by_id[variant.variant_id]["model_status"] == "NOT_RUN"
    ]
    durations: dict[str, float] = {
        str(row["variant_id"]): float(row["duration_seconds"] or 0.0)
        for row in store.rows
    }
    metric_inputs: list[dict[str, str]] = [
        {
            "variant_id": variant.variant_id,
            "path": str(
                (store.output / str(row_by_id[variant.variant_id]["structure_path"])).resolve()
            ),
        }
        for variant in pending_variants
        if row_by_id[variant.variant_id]["model_status"] in {"OK", "NOT_REQUESTED"}
        and row_by_id[variant.variant_id]["metrics_status"] == "NOT_RUN"
    ]
    checks = config["quality_control"]["required_selections"]
    if modeling_enabled and model_variants:
        assert schrodinger_probe is not None and schrodinger_probe.location is not None
        by_id = {variant.variant_id: variant for variant in report.variants}
        completed: set[str] = set()
        with tempfile.TemporaryDirectory(prefix="mutflow_model_") as temporary:
            staging = Path(temporary).resolve()
            try:
                results = iter_schrodinger_models(
                    schrodinger_probe.location,
                    config,
                    model_variants,
                    report.input_path,
                    staging,
                )
                for result in results:
                    variant_id = str(result["variant_id"])
                    variant = by_id[variant_id]
                    durations[variant_id] = float(result.get("duration_seconds", 0.0))
                    if result.get("status") == "FAILED":
                        store.update_row(
                            variant_id,
                            model_status="FAILED",
                            mutation_readback_status=result.get(
                                "mutation_readback_status", "FAILED"
                            ),
                            metrics_status="SKIPPED",
                            qc_status="SKIPPED",
                            overall_status="FAILED",
                            duration_seconds=round(durations[variant_id], 6),
                            issue_codes=";".join(result.get("issue_codes", [])),
                            notes=str(result.get("notes", "modeling failed"))
                            .replace("\r", " ")
                            .replace("\n", " "),
                        )
                        store.append_log(
                            "ERROR",
                            "variant_failed",
                            f"variant_id={variant_id} stage=schrodinger",
                        )
                        completed.add(variant_id)
                        continue

                    staged_outputs = result.get("outputs", {})
                    expected_formats = config["output"]["structure_formats"]
                    promoted: list[Path] = []
                    try:
                        sources: dict[str, Path] = {}
                        for output_format in expected_formats:
                            source = Path(staged_outputs[output_format]).resolve()
                            expected = (staging / f"{variant_id}.{output_format}").resolve()
                            if source != expected or not source.is_file():
                                raise ValueError(
                                    f"unexpected staged {output_format} output for {variant_id}"
                                )
                            sources[output_format] = source
                        for output_format in expected_formats:
                            target = store.structures / f"{variant_id}.{output_format}"
                            atomic_copy(sources[output_format], target)
                            promoted.append(target)
                    except Exception:
                        for path in promoted:
                            path.unlink(missing_ok=True)
                        raise

                    primary_format = "pdb" if "pdb" in expected_formats else expected_formats[0]
                    primary = store.structures / f"{variant_id}.{primary_format}"
                    readback_status = str(result["mutation_readback_status"])
                    issue_codes = ";".join(result.get("issue_codes", []))
                    details = result.get("details", {})
                    detail_note = ""
                    if details:
                        detail_note = (
                            f"selected_residues={details.get('selected_residue_count')}, "
                            f"selected_atoms={details.get('selected_atom_count')}, "
                            f"selected_fraction={details.get('selected_fraction')}, "
                            "fixed_heavy_max_displacement_angstrom="
                            f"{details.get('fixed_heavy_atom_max_displacement_angstrom')}"
                        )
                    notes = "; ".join(
                        part for part in (str(result.get("notes", "")), detail_note) if part
                    )
                    overall = "CHECK" if readback_status == "CHECK" else "OK"
                    waiting_for_pymol = pymol_requested
                    store.update_row(
                        variant_id,
                        structure_path=primary.relative_to(store.output).as_posix(),
                        model_status=("NOT_REQUESTED" if variant.mutation_count == 0 else "OK"),
                        mutation_readback_status=readback_status,
                        metrics_status="NOT_RUN" if waiting_for_pymol else "NOT_REQUESTED",
                        qc_status="NOT_RUN" if checks else "OK",
                        overall_status="NOT_RUN" if waiting_for_pymol else overall,
                        duration_seconds=round(durations[variant_id], 6),
                        issue_codes=issue_codes,
                        notes=notes,
                    )
                    if waiting_for_pymol:
                        pdb_path = store.structures / f"{variant_id}.pdb"
                        metric_inputs.append(
                            {"variant_id": variant_id, "path": str(pdb_path.resolve())}
                        )
                    else:
                        store.append_log(
                            "INFO", "variant_completed", f"variant_id={variant_id}"
                        )
                    completed.add(variant_id)
            except KeyboardInterrupt:
                store.mark_interrupted()
                raise
            except (ModelingExecutionError, KeyError, OSError, ValueError) as exc:
                for variant in model_variants:
                    if variant.variant_id in completed:
                        continue
                    store.update_row(
                        variant.variant_id,
                        model_status="FAILED",
                        mutation_readback_status="SKIPPED",
                        metrics_status="SKIPPED",
                        qc_status="SKIPPED",
                        overall_status="FAILED",
                        duration_seconds=None,
                        issue_codes="SCHRODINGER_BATCH_FAILED",
                        notes=str(exc).replace("\r", " ").replace("\n", " "),
                    )
                    store.append_log(
                        "ERROR",
                        "variant_failed",
                        f"variant_id={variant.variant_id} stage=schrodinger_batch",
                    )
    elif not modeling_enabled:
        for variant in model_variants:
            started = time.perf_counter()
            try:
                source = _source_for_variant(variant, report.input_path)
                suffix = source.suffix.lower() or ".pdb"
                target = store.structures / f"{variant.variant_id}{suffix}"
                atomic_copy(source, target)
                durations[variant.variant_id] = time.perf_counter() - started
                waiting_for_pymol = pymol_requested
                store.update_row(
                    variant.variant_id,
                    structure_path=target.relative_to(store.output).as_posix(),
                    model_status="NOT_REQUESTED",
                    mutation_readback_status=(
                        "NOT_REQUESTED" if variant.mutation_count == 0 else "OK"
                    ),
                    metrics_status="NOT_RUN" if waiting_for_pymol else "NOT_REQUESTED",
                    qc_status="NOT_RUN" if checks else "OK",
                    overall_status="NOT_RUN" if waiting_for_pymol else "OK",
                    duration_seconds=round(durations[variant.variant_id], 6),
                    issue_codes="",
                    notes="existing structure validated and copied",
                )
                if waiting_for_pymol:
                    metric_inputs.append(
                        {"variant_id": variant.variant_id, "path": str(target.resolve())}
                    )
                else:
                    store.append_log(
                        "INFO", "variant_completed", f"variant_id={variant.variant_id}"
                    )
            except KeyboardInterrupt:
                store.mark_interrupted()
                raise
            except Exception as exc:
                durations[variant.variant_id] = time.perf_counter() - started
                store.update_row(
                    variant.variant_id,
                    model_status="NOT_REQUESTED",
                    mutation_readback_status=(
                        "NOT_REQUESTED" if variant.mutation_count == 0 else "FAILED"
                    ),
                    metrics_status="SKIPPED",
                    qc_status="FAILED",
                    overall_status="FAILED",
                    duration_seconds=round(durations[variant.variant_id], 6),
                    issue_codes="EXISTING_STRUCTURE_COPY_FAILED",
                    notes=str(exc).replace("\r", " ").replace("\n", " "),
                )
                store.append_log(
                    "ERROR", "variant_failed", f"variant_id={variant.variant_id} error={exc}"
                )

    if pymol_requested and metric_inputs:
        assert pymol_probe is not None and pymol_probe.location is not None
        metric_run_inputs = list(metric_inputs)
        wants_wt_delta = any(
            metric.get("calculate_wt_delta", False)
            for section in ("sasa", "minimum_distance")
            for metric in config["metrics"][section]
        )
        if wants_wt_delta and not any(
            item["variant_id"] == "WT" for item in metric_run_inputs
        ):
            wt_row = row_by_id.get("WT")
            if wt_row is not None and wt_row.get("structure_path"):
                metric_run_inputs.insert(
                    0,
                    {
                        "variant_id": "WT",
                        "path": str(
                            (store.output / str(wt_row["structure_path"])).resolve()
                        ),
                    },
                )
        batch_started = time.perf_counter()
        try:
            metric_results = run_pymol_metrics(
                pymol_probe.location,
                config,
                metric_run_inputs,
            )
        except KeyboardInterrupt:
            store.mark_interrupted()
            raise
        except MetricExecutionError as exc:
            metric_results = {}
            for item in metric_inputs:
                variant_id = item["variant_id"]
                current = next(row for row in store.rows if row["variant_id"] == variant_id)
                store.update_row(
                    variant_id,
                    metrics_status="FAILED",
                    qc_status="SKIPPED",
                    overall_status="FAILED",
                    issue_codes=";".join(
                        part for part in (current["issue_codes"], "PYMOL_BATCH_FAILED") if part
                    ),
                    notes="; ".join(
                        part
                        for part in (
                            current["notes"],
                            str(exc).replace("\r", " ").replace("\n", " "),
                        )
                        if part
                    ),
                )
                store.append_log(
                    "ERROR", "variant_failed", f"variant_id={variant_id} error={exc}"
                )
        else:
            batch_duration = time.perf_counter() - batch_started
            share = batch_duration / len(metric_inputs)
            wt_result = metric_results.get("WT")
            for item in metric_inputs:
                variant_id = item["variant_id"]
                result = metric_results[variant_id]
                current = next(row for row in store.rows if row["variant_id"] == variant_id)
                durations[variant_id] = durations.get(variant_id, 0.0) + share
                if result.get("status") != "OK":
                    store.update_row(
                        variant_id,
                        metrics_status="FAILED",
                        qc_status="SKIPPED",
                        overall_status="FAILED",
                        duration_seconds=round(durations[variant_id], 6),
                        issue_codes=";".join(
                            part
                            for part in (current["issue_codes"], "PYMOL_VARIANT_FAILED")
                            if part
                        ),
                        notes="; ".join(
                            part
                            for part in (
                                current["notes"],
                                str(result.get("error", "PyMOL metric failure")),
                            )
                            if part
                        ),
                    )
                    store.append_log(
                        "ERROR",
                        "variant_failed",
                        f"variant_id={variant_id} error={result.get('error', '')}",
                    )
                    continue

                updates: dict[str, Any] = {}
                metric_error = ""
                for metric in config["metrics"]["sasa"]:
                    name = metric["name"]
                    slug = normalized_metric_name(name)
                    value = float(result["metrics"][name])
                    updates[f"metric__{slug}__angstrom2"] = round(value, 6)
                    if metric["calculate_wt_delta"]:
                        if wt_result is None or wt_result.get("status") != "OK":
                            metric_error = "WT baseline is unavailable for requested SASA delta"
                            break
                        wt_value = float(wt_result["metrics"][name])
                        updates[f"metric__{slug}__delta_angstrom2"] = round(
                            value - wt_value, 6
                        )
                for metric in config["metrics"]["minimum_distance"]:
                    name = metric["name"]
                    slug = normalized_metric_name(name)
                    value = float(result["metrics"][name])
                    updates[f"metric__{slug}__angstrom"] = round(value, 6)
                    if metric["calculate_wt_delta"]:
                        if wt_result is None or wt_result.get("status") != "OK":
                            metric_error = "WT baseline is unavailable for requested distance delta"
                            break
                        wt_value = float(wt_result["metrics"][name])
                        updates[f"metric__{slug}__delta_angstrom"] = round(
                            value - wt_value, 6
                        )
                for metric in config["metrics"]["clashes"]:
                    name = metric["name"]
                    slug = normalized_metric_name(name)
                    values = result["metrics"][name]
                    updates[f"metric__{slug}__clash_count"] = int(values["clash_count"])
                    updates[f"metric__{slug}__severe_clash_count"] = int(
                        values["severe_clash_count"]
                    )
                    max_overlap = values["max_overlap_angstrom"]
                    minimum_contact = values["min_contact_distance_angstrom"]
                    updates[f"metric__{slug}__max_overlap_angstrom"] = (
                        None if max_overlap is None else round(float(max_overlap), 6)
                    )
                    updates[f"metric__{slug}__min_contact_distance_angstrom"] = round(
                        float(minimum_contact), 6
                    )
                    updates[f"metric__{slug}__worst_pair"] = values["worst_pair"]
                if metric_error:
                    store.update_row(
                        variant_id,
                        metrics_status="FAILED",
                        qc_status="SKIPPED",
                        overall_status="FAILED",
                        duration_seconds=round(durations[variant_id], 6),
                        issue_codes=";".join(
                            part
                            for part in (current["issue_codes"], "WT_BASELINE_UNAVAILABLE")
                            if part
                        ),
                        notes="; ".join(
                            part for part in (current["notes"], metric_error) if part
                        ),
                    )
                    continue

                qc_status, qc_codes, qc_notes = _selection_qc(
                    result.get("selection_counts", {}), checks
                )
                prior_check = current["model_status"] == "CHECK" or current[
                    "mutation_readback_status"
                ] == "CHECK"
                overall = "FAILED" if qc_status == "FAILED" else (
                    "CHECK" if qc_status == "CHECK" or prior_check else "OK"
                )
                updates.update(
                    {
                        "metrics_status": "OK" if report.metric_names else "NOT_REQUESTED",
                        "qc_status": qc_status,
                        "overall_status": overall,
                        "duration_seconds": round(durations[variant_id], 6),
                        "issue_codes": ";".join(
                            part
                            for part in (current["issue_codes"], *qc_codes)
                            if part
                        ),
                        "notes": "; ".join(
                            part for part in (current["notes"], *qc_notes) if part
                        ),
                    }
                )
                store.update_row(variant_id, **updates)
                store.append_log(
                    "INFO", "variant_completed", f"variant_id={variant_id}"
                )
            store.append_log(
                "INFO",
                "pymol_batch_completed",
                (
                    f"variant_count={len(metric_inputs)} "
                    f"reference_count={len(metric_run_inputs) - len(metric_inputs)}"
                ),
            )

    store.finalize()
    return store
