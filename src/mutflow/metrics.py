from __future__ import annotations

import json
import os
import re
import subprocess
from pathlib import Path
from typing import Any

from mutflow.backends import BackendLocation


_RESULT_MARKER = "MUTFLOW_RESULT_JSON="


class MetricExecutionError(RuntimeError):
    """Raised when an external metric worker cannot return usable results."""


def normalized_metric_name(name: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    if not normalized:
        raise ValueError(f"metric name cannot form a result column: {name!r}")
    return normalized


def metric_result_columns(config: dict[str, Any]) -> list[str]:
    columns: list[str] = []
    for metric in config["metrics"]["sasa"]:
        slug = normalized_metric_name(metric["name"])
        columns.append(f"metric__{slug}__angstrom2")
        if metric["calculate_wt_delta"]:
            columns.append(f"metric__{slug}__delta_angstrom2")
    for metric in config["metrics"]["minimum_distance"]:
        slug = normalized_metric_name(metric["name"])
        columns.append(f"metric__{slug}__angstrom")
        if metric["calculate_wt_delta"]:
            columns.append(f"metric__{slug}__delta_angstrom")
    for metric in config["metrics"]["clashes"]:
        slug = normalized_metric_name(metric["name"])
        columns.extend(
            [
                f"metric__{slug}__clash_count",
                f"metric__{slug}__severe_clash_count",
                f"metric__{slug}__max_overlap_angstrom",
                f"metric__{slug}__min_contact_distance_angstrom",
                f"metric__{slug}__worst_pair",
            ]
        )
    return columns


def metric_definitions(config: dict[str, Any]) -> list[dict[str, Any]]:
    definitions: list[dict[str, Any]] = []
    for metric in config["metrics"]["sasa"]:
        slug = normalized_metric_name(metric["name"])
        result_columns = [
            {"name": f"metric__{slug}__angstrom2", "unit": "angstrom^2"}
        ]
        if metric["calculate_wt_delta"]:
            result_columns.append(
                {"name": f"metric__{slug}__delta_angstrom2", "unit": "angstrom^2"}
            )
        definitions.append(
            {
                "name": metric["name"],
                "type": "sasa",
                "selections": {"selection": metric["selection"]},
                "parameters": {
                    "dot_solvent": metric["dot_solvent"],
                    "solvent_radius_angstrom": metric["solvent_radius_angstrom"],
                    "dot_density": metric["dot_density"],
                    "calculate_wt_delta": metric["calculate_wt_delta"],
                },
                "result_columns": result_columns,
            }
        )
    for metric in config["metrics"]["minimum_distance"]:
        slug = normalized_metric_name(metric["name"])
        result_columns = [
            {"name": f"metric__{slug}__angstrom", "unit": "angstrom"}
        ]
        if metric["calculate_wt_delta"]:
            result_columns.append(
                {"name": f"metric__{slug}__delta_angstrom", "unit": "angstrom"}
            )
        definitions.append(
            {
                "name": metric["name"],
                "type": "minimum_distance",
                "selections": {
                    "selection_a": metric["selection_a"],
                    "selection_b": metric["selection_b"],
                },
                "parameters": {
                    "aggregation": "minimum",
                    "heavy_atoms_only": True,
                    "calculate_wt_delta": metric["calculate_wt_delta"],
                },
                "result_columns": result_columns,
            }
        )
    for metric in config["metrics"]["clashes"]:
        slug = normalized_metric_name(metric["name"])
        definitions.append(
            {
                "name": metric["name"],
                "type": "clash",
                "selections": {
                    "selection_a": metric["selection_a"],
                    "selection_b": metric["selection_b"],
                },
                "parameters": {
                    "clash_overlap_angstrom": metric["clash_overlap_angstrom"],
                    "severe_overlap_angstrom": metric["severe_overlap_angstrom"],
                    "heavy_atoms_only": True,
                    "vdw_radii": metric["vdw_radii"],
                    "exclude_bonded_pairs": metric["exclude_bonded_pairs"],
                },
                "result_columns": [
                    {"name": f"metric__{slug}__clash_count", "unit": "count"},
                    {"name": f"metric__{slug}__severe_clash_count", "unit": "count"},
                    {"name": f"metric__{slug}__max_overlap_angstrom", "unit": "angstrom"},
                    {
                        "name": f"metric__{slug}__min_contact_distance_angstrom",
                        "unit": "angstrom",
                    },
                    {"name": f"metric__{slug}__worst_pair", "unit": None},
                ],
            }
        )
    return definitions


def run_pymol_metrics(
    location: BackendLocation,
    config: dict[str, Any],
    structures: list[dict[str, str]],
    *,
    timeout_seconds: float | None = None,
) -> dict[str, dict[str, Any]]:
    worker = Path(__file__).with_name("pymol_worker.py").resolve()
    request = {
        "structures": structures,
        "sasa": config["metrics"]["sasa"],
        "minimum_distance": config["metrics"]["minimum_distance"],
        "clashes": config["metrics"]["clashes"],
        "required_selections": config["quality_control"]["required_selections"],
    }
    timeout = timeout_seconds or max(30.0, 10.0 * len(structures))
    creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    try:
        completed = subprocess.run(
            [str(location.launcher), "-B", str(worker)],
            input=json.dumps(request),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            check=False,
            creationflags=creationflags,
        )
    except subprocess.TimeoutExpired as exc:
        raise MetricExecutionError(f"PyMOL metric batch timed out after {timeout:g} seconds") from exc
    except OSError as exc:
        raise MetricExecutionError(f"PyMOL metric worker could not start: {exc}") from exc
    combined = "\n".join(part for part in (completed.stdout, completed.stderr) if part)
    if completed.returncode != 0:
        tail = combined.strip().splitlines()[-1] if combined.strip() else "no diagnostic output"
        raise MetricExecutionError(
            f"PyMOL metric worker exited {completed.returncode}: {tail}"
        )
    payload: dict[str, Any] | None = None
    for line in completed.stdout.splitlines():
        if line.startswith(_RESULT_MARKER):
            try:
                payload = json.loads(line[len(_RESULT_MARKER):])
            except json.JSONDecodeError as exc:
                raise MetricExecutionError(f"PyMOL worker returned invalid JSON: {exc}") from exc
            break
    if payload is None:
        raise MetricExecutionError("PyMOL worker returned no result marker")
    results = payload.get("results")
    if not isinstance(results, list):
        raise MetricExecutionError("PyMOL worker result payload has no results list")
    by_id: dict[str, dict[str, Any]] = {}
    for result in results:
        variant_id = str(result.get("variant_id", ""))
        if not variant_id or variant_id in by_id:
            raise MetricExecutionError("PyMOL worker returned missing or duplicate variant_id")
        by_id[variant_id] = result
    requested = {item["variant_id"] for item in structures}
    if set(by_id) != requested:
        raise MetricExecutionError("PyMOL worker result identifiers do not match request")
    return by_id
