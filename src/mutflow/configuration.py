from __future__ import annotations

import copy
import json
from importlib import resources
from pathlib import Path
from typing import Any

import yaml
from jsonschema import Draft202012Validator


class ConfigurationError(ValueError):
    """Raised when workflow configuration cannot be loaded or validated."""


def _workflow_schema() -> dict[str, Any]:
    schema_ref = resources.files("mutflow.schemas").joinpath("workflow.schema.json")
    return json.loads(schema_ref.read_text(encoding="utf-8"))


def load_workflow(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"workflow file not found: {resolved}")
    try:
        loaded = yaml.safe_load(resolved.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise ConfigurationError(f"cannot read workflow YAML: {exc}") from exc
    if not isinstance(loaded, dict):
        raise ConfigurationError("workflow root must be a mapping")

    return validate_workflow_data(loaded)


def validate_workflow_data(raw: dict[str, Any]) -> dict[str, Any]:
    """Validate an in-memory workflow and return its normalized copy."""

    validator = Draft202012Validator(_workflow_schema())
    errors = sorted(validator.iter_errors(raw), key=lambda item: list(item.absolute_path))
    if errors:
        messages: list[str] = []
        for error in errors[:10]:
            location = ".".join(str(part) for part in error.absolute_path) or "<root>"
            messages.append(f"{location}: {error.message}")
        if len(errors) > 10:
            messages.append(f"... and {len(errors) - 10} more schema error(s)")
        raise ConfigurationError("; ".join(messages))
    return normalize_workflow(raw)


def normalize_workflow(raw: dict[str, Any]) -> dict[str, Any]:
    config = copy.deepcopy(raw)
    config["input"].setdefault("format", "auto")
    config["variants"].setdefault("include_wt", True)
    config["modeling"].setdefault("strategy", "sequential_from_wt")
    config["metrics"] = config.get("metrics") or {}
    config["metrics"].setdefault("sasa", [])
    config["metrics"].setdefault("minimum_distance", [])
    config["metrics"].setdefault("clashes", [])
    for metric in config["metrics"]["sasa"]:
        metric.setdefault("calculate_wt_delta", False)
        metric.setdefault("dot_solvent", True)
        metric.setdefault("solvent_radius_angstrom", 1.4)
        metric.setdefault("dot_density", 3)
    for metric in config["metrics"]["minimum_distance"]:
        metric.setdefault("calculate_wt_delta", False)
    for metric in config["metrics"]["clashes"]:
        metric.setdefault("clash_overlap_angstrom", 0.4)
        metric.setdefault("severe_overlap_angstrom", 0.8)
        metric.setdefault(
            "vdw_radii",
            {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80},
        )
        metric.setdefault("exclude_bonded_pairs", False)
    config.setdefault("quality_control", {})
    config["quality_control"].setdefault("required_selections", [])
    config["output"].setdefault("structure_formats", ["pdb"])
    config["output"].setdefault("keep_intermediates", False)
    config.setdefault("execution", {})
    config["execution"].setdefault("max_variants", 500)
    config.setdefault("backends", {})
    config.setdefault("extensions", {})
    return config


def resolve_workflow_paths(config: dict[str, Any], workflow_path: Path) -> dict[str, Any]:
    resolved = copy.deepcopy(config)
    base = workflow_path.expanduser().resolve().parent

    def resolve_value(value: str) -> str:
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = base / path
        return str(path.resolve())

    resolved["input"]["structure"] = resolve_value(resolved["input"]["structure"])
    if resolved["variants"]["mode"] == "explicit":
        resolved["variants"]["file"] = resolve_value(resolved["variants"]["file"])
    resolved["output"]["directory"] = resolve_value(resolved["output"]["directory"])
    if "schrodinger" in resolved["backends"]:
        resolved["backends"]["schrodinger"]["home"] = resolve_value(
            resolved["backends"]["schrodinger"]["home"]
        )
    if "pymol" in resolved["backends"]:
        resolved["backends"]["pymol"]["launcher"] = resolve_value(
            resolved["backends"]["pymol"]["launcher"]
        )
    return resolved
