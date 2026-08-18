from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Callable

import yaml

from mutflow.configuration import ConfigurationError, validate_workflow_data


class InitError(ValueError):
    """Raised when deterministic workflow initialization cannot complete."""


InputFn = Callable[[str], str]
OutputFn = Callable[[str], None]


def _ask(
    prompt: str,
    *,
    input_fn: InputFn,
    output_fn: OutputFn,
    default: str | None = None,
) -> str:
    suffix = f" [{default}]" if default is not None else ""
    while True:
        value = input_fn(f"{prompt}{suffix}: ").strip()
        if value:
            return value
        if default is not None:
            return default
        output_fn("A value is required.")


def _choice(
    prompt: str,
    choices: tuple[str, ...],
    *,
    input_fn: InputFn,
    output_fn: OutputFn,
    default: str,
) -> str:
    allowed = {choice.lower(): choice for choice in choices}
    while True:
        value = _ask(
            f"{prompt} ({'/'.join(choices)})",
            input_fn=input_fn,
            output_fn=output_fn,
            default=default,
        ).lower()
        if value in allowed:
            return allowed[value]
        output_fn(f"Choose one of: {', '.join(choices)}.")


def _yes_no(
    prompt: str,
    *,
    input_fn: InputFn,
    output_fn: OutputFn,
    default: bool,
) -> bool:
    default_text = "yes" if default else "no"
    while True:
        value = _ask(
            f"{prompt} (yes/no)",
            input_fn=input_fn,
            output_fn=output_fn,
            default=default_text,
        ).lower()
        if value in {"yes", "y"}:
            return True
        if value in {"no", "n"}:
            return False
        output_fn("Enter yes or no.")


def _number(
    prompt: str,
    *,
    input_fn: InputFn,
    output_fn: OutputFn,
    default: float,
    minimum: float = 0.0,
) -> float:
    while True:
        text = _ask(
            prompt,
            input_fn=input_fn,
            output_fn=output_fn,
            default=str(default),
        )
        try:
            value = float(text)
        except ValueError:
            output_fn("Enter a number.")
            continue
        if value <= minimum:
            output_fn(f"Enter a number greater than {minimum}.")
            continue
        return value


def _nonnegative_integer(
    prompt: str,
    *,
    input_fn: InputFn,
    output_fn: OutputFn,
) -> int:
    while True:
        text = _ask(prompt, input_fn=input_fn, output_fn=output_fn)
        try:
            value = int(text)
        except ValueError:
            output_fn("Enter a whole number.")
            continue
        if value < 0:
            output_fn("Enter zero or a positive whole number.")
            continue
        return value


_SITE_PATTERN = re.compile(
    r"^(?P<chain>[^:\s]+):(?P<number>-?\d+)(?P<icode>[A-Za-z]?):"
    r"(?P<wt>[ACDEFGHIKLMNPQRSTVWY])$"
)


def _parse_sites(text: str) -> list[dict[str, object]]:
    sites: list[dict[str, object]] = []
    seen: set[tuple[str, int, str]] = set()
    for token in (part.strip() for part in text.split(",")):
        match = _SITE_PATTERN.fullmatch(token)
        if match is None:
            raise InitError(
                f"invalid saturation site {token!r}; use CHAIN:NUMBER:WT, for example A:10:G"
            )
        site = {
            "chain": match.group("chain"),
            "residue_number": int(match.group("number")),
            "insertion_code": match.group("icode").upper(),
            "expected_wt": match.group("wt"),
        }
        key = (str(site["chain"]), int(site["residue_number"]), str(site["insertion_code"]))
        if key in seen:
            raise InitError(f"duplicate saturation site: {token}")
        seen.add(key)
        sites.append(site)
    if not sites:
        raise InitError("at least one saturation site is required")
    return sites


def _collect_metrics(input_fn: InputFn, output_fn: OutputFn) -> dict[str, list[dict]]:
    metrics: dict[str, list[dict]] = {
        "sasa": [],
        "minimum_distance": [],
        "clashes": [],
    }
    used_names: set[str] = set()
    while True:
        kind = _choice(
            "Add a PyMOL metric",
            ("none", "sasa", "distance", "clash"),
            input_fn=input_fn,
            output_fn=output_fn,
            default="none",
        )
        if kind == "none":
            return metrics
        name = _ask(
            "Metric name",
            input_fn=input_fn,
            output_fn=output_fn,
        )
        if name.lower() in used_names:
            output_fn("Metric names must be unique.")
            continue
        used_names.add(name.lower())
        if kind == "sasa":
            metrics["sasa"].append(
                {
                    "name": name,
                    "selection": _ask(
                        "PyMOL selection",
                        input_fn=input_fn,
                        output_fn=output_fn,
                    ),
                    "calculate_wt_delta": _yes_no(
                        "Calculate mutant-minus-WT delta",
                        input_fn=input_fn,
                        output_fn=output_fn,
                        default=False,
                    ),
                }
            )
        elif kind == "distance":
            metrics["minimum_distance"].append(
                {
                    "name": name,
                    "selection_a": _ask(
                        "PyMOL selection A",
                        input_fn=input_fn,
                        output_fn=output_fn,
                    ),
                    "selection_b": _ask(
                        "PyMOL selection B",
                        input_fn=input_fn,
                        output_fn=output_fn,
                    ),
                    "calculate_wt_delta": _yes_no(
                        "Calculate mutant-minus-WT delta",
                        input_fn=input_fn,
                        output_fn=output_fn,
                        default=False,
                    ),
                }
            )
        else:
            clash = _number(
                "Clash overlap threshold in angstrom",
                input_fn=input_fn,
                output_fn=output_fn,
                default=0.4,
                minimum=-1e-12,
            )
            severe = _number(
                "Severe clash overlap threshold in angstrom",
                input_fn=input_fn,
                output_fn=output_fn,
                default=0.8,
                minimum=-1e-12,
            )
            if severe < clash:
                raise InitError("severe clash threshold must be at least the clash threshold")
            metrics["clashes"].append(
                {
                    "name": name,
                    "selection_a": _ask(
                        "PyMOL selection A",
                        input_fn=input_fn,
                        output_fn=output_fn,
                    ),
                    "selection_b": _ask(
                        "PyMOL selection B",
                        input_fn=input_fn,
                        output_fn=output_fn,
                    ),
                    "clash_overlap_angstrom": clash,
                    "severe_overlap_angstrom": severe,
                }
            )


def _collect_qc(input_fn: InputFn, output_fn: OutputFn) -> list[dict]:
    checks: list[dict] = []
    used_names: set[str] = set()
    while _yes_no(
        "Add a required-selection QC check",
        input_fn=input_fn,
        output_fn=output_fn,
        default=False,
    ):
        mode = _choice(
            "Count rule",
            ("expected", "minimum", "maximum"),
            input_fn=input_fn,
            output_fn=output_fn,
            default="minimum",
        )
        name = _ask("QC name", input_fn=input_fn, output_fn=output_fn)
        if name.lower() in used_names:
            output_fn("QC names must be unique.")
            continue
        used_names.add(name.lower())
        checks.append(
            {
                "name": name,
                "selection": _ask(
                    "PyMOL selection",
                    input_fn=input_fn,
                    output_fn=output_fn,
                ),
                f"{mode}_count": _nonnegative_integer(
                    "Required atom count",
                    input_fn=input_fn,
                    output_fn=output_fn,
                ),
                "severity": _choice(
                    "Failure handling",
                    ("fail", "check"),
                    input_fn=input_fn,
                    output_fn=output_fn,
                    default="fail",
                ),
            }
        )
    return checks


def _write_new_file(destination: Path, text: str) -> None:
    if not destination.parent.is_dir():
        raise InitError(f"destination directory does not exist: {destination.parent}")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    try:
        descriptor = os.open(destination, flags, 0o666)
    except FileExistsError as exc:
        raise InitError(f"refusing to overwrite existing workflow: {destination}") from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
    except Exception:
        destination.unlink(missing_ok=True)
        raise


def run_init(
    destination: Path,
    *,
    input_fn: InputFn = input,
    output_fn: OutputFn = print,
) -> Path:
    """Interactively create one workflow YAML without invoking scientific backends."""
    destination = destination.expanduser().resolve()
    if destination.exists():
        raise InitError(f"refusing to overwrite existing workflow: {destination}")

    output_fn("MUTFLOW INIT")
    output_fn("This wizard records your choices; it does not choose scientific targets.")
    structure = _ask(
        "WT PDB path",
        input_fn=input_fn,
        output_fn=output_fn,
    )
    structure_format = _choice(
        "Input format",
        ("pdb",),
        input_fn=input_fn,
        output_fn=output_fn,
        default="pdb",
    )
    mode = _choice(
        "Mutation definition mode",
        ("explicit", "saturation"),
        input_fn=input_fn,
        output_fn=output_fn,
        default="explicit",
    )
    if mode == "explicit":
        variants = {
            "mode": "explicit",
            "include_wt": True,
            "file": _ask(
                "Mutation CSV path",
                input_fn=input_fn,
                output_fn=output_fn,
            ),
        }
        modeling_enabled = _yes_no(
            "Generate mutant structures with Schrödinger",
            input_fn=input_fn,
            output_fn=output_fn,
            default=True,
        )
    else:
        site_text = _ask(
            "Saturation sites as CHAIN:NUMBER:WT, comma separated",
            input_fn=input_fn,
            output_fn=output_fn,
        )
        variants = {
            "mode": "saturation_single",
            "include_wt": True,
            "sites": _parse_sites(site_text),
        }
        modeling_enabled = True
        output_fn("Saturation mode requires Schrödinger modeling.")

    if modeling_enabled:
        output_fn("Schrödinger modeling requires a chemically prepared WT PDB.")
        modeling = {
            "enabled": True,
            "backend": "schrodinger",
            "strategy": "sequential_from_wt",
            "local_minimization": {
                "radius_angstrom": _number(
                    "Local-minimization radius in angstrom",
                    input_fn=input_fn,
                    output_fn=output_fn,
                    default=5.0,
                ),
                "include_nearby_nonprotein": _yes_no(
                    "Include nearby nonprotein residues",
                    input_fn=input_fn,
                    output_fn=output_fn,
                    default=True,
                ),
            },
        }
    else:
        modeling = {"enabled": False}

    metrics = _collect_metrics(input_fn, output_fn)
    checks = _collect_qc(input_fn, output_fn)
    output_directory = _ask(
        "Output directory",
        input_fn=input_fn,
        output_fn=output_fn,
        default="./output",
    )

    backends: dict[str, dict[str, str]] = {}
    if modeling_enabled:
        schrodinger_home = _ask(
            "Schrödinger installation root (leave blank for auto-discovery)",
            input_fn=input_fn,
            output_fn=output_fn,
            default="",
        )
        if schrodinger_home:
            backends["schrodinger"] = {"home": schrodinger_home}
    if any(metrics.values()) or checks:
        pymol_launcher = _ask(
            "PyMOL launcher/environment path (leave blank for auto-discovery)",
            input_fn=input_fn,
            output_fn=output_fn,
            default="",
        )
        if pymol_launcher:
            backends["pymol"] = {"launcher": pymol_launcher}

    workflow: dict[str, object] = {
        "schema_version": "1.0",
        "input": {"structure": structure, "format": structure_format},
        "variants": variants,
        "modeling": modeling,
        "metrics": metrics,
        "quality_control": {"required_selections": checks},
        "output": {
            "directory": output_directory,
            "structure_formats": ["pdb"],
            "keep_intermediates": False,
        },
    }
    if backends:
        workflow["backends"] = backends
    try:
        validate_workflow_data(workflow)
    except ConfigurationError as exc:
        raise InitError(f"generated workflow is invalid: {exc}") from exc

    rendered = yaml.safe_dump(
        workflow,
        sort_keys=False,
        allow_unicode=True,
        default_flow_style=False,
    )
    try:
        _write_new_file(destination, rendered)
    except InitError:
        raise
    except OSError as exc:
        raise InitError(f"cannot create workflow {destination}: {exc}") from exc
    output_fn(f"created: {destination}")
    output_fn(f"next: mutflow preflight {destination}")
    return destination
