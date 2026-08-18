"""Minimal stdin/stdout worker executed by the external PyMOL Python runtime."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any


MARKER = "MUTFLOW_RESULT_JSON="


def _heavy_atoms(cmd: Any, selection: str) -> list[Any]:
    atoms = list(cmd.get_model(selection).atom)
    return [atom for atom in atoms if str(atom.symbol).strip().upper() != "H"]


def _atom_label(atom: Any) -> str:
    return "/".join(
        value
        for value in (
            str(atom.chain).strip(),
            str(atom.resn).strip(),
            str(atom.resi).strip(),
            str(atom.name).strip(),
        )
        if value
    )


def _element(atom: Any) -> str:
    symbol = str(atom.symbol).strip()
    if symbol:
        return symbol[0].upper() + symbol[1:].lower()
    name = str(atom.name).strip().lstrip("0123456789")
    return name[:1].upper()


def _minimum_distance(cmd: Any, metric: dict[str, Any]) -> float:
    atoms_a = _heavy_atoms(cmd, metric["selection_a"])
    atoms_b = _heavy_atoms(cmd, metric["selection_b"])
    if not atoms_a or not atoms_b:
        raise ValueError(f"distance metric {metric['name']!r} has an empty selection")
    distances = [
        math.dist(atom_a.coord, atom_b.coord)
        for atom_a in atoms_a
        for atom_b in atoms_b
        if atom_a.index != atom_b.index
    ]
    if not distances:
        raise ValueError(f"distance metric {metric['name']!r} has no distinct atom pair")
    value = min(distances)
    if not math.isfinite(value):
        raise ValueError(f"distance metric {metric['name']!r} is not finite")
    return value


def _clash_metrics(cmd: Any, metric: dict[str, Any]) -> dict[str, Any]:
    atoms_a = _heavy_atoms(cmd, metric["selection_a"])
    atoms_b = _heavy_atoms(cmd, metric["selection_b"])
    if not atoms_a or not atoms_b:
        raise ValueError(f"clash metric {metric['name']!r} has an empty selection")
    radii = {str(key): float(value) for key, value in metric["vdw_radii"].items()}
    clash_threshold = float(metric["clash_overlap_angstrom"])
    severe_threshold = float(metric["severe_overlap_angstrom"])
    seen: set[tuple[int, int]] = set()
    clash_count = 0
    severe_count = 0
    maximum_overlap: float | None = None
    minimum_contact: float | None = None
    worst_pair: str | None = None
    for atom_a in atoms_a:
        for atom_b in atoms_b:
            if atom_a.index == atom_b.index:
                continue
            pair_key = tuple(sorted((int(atom_a.index), int(atom_b.index))))
            if pair_key in seen:
                continue
            seen.add(pair_key)
            element_a = _element(atom_a)
            element_b = _element(atom_b)
            if element_a not in radii or element_b not in radii:
                raise ValueError(
                    f"clash metric {metric['name']!r} has no VDW radius for "
                    f"{element_a or '?'} or {element_b or '?'}"
                )
            distance = math.dist(atom_a.coord, atom_b.coord)
            minimum_contact = (
                distance if minimum_contact is None else min(minimum_contact, distance)
            )
            overlap = radii[element_a] + radii[element_b] - distance
            if overlap > clash_threshold:
                clash_count += 1
                if maximum_overlap is None or overlap > maximum_overlap:
                    maximum_overlap = overlap
                    worst_pair = f"{_atom_label(atom_a)}--{_atom_label(atom_b)}"
            if overlap > severe_threshold:
                severe_count += 1
    if minimum_contact is None:
        raise ValueError(f"clash metric {metric['name']!r} has no distinct atom pair")
    return {
        "clash_count": clash_count,
        "severe_clash_count": severe_count,
        "max_overlap_angstrom": maximum_overlap,
        "min_contact_distance_angstrom": minimum_contact,
        "worst_pair": worst_pair,
    }


def _measure_structure(cmd: Any, item: dict[str, Any], request: dict[str, Any]) -> dict[str, Any]:
    variant_id = item["variant_id"]
    path = Path(item["path"])
    output: dict[str, Any] = {
        "variant_id": variant_id,
        "status": "OK",
        "metrics": {},
        "selection_counts": {},
        "error": "",
    }
    try:
        cmd.delete("all")
        cmd.load(str(path), "structure")
        cmd.do("flag ignore, all, clear")
        for metric in request.get("sasa", []):
            cmd.set("dot_solvent", 1 if metric["dot_solvent"] else 0)
            cmd.set("solvent_radius", float(metric["solvent_radius_angstrom"]))
            cmd.set("dot_density", int(metric["dot_density"]))
            cmd.rebuild()
            selection = metric["selection"]
            count = int(cmd.count_atoms(selection))
            if count < 1:
                raise ValueError(
                    f"SASA selection {metric['name']!r} matched zero atoms"
                )
            value = float(cmd.get_area(selection))
            if not math.isfinite(value):
                raise ValueError(f"SASA metric {metric['name']!r} is not finite")
            output["metrics"][metric["name"]] = value
        for metric in request.get("minimum_distance", []):
            output["metrics"][metric["name"]] = _minimum_distance(cmd, metric)
        for metric in request.get("clashes", []):
            output["metrics"][metric["name"]] = _clash_metrics(cmd, metric)
        for check in request.get("required_selections", []):
            output["selection_counts"][check["name"]] = int(
                cmd.count_atoms(check["selection"])
            )
    except Exception as exc:
        output["status"] = "FAILED"
        output["error"] = str(exc).replace("\r", " ").replace("\n", " ")
    finally:
        cmd.delete("all")
    return output


def main() -> int:
    request = json.load(sys.stdin)
    from pymol import cmd

    results = [
        _measure_structure(cmd, item, request)
        for item in request.get("structures", [])
    ]
    print(MARKER + json.dumps({"results": results}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
