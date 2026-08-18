from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Any


_RESULT_MARKER = "MUTFLOW_MODEL_JSON="

_AA_1_TO_3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIE",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}
_PROTEIN_RESIDUES = frozenset(_AA_1_TO_3.values()) | {"HIS", "HID", "HIP"}
_HISTIDINE_ALIASES = frozenset({"HIS", "HID", "HIP"})


def _load_modules() -> tuple[Any, Any, Any, Any]:
    from schrodinger import structure
    from schrodinger.application.bioluminate.protein import get_residues_within
    from schrodinger.application.bioluminate.protein.mutator import ProteinMutator
    from schrodinger.forcefield import minimizer

    return structure, get_residues_within, ProteinMutator, minimizer


def _emit(payload: dict[str, Any]) -> None:
    print(_RESULT_MARKER + json.dumps(payload, sort_keys=True), flush=True)


def _residue_name(residue: Any) -> str:
    return str(getattr(residue, "pdbres", "")).strip().upper()


def _normalize_inscode(value: Any) -> str:
    return "" if value is None else str(value).strip().upper()


def _residue_key(residue: Any) -> tuple[str, int, str]:
    return (
        str(residue.chain).strip(),
        int(residue.resnum),
        _normalize_inscode(getattr(residue, "inscode", "")),
    )


def _find_residue(structure: Any, chain: str, residue_number: int, insertion_code: str) -> Any:
    wanted = (chain.strip(), int(residue_number), _normalize_inscode(insertion_code))
    matches = [residue for residue in structure.residue if _residue_key(residue) == wanted]
    if len(matches) != 1:
        raise ValueError(
            f"residue {wanted[0]!r}:{wanted[1]}{wanted[2]} matched {len(matches)} entries"
        )
    return matches[0]


def _read_structure(structure_module: Any, path: Path) -> Any:
    with structure_module.StructureReader(str(path)) as reader:
        return next(reader)


def _classify_readback(target: str, observed: str) -> str:
    if observed == target:
        return "OK"
    if target == "HIE" and observed in _HISTIDINE_ALIASES:
        return "CHECK"
    return "FAILED"


def _verify_mutations(structure: Any, mutations: list[dict[str, Any]]) -> str:
    status = "OK"
    for mutation in mutations:
        residue = _find_residue(
            structure,
            mutation["chain"],
            mutation["residue_number"],
            mutation.get("insertion_code", ""),
        )
        observed = _residue_name(residue)
        target = _AA_1_TO_3[mutation["target"]]
        readback = _classify_readback(target, observed)
        if readback == "FAILED":
            raise ValueError(
                f"mutation readback failed at {mutation['chain']}:{mutation['residue_number']}"
                f"{mutation.get('insertion_code', '')}: expected {target}, observed {observed}"
            )
        if readback == "CHECK":
            status = "CHECK"
    return status


def _atom_xyz(atom: Any) -> tuple[float, float, float]:
    return float(atom.x), float(atom.y), float(atom.z)


def _squared_distance(first: Any, second: Any) -> float:
    x1, y1, z1 = _atom_xyz(first)
    x2, y2, z2 = _atom_xyz(second)
    return (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2


def _nearby_nonprotein(structure: Any, targets: list[Any], radius: float) -> list[Any]:
    cutoff_squared = radius * radius
    target_keys = {_residue_key(residue) for residue in targets}
    target_atoms = [atom for residue in targets for atom in residue.atom]
    nearby: list[Any] = []
    for residue in structure.residue:
        if _residue_key(residue) in target_keys or _residue_name(residue) in _PROTEIN_RESIDUES:
            continue
        if any(
            _squared_distance(atom, target_atom) <= cutoff_squared
            for atom in residue.atom
            for target_atom in target_atoms
        ):
            nearby.append(residue)
    return nearby


def _atom_indices(residues: list[Any]) -> set[int]:
    return {int(atom.index) for residue in residues for atom in residue.atom}


def _fixed_heavy_max_displacement(
    structure: Any,
    before: dict[int, tuple[float, float, float]],
    selected: set[int],
) -> float | None:
    displacements: list[float] = []
    for atom in structure.atom:
        index = int(atom.index)
        if index in selected or str(getattr(atom, "element", "")).strip().upper() == "H":
            continue
        if index not in before:
            continue
        x0, y0, z0 = before[index]
        x1, y1, z1 = _atom_xyz(atom)
        displacements.append(math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2))
    return max(displacements) if displacements else None


def _sequential_mutation(
    structure: Any,
    mutations: list[dict[str, Any]],
    protein_mutator: Any,
) -> tuple[Any, str]:
    current = structure
    processed: list[dict[str, Any]] = []
    for mutation in mutations:
        residue = _find_residue(
            current,
            mutation["chain"],
            mutation["residue_number"],
            mutation.get("insertion_code", ""),
        )
        expected = _AA_1_TO_3[mutation["wt"]]
        observed = _residue_name(residue)
        if observed != expected:
            raise ValueError(
                f"WT mismatch at {mutation['chain']}:{mutation['residue_number']}"
                f"{mutation.get('insertion_code', '')}: expected {expected}, observed {observed}"
            )
        mutation_tuple = (
            mutation["chain"],
            int(mutation["residue_number"]),
            str(getattr(residue, "inscode", "")),
            _AA_1_TO_3[mutation["target"]],
        )
        mutator = protein_mutator(
            current,
            [mutation_tuple],
            concurrent=1,
            sequential=False,
            idealize=False,
        )
        current = next(mutator.generate()).struct
        processed.append(mutation)
        _verify_mutations(current, processed)
    return current, _verify_mutations(current, mutations)


def _local_minimize(
    structure: Any,
    mutations: list[dict[str, Any]],
    local_config: dict[str, Any],
    get_residues_within: Any,
    minimizer: Any,
) -> tuple[Any, dict[str, Any]]:
    radius = float(local_config["radius_angstrom"])
    targets = [
        _find_residue(
            structure,
            mutation["chain"],
            mutation["residue_number"],
            mutation.get("insertion_code", ""),
        )
        for mutation in mutations
    ]
    selected: dict[tuple[str, int, str], Any] = {
        _residue_key(residue): residue for residue in targets
    }
    for target in targets:
        for residue in get_residues_within(structure, [target], within=radius, ca=False):
            if _residue_name(residue) in _PROTEIN_RESIDUES:
                selected[_residue_key(residue)] = residue
    if local_config["include_nearby_nonprotein"]:
        for residue in _nearby_nonprotein(structure, targets, radius):
            selected[_residue_key(residue)] = residue

    selected_residues = list(selected.values())
    selected_atoms = _atom_indices(selected_residues)
    target_atoms = _atom_indices(targets)
    if not selected_atoms or not target_atoms.issubset(selected_atoms):
        raise ValueError("local-minimization selection does not contain every mutated residue")
    total_atoms = int(getattr(structure, "atom_total", len(list(structure.atom))))
    before = {int(atom.index): _atom_xyz(atom) for atom in structure.atom}
    minimized = minimizer.minimize_substructure(structure, selected_atoms)
    if minimized is None:
        localmin = structure
        return_status = "returned_none_used_input_structure"
    elif hasattr(minimized, "atom") and hasattr(minimized, "residue"):
        localmin = minimized
        return_status = "returned_structure_used_return"
    else:
        raise TypeError(
            f"unexpected minimize_substructure return type: {type(minimized).__name__}"
        )
    details = {
        "radius_angstrom": radius,
        "include_nearby_nonprotein": bool(local_config["include_nearby_nonprotein"]),
        "selected_residue_count": len(selected_residues),
        "selected_atom_count": len(selected_atoms),
        "total_atom_count": total_atoms,
        "selected_fraction": (len(selected_atoms) / total_atoms if total_atoms else None),
        "fixed_heavy_atom_max_displacement_angstrom": _fixed_heavy_max_displacement(
            localmin, before, selected_atoms
        ),
        "minimizer_return": return_status,
    }
    return localmin, details


def _write_outputs(
    structure_module: Any,
    modeled: Any,
    staging: Path,
    variant_id: str,
    formats: list[str],
) -> dict[str, str]:
    outputs: dict[str, str] = {}
    for output_format in formats:
        path = staging / f"{variant_id}.{output_format}"
        writer_format = "maestro" if output_format == "maegz" else output_format
        modeled.write(str(path), format=writer_format)
        _read_structure(structure_module, path)
        outputs[output_format] = str(path.resolve())
    return outputs


def _run_variant(
    request: dict[str, Any],
    variant: dict[str, Any],
    modules: tuple[Any, Any, Any, Any],
) -> dict[str, Any]:
    structure_module, get_residues_within, protein_mutator, minimizer = modules
    started = time.perf_counter()
    mutations = list(variant["mutations"])
    try:
        model = _read_structure(structure_module, Path(request["input_structure"]))
        if mutations:
            model, raw_readback = _sequential_mutation(model, mutations, protein_mutator)
            model, details = _local_minimize(
                model,
                mutations,
                request["local_minimization"],
                get_residues_within,
                minimizer,
            )
            local_readback = _verify_mutations(model, mutations)
            readback = "CHECK" if "CHECK" in {raw_readback, local_readback} else "OK"
        else:
            details = {}
            readback = "NOT_REQUESTED"
        outputs = _write_outputs(
            structure_module,
            model,
            Path(request["staging_directory"]),
            variant["variant_id"],
            request["structure_formats"],
        )
        issue_codes = ["MUTATION_READBACK_HIS_ALIAS"] if readback == "CHECK" else []
        return {
            "variant_id": variant["variant_id"],
            "status": "CHECK" if readback == "CHECK" else "OK",
            "mutation_readback_status": readback,
            "issue_codes": issue_codes,
            "outputs": outputs,
            "details": details,
            "duration_seconds": round(time.perf_counter() - started, 6),
            "notes": "sequential mutation from fresh WT and union local minimization completed"
            if mutations
            else "WT structure exported without mutation or minimization",
        }
    except Exception as exc:
        return {
            "variant_id": variant["variant_id"],
            "status": "FAILED",
            "mutation_readback_status": "FAILED" if mutations else "NOT_REQUESTED",
            "issue_codes": ["SCHRODINGER_MODEL_FAILED"],
            "outputs": {},
            "details": {},
            "duration_seconds": round(time.perf_counter() - started, 6),
            "notes": f"{type(exc).__name__}: {exc}",
        }


def main() -> int:
    try:
        request = json.loads(sys.stdin.read())
        staging = Path(request["staging_directory"]).resolve()
        staging.mkdir(parents=True, exist_ok=True)
        modules = _load_modules()
        for variant in request["variants"]:
            _emit(_run_variant(request, variant, modules))
        return 0
    except Exception as exc:
        print(f"MutFlow Schrödinger worker failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
