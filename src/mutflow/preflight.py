from __future__ import annotations

import csv
import hashlib
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mutflow.backends import BackendProbe, inspect_pymol, inspect_schrodinger
from mutflow.configuration import ConfigurationError, load_workflow
from mutflow.mutations import (
    Mutation,
    STANDARD_AMINO_ACIDS,
    canonical_variant,
    parse_variant,
    variant_id_from_name,
)
from mutflow.pdbio import PDBError, ResidueIdentity, read_protein_residues


class PreflightError(ValueError):
    """Raised when preflight cannot even construct a report."""


@dataclass(frozen=True)
class Issue:
    level: str
    code: str
    message: str


@dataclass(frozen=True)
class PlannedVariant:
    variant_id: str
    mutation_name: str
    mutations: tuple[Mutation, ...]
    source_structure: Path | None = None

    @property
    def mutation_count(self) -> int:
        return len(self.mutations)


@dataclass
class PreflightReport:
    workflow_path: Path
    input_path: Path
    input_sha256: str
    variant_mode: str
    variants: list[PlannedVariant]
    modeling_enabled: bool
    metric_names: list[str]
    output_path: Path
    structure_formats: list[str]
    max_variants: int
    backend_summaries: list[str] = field(default_factory=list)
    issues: list[Issue] = field(default_factory=list)

    @property
    def status(self) -> str:
        if any(issue.level == "ERROR" for issue in self.issues):
            return "REFUSED"
        if any(issue.level == "CHECK" for issue in self.issues):
            return "CHECK"
        return "READY"

    def render(self) -> str:
        by_count = Counter(variant.mutation_count for variant in self.variants)
        counts = ", ".join(f"{count}:{by_count[count]}" for count in sorted(by_count)) or "none"
        metrics = ", ".join(self.metric_names) if self.metric_names else "none"
        final_file_count = len(self.variants) * len(self.structure_formats) + 3
        lines = [
            "MUTFLOW PREFLIGHT",
            f"workflow: {self.workflow_path}",
            f"input: {self.input_path}",
            f"input_sha256: {self.input_sha256}",
            f"variant_mode: {self.variant_mode}",
            f"variant_count: {len(self.variants)}",
            f"variant_limit: {self.max_variants}",
            f"variants_by_mutation_count: {counts}",
            f"modeling: {'enabled' if self.modeling_enabled else 'disabled'}",
            f"metrics: {metrics}",
            f"output: {self.output_path}",
            f"planned_final_file_count: {final_file_count}",
            f"core_python: {sys.version.split()[0]}",
        ]
        if self.backend_summaries:
            lines.append("backends:")
            lines.extend(f"  - {summary}" for summary in self.backend_summaries)
        else:
            lines.append("backends: none requested")
        if self.issues:
            lines.append("issues:")
            lines.extend(
                f"  - [{issue.level}] {issue.code}: {issue.message}" for issue in self.issues
            )
        else:
            lines.append("issues: none")
        lines.append(f"status: {self.status}")
        lines.append("files_written: 0")
        lines.append("scientific_backends_invoked: 0")
        return "\n".join(lines)


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


def _resolve_mutation(
    mutation: Mutation,
    residues: dict[tuple[str, int, str], ResidueIdentity],
) -> Mutation:
    if mutation.chain is not None:
        resolved = mutation
        identity = residues.get(
            (mutation.chain, mutation.residue_number, mutation.insertion_code)
        )
        if identity is None:
            raise ValueError(
                f"residue not found for mutation {mutation.canonical}: "
                f"chain {mutation.chain!r} residue {mutation.residue_number}{mutation.insertion_code}"
            )
    else:
        candidates = [
            identity
            for identity in residues.values()
            if identity.residue_number == mutation.residue_number
            and identity.insertion_code == mutation.insertion_code
            and identity.one_letter == mutation.wt
        ]
        if len(candidates) != 1:
            raise ValueError(
                f"mutation {mutation.canonical} has no chain and resolves to "
                f"{len(candidates)} matching residues; specify the chain explicitly"
            )
        identity = candidates[0]
        resolved = mutation.with_chain(identity.chain)
    if identity.one_letter != mutation.wt:
        raise ValueError(
            f"WT mismatch for {resolved.canonical}: structure contains "
            f"{identity.one_letter} ({identity.residue_name})"
        )
    return resolved


def _expand_saturation(
    config: dict[str, Any],
    residues: dict[tuple[str, int, str], ResidueIdentity],
) -> list[PlannedVariant]:
    variants: list[PlannedVariant] = []
    if config["variants"]["include_wt"]:
        variants.append(PlannedVariant("WT", "WT", ()))
    seen_sites: set[tuple[str, int, str]] = set()
    for site in config["variants"]["sites"]:
        chain = site["chain"]
        residue_number = site["residue_number"]
        insertion_code = site.get("insertion_code", "").upper()
        key = (chain, residue_number, insertion_code)
        if key in seen_sites:
            raise ValueError(
                f"duplicate saturation site: chain {chain!r} residue {residue_number}{insertion_code}"
            )
        seen_sites.add(key)
        identity = residues.get(key)
        if identity is None:
            raise ValueError(
                f"saturation site not found: chain {chain!r} residue "
                f"{residue_number}{insertion_code}"
            )
        expected = site.get("expected_wt")
        if expected and identity.one_letter != expected:
            raise ValueError(
                f"expected WT mismatch at chain {chain!r} residue {residue_number}{insertion_code}: "
                f"configured {expected}, structure {identity.one_letter}"
            )
        for target in STANDARD_AMINO_ACIDS:
            if target == identity.one_letter:
                continue
            mutation = Mutation(
                chain=chain,
                wt=identity.one_letter,
                residue_number=residue_number,
                insertion_code=insertion_code,
                target=target,
            )
            name = mutation.canonical
            variants.append(
                PlannedVariant(variant_id_from_name(name), name, (mutation,))
            )
    return variants


def _read_explicit_variants(
    config: dict[str, Any],
    base: Path,
    wt_residues: dict[tuple[str, int, str], ResidueIdentity],
) -> list[PlannedVariant]:
    csv_path = _resolve(base, config["variants"]["file"])
    if not csv_path.is_file():
        raise ValueError(f"explicit mutation file not found: {csv_path}")
    try:
        with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            headers = reader.fieldnames or []
            if not {"variant_id", "mutations"}.issubset(headers):
                raise ValueError("explicit CSV requires variant_id and mutations columns")
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValueError(f"cannot read explicit mutation CSV: {exc}") from exc
    if not rows:
        raise ValueError("explicit mutation CSV has no rows")

    modeling_enabled = config["modeling"]["enabled"]
    if modeling_enabled and "structure_path" in headers and any(
        (row.get("structure_path") or "").strip() for row in rows
    ):
        raise ValueError("structure_path must be empty when modeling is enabled")
    if not modeling_enabled and "structure_path" not in headers:
        raise ValueError("metrics-only explicit CSV requires a structure_path column")

    variants: list[PlannedVariant] = []
    if config["variants"]["include_wt"]:
        variants.append(PlannedVariant("WT", "WT", ()))
    seen_ids = {"WT"} if config["variants"]["include_wt"] else set()
    seen_names: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        variant_id = (row.get("variant_id") or "").strip()
        mutation_text = (row.get("mutations") or "").strip()
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", variant_id):
            raise ValueError(f"row {row_number}: unsafe or empty variant_id {variant_id!r}")
        if variant_id in seen_ids:
            raise ValueError(f"row {row_number}: duplicate variant_id {variant_id!r}")
        seen_ids.add(variant_id)
        try:
            mutations = tuple(
                _resolve_mutation(mutation, wt_residues)
                for mutation in parse_variant(mutation_text)
            )
        except ValueError as exc:
            raise ValueError(f"row {row_number}: {exc}") from exc
        mutation_name = canonical_variant(mutations)
        if mutation_name in seen_names:
            raise ValueError(f"row {row_number}: duplicate normalized variant {mutation_name}")
        seen_names.add(mutation_name)

        source_structure: Path | None = None
        if not modeling_enabled:
            structure_value = (row.get("structure_path") or "").strip()
            if not structure_value:
                raise ValueError(f"row {row_number}: missing structure_path")
            source_structure = _resolve(base, structure_value)
            if not source_structure.is_file():
                raise ValueError(
                    f"row {row_number}: existing structure not found: {source_structure}"
                )
            if source_structure.suffix.lower() != ".pdb":
                raise ValueError(
                    f"row {row_number}: metrics-only vertical slice currently requires PDB"
                )
            observed = read_protein_residues(source_structure)
            for mutation in mutations:
                identity = observed.get(
                    (mutation.chain or "", mutation.residue_number, mutation.insertion_code)
                )
                if identity is None or identity.one_letter != mutation.target:
                    found = identity.one_letter if identity else "missing"
                    raise ValueError(
                        f"row {row_number}: mutation readback mismatch for "
                        f"{mutation.canonical}; observed {found}"
                    )
        variants.append(
            PlannedVariant(variant_id, mutation_name, mutations, source_structure)
        )
    return variants


def _metric_names(config: dict[str, Any]) -> list[str]:
    metrics = config["metrics"]
    for metric in metrics["clashes"]:
        if metric["severe_overlap_angstrom"] < metric["clash_overlap_angstrom"]:
            raise ValueError(
                f"clash metric {metric['name']!r} has severe threshold below clash threshold"
            )
    names = [
        item["name"]
        for section in ("sasa", "minimum_distance", "clashes")
        for item in metrics[section]
    ]
    normalized = [re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_") for name in names]
    if len(normalized) != len(set(normalized)):
        raise ValueError("configured metric names collide after column normalization")
    qc_names = [item["name"] for item in config["quality_control"]["required_selections"]]
    if len(qc_names) != len(set(qc_names)):
        raise ValueError("duplicate required-selection QC name")
    return names


def _probe_issue(probe: BackendProbe) -> Issue | None:
    if probe.status == "OK":
        return None
    if probe.status == "NOT_PROBED":
        return Issue("CHECK", f"{probe.name.upper()}_NOT_PROBED", probe.message)
    return Issue("ERROR", f"{probe.name.upper()}_{probe.status}", probe.message)


def _backend_inspection(
    config: dict[str, Any],
    base: Path,
    metrics_requested: bool,
    *,
    run_probes: bool,
    auto_discover: bool,
) -> tuple[list[Issue], list[str]]:
    issues: list[Issue] = []
    summaries: list[str] = []
    backends = config["backends"]
    if config["modeling"]["enabled"]:
        configured = backends.get("schrodinger", {}).get("home")
        if configured and not Path(configured).expanduser().is_absolute():
            configured = str(_resolve(base, configured))
        probe = inspect_schrodinger(
            configured, run_probe=run_probes, auto_discover=auto_discover
        )
        summaries.append(probe.summary)
        if issue := _probe_issue(probe):
            issues.append(issue)

    if metrics_requested:
        configured = backends.get("pymol", {}).get("launcher")
        if configured and not Path(configured).expanduser().is_absolute():
            configured = str(_resolve(base, configured))
        probe = inspect_pymol(
            configured, run_probe=run_probes, auto_discover=auto_discover
        )
        summaries.append(probe.summary)
        if issue := _probe_issue(probe):
            issues.append(issue)
    return issues, summaries


def run_preflight(
    workflow_path: Path,
    *,
    run_backend_probes: bool = True,
    auto_discover_backends: bool = True,
    resume: bool = False,
) -> PreflightReport:
    try:
        workflow_path = workflow_path.expanduser().resolve()
        config = load_workflow(workflow_path)
    except ConfigurationError as exc:
        raise PreflightError(str(exc)) from exc
    base = workflow_path.parent
    input_path = _resolve(base, config["input"]["structure"])
    if not input_path.is_file():
        raise PreflightError(f"input structure not found: {input_path}")

    configured_format = config["input"]["format"]
    structure_format = input_path.suffix.lower().lstrip(".") if configured_format == "auto" else configured_format
    if structure_format != "pdb":
        raise PreflightError(
            "this first preflight slice can inspect PDB only; MAEGZ support belongs in the Schrödinger adapter"
        )
    try:
        residues = read_protein_residues(input_path)
    except PDBError as exc:
        raise PreflightError(str(exc)) from exc

    try:
        if config["variants"]["mode"] == "saturation_single":
            if not config["modeling"]["enabled"]:
                raise ValueError("saturation mode requires modeling.enabled: true")
            variants = _expand_saturation(config, residues)
        else:
            variants = _read_explicit_variants(config, base, residues)
        metric_names = _metric_names(config)
    except (ValueError, PDBError) as exc:
        raise PreflightError(str(exc)) from exc

    wants_delta = any(
        metric.get("calculate_wt_delta", False)
        for section in ("sasa", "minimum_distance")
        for metric in config["metrics"][section]
    )
    issues: list[Issue] = []
    max_variants = int(config["execution"]["max_variants"])
    if len(variants) > max_variants:
        issues.append(
            Issue(
                "ERROR",
                "VARIANT_LIMIT_EXCEEDED",
                (
                    f"expanded variant count {len(variants)} exceeds configured limit "
                    f"{max_variants}; increase execution.max_variants only after "
                    "reviewing the expanded batch"
                ),
            )
        )
    if wants_delta and not config["variants"]["include_wt"]:
        issues.append(
            Issue("ERROR", "WT_BASELINE_REQUIRED", "a WT-relative metric was requested without include_wt")
        )

    output_path = _resolve(base, config["output"]["directory"])
    if resume:
        if not output_path.exists():
            issues.append(
                Issue(
                    "ERROR",
                    "RESUME_OUTPUT_MISSING",
                    f"resume output does not exist: {output_path}",
                )
            )
        elif not output_path.is_dir():
            issues.append(
                Issue(
                    "ERROR",
                    "RESUME_OUTPUT_INVALID",
                    f"resume output is not a directory: {output_path}",
                )
            )
    elif output_path.exists():
        nonempty = output_path.is_file() or any(output_path.iterdir())
        if nonempty:
            issues.append(
                Issue("ERROR", "OUTPUT_NOT_EMPTY", f"existing output would be reused: {output_path}")
            )
    if not config["modeling"]["enabled"] and "maegz" in config["output"]["structure_formats"]:
        issues.append(
            Issue("ERROR", "MAEGZ_REQUIRES_SCHRODINGER", "metrics-only mode cannot create MAEGZ output")
        )

    metrics_requested = bool(metric_names or config["quality_control"]["required_selections"])
    if (
        config["modeling"]["enabled"]
        and metrics_requested
        and "pdb" not in config["output"]["structure_formats"]
    ):
        issues.append(
            Issue(
                "ERROR",
                "PYMOL_METRICS_REQUIRE_PDB_OUTPUT",
                "a modeling workflow with PyMOL metrics must retain PDB output",
            )
        )
    backend_issues, backend_summaries = _backend_inspection(
        config,
        base,
        metrics_requested,
        run_probes=run_backend_probes,
        auto_discover=auto_discover_backends,
    )
    issues.extend(backend_issues)
    return PreflightReport(
        workflow_path=workflow_path,
        input_path=input_path,
        input_sha256=_sha256(input_path),
        variant_mode=config["variants"]["mode"],
        variants=variants,
        modeling_enabled=config["modeling"]["enabled"],
        metric_names=metric_names,
        output_path=output_path,
        structure_formats=list(config["output"]["structure_formats"]),
        max_variants=max_variants,
        backend_summaries=backend_summaries,
        issues=issues,
    )
