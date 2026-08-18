from __future__ import annotations

import json
import os
import subprocess
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from mutflow.backends import BackendLocation
from mutflow.preflight import PlannedVariant


_RESULT_MARKER = "MUTFLOW_MODEL_JSON="


class ModelingExecutionError(RuntimeError):
    """Raised when the external Schrödinger worker cannot complete its protocol."""


def _request_variant(variant: PlannedVariant) -> dict[str, Any]:
    return {
        "variant_id": variant.variant_id,
        "mutations": [
            {
                "chain": mutation.chain or "",
                "wt": mutation.wt,
                "residue_number": mutation.residue_number,
                "insertion_code": mutation.insertion_code,
                "target": mutation.target,
            }
            for mutation in variant.mutations
        ],
    }


def iter_schrodinger_models(
    location: BackendLocation,
    config: dict[str, Any],
    variants: list[PlannedVariant],
    input_structure: Path,
    staging_directory: Path,
) -> Iterator[dict[str, Any]]:
    worker = Path(__file__).with_name("schrodinger_worker.py").resolve()
    request = {
        "input_structure": str(input_structure.resolve()),
        "staging_directory": str(staging_directory.resolve()),
        "structure_formats": list(config["output"]["structure_formats"]),
        "local_minimization": dict(config["modeling"]["local_minimization"]),
        "variants": [_request_variant(variant) for variant in variants],
    }
    requested = {variant.variant_id for variant in variants}
    seen: set[str] = set()
    diagnostic_tail: list[str] = []
    creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    try:
        process = subprocess.Popen(
            [str(location.launcher), "python3", "-B", str(worker)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
            creationflags=creationflags,
        )
    except OSError as exc:
        raise ModelingExecutionError(f"Schrödinger model worker could not start: {exc}") from exc

    try:
        assert process.stdin is not None
        assert process.stdout is not None
        process.stdin.write(json.dumps(request))
        process.stdin.close()
        for line in process.stdout:
            text = line.rstrip("\r\n")
            if not text.startswith(_RESULT_MARKER):
                if text:
                    diagnostic_tail.append(text)
                    diagnostic_tail = diagnostic_tail[-10:]
                continue
            try:
                result = json.loads(text[len(_RESULT_MARKER) :])
            except json.JSONDecodeError as exc:
                raise ModelingExecutionError(
                    f"Schrödinger model worker returned invalid JSON: {exc}"
                ) from exc
            variant_id = str(result.get("variant_id", ""))
            if variant_id not in requested or variant_id in seen:
                raise ModelingExecutionError(
                    f"Schrödinger model worker returned unexpected or duplicate variant_id: {variant_id!r}"
                )
            seen.add(variant_id)
            yield result
        return_code = process.wait()
        if return_code != 0:
            tail = diagnostic_tail[-1] if diagnostic_tail else "no diagnostic output"
            raise ModelingExecutionError(
                f"Schrödinger model worker exited {return_code}: {tail}"
            )
        missing = requested - seen
        if missing:
            raise ModelingExecutionError(
                f"Schrödinger model worker returned no result for: {', '.join(sorted(missing))}"
            )
    finally:
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
