from __future__ import annotations

import csv
import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


RESULT_FIELDS = [
    "schema_version",
    "variant_id",
    "mutation_name",
    "mutation_count",
    "structure_path",
    "model_status",
    "mutation_readback_status",
    "metrics_status",
    "qc_status",
    "overall_status",
    "duration_seconds",
    "issue_codes",
    "notes",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _atomic_target_temp(target: Path) -> tuple[int, Path]:
    target.parent.mkdir(parents=True, exist_ok=True)
    temp_parent = (
        target.parent.parent.parent
        if target.parent.name == "structures"
        else target.parent.parent
    )
    temp_parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=temp_parent)
    return fd, Path(name)


def atomic_write_text(target: Path, text: str) -> None:
    fd, temporary = _atomic_target_temp(target)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(target)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def atomic_write_json(target: Path, payload: dict[str, Any]) -> None:
    atomic_write_text(target, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def atomic_write_results(
    target: Path,
    rows: list[dict[str, Any]],
    fieldnames: list[str] | None = None,
) -> None:
    fd, temporary = _atomic_target_temp(target)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames or RESULT_FIELDS,
                extrasaction="raise",
            )
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(target)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def atomic_copy(source: Path, target: Path) -> None:
    fd, temporary = _atomic_target_temp(target)
    os.close(fd)
    try:
        shutil.copyfile(source, temporary)
        temporary.replace(target)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


@dataclass
class RunStore:
    output: Path
    rows: list[dict[str, Any]]
    run_document: dict[str, Any]
    fieldnames: list[str]

    @property
    def structures(self) -> Path:
        return self.output / "structures"

    @property
    def results_path(self) -> Path:
        return self.output / "results.csv"

    @property
    def run_path(self) -> Path:
        return self.output / "run.json"

    @property
    def log_path(self) -> Path:
        return self.output / "run.log"

    @classmethod
    def create(
        cls,
        output: Path,
        rows: list[dict[str, Any]],
        run_document: dict[str, Any],
        fieldnames: list[str] | None = None,
    ) -> "RunStore":
        if output.exists():
            raise FileExistsError(f"output already exists: {output}")
        output.mkdir(parents=True, exist_ok=False)
        store = cls(output, rows, run_document, fieldnames or list(RESULT_FIELDS))
        store.structures.mkdir()
        try:
            atomic_write_results(store.results_path, store.rows, store.fieldnames)
            atomic_write_json(store.run_path, store.run_document)
            atomic_write_text(
                store.log_path,
                f"{utc_now()} INFO run_initialized run_id={run_document['run_id']}\n",
            )
        except Exception:
            # The newly created directory contains only MutFlow-owned partial initialization.
            shutil.rmtree(output, ignore_errors=True)
            raise
        return store

    @classmethod
    def open_existing(cls, output: Path) -> "RunStore":
        expected = {"structures", "results.csv", "run.json", "run.log"}
        if not output.is_dir():
            raise FileNotFoundError(f"resume output directory not found: {output}")
        actual = {item.name for item in output.iterdir()}
        if actual != expected:
            raise ValueError(
                f"resume output tree mismatch: expected {sorted(expected)}, got {sorted(actual)}"
            )
        structures = output / "structures"
        if not structures.is_dir():
            raise ValueError("resume structures path is not a directory")
        for name in ("results.csv", "run.json", "run.log"):
            if not (output / name).is_file():
                raise ValueError(f"resume {name} is not a file")
        try:
            run_document = json.loads((output / "run.json").read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise ValueError(f"cannot read resume run.json: {exc}") from exc
        try:
            with (output / "results.csv").open(
                "r", encoding="utf-8-sig", newline=""
            ) as handle:
                reader = csv.DictReader(handle)
                fieldnames = list(reader.fieldnames or [])
                rows = list(reader)
        except (OSError, UnicodeError, csv.Error) as exc:
            raise ValueError(f"cannot read resume results.csv: {exc}") from exc
        if not fieldnames:
            raise ValueError("resume results.csv has no header")
        for row in rows:
            try:
                row["mutation_count"] = int(row["mutation_count"])
                row["duration_seconds"] = (
                    None if row["duration_seconds"] == "" else float(row["duration_seconds"])
                )
                for name in fieldnames:
                    if not name.startswith("metric__") or row[name] == "":
                        continue
                    if name.endswith("__clash_count") or name.endswith(
                        "__severe_clash_count"
                    ):
                        row[name] = int(row[name])
                    elif not name.endswith("__worst_pair"):
                        row[name] = float(row[name])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"cannot parse resume results row: {exc}") from exc
        return cls(output.resolve(), rows, run_document, fieldnames)

    def begin_resume(self, command: str) -> None:
        execution = self.run_document.setdefault("execution", {})
        execution["new_or_resumed"] = "resumed"
        execution["resume_count"] = int(execution.get("resume_count", 0)) + 1
        execution["last_resume_command"] = command
        self.run_document["state"] = "RUNNING"
        self.run_document["finished_at"] = None
        atomic_write_json(self.run_path, self.run_document)
        self.append_log(
            "INFO", "run_resumed", f"resume_count={execution['resume_count']}"
        )

    def append_log(self, level: str, event: str, detail: str = "") -> None:
        safe_detail = detail.replace("\r", " ").replace("\n", " ").strip()
        suffix = f" {safe_detail}" if safe_detail else ""
        with self.log_path.open("a", encoding="utf-8", newline="") as handle:
            handle.write(f"{utc_now()} {level} {event}{suffix}\n")

    def update_row(self, variant_id: str, **updates: Any) -> None:
        matches = [row for row in self.rows if row["variant_id"] == variant_id]
        if len(matches) != 1:
            raise KeyError(f"expected one results row for {variant_id!r}, got {len(matches)}")
        unknown = set(updates) - set(self.fieldnames)
        if unknown:
            raise KeyError(f"unknown results fields: {sorted(unknown)}")
        matches[0].update(updates)
        atomic_write_results(self.results_path, self.rows, self.fieldnames)

    def finalize(self) -> None:
        counts = {"ok": 0, "check": 0, "failed": 0, "not_run": 0}
        for row in self.rows:
            status = row["overall_status"]
            key = {"OK": "ok", "CHECK": "check", "FAILED": "failed", "NOT_RUN": "not_run"}[status]
            counts[key] += 1
        counts["requested"] = len(self.rows)
        self.run_document["counts"] = {
            "requested": counts["requested"],
            "ok": counts["ok"],
            "check": counts["check"],
            "failed": counts["failed"],
            "not_run": counts["not_run"],
        }
        if counts["not_run"]:
            state = "INTERRUPTED"
        elif counts["failed"]:
            state = "FAILED"
        elif counts["check"]:
            state = "COMPLETED_WITH_CHECKS"
        else:
            state = "COMPLETED"
        self.run_document["state"] = state
        self.run_document["finished_at"] = utc_now()
        atomic_write_json(self.run_path, self.run_document)
        self.append_log("INFO", "run_finalized", f"state={state}")

    def mark_interrupted(self) -> None:
        counts = {"ok": 0, "check": 0, "failed": 0, "not_run": 0}
        for row in self.rows:
            key = {
                "OK": "ok",
                "CHECK": "check",
                "FAILED": "failed",
                "NOT_RUN": "not_run",
            }[row["overall_status"]]
            counts[key] += 1
        self.run_document["counts"] = {
            "requested": len(self.rows),
            **counts,
        }
        self.run_document["state"] = "INTERRUPTED"
        self.run_document["finished_at"] = utc_now()
        atomic_write_json(self.run_path, self.run_document)
        self.append_log("WARNING", "run_interrupted", f"not_run={counts['not_run']}")
