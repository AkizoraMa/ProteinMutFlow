from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


_MARKER = "MUTFLOW_PROBE_JSON="


@dataclass(frozen=True)
class BackendLocation:
    name: str
    launcher: Path
    root: Path
    source: str


@dataclass(frozen=True)
class BackendProbe:
    name: str
    status: str
    location: BackendLocation | None
    versions: dict[str, str] = field(default_factory=dict)
    message: str = ""

    @property
    def summary(self) -> str:
        if self.location is None:
            return f"{self.name}: {self.status} ({self.message})"
        version_text = ", ".join(f"{key}={value}" for key, value in self.versions.items())
        suffix = f"; {version_text}" if version_text else ""
        return (
            f"{self.name}: {self.status}; launcher={self.location.launcher}; "
            f"source={self.location.source}{suffix}"
        )


def _schrodinger_from_root(root: str | Path, source: str) -> BackendLocation | None:
    root_path = Path(root).expanduser()
    launcher = root_path / ("run.exe" if os.name == "nt" else "run")
    if launcher.is_file():
        return BackendLocation("schrodinger", launcher.resolve(), root_path.resolve(), source)
    return None


def _registry_schrodinger_roots() -> list[Path]:
    if os.name != "nt":
        return []
    try:
        import winreg
    except ImportError:
        return []

    roots: list[Path] = []
    uninstall = r"SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall"
    views = (winreg.KEY_WOW64_64KEY, winreg.KEY_WOW64_32KEY)
    for hive in (winreg.HKEY_LOCAL_MACHINE, winreg.HKEY_CURRENT_USER):
        for view in views:
            try:
                parent = winreg.OpenKey(hive, uninstall, 0, winreg.KEY_READ | view)
            except OSError:
                continue
            with parent:
                for index in range(winreg.QueryInfoKey(parent)[0]):
                    try:
                        name = winreg.EnumKey(parent, index)
                        with winreg.OpenKey(parent, name) as entry:
                            display, _ = winreg.QueryValueEx(entry, "DisplayName")
                            location, _ = winreg.QueryValueEx(entry, "InstallLocation")
                    except OSError:
                        continue
                    if "schrodinger" in str(display).lower() and str(location).strip():
                        roots.append(Path(str(location)))
    return roots


def discover_schrodinger(explicit_home: str | None = None) -> BackendLocation | None:
    if explicit_home:
        return _schrodinger_from_root(explicit_home, "workflow")
    for variable in ("SCHRODINGER", "SCHRODINGER_HOME"):
        if value := os.environ.get(variable):
            if found := _schrodinger_from_root(value, f"environment:{variable}"):
                return found
    if launcher := shutil.which("run.exe" if os.name == "nt" else "run"):
        path = Path(launcher).resolve()
        return BackendLocation("schrodinger", path, path.parent, "PATH")
    seen: set[Path] = set()
    for root in _registry_schrodinger_roots():
        resolved = root.expanduser().resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if found := _schrodinger_from_root(resolved, "windows_registry"):
            return found
    return None


def _pymol_from_candidate(candidate: str | Path, source: str) -> BackendLocation | None:
    path = Path(candidate).expanduser()
    if path.name.lower() == "pymolwin.exe":
        return BackendLocation("pymol", path.resolve(), path.parent.resolve(), source)
    if path.is_file() and path.name.lower() in {"python", "python.exe"}:
        root = path.parent
        if any((root / "conda-meta").glob("pymol*.json")):
            return BackendLocation("pymol", path.resolve(), root.resolve(), source)
    if path.is_file() and path.name.lower() in {"pymol", "pymol.exe"}:
        root = path.parent.parent if path.parent.name.lower() in {"scripts", "bin"} else path.parent
        python = root / ("python.exe" if os.name == "nt" else "bin/python")
        if python.is_file():
            return BackendLocation("pymol", python.resolve(), root.resolve(), source)
    if path.is_dir():
        python = path / ("python.exe" if os.name == "nt" else "bin/python")
        if python.is_file() and any((path / "conda-meta").glob("pymol*.json")):
            return BackendLocation("pymol", python.resolve(), path.resolve(), source)
    return None


def _conda_environment_roots() -> list[Path]:
    roots: list[Path] = []
    if prefix := os.environ.get("CONDA_PREFIX"):
        roots.append(Path(prefix))
    environments_file = Path.home() / ".conda" / "environments.txt"
    if environments_file.is_file():
        try:
            roots.extend(
                Path(line.strip())
                for line in environments_file.read_text(encoding="utf-8").splitlines()
                if line.strip()
            )
        except (OSError, UnicodeError):
            pass
    return roots


def discover_pymol(explicit_launcher: str | None = None) -> BackendLocation | None:
    if explicit_launcher:
        return _pymol_from_candidate(explicit_launcher, "workflow")
    if value := os.environ.get("MUTFLOW_PYMOL_PYTHON"):
        if found := _pymol_from_candidate(value, "environment:MUTFLOW_PYMOL_PYTHON"):
            return found
    for executable in ("pymol.exe", "pymol"):
        if launcher := shutil.which(executable):
            if found := _pymol_from_candidate(launcher, "PATH"):
                return found
    seen: set[Path] = set()
    for root in _conda_environment_roots():
        resolved = root.expanduser().resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if found := _pymol_from_candidate(resolved, "conda_environment_registry"):
            return found
    return None


def _run_probe(command: list[str], timeout_seconds: float) -> tuple[dict[str, str], str]:
    creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout_seconds,
            check=False,
            creationflags=creationflags,
        )
    except subprocess.TimeoutExpired:
        return {}, f"probe timed out after {timeout_seconds:g} seconds"
    except OSError as exc:
        return {}, f"probe could not start: {exc}"
    combined = "\n".join(part for part in (completed.stdout, completed.stderr) if part)
    if completed.returncode != 0:
        tail = combined.strip().splitlines()[-1] if combined.strip() else "no diagnostic output"
        return {}, f"probe exited {completed.returncode}: {tail}"
    for line in combined.splitlines():
        if line.startswith(_MARKER):
            try:
                payload = json.loads(line[len(_MARKER):])
            except json.JSONDecodeError as exc:
                return {}, f"probe returned invalid JSON: {exc}"
            return {str(key): str(value) for key, value in payload.items()}, ""
    return {}, "probe succeeded but did not return the expected marker"


def probe_schrodinger(location: BackendLocation, timeout_seconds: float = 30) -> BackendProbe:
    version_file = location.root / "version.txt"
    suite_version = ""
    if version_file.is_file():
        try:
            suite_version = version_file.read_text(encoding="utf-8", errors="replace").strip()
        except OSError:
            pass
    code = (
        "import json,sys; "
        "from schrodinger import structure; "
        "from schrodinger.application.bioluminate.protein import get_residues_within; "
        "from schrodinger.application.bioluminate.protein.mutator import ProteinMutator; "
        "from schrodinger.forcefield import minimizer; "
        "print('MUTFLOW_PROBE_JSON='+json.dumps({"
        "'python':sys.version.split()[0],"
        "'structure':bool(structure),"
        "'protein_mutator':bool(ProteinMutator),"
        "'get_residues_within':bool(get_residues_within),"
        "'minimize_substructure':hasattr(minimizer,'minimize_substructure')},sort_keys=True))"
    )
    versions, message = _run_probe(
        [str(location.launcher), "python3", "-B", "-c", code], timeout_seconds
    )
    if suite_version:
        versions = {"suite": suite_version, **versions}
    required_true = (
        versions.get("structure") == "True"
        and versions.get("protein_mutator") == "True"
        and versions.get("get_residues_within") == "True"
        and versions.get("minimize_substructure") == "True"
    )
    status = "OK" if not message and required_true else "FAILED"
    if not message and not required_true:
        message = "one or more required Schrödinger APIs are unavailable"
    return BackendProbe("schrodinger", status, location, versions, message)


def probe_pymol(location: BackendLocation, timeout_seconds: float = 30) -> BackendProbe:
    if location.launcher.name.lower() == "pymolwin.exe":
        return BackendProbe(
            "pymol",
            "FAILED",
            location,
            message="PyMOLWin.exe is not an environment-aware headless launcher",
        )
    code = (
        "import json,sys; import pymol; from pymol import cmd; "
        "version=cmd.get_version(); "
        "print('MUTFLOW_PROBE_JSON='+json.dumps({"
        "'python':sys.version.split()[0],'pymol':str(version[0])},sort_keys=True))"
    )
    versions, message = _run_probe(
        [str(location.launcher), "-B", "-c", code], timeout_seconds
    )
    status = "OK" if not message and versions.get("pymol") else "FAILED"
    return BackendProbe("pymol", status, location, versions, message)


def inspect_schrodinger(
    explicit_home: str | None = None,
    *,
    run_probe: bool = True,
    auto_discover: bool = True,
) -> BackendProbe:
    if not auto_discover and not explicit_home:
        return BackendProbe(
            "schrodinger", "NOT_PROBED", None, message="discovery and runtime probe disabled"
        )
    location = discover_schrodinger(explicit_home) if auto_discover or explicit_home else None
    if location is None:
        return BackendProbe("schrodinger", "NOT_FOUND", None, message="installation not found")
    if not run_probe:
        return BackendProbe("schrodinger", "NOT_PROBED", location, message="runtime probe disabled")
    return probe_schrodinger(location)


def inspect_pymol(
    explicit_launcher: str | None = None,
    *,
    run_probe: bool = True,
    auto_discover: bool = True,
) -> BackendProbe:
    if not auto_discover and not explicit_launcher:
        return BackendProbe(
            "pymol", "NOT_PROBED", None, message="discovery and runtime probe disabled"
        )
    location = discover_pymol(explicit_launcher) if auto_discover or explicit_launcher else None
    if location is None:
        return BackendProbe("pymol", "NOT_FOUND", None, message="environment not found")
    if location.launcher.name.lower() == "pymolwin.exe":
        return BackendProbe(
            "pymol", "FAILED", location, message="unsafe direct PyMOLWin.exe launcher"
        )
    if not run_probe:
        return BackendProbe("pymol", "NOT_PROBED", location, message="runtime probe disabled")
    return probe_pymol(location)
