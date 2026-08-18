# Installation

## Core Environment

ProteinMutFlow currently targets CPython `>=3.11,<3.13` in an isolated virtual environment.

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install .
.venv\Scripts\mutflow --version
```

The core environment contains the CLI, YAML/schema validation, orchestration, and state handling. Do not install these packages into vendor-managed Schrödinger or PyMOL Python environments.

## Capability-Based Dependencies

Schrödinger is required only for mutation modeling and local minimization. PyMOL is required only for configured SASA, distance, clash, and selection-count work. `mutflow init` remains usable without either backend; `preflight` refuses a requested capability when its backend is unavailable.

The verified reference environment is:

- Windows;
- core Python 3.12.6;
- Schrödinger Suite 2025-1, Build 129, with Python 3.11.4;
- PyMOL 3.0.3 with Python 3.10.14.

These are reference versions, not minimum compatibility claims.

## Backend Discovery

ProteinMutFlow can use explicit backend paths in `workflow.yaml` or discover common Windows/Conda installations. PyMOL must be launched through its environment-aware Python runtime; direct `PyMOLWin.exe` launching is rejected.

Run preflight before scientific work:

```powershell
mutflow preflight workflow.yaml
```

Preflight may perform short import/version probes, but it writes no file and performs no mutation, minimization, or metric calculation.

## Prepared-Structure Requirement

Schrödinger force-field operations require an appropriately prepared input. PDB readability is not evidence of chemical preparation. ProteinMutFlow alpha does not silently run Protein Preparation Wizard.
