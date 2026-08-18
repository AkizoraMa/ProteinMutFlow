# Contributing

ProteinMutFlow welcomes narrowly scoped fixes, tests, documentation improvements, and backend compatibility reports.

## Development Setup

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install .
.venv\Scripts\python -m unittest discover -s tests -v
```

## Contribution Rules

- Keep protein sites, chains, ligands, selections, and paths configurable; do not add project-specific constants to `src/`.
- Add or update tests for behavior changes.
- Preserve the compact default output contract.
- Do not commit structures, result tables, license material, credentials, or private project data without explicit redistribution approval.
- Do not claim a backend or operating-system combination is supported without reproducible evidence.
- Separate structural workflow evidence from biological or experimental interpretation.

Licensed Schrödinger jobs and local PyMOL compatibility runs are not expected in public CI. A compatibility contribution should include versions, a public input, the exact configuration, and a data-minimized acceptance summary.
