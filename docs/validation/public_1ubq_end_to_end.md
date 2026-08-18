# Public 1UBQ End-to-End Acceptance

Status: `ACCEPTED`.

On 2026-08-18, an installed `protein-mutflow` alpha package completed the documented public workflow through real Schrödinger modeling and headless PyMOL metrics/QC.

Reference environment:

- core Python 3.12.6;
- Schrödinger Suite 2025-1, Build 129, Python 3.11.4;
- PyMOL 3.0.3, Python 3.10.14;
- Windows.

Accepted objects:

| Variant | Mutation | Overall |
|---|---|---|
| `WT` | `WT` | `OK` |
| `M001` | `A:G10A` | `OK` |
| `M002` | `A:G10A_A:K11A` | `OK` |

Mutation readback, finite configured metrics, selection-count QC, schemas, counts, and the compact output tree passed. A completed-run `--resume` attempt was refused without changing or adding output files.

Downloaded/prepared PDBs and generated metric values are not distributed. This acceptance demonstrates workflow integration, not biological or experimental validity.
