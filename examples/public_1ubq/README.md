# Public 1UBQ End-to-End Example

- Example status: `END_TO_END_ACCEPTED`
- Real modeling-only smoke: `ACCEPTED`
- Real modeling-plus-PyMOL run: `ACCEPTED`
- Data source: public RCSB 1UBQ; no ecHint structure, mutation table, or result is included

This example is the public, non-biological demonstration of ProteinMutFlow's ordinary path. It requests WT, one single substitution, and one double substitution, then measures generic whole-protein and mutation-region geometry. The metrics demonstrate orchestration and traceability; they are not a score, activity prediction, or biological recommendation.

## 1. Obtain The Public Structure

Download the PDB-format 1UBQ structure from RCSB:

```text
https://files.rcsb.org/download/1UBQ.pdb
```

Place it in this directory as `1UBQ.pdb`. The accepted download had SHA256:

```text
d4a6812d8951cf6594e6a0763f089e35f5a80b62acb3c117b2c5565228a7b161
```

On PowerShell, inspect the downloaded file with:

```powershell
Get-FileHash -Algorithm SHA256 .\1UBQ.pdb
```

## 2. Prepare The Structure

Schrödinger mutation/minimization requires a chemically prepared input. A readable raw crystal PDB is not sufficient, and ProteinMutFlow alpha does not run Protein Preparation Wizard automatically.

The accepted adapter smoke used Schrödinger 2025-1 and the deliberately limited reproducibility command below:

```powershell
prepwizard -noepik -noprotassign -noimpref .\1UBQ.pdb .\1UBQ_prepared.pdb
```

Run it from a Schrödinger-aware shell or replace `prepwizard` with the absolute launcher path for the local installation. These settings reproduce the adapter smoke protocol; they are not a universal preparation recommendation for scientific studies.

The accepted prepared file had SHA256:

```text
c0afbc9ecbed64ceb113cb2550873a617460fc1f468fc76f581276ce8731ddde
```

Preparation output may differ across Schrödinger releases. A different hash therefore requires a new compatibility record rather than being silently presented as the accepted reference input.

## 3. Requested Objects

`mutations.csv` defines:

| Variant | Mutation |
|---|---|
| `WT` | automatically included baseline |
| `M001` | `A:G10A` |
| `M002` | `A:G10A_A:K11A` |

Every modeled variant starts from a fresh prepared WT. The double variant is applied sequentially, and local minimization uses the union of residues within 5 angstrom of the mutated sites.

## 4. Demonstration Metrics

`workflow.yaml` requests:

- whole-protein SASA and mutant-minus-WT delta;
- residue 10/11 region SASA and mutant-minus-WT delta;
- minimum heavy-atom distance from residue 10/11 sidechains to the distal protein, excluding residues 9-12 from the opposite selection;
- transparent VDW overlap counts for the same non-neighbor selections;
- exact residue 10/11 backbone atom-count QC;
- non-empty residue 10/11 heavy-sidechain QC.

The distance/clash selections exclude adjacent residues to avoid treating peptide-bond neighbors as the intended demonstration contact. `exclude_bonded_pairs` remains `false`, and that method choice is recorded in `run.json`.

## 5. Run

From this directory:

```powershell
mutflow preflight .\workflow.yaml
mutflow run .\workflow.yaml
```

Preflight may perform short headless import/version probes, but it writes no file and performs no scientific calculation. `run` invokes the separately installed and licensed Schrödinger backend and the configured/discovered headless PyMOL environment.

The successful default tree must remain exactly:

```text
output/
├── structures/
├── results.csv
├── run.json
└── run.log
```

If execution is interrupted after the output is safely initialized:

```powershell
mutflow run --resume .\workflow.yaml
```

Resume is compatibility-gated and processes only `NOT_RUN` stages. It does not automatically retry terminal failures.

## 6. Acceptance Boundary

An explicitly authorized installed-package run on 2026-08-18 completed the documented path through Schrödinger modeling and headless PyMOL metrics/QC with WT, `M001`, and `M002` all `OK`. Mutation readback, finite numeric metrics, selection QC, schema/count reconciliation, and the exact compact output tree passed on the recorded Windows reference environment: Schrödinger 2025-1 with Python 3.11.4 and PyMOL 3.0.3 with Python 3.10.14.

This is a workflow-integration acceptance, not evidence of biological quality, activity, experimental validation, performance, automatic structure preparation, or compatibility with other software versions and operating systems. Downloaded/prepared structures and generated results are intentionally not committed.
