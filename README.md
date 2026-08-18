# ProteinMutFlow

ProteinMutFlow is a deterministic, configuration-driven toolkit for batch single-site and multi-site protein mutation modeling, local minimization, structural metric extraction, and quality control. The installed command remains `mutflow`.

Status: `0.1.0-alpha.1`. This is a Windows-first alpha release validated on one recorded Schrödinger/PyMOL workstation combination.

## Why It Exists

Modeling one mutation by hand is manageable. Repeating structure loading, mutation, local minimization, measurement, naming, failure inspection, and result aggregation across dozens or hundreds of variants is slow and error-prone. ProteinMutFlow turns those repeated operations into one auditable path:

```text
workflow.yaml + prepared WT structure
  -> deterministic preflight
  -> Schrödinger mutation and local minimization
  -> headless PyMOL measurements and QC
  -> structures/ + results.csv + run.json + run.log
```

The toolkit automates choices supplied by the user. It does not select biologically meaningful sites, invent PyMOL selections, score activity, calculate DeltaDeltaG, or make Wet Lab recommendations.

## Quick Start

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install .
.venv\Scripts\mutflow init
.venv\Scripts\mutflow preflight workflow.yaml
.venv\Scripts\mutflow run workflow.yaml
```

`mutflow init` is optional. It writes one schema-valid YAML file from explicit user choices, refuses to overwrite an existing file, and invokes no scientific backend.

See [installation](docs/installation.md) and [configuration](docs/configuration.md) for the dependency and input contracts.

## Capabilities

- explicit single-site and multi-site mutation tables;
- deterministic single-site saturation expansion;
- sequential multi-site modeling from a fresh WT structure;
- configurable-radius local minimization through Schrödinger;
- configurable PyMOL SASA, minimum-distance, and VDW-overlap measurements;
- mutation readback and required-selection/atom-count QC;
- per-variant failure isolation and compact provenance;
- compatibility-gated resume without extra state files;
- a default 500-variant safety ceiling that can be deliberately raised in YAML.

## External Dependencies

The open-source core uses Python `>=3.11,<3.13`, PyYAML, and jsonschema. Scientific capabilities are enabled only when their external runtimes are available:

| Requested work | Schrödinger | PyMOL |
|---|---:|---:|
| configuration and non-scientific checks | no | no |
| mutation modeling and local minimization | yes | no |
| configured static metrics on existing PDBs | no | yes |
| full modeling-plus-metrics workflow | yes | yes |

Schrödinger remains separately installed and licensed. PyMOL is also installed separately. Neither application is bundled or licensed by this repository.

## Public Example And Validation

[`examples/public_1ubq/`](examples/public_1ubq/) defines the public WT/single/double modeling-plus-metrics example. The repository contains its mutation table, workflow, preparation instructions, checksums, and validation boundary; it excludes the downloaded/prepared PDB and generated results.

On the recorded reference workstation, WT, `A:G10A`, and `A:G10A_A:K11A` completed the installed-package Schrödinger-plus-PyMOL path with successful readback, configured metrics/QC, and overall status `OK`. See [validation evidence](docs/validation/public_1ubq_end_to_end.md).

The underlying automation approach was also exercised in the private ecHint iGEM project on 57 single mutants plus WT and 136 explicitly defined double/triple variants. These counts demonstrate workflow scale and operational validation, not experimental or predictive accuracy. Per-variant ecHint data remain private.

## Default Outputs

```text
output/
├── structures/
├── results.csv
├── run.json
└── run.log
```

PDB is the default portable structure format. MAEGZ is opt-in and requires Schrödinger. Successful default runs retain no raw mutations or stage-specific intermediate tables.

## Development

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install .
.venv\Scripts\python -m unittest discover -s tests -v
```

Continuous integration runs only the core, schema, synthetic-fixture, and mocked-backend tests. Licensed scientific jobs are documented as acceptance evidence rather than executed in public CI.

## Scientific Boundaries

Read [scientific boundaries](docs/scientific_boundaries.md) before interpreting outputs. ProteinMutFlow reports modeled structures and configured static measurements; it does not establish biological quality, activity, experimental success, or compatibility outside the recorded reference environment.

## License And Citation

ProteinMutFlow is released under the [MIT License](LICENSE). External tools and dependencies retain their own licenses; see [third-party notices](THIRD_PARTY_NOTICES.md).

Citation metadata are provided in [CITATION.cff](CITATION.cff).
