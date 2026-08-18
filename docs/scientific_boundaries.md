# Scientific Boundaries

ProteinMutFlow is workflow engineering around user-selected scientific operations.

It can:

- model explicitly defined mutations through a configured Schrödinger backend;
- locally minimize configured neighborhoods;
- measure configured PyMOL SASA, distances, VDW overlaps, and selection counts;
- verify mutation readback, isolate failures, and record versions and hashes.

It does not:

- select biologically meaningful sites;
- invent ligand, pocket, or atom selections;
- run a universal protein-preparation protocol;
- calculate DeltaDeltaG through PyMOL;
- predict activity, fitness, or experimental success;
- perform molecular dynamics;
- generate a Wet Lab shortlist or biological recommendation.

Static structures and measurements are computational evidence, not experimental validation. Backend success on one public structure does not prove compatibility with arbitrary structures, software versions, or operating systems.
