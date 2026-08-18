# Configuration

The normal path uses one `workflow.yaml` plus a WT structure. Experienced users may write YAML directly; `mutflow init` provides an optional deterministic wizard.

## Mutation Inputs

- `saturation_single`: the user supplies chain/residue sites and ProteinMutFlow expands the 19 non-WT standard amino-acid substitutions per site.
- `explicit`: a CSV supplies stable `variant_id` and `mutations` fields and may contain single-site or multi-site variants.
- Metrics-only explicit mode adds `structure_path` for each existing mutant PDB.

Mutation-site selection and combinatorial library design remain user responsibilities.

## Batch Ceiling

The normalized default is:

```yaml
execution:
  max_variants: 500
```

If deterministic expansion exceeds this value, preflight returns `REFUSED`. Raising the value in YAML is the explicit acknowledgement for a larger batch, and the normalized value is retained in `run.json`.

## Outputs

PDB is the default structure format. MAEGZ is opt-in and requires Schrödinger. A successful default run contains only `structures/`, `results.csv`, `run.json`, and `run.log`.

The authoritative machine-readable contract is [`src/mutflow/schemas/workflow.schema.json`](../src/mutflow/schemas/workflow.schema.json). Working examples are in [`examples/`](../examples/).
