# Schrödinger Adapter Smoke

Status: `ACCEPTED` with a prepared-input precondition.

On 2026-08-18, the installed adapter processed public prepared 1UBQ as WT, `A:G10A`, and `A:G10A_A:K11A` through the ordinary `mutflow run` path on Schrödinger Suite 2025-1. All three final PDBs were readable; both mutants passed residue-identity readback and local-minimization boundary checks.

An unprepared 1UBQ negative control preserved real atom-typing failures instead of labeling them successful. This established the explicit requirement that modeling input be chemically prepared for the requested force-field operation.

The test establishes adapter integration on the recorded workstation. It does not establish automatic protein preparation, biological quality, performance, or broad compatibility.
