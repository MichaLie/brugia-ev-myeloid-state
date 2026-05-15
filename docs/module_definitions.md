# Module Definitions

The deterministic module labels were assigned after the fixed 76-gene Core EV-response signature had been defined. They were not used to select genes for the signature.

## Primary split

- The fixed signature contains 56 down-regulated genes and 20 up-regulated genes.
- All 20 up-regulated genes were assigned to Module B, the MYC-associated proliferative-metabolic arm.
- The 56 down-regulated genes were split secondarily into Modules A and C.

## Down-arm split

The down-arm split used a prespecified external-projection rule based only on available construction-stage resting contexts:

- Module A: down-regulated genes with mean resting-context concordance of at least 2/3 across available construction-stage resting contexts.
- Module C: remaining down-regulated genes with mean resting-context concordance below 2/3 across available construction-stage resting contexts.

The `resting_contexts_available` column records the denominator used for each gene because not every gene was evaluable in all construction contexts.

This yielded:

| Module | Count | Direction | Assignment rule |
|---|---:|---|---|
| Module A, myeloid-lineage state | 34 | Down | Mean resting-context concordance >= 2/3 across available construction-only resting external contexts. |
| Module B, proliferative-metabolic | 20 | Up | All up-regulated genes in the fixed signature. |
| Module C, context-sensitive residual | 22 | Down | Mean resting-context concordance < 2/3 across available construction-only resting external contexts. |

## Interpretation guardrails

- Module A and Module C are descriptive substructure within the down arm.
- Module B is the up arm and is interpreted primarily as the MYC-associated proliferative-metabolic arm.
- The module labels support arm-level interpretation, figure organization, and external-projection summaries.
- The labels are not independent discovery claims and should not be described as causal mechanisms without perturbation experiments.

The row-level assignments are in `data/evidence/core/module_membership.csv`.
