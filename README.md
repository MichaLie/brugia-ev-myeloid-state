# Brugia EV Myeloid State

This repository is a clean, journal-agnostic reproducibility package for the reconstructed *Brugia malayi* extracellular-vesicle-induced composite myeloid state. It contains the frozen 76-gene signature, derived evidence tables, cached public GEO inputs used for orthogonal reruns, figure assets, and the scripts used to rebuild the main derived outputs.

## Repository scope

- `data/evidence/`: frozen signature tables and derived evidence used across the repository.
- `data/public_inputs/geo/`: public GEO input files cached locally for reruns.
- `results/tables/`: derived summary tables.
- `results/validation/`: orthogonal perturbation outputs derived from public GEO inputs.
- `results/figures/`: final figure assets and legends.
- `scripts/`: analysis and figure-building scripts rewritten to run from this repository layout.
- `docs/`: repository notes, data manifest, and derived-table summaries.
- `environment/`: lightweight dependency listing.

## Intentional exclusions

- No manuscript draft files.
- No internal review notes, revision logs, or editorial checklists.
- No copyrighted source-paper PDFs.
- No journal-specific submission forms or cover-letter material.

## Quick start

1. Create a Python environment and install `environment/requirements.txt`.
2. Rebuild orthogonal perturbation outputs as needed with:
   - `python3 scripts/run_gse187403_validation.py`
   - `python3 scripts/run_gse143170_temsirolimus_sensitivity.py`
   - `python3 scripts/run_gse143170_paired_deseq2.py`
3. Rebuild derived summary tables with:
   - `python3 scripts/build_manuscript_tables.py`
   - `python3 scripts/build_claim_support_summary.py`
4. Rebuild figures with:
   - `python3 scripts/build_figures.py`

## Notes

- The GEO input files in `data/public_inputs/geo/` are public-source inputs mirrored locally for convenience.
- Repository result files intentionally use descriptive names rather than journal-specific supplementary-table numbering so the same outputs remain stable across manuscripts and resubmissions.
- Before an archival DOI release, choose an explicit repository license and archive a tagged release to a DOI-minting service such as Zenodo.
