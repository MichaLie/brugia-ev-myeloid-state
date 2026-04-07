# *Brugia malayi* EV-induced composite myeloid state

This repository contains the data products and analysis code for a reproducibility-first reconstruction of the host transcriptional response to *Brugia malayi* extracellular vesicles. It includes the frozen 76-gene signature, derived evidence tables, orthogonal perturbation outputs, final figure assets, and the scripts used to regenerate the main repository outputs.

## Repository layout

- `data/evidence/`
  Frozen signature tables and supporting derived evidence.
- `data/public_inputs/geo/`
  Cached public GEO files used for orthogonal reruns.
- `results/tables/`
  Derived summary tables.
- `results/validation/`
  Outputs from the orthogonal perturbation analyses.
- `results/figures/`
  Final figure assets and figure legends.
- `scripts/`
  Analysis and figure-building scripts.
- `docs/`
  Data manifest and summary documentation.
- `environment/`
  Python dependency listing.

## Getting started

1. Create a Python environment.
2. Install the dependencies listed in `environment/requirements.txt`.
3. Rebuild orthogonal perturbation outputs:
   - `python3 scripts/run_gse187403_validation.py`
   - `python3 scripts/run_gse143170_temsirolimus_sensitivity.py`
   - `python3 scripts/run_gse143170_paired_deseq2.py`
4. Rebuild derived summary tables:
   - `python3 scripts/build_manuscript_tables.py`
   - `python3 scripts/build_claim_support_summary.py`
5. Rebuild figures:
   - `python3 scripts/build_figures.py`

## Data provenance

- The files in `data/public_inputs/geo/` are cached copies of public GEO resources used in the orthogonal reruns.
- Derived outputs in `results/` are generated from repository inputs using the scripts in `scripts/`.
- Additional file-level summaries are provided in `docs/public_data_manifest.md` and `docs/derived_table_manifest.md`.
