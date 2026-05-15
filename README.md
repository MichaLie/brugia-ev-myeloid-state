# *Brugia malayi* EV-induced composite myeloid state

This repository contains the reproducibility package for a transcriptomic reconstruction of the host response to *Brugia malayi* extracellular vesicles in human myeloid cells. It includes the fixed 76-gene Core EV-response signature, derived evidence tables, selected public-input files for local reruns, orthogonal perturbation outputs, final figure assets, and scripts for regenerating publication-ready repository outputs.

The package is intentionally evidence-bounded. Some analyses are rerunnable from cached public input files in this repository, whereas other layers are represented as fixed derived evidence tables produced from public datasets and preserved here for inspection, figure regeneration, and claim verification. The repository does not contain drafts, publication-package files, literature PDFs, or non-public data.

## Repository layout

- `data/evidence/`
  Fixed signature tables and supporting derived evidence.
- `data/public_inputs/geo/`
  Cached public GEO files used for selected local reruns.
- `results/tables/`
  Derived summary tables.
- `results/validation/`
  Outputs from orthogonal perturbation analyses.
- `results/figures/`
  Final figure assets and figure legends.
- `scripts/`
  Analysis and figure-building scripts.
- `docs/`
  Data manifest and summary documentation.
- `environment/`
  Python dependency listing and local verification version record.

## What can be rerun directly

1. Create a Python environment and install the dependencies listed in `environment/requirements.txt`.
2. Rebuild selected orthogonal perturbation outputs from cached public inputs:
   - `python3 scripts/run_gse187403_validation.py`
   - `python3 scripts/run_gse143170_temsirolimus_sensitivity.py`
   - `python3 scripts/run_gse143170_paired_deseq2.py`
3. Rebuild derived summary tables:
   - `python3 scripts/build_summary_tables.py`
   - `python3 scripts/build_claim_support_summary.py`
4. Rebuild the figure assets from the fixed evidence tables and regenerated summary outputs:
   - `python3 scripts/build_figures.py`

The derivation signature, transcription-factor enrichment, miRNA-target/motif tests, TLR projection, external projection, and comparator layers are provided as fixed evidence tables under `data/evidence/`. Their dataset provenance and repository representation are documented in `docs/public_data_manifest.md`.

The local package versions used for verification are recorded in `environment/software_versions_verified.txt`.

## Data provenance

- All analyzed datasets are public GEO resources; no new human participant data or non-public datasets are included.
- The files in `data/public_inputs/geo/` are cached copies of public GEO resources used in selected local reruns.
- Derived evidence tables under `data/evidence/` support the fixed signature, regulatory analyses, miRNA-target/motif tests, TLR projection, external projection, tolerance comparison, and comparator-state analyses.
- Derived outputs in `results/` can be regenerated from repository inputs using the scripts in `scripts/`.
- File-level summaries are provided in `docs/public_data_manifest.md`, `docs/derived_table_manifest.md`, and `docs/module_definitions.md`.

## Citation and release

When citing this repository, include the exact Git commit hash used for the analysis package. Before archival release, create a GitHub release and archive it with a persistent DOI provider such as Zenodo or Figshare. A release checklist is provided in `docs/release_checklist.md`.

## License

Code is licensed under the MIT License. Documentation, figure assets, legends, manifests, and derived data tables are licensed under CC BY 4.0. See `LICENSE.md`.
