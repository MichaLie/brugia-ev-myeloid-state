# Public Data Manifest

This repository uses public GEO resources only. Cached public inputs are included for selected local reruns; other public datasets are represented by fixed derived evidence tables that support claim inspection and figure rebuilding.

## Cached public inputs

| Local file | Upstream source | Local use |
|---|---|---|
| `data/public_inputs/geo/GSE143170_count_matrix_RNAseq.txt.gz` | GEO `GSE143170` | Count matrix used by `scripts/run_gse143170_temsirolimus_sensitivity.py` and `scripts/run_gse143170_paired_deseq2.py`. |
| `data/public_inputs/geo/GSE143170_series_matrix.txt.gz` | GEO `GSE143170` | Series matrix retained for sample-level reference. |
| `data/public_inputs/geo/GSE187403_M_vs_L_Differential_Expression.txt.gz` | GEO `GSE187403` | Public differential-expression table for `LPS vs media`. |
| `data/public_inputs/geo/GSE187403_L_vs_AL_Differential_Expression.txt.gz` | GEO `GSE187403` | Public differential-expression table for `AZD2014+LPS vs LPS`. |
| `data/public_inputs/geo/GSE187403_M_vs_AL_Differential_Expression.txt.gz` | GEO `GSE187403` | Public differential-expression table for `AZD2014+LPS vs media`. |

## Dataset-level provenance

| Dataset | Analysis role | Repository representation |
|---|---|---|
| `GSE263693` | EV-exposed human monocyte derivation dataset. | Fixed signature and enrichment evidence under `data/evidence/core/`. Raw count input is not cached in this repository. |
| `GSE263690` | EV-exposed human moDC derivation dataset. | Fixed signature and enrichment evidence under `data/evidence/core/`. Raw count input is not cached in this repository. |
| `GSE107011` | Baseline myeloid-state alignment. | `data/evidence/core/baseline_myeloid_state_correlations.csv`. |
| `GSE147310` / `GSE147314` | TLR-stimulation projection across human myeloid cell types. | `data/evidence/tlr_projection/tlr78_projection_all_genes.csv` and `data/evidence/tlr_projection/tlr78_concordance_summary.csv`. |
| `GSE187403` | Primary-monocyte LPS and AZD2014 orthogonal perturbation testing. | Cached public differential-expression inputs plus outputs under `results/validation/`; rerunnable with `scripts/run_gse187403_validation.py`. |
| `GSE143170` | Paired moDC temsirolimus and acidic-pH orthogonal perturbation testing. | Cached count/series inputs plus outputs under `results/validation/` and `results/tables/`; rerunnable with `scripts/run_gse143170_temsirolimus_sensitivity.py` and `scripts/run_gse143170_paired_deseq2.py`. |
| `GSE360` | Direct human *Brugia* exposure external projection. | Fixed evidence in `data/evidence/external/external_contrast_overview.csv` and derived summaries in `results/tables/`. |
| `GSE2135` | Natural human filarial infection external projection. | Fixed evidence in `data/evidence/external/external_contrast_overview.csv` and derived summaries in `results/tables/`. |
| `GSE42694` | Same-parasite L3 Langerhans-cell exposure stress test. | `data/evidence/external/gse42694_contrast_results.csv`. |
| `GSE250463` | Cross-helminth moDC perturbation external projection. | `data/evidence/external/gse250463_contrast_results.csv`. |
| `GSE155698` | PDAC myeloid comparator-state analysis. | Fixed comparator evidence under `data/evidence/comparators/`. |
| `GSE111672` | PDAC myeloid comparator-state analysis. | Fixed comparator evidence under `data/evidence/comparators/`. |
| `GSE189847` | Canonical macrophage-polarization comparator analysis. | Fixed comparator evidence under `data/evidence/comparators/`. |

## Derived evidence bundle

- `data/evidence/core/`: fixed 76-gene signature, cross-dataset overlap, Hallmark enrichment, transcription-factor enrichment, module membership, motif tests, miRNA-target rerun summaries, and MTOR positive-control traceability.
- `data/evidence/regulon/`: SPI1/RUNX1 regulon-coverage summaries.
- `data/evidence/tlr_projection/`: fixed TLR-context gene-level and module-level projection summaries.
- `data/evidence/external/`: fixed external-projection evidence tables.
- `data/evidence/tolerance/`: tolerance-marker overlap and external-behavior summaries.
- `data/evidence/comparators/`: PDAC myeloid and macrophage-polarization comparator summaries.
- `results/validation/`: outputs from local orthogonal perturbation scripts.
- `results/tables/`: derived summaries rebuilt from the evidence and perturbation-output layers.

## Scope note

The package is not a full raw-data mirror for every GEO dataset. It provides fixed evidence for claim inspection and cached public inputs for selected perturbation reruns.

The fixed external-projection support flag required total directional concordance >= 0.60 and empirical direction-test P < 0.05. No external contrast met both criteria.
