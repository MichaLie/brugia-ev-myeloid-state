# Figure Source Manifest

## Build

- Build script: `scripts/build_figures.py`
- Render command: `python3 scripts/build_figures.py`
- Output formats: publication PDF and 600 dpi PNG for all main and supplemental figures

## Main figures

| Figure | Output directory | Primary source files |
|---|---|---|
| Figure 1 | `results/figures/main/` | `data/evidence/core/frozen_signature_cross_dataset_overlap.csv`; `data/evidence/core/frozen_76_gene_signature.csv`; `data/evidence/core/frozen_signature_hallmark_enrichment.csv` |
| Figure 2 | `results/figures/main/` | `data/evidence/core/preregistered_tf_enrichment_primary.csv`; `data/evidence/regulon/regulon_coverage_summary.csv`; `data/evidence/core/module_membership.csv`; `data/evidence/core/baseline_myeloid_state_correlations.csv` |
| Figure 3 | `results/figures/main/` | `data/evidence/core/regulator_panel_target_rerun_results.csv`; `data/evidence/core/frozen_signature_6mer_enrichment.csv`; `data/evidence/core/frozen_signature_7mer_enrichment.csv` |
| Figure 4 | `results/figures/main/` | `results/tables/orthogonal_validation_summary.csv`; `results/tables/external_transport_context_summary.csv`; `results/tables/gse143170_paired_deseq2_concordance_summary.csv`; `results/validation/gse187403_paired_score_tests.csv`; `results/validation/gse143170_paired_score_tests.csv` |
| Figure 5 | `results/figures/main/` | `data/evidence/external/external_arm_transport_contrast_summary.csv`; `results/tables/external_transport_context_summary.csv` |

## Supplemental figures

| Figure | Output directory | Primary source files |
|---|---|---|
| Figure S1 | `results/figures/supplementary/` | `data/evidence/tlr_projection/tlr78_projection_all_genes.csv`; `data/evidence/tlr_projection/tlr78_concordance_summary.csv` |
| Figure S2 | `results/figures/supplementary/` | `data/evidence/core/preregistered_tf_enrichment_primary.csv`; `data/evidence/core/module_gene_level_concordance.csv` |
| Figure S3 | `results/figures/supplementary/` | `results/validation/gse187403_gene_level_validation.csv`; `results/validation/gse143170_gene_level_sensitivity.csv` |
| Figure S4 | `results/figures/supplementary/` | `data/evidence/external/external_arm_transport_contrast_summary.csv`; `data/evidence/external/external_standardized_module_contrast_summary.csv` |
| Figure S5 | `results/figures/supplementary/` | `data/evidence/tolerance/tolerance_signature_overlap.csv`; `data/evidence/tolerance/tolerance_markers_in_external.csv`; `data/evidence/comparators/pdac_cross_cohort_state_geometry.csv`; `data/evidence/comparators/macrophage_core_by_polarization.csv` |
