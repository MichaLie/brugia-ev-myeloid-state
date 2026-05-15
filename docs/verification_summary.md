# Verification Summary

Verification date: 2026-05-15

## Reruns Completed

- `python3 scripts/run_gse187403_validation.py`
- `python3 scripts/run_gse143170_temsirolimus_sensitivity.py`
- `python3 scripts/run_gse143170_paired_deseq2.py`
- `python3 scripts/build_summary_tables.py`
- `python3 scripts/build_claim_support_summary.py`
- `python3 scripts/build_figures.py`
- `python3 -m py_compile scripts/*.py`

## Integrity Checks

- Required package files were present after the rerun.
- The fixed signature contains 76 genes: 56 down-regulated and 20 up-regulated.
- Deterministic module counts match the analysis framing: Module A = 34, Module B = 20, Module C = 22.
- Module-assignment row text refers to available resting construction contexts, and `resting_contexts_available` records the per-gene denominator.
- The preserved 17-gene regulator panel remains negative for direct hsa-miR-100-5p hits in the rerun resources.
- The MTOR positive control is represented as `data/evidence/core/mir100_mtor_positive_control.csv` and is read by the Figure 3 build.
- `docs/public_data_manifest.md` lists all GEO accessions used by the analysis.
- External-projection evidence tables record the fixed support rule: total directional concordance >= 0.60 and empirical direction-test P < 0.05.
- `environment/requirements.txt` includes `gseapy` and `scikit-learn`; local verification versions are recorded in `environment/software_versions_verified.txt`.
- The corrected `GSE143170` paired DESeq2 script writes the downstream-expected `gse143170_paired_deseq2_*` outputs.
- Figure S3 uses the paired DESeq2 GSE143170 gene-level table.
- Main and supplemental figure PNGs passed a nonblank image check after rebuild.

## Notes

The figure build emitted matplotlib `constrained_layout` warnings during export, but the build completed and checked PNG outputs had valid dimensions and nonblank pixel variation.

`LICENSE.md` records MIT licensing for code and CC BY 4.0 licensing for documentation, figure assets, legends, manifests, and derived data tables.
