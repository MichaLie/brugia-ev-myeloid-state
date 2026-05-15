# Derived Table Manifest

## Included tables
- `claim_support_stats.csv`
  Suggested caption: Claim-facing statistical anchors for the main text and extended materials.
- `orthogonal_validation_summary.csv`
  Suggested caption: Summary of TLR-context projection and orthogonal perturbation testing across GSE147310, GSE187403, and GSE143170.
- `gse143170_paired_deseq2_concordance_summary.csv`
  Suggested caption: Paired DESeq2 concordance summary for GSE143170 temsirolimus and acidic-pH contrasts.
- `gse143170_paired_deseq2_crosscheck.csv`
  Suggested caption: Crosscheck between the paired DESeq2 rerun and the normalized-count sensitivity rerun for GSE143170.
- `gse143170_paired_deseq2_gene_level.csv`
  Suggested caption: Gene-level paired DESeq2 output for fixed-signature genes in GSE143170.
- `external_transport_contrast_summary.csv`
  Suggested caption: Contrast-level external projection summary across direct Brugia exposure, natural infection, same-parasite stress-test, and cross-helminth datasets.
- `external_transport_context_summary.csv`
  Suggested caption: Context-level external-projection summary showing differential down-arm and up-arm behavior across external settings.
- `comparator_state_summary.csv`
  Suggested caption: Comparator-state summary for PDAC myeloid cohorts and canonical macrophage polarization, emphasizing pathway-level neighborhood rather than shared mechanism.

## External support rule
- The preserved `overall_supportive` flag for external contrasts requires total directional concordance >= 0.60 and empirical direction-test P < 0.05.
- The context-level `n_supportive_contrasts` field counts contrasts satisfying that fixed rule.

## Associated summaries
- `gse143170_paired_deseq2_summary.md`
  Suggested use: human-readable summary of the paired DESeq2 rerun.

## Repository use
- These tables are derived summaries built from the repository evidence and perturbation-output layers.
- They are suitable for supplements, appendices, or public repository browsing without publication-package files.
