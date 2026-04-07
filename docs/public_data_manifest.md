# Public Data Manifest

## Cached public inputs

| Local file | Upstream source | Notes |
|---|---|---|
| `data/public_inputs/geo/GSE143170_count_matrix_RNAseq.txt.gz` | GEO `GSE143170` | Count matrix used for the paired moDC reruns. |
| `data/public_inputs/geo/GSE143170_series_matrix.txt.gz` | GEO `GSE143170` | Series matrix retained for sample-level reference. |
| `data/public_inputs/geo/GSE187403_M_vs_L_Differential_Expression.txt.gz` | GEO `GSE187403` | Differential-expression table for `LPS vs media`. |
| `data/public_inputs/geo/GSE187403_L_vs_AL_Differential_Expression.txt.gz` | GEO `GSE187403` | Differential-expression table for `AZD2014+LPS vs LPS`. |
| `data/public_inputs/geo/GSE187403_M_vs_AL_Differential_Expression.txt.gz` | GEO `GSE187403` | Differential-expression table for `AZD2014+LPS vs media`. |

## Derived evidence bundle

- `data/evidence/core/frozen_76_gene_signature.csv`: frozen 76-gene signature used throughout the repository.
- `data/evidence/tlr_projection/`: TLR-context projection results derived from `GSE147310`.
- `results/validation/`: orthogonal perturbation outputs derived from `GSE187403` and `GSE143170`.
- `results/tables/`: derived summary tables built from the evidence and orthogonal-output layers.
