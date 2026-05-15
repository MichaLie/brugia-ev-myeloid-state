# GSE143170 Paired DESeq2 Perturbation Summary

## Design
- Dataset: GSE143170
- Samples used: paired donor design within monocyte-derived dendritic cell differentiation
- Contrasts rerun locally: `Tems vs pH73` and `pH65 vs pH73`
- Differential framework: `pydeseq2` with design `~ donor + condition`

## Key findings
- `Tems vs pH73` reproduced 62/73 evaluable signature genes in the EV-consistent direction (0.849; one-sided binomial p = 4.545e-10).
- Module B / up-arm concordance under temsirolimus was 17/19 (0.895; p = 0.0003643).
- The acidic-pH control did not reproduce Module B: 8/19 (0.421; p = 0.8204).
- Existing repository-local score tests remain directionally aligned with the DESeq2 rerun: core EV-like score delta = 1.494 (sign-flip p = 0.03125); Module B score delta = 0.738 (p = 0.03125).
- DESeq2 and the earlier normalized-count sensitivity rerun were tightly aligned for `Tems vs pH73` (Pearson r = 0.995; Spearman rho = 0.995).

## Interpretation
- This rerun supports treating the temsirolimus result as a stable orthogonal perturbation result rather than only a sensitivity note.
- The analysis-level interpretation remains context dependence: temsirolimus during moDC differentiation reproduces Module B, whereas acidic pH does not.
- Negative-control score behavior remained non-supportive for `pH65 vs pH73` (Module B score delta = -0.084; sign-flip p = 0.31250).
