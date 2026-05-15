# GSE187403 Orthogonal Perturbation Summary

## Dataset
- GEO series: GSE187403
- Title: Transcriptional effects of catalytic mTOR inhibition among LPS-stimulated primary human monocytes
- Design: six donors; media, LPS, and AZD2014+LPS conditions

## Key findings
- LPS alone reproduced the EV down arm better than the up arm (44/56, 0.786 vs 8/20, 0.400).
- AZD2014 within LPS did not preferentially reconstruct the EV up arm (10/20, 0.500) and remained less aligned than the down arm (37/56, 0.661).
- The direct combined perturbation (AZD2014+LPS vs media) achieved 54/76 full-signature concordance (0.711).
- In leave-one-gene-out prediction, the best reconstruction model was `lps_only` with Pearson r = 0.469, Spearman rho = 0.541, and R² = 0.082.
- For comparison, direct combined contrast alone reached LOOCV Pearson r = 0.432; LPS alone = 0.469; AZD2014-within-LPS alone = 0.298.

## Sample-score behavior
- Core EV-like score: LPS vs media mean delta = 0.840 (sign-flip p = 0.031); AZD2014+LPS vs LPS = 0.211 (p = 0.344); AZD2014+LPS vs media = 1.051 (p = 0.031).
- Module B EV-like score decreased under AZD2014 within LPS (mean delta = -0.355, p = 0.031), which argues against a simple claim that acute catalytic mTOR inhibition in this monocyte context reproduces the MYC-associated arm.

## Interpretation
- This dataset is useful because it supplies the missing mTOR-inhibition side of the two-component model in primary human monocytes. The most important question is not whether it perfectly reproduces the EV state, but whether LPS and mTOR inhibition contribute different arms.
- Here, the strongest reproducible signal is still the LPS-like down arm. The mTOR-inhibitor step adds only modest full-signature support and does not independently recover Module B in this acute monocyte setting.
