# GSE143170 Temsirolimus Sensitivity Summary

## Dataset
- GEO series: GSE143170
- System: human monocyte-derived dendritic cells differentiated for 5 days in vitro
- Conditions used here: temsirolimus vs pH 7.3 control (paired by donor, n = 6)
- Note: this is a sensitivity analysis from normalized count data, not a formal DESeq2 rerun

## Key findings
- Temsirolimus showed 61/73 full-signature concordance (0.836).
- Concordance was strong for both the down arm (44/54, 0.815) and the up arm (17/19, 0.895).
- Temsirolimus increased the core EV-like score by 1.494 (sign-flip p = 0.031) and the Module B EV-like score by 0.738 (p = 0.031).

## Interpretation
- Unlike the acute AZD2014-within-LPS monocyte dataset, this moDC differentiation context shows a strong EV-like response under mTORC1 inhibition.
- This supports the interpretive caution that the MYC-associated arm is compatible with mTOR-linked control, but context dependent rather than universally reproduced by every mTOR perturbation.
