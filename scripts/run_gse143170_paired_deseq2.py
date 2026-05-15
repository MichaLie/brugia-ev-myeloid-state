#!/usr/bin/env python3

from __future__ import annotations

import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import binomtest, pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
RAW_DIR = DATA_DIR / "public_inputs" / "geo"
LEGACY_OUT_DIR = RESULTS_DIR / "validation"
OUT_DIR = RESULTS_DIR / "tables"

SIG_PATH = DATA_DIR / "evidence" / "core" / "frozen_76_gene_signature.csv"
MODULE_PATH = DATA_DIR / "evidence" / "core" / "module_membership.csv"
COUNT_PATH = RAW_DIR / "GSE143170_count_matrix_RNAseq.txt.gz"

warnings.filterwarnings("ignore", category=RuntimeWarning)


def safe_corr(func, x: pd.Series, y: pd.Series) -> float:
    valid = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(valid) < 3:
        return math.nan
    if valid["x"].nunique() < 2 or valid["y"].nunique() < 2:
        return math.nan
    return float(func(valid["x"], valid["y"])[0])


def short_module_name(module: str) -> str:
    if module.startswith("A_"):
        return "module_a"
    if module.startswith("B_"):
        return "module_b"
    if module.startswith("C_"):
        return "module_c"
    return module.lower()


def load_signature() -> pd.DataFrame:
    sig = pd.read_csv(SIG_PATH)
    modules = pd.read_csv(MODULE_PATH)[["gene", "module"]]
    sig = sig.merge(modules, on="gene", how="left")
    sig["ev_mean_log2fc"] = sig[["log2FoldChange_monocytes", "log2FoldChange_moDCs"]].mean(axis=1)
    sig["arm"] = np.where(sig["direction"].eq("up"), "up_arm", "down_arm")
    sig["module_group"] = sig["module"].map(short_module_name)
    return sig[["gene", "direction", "arm", "module_group", "ev_mean_log2fc"]].copy()


def load_counts() -> pd.DataFrame:
    return pd.read_csv(COUNT_PATH, sep="\t", index_col=0, compression="gzip")


def run_deseq_contrast(
    counts: pd.DataFrame,
    cond_a: str,
    cond_b: str,
    contrast_name: str,
) -> pd.DataFrame:
    sample_prefixes = (f"{cond_a}_", f"{cond_b}_")
    samples = [col for col in counts.columns if col.startswith(sample_prefixes)]
    subset = counts[samples].copy()
    subset = subset.loc[subset.sum(axis=1) >= 10].copy()

    metadata = pd.DataFrame(
        {
            "condition": [cond_a if sample.startswith(f"{cond_a}_") else cond_b for sample in samples],
            "donor": [sample.split("_")[1] for sample in samples],
        },
        index=samples,
    )

    dds = DeseqDataSet(
        counts=subset.T.astype(int),
        metadata=metadata,
        design="~ donor + condition",
        quiet=True,
        min_replicates=2,
        refit_cooks=False,
        n_cpus=1,
    )
    dds.deseq2()

    stats = DeseqStats(dds, contrast=["condition", cond_b, cond_a], quiet=True, n_cpus=1)
    stats.summary()
    result = stats.results_df.reset_index().rename(columns={"index": "gene"})
    return result.rename(
        columns={
            "baseMean": f"{contrast_name}_base_mean",
            "log2FoldChange": f"{contrast_name}_log2fc",
            "lfcSE": f"{contrast_name}_lfc_se",
            "stat": f"{contrast_name}_stat",
            "pvalue": f"{contrast_name}_pvalue",
            "padj": f"{contrast_name}_padj",
        }
    )


def summarize_concordance(merged: pd.DataFrame, contrast_name: str) -> pd.DataFrame:
    lfc_col = f"{contrast_name}_log2fc"
    padj_col = f"{contrast_name}_padj"
    merged = merged.copy()
    merged[f"{contrast_name}_concordant"] = np.sign(merged[lfc_col]) == np.sign(merged["ev_mean_log2fc"])

    groups = {
        "full_signature": merged[lfc_col].notna(),
        "down_arm": merged[lfc_col].notna() & merged["arm"].eq("down_arm"),
        "up_arm": merged[lfc_col].notna() & merged["arm"].eq("up_arm"),
        "module_a": merged[lfc_col].notna() & merged["module_group"].eq("module_a"),
        "module_b": merged[lfc_col].notna() & merged["module_group"].eq("module_b"),
        "module_c": merged[lfc_col].notna() & merged["module_group"].eq("module_c"),
    }

    rows = []
    for group_name, mask in groups.items():
        subset = merged.loc[mask].copy()
        n = len(subset)
        n_concordant = int(subset[f"{contrast_name}_concordant"].sum()) if n else 0
        rows.append(
            {
                "contrast": contrast_name,
                "group": group_name,
                "n_available": n,
                "n_concordant": n_concordant,
                "concordance_fraction": float(n_concordant / n) if n else math.nan,
                "binom_p_greater": float(binomtest(n_concordant, n, 0.5, alternative="greater").pvalue) if n else math.nan,
                "pearson_r": safe_corr(pearsonr, subset["ev_mean_log2fc"], subset[lfc_col]),
                "spearman_rho": safe_corr(spearmanr, subset["ev_mean_log2fc"], subset[lfc_col]),
                "median_abs_lfc": float(subset[lfc_col].abs().median()) if n else math.nan,
                "median_padj": float(subset[padj_col].median()) if n else math.nan,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    signature = load_signature()
    counts = load_counts()

    tems = run_deseq_contrast(counts, cond_a="pH73", cond_b="Tems", contrast_name="tems_vs_ph73")
    ph65 = run_deseq_contrast(counts, cond_a="pH73", cond_b="pH65", contrast_name="ph65_vs_ph73")

    merged = signature.merge(
        tems[
            [
                "gene",
                "tems_vs_ph73_base_mean",
                "tems_vs_ph73_log2fc",
                "tems_vs_ph73_lfc_se",
                "tems_vs_ph73_stat",
                "tems_vs_ph73_pvalue",
                "tems_vs_ph73_padj",
            ]
        ],
        on="gene",
        how="left",
    ).merge(
        ph65[
            [
                "gene",
                "ph65_vs_ph73_base_mean",
                "ph65_vs_ph73_log2fc",
                "ph65_vs_ph73_lfc_se",
                "ph65_vs_ph73_stat",
                "ph65_vs_ph73_pvalue",
                "ph65_vs_ph73_padj",
            ]
        ],
        on="gene",
        how="left",
    )

    for contrast_name in ("tems_vs_ph73", "ph65_vs_ph73"):
        merged[f"{contrast_name}_concordant"] = np.sign(merged[f"{contrast_name}_log2fc"]) == np.sign(
            merged["ev_mean_log2fc"]
        )

    summary = pd.concat(
        [
            summarize_concordance(merged, "tems_vs_ph73"),
            summarize_concordance(merged, "ph65_vs_ph73"),
        ],
        ignore_index=True,
    )

    legacy_gene = pd.read_csv(LEGACY_OUT_DIR / "gse143170_gene_level_sensitivity.csv")
    crosscheck_rows = []
    for contrast_name, legacy_col in [
        ("tems_vs_ph73", "tems_vs_ph73_log2fc"),
        ("ph65_vs_ph73", "ph65_vs_ph73_log2fc"),
    ]:
        comp = merged[["gene", f"{contrast_name}_log2fc"]].merge(
            legacy_gene[["gene", legacy_col]].rename(columns={legacy_col: f"legacy_{legacy_col}"}),
            on="gene",
            how="inner",
        )
        crosscheck_rows.append(
            {
                "contrast": contrast_name,
                "n_shared_genes": int(comp.dropna().shape[0]),
                "pearson_r": safe_corr(
                    pearsonr,
                    comp[f"{contrast_name}_log2fc"],
                    comp[f"legacy_{legacy_col}"],
                ),
                "spearman_rho": safe_corr(
                    spearmanr,
                    comp[f"{contrast_name}_log2fc"],
                    comp[f"legacy_{legacy_col}"],
                ),
            }
        )
    crosscheck = pd.DataFrame(crosscheck_rows)

    paired = pd.read_csv(LEGACY_OUT_DIR / "gse143170_paired_score_tests.csv")

    merged.sort_values(["direction", "module_group", "gene"]).to_csv(
        OUT_DIR / "gse143170_paired_deseq2_gene_level.csv",
        index=False,
    )
    summary.to_csv(OUT_DIR / "gse143170_paired_deseq2_concordance_summary.csv", index=False)
    crosscheck.to_csv(OUT_DIR / "gse143170_paired_deseq2_crosscheck.csv", index=False)

    tems_full = summary[(summary["contrast"].eq("tems_vs_ph73")) & (summary["group"].eq("full_signature"))].iloc[0]
    tems_modb = summary[(summary["contrast"].eq("tems_vs_ph73")) & (summary["group"].eq("module_b"))].iloc[0]
    ph65_modb = summary[(summary["contrast"].eq("ph65_vs_ph73")) & (summary["group"].eq("module_b"))].iloc[0]
    tems_core = paired[(paired["contrast"].eq("tems_vs_ph73")) & (paired["score_name"].eq("core_ev_score"))].iloc[0]
    tems_modb_score = paired[(paired["contrast"].eq("tems_vs_ph73")) & (paired["score_name"].eq("module_b_ev_score"))].iloc[0]
    ph65_modb_score = paired[(paired["contrast"].eq("ph65_vs_ph73")) & (paired["score_name"].eq("module_b_ev_score"))].iloc[0]
    cross = crosscheck[crosscheck["contrast"].eq("tems_vs_ph73")].iloc[0]

    lines = [
        "# GSE143170 Paired DESeq2 Perturbation Summary",
        "",
        "## Design",
        "- Dataset: GSE143170",
        "- Samples used: paired donor design within monocyte-derived dendritic cell differentiation",
        "- Contrasts rerun locally: `Tems vs pH73` and `pH65 vs pH73`",
        "- Differential framework: `pydeseq2` with design `~ donor + condition`",
        "",
        "## Key findings",
        (
            f"- `Tems vs pH73` reproduced {int(tems_full['n_concordant'])}/{int(tems_full['n_available'])} "
            f"evaluable signature genes in the EV-consistent direction ({tems_full['concordance_fraction']:.3f}; "
            f"one-sided binomial p = {tems_full['binom_p_greater']:.4g})."
        ),
        (
            f"- Module B / up-arm concordance under temsirolimus was "
            f"{int(tems_modb['n_concordant'])}/{int(tems_modb['n_available'])} "
            f"({tems_modb['concordance_fraction']:.3f}; p = {tems_modb['binom_p_greater']:.4g})."
        ),
        (
            f"- The acidic-pH control did not reproduce Module B: "
            f"{int(ph65_modb['n_concordant'])}/{int(ph65_modb['n_available'])} "
            f"({ph65_modb['concordance_fraction']:.3f}; p = {ph65_modb['binom_p_greater']:.4g})."
        ),
        (
            f"- Existing repository-local score tests remain directionally aligned with the DESeq2 rerun: "
            f"core EV-like score delta = {tems_core['mean_delta']:.3f} (sign-flip p = {tems_core['signflip_p']:.5f}); "
            f"Module B score delta = {tems_modb_score['mean_delta']:.3f} "
            f"(p = {tems_modb_score['signflip_p']:.5f})."
        ),
        (
            f"- DESeq2 and the earlier normalized-count sensitivity rerun were tightly aligned for `Tems vs pH73` "
            f"(Pearson r = {cross['pearson_r']:.3f}; Spearman rho = {cross['spearman_rho']:.3f})."
        ),
        "",
        "## Interpretation",
        "- This rerun supports treating the temsirolimus result as a stable orthogonal perturbation result rather than only a sensitivity note.",
        "- The analysis-level interpretation remains context dependence: temsirolimus during moDC differentiation reproduces Module B, whereas acidic pH does not.",
        (
            f"- Negative-control score behavior remained non-supportive for `pH65 vs pH73` "
            f"(Module B score delta = {ph65_modb_score['mean_delta']:.3f}; "
            f"sign-flip p = {ph65_modb_score['signflip_p']:.5f})."
        ),
    ]

    (OUT_DIR / "gse143170_paired_deseq2_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
