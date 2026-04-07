#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
RAW_DIR = DATA_DIR / "public_inputs" / "geo"
OUT_DIR = RESULTS_DIR / "validation"

SIG_PATH = DATA_DIR / "evidence" / "core" / "frozen_76_gene_signature.csv"
MODULE_PATH = DATA_DIR / "evidence" / "core" / "module_membership.csv"

CONTRAST_PATHS = {
    "lps_vs_media": RAW_DIR / "GSE187403_M_vs_L_Differential_Expression.txt.gz",
    "mtori_lps_vs_lps": RAW_DIR / "GSE187403_L_vs_AL_Differential_Expression.txt.gz",
    "mtori_lps_vs_media": RAW_DIR / "GSE187403_M_vs_AL_Differential_Expression.txt.gz",
}

CONTRAST_LABELS = {
    "lps_vs_media": "LPS vs media",
    "mtori_lps_vs_lps": "AZD2014+LPS vs LPS",
    "mtori_lps_vs_media": "AZD2014+LPS vs media",
}

SCORE_CONTRASTS = [
    ("M", "L", "lps_vs_media"),
    ("L", "AL", "mtori_lps_vs_lps"),
    ("M", "AL", "mtori_lps_vs_media"),
]

ANCHOR_GENES = [
    "MYC",
    "MAX",
    "SPI1",
    "RUNX1",
    "MTOR",
    "IRAK1",
    "SLC7A11",
    "NCL",
    "LARS",
    "NOP58",
    "CSF3R",
    "PLD3",
    "PLA2G15",
    "METTL7B",
]


@dataclass
class ContrastData:
    name: str
    label: str
    gene_table: pd.DataFrame
    expr_table: pd.DataFrame
    sample_columns: list[str]


def safe_corr(func, x: pd.Series, y: pd.Series) -> float:
    valid = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(valid) < 3:
        return math.nan
    if valid["x"].nunique() < 2 or valid["y"].nunique() < 2:
        return math.nan
    return float(func(valid["x"], valid["y"])[0])


def cosine_similarity(x: pd.Series, y: pd.Series) -> float:
    valid = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(valid) < 3:
        return math.nan
    xvals = valid["x"].to_numpy(dtype=float)
    yvals = valid["y"].to_numpy(dtype=float)
    denom = np.linalg.norm(xvals) * np.linalg.norm(yvals)
    if denom == 0:
        return math.nan
    return float(np.dot(xvals, yvals) / denom)


def exact_paired_signflip_p(diffs: pd.Series) -> float:
    diffs = diffs.dropna().astype(float)
    if diffs.empty:
        return math.nan
    arr = diffs.to_numpy()
    observed = abs(arr.mean())
    total = 0
    extreme = 0
    for mask in range(1 << len(arr)):
        signs = np.array([1 if (mask >> i) & 1 else -1 for i in range(len(arr))], dtype=float)
        stat = abs((arr * signs).mean())
        total += 1
        if stat >= observed - 1e-12:
            extreme += 1
    return float(extreme / total)


def fit_no_intercept(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return beta, x @ beta


def loocv_predictions(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    preds = np.full_like(y, fill_value=np.nan, dtype=float)
    for idx in range(len(y)):
        train_mask = np.ones(len(y), dtype=bool)
        train_mask[idx] = False
        beta, _ = fit_no_intercept(x[train_mask], y[train_mask])
        preds[idx] = float(x[idx] @ beta)
    return preds


def regression_metrics(y_true: pd.Series, y_pred: pd.Series) -> dict[str, float]:
    valid = pd.DataFrame({"y_true": y_true, "y_pred": y_pred}).dropna()
    if len(valid) < 3:
        return {
            "pearson_r": math.nan,
            "spearman_rho": math.nan,
            "cosine_similarity": math.nan,
            "r2": math.nan,
            "mae": math.nan,
        }
    sse = float(((valid["y_true"] - valid["y_pred"]) ** 2).sum())
    sst = float(((valid["y_true"] - valid["y_true"].mean()) ** 2).sum())
    return {
        "pearson_r": safe_corr(pearsonr, valid["y_true"], valid["y_pred"]),
        "spearman_rho": safe_corr(spearmanr, valid["y_true"], valid["y_pred"]),
        "cosine_similarity": cosine_similarity(valid["y_true"], valid["y_pred"]),
        "r2": float(1 - (sse / sst)) if sst > 0 else math.nan,
        "mae": float((valid["y_true"] - valid["y_pred"]).abs().mean()),
    }


def ensure_dirs() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)


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
    modules = pd.read_csv(MODULE_PATH)[["gene", "module"]].copy()
    sig = sig.merge(modules, on="gene", how="left")
    sig["ev_mean_log2fc"] = sig[["log2FoldChange_monocytes", "log2FoldChange_moDCs"]].mean(axis=1)
    sig["module_group"] = sig["module"].map(short_module_name)
    sig["arm"] = np.where(sig["direction"] == "up", "up_arm", "down_arm")
    return sig[
        [
            "gene",
            "direction",
            "arm",
            "module",
            "module_group",
            "log2FoldChange_monocytes",
            "log2FoldChange_moDCs",
            "ev_mean_log2fc",
        ]
    ].copy()


def load_contrast(name: str, label: str, path: Path) -> ContrastData:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    df = df.rename(columns={"Gene Symbol": "gene", "Log2 Fold Change": "log2fc", "FDR Adj p Value": "fdr"})
    df["gene"] = df["gene"].astype(str).str.strip()
    df = df[df["gene"].ne("")].copy()

    sample_columns = [col for col in df.columns if re.match(r"^S\d+_(M|L|AL)_", str(col))]
    numeric_cols = ["log2fc", "fdr", "Base Mean", "M", "L", "AL", *sample_columns]
    present_numeric = [col for col in numeric_cols if col in df.columns]

    grouped = df.groupby("gene", as_index=False)[present_numeric].mean()
    gene_table = grouped.rename(
        columns={
            "log2fc": f"{name}_log2fc",
            "fdr": f"{name}_fdr",
            "Base Mean": f"{name}_base_mean",
        }
    )
    expr_table = grouped.set_index("gene")[sample_columns].sort_index()
    return ContrastData(name=name, label=label, gene_table=gene_table, expr_table=expr_table, sample_columns=sample_columns)


def combine_expression_tables(contrasts: list[ContrastData]) -> pd.DataFrame:
    merged = pd.concat([contrast.expr_table for contrast in contrasts], axis=1)
    merged = merged.loc[:, ~merged.columns.duplicated()].copy()
    merged = merged.groupby(level=0).mean()
    return merged


def summarize_concordance(signature: pd.DataFrame, contrast: ContrastData) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = signature.merge(contrast.gene_table[["gene", f"{contrast.name}_log2fc", f"{contrast.name}_fdr"]], on="gene", how="left")
    lfc_col = f"{contrast.name}_log2fc"
    fdr_col = f"{contrast.name}_fdr"
    merged["available"] = merged[lfc_col].notna()
    merged["concordant"] = np.sign(merged[lfc_col]) == np.sign(merged["ev_mean_log2fc"])

    group_defs = {
        "full_signature": merged["available"],
        "down_arm": merged["available"] & merged["arm"].eq("down_arm"),
        "up_arm": merged["available"] & merged["arm"].eq("up_arm"),
        "module_a": merged["available"] & merged["module_group"].eq("module_a"),
        "module_b": merged["available"] & merged["module_group"].eq("module_b"),
        "module_c": merged["available"] & merged["module_group"].eq("module_c"),
    }

    rows = []
    for group_name, mask in group_defs.items():
        subset = merged.loc[mask].copy()
        n = len(subset)
        concordant = int(subset["concordant"].sum()) if n else 0
        rows.append(
            {
                "contrast": contrast.name,
                "contrast_label": contrast.label,
                "group": group_name,
                "n_available": n,
                "n_concordant": concordant,
                "concordance_fraction": float(concordant / n) if n else math.nan,
                "pearson_r": safe_corr(pearsonr, subset["ev_mean_log2fc"], subset[lfc_col]),
                "spearman_rho": safe_corr(spearmanr, subset["ev_mean_log2fc"], subset[lfc_col]),
                "cosine_similarity": cosine_similarity(subset["ev_mean_log2fc"], subset[lfc_col]),
                "median_abs_lfc": float(subset[lfc_col].abs().median()) if n else math.nan,
                "median_fdr": float(subset[fdr_col].median()) if n else math.nan,
            }
        )
    return pd.DataFrame(rows), merged


def build_sample_metadata(expression: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample in expression.columns:
        match = re.match(r"^(S\d+)_(M|L|AL)_", sample)
        if not match:
            continue
        rows.append({"sample": sample, "donor": match.group(1), "condition": match.group(2)})
    return pd.DataFrame(rows).sort_values(["donor", "condition"])


def compute_sample_scores(signature: pd.DataFrame, expression: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = [gene for gene in signature["gene"] if gene in expression.index]
    expr = expression.loc[genes].copy()
    gene_mean = expr.mean(axis=1)
    gene_std = expr.std(axis=1, ddof=0).replace(0, np.nan)
    zexpr = expr.sub(gene_mean, axis=0).div(gene_std, axis=0)

    gene_sets = {
        "core_ev_score": (
            signature.loc[signature["arm"].eq("up_arm"), "gene"].tolist(),
            signature.loc[signature["arm"].eq("down_arm"), "gene"].tolist(),
        ),
        "up_arm_ev_score": (signature.loc[signature["arm"].eq("up_arm"), "gene"].tolist(), []),
        "down_arm_ev_score": ([], signature.loc[signature["arm"].eq("down_arm"), "gene"].tolist()),
        "module_a_ev_score": ([], signature.loc[signature["module_group"].eq("module_a"), "gene"].tolist()),
        "module_b_ev_score": (signature.loc[signature["module_group"].eq("module_b"), "gene"].tolist(), []),
        "module_c_ev_score": ([], signature.loc[signature["module_group"].eq("module_c"), "gene"].tolist()),
    }

    scores = {}
    for score_name, (positive_genes, negative_genes) in gene_sets.items():
        pos_genes = [gene for gene in positive_genes if gene in zexpr.index]
        neg_genes = [gene for gene in negative_genes if gene in zexpr.index]
        pos_term = zexpr.loc[pos_genes].mean(axis=0) if pos_genes else 0.0
        neg_term = zexpr.loc[neg_genes].mean(axis=0) if neg_genes else 0.0
        if isinstance(pos_term, float):
            pos_term = pd.Series(pos_term, index=zexpr.columns)
        if isinstance(neg_term, float):
            neg_term = pd.Series(neg_term, index=zexpr.columns)
        scores[score_name] = pos_term - neg_term

    score_df = pd.DataFrame(scores)
    meta = build_sample_metadata(expression)
    long_scores = (
        score_df.reset_index(names="sample")
        .melt(id_vars="sample", var_name="score_name", value_name="score_value")
        .merge(meta, on="sample", how="left")
    )

    paired_rows = []
    for score_name in sorted(long_scores["score_name"].unique()):
        subset = long_scores[long_scores["score_name"].eq(score_name)].copy()
        wide = subset.pivot(index="donor", columns="condition", values="score_value")
        for cond_a, cond_b, contrast_name in SCORE_CONTRASTS:
            if cond_a not in wide.columns or cond_b not in wide.columns:
                continue
            diffs = wide[cond_b] - wide[cond_a]
            paired_rows.append(
                {
                    "score_name": score_name,
                    "contrast": contrast_name,
                    "contrast_label": CONTRAST_LABELS[contrast_name],
                    "n_donors": int(diffs.notna().sum()),
                    "mean_delta": float(diffs.mean()),
                    "median_delta": float(diffs.median()),
                    "signflip_p": exact_paired_signflip_p(diffs),
                }
            )

    return long_scores, pd.DataFrame(paired_rows)


def additive_model_summary(signature: pd.DataFrame, contrasts: list[ContrastData]) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = signature[["gene", "ev_mean_log2fc", "arm", "module_group"]].copy()
    for contrast in contrasts:
        merged = merged.merge(contrast.gene_table[["gene", f"{contrast.name}_log2fc"]], on="gene", how="left")

    required_cols = [
        "lps_vs_media_log2fc",
        "mtori_lps_vs_lps_log2fc",
        "mtori_lps_vs_media_log2fc",
    ]
    complete = merged.dropna(subset=required_cols).copy()

    model_specs = {
        "lps_only": ["lps_vs_media_log2fc"],
        "mtori_only": ["mtori_lps_vs_lps_log2fc"],
        "combined_direct_only": ["mtori_lps_vs_media_log2fc"],
        "additive_two_component": ["lps_vs_media_log2fc", "mtori_lps_vs_lps_log2fc"],
    }

    rows = []
    predictions = complete[["gene", "ev_mean_log2fc", "arm", "module_group"]].copy()
    y = complete["ev_mean_log2fc"].to_numpy(dtype=float)

    for model_name, columns in model_specs.items():
        x = complete[columns].to_numpy(dtype=float)
        beta, fitted = fit_no_intercept(x, y)
        loocv = loocv_predictions(x, y)
        fit_stats = regression_metrics(pd.Series(y), pd.Series(fitted))
        loocv_stats = regression_metrics(pd.Series(y), pd.Series(loocv))

        row = {
            "model": model_name,
            "predictors": ",".join(columns),
            "n_genes": len(complete),
            "beta_json": json.dumps({col: float(val) for col, val in zip(columns, beta)}),
        }
        row.update({f"fit_{k}": v for k, v in fit_stats.items()})
        row.update({f"loocv_{k}": v for k, v in loocv_stats.items()})
        rows.append(row)

        predictions[f"{model_name}_fitted"] = fitted
        predictions[f"{model_name}_loocv"] = loocv

    return pd.DataFrame(rows), predictions


def anchor_gene_table(contrasts: list[ContrastData]) -> pd.DataFrame:
    merged = None
    for contrast in contrasts:
        cols = [col for col in contrast.gene_table.columns if col.endswith("_log2fc") or col.endswith("_fdr")]
        frame = contrast.gene_table[["gene", *cols]].copy()
        merged = frame if merged is None else merged.merge(frame, on="gene", how="outer")
    anchor_df = merged[merged["gene"].isin(ANCHOR_GENES)].copy()
    return anchor_df.sort_values("gene")


def plot_summary(concordance_df: pd.DataFrame, sample_scores: pd.DataFrame, additive_df: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="talk")

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    plot_groups = ["full_signature", "down_arm", "up_arm", "module_b"]
    plot_data = concordance_df[concordance_df["group"].isin(plot_groups)].copy()
    plot_data["group"] = plot_data["group"].map(
        {
            "full_signature": "Full signature",
            "down_arm": "Down arm",
            "up_arm": "Up arm",
            "module_b": "Module B",
        }
    )
    sns.barplot(
        data=plot_data,
        x="group",
        y="concordance_fraction",
        hue="contrast_label",
        palette="Set2",
        ax=axes[0, 0],
    )
    axes[0, 0].set_ylim(0, 1)
    axes[0, 0].set_xlabel("")
    axes[0, 0].set_ylabel("Concordance fraction")
    axes[0, 0].set_title("Gene-level concordance with EV signature")
    axes[0, 0].legend(frameon=False, title="")

    corr_models = additive_df[["model", "loocv_pearson_r", "loocv_spearman_rho", "loocv_r2"]].copy()
    corr_models = corr_models.melt(id_vars="model", var_name="metric", value_name="value")
    corr_models["metric"] = corr_models["metric"].map(
        {
            "loocv_pearson_r": "LOOCV Pearson r",
            "loocv_spearman_rho": "LOOCV Spearman rho",
            "loocv_r2": "LOOCV R²",
        }
    )
    sns.barplot(data=corr_models, x="model", y="value", hue="metric", palette="muted", ax=axes[0, 1])
    axes[0, 1].set_xlabel("")
    axes[0, 1].set_ylabel("Model fit")
    axes[0, 1].set_title("Cross-validated reconstruction of EV effect")
    axes[0, 1].tick_params(axis="x", rotation=20)
    axes[0, 1].legend(frameon=False, title="")

    for ax, score_name, title in [
        (axes[1, 0], "down_arm_ev_score", "Down-arm EV-like score"),
        (axes[1, 1], "module_b_ev_score", "Module B EV-like score"),
    ]:
        subset = sample_scores[sample_scores["score_name"].eq(score_name)].copy()
        order = ["M", "L", "AL"]
        for donor, donor_df in subset.groupby("donor"):
            donor_df = donor_df.set_index("condition").reindex(order).reset_index()
            ax.plot(donor_df["condition"], donor_df["score_value"], marker="o", linewidth=1.2, alpha=0.6)
        mean_df = subset.groupby("condition", as_index=False)["score_value"].mean()
        mean_df["condition"] = pd.Categorical(mean_df["condition"], categories=order, ordered=True)
        mean_df = mean_df.sort_values("condition")
        ax.plot(mean_df["condition"], mean_df["score_value"], color="black", marker="o", linewidth=3)
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("Higher = more EV-like")

    fig.suptitle("GSE187403 orthogonal validation of the EV composite-state model", y=1.02, fontsize=18)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "gse187403_validation_summary.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_summary(
    concordance_df: pd.DataFrame,
    paired_scores: pd.DataFrame,
    additive_df: pd.DataFrame,
) -> None:
    def metric_row(contrast: str, group: str) -> pd.Series:
        return concordance_df[
            concordance_df["contrast"].eq(contrast) & concordance_df["group"].eq(group)
        ].iloc[0]

    lps_down = metric_row("lps_vs_media", "down_arm")
    lps_up = metric_row("lps_vs_media", "up_arm")
    mtori_down = metric_row("mtori_lps_vs_lps", "down_arm")
    mtori_up = metric_row("mtori_lps_vs_lps", "up_arm")
    combined_full = metric_row("mtori_lps_vs_media", "full_signature")

    additive_best = additive_df.sort_values("loocv_pearson_r", ascending=False).iloc[0]
    direct_model = additive_df[additive_df["model"].eq("combined_direct_only")].iloc[0]
    lps_model = additive_df[additive_df["model"].eq("lps_only")].iloc[0]
    mtori_model = additive_df[additive_df["model"].eq("mtori_only")].iloc[0]

    core_scores = paired_scores[paired_scores["score_name"].eq("core_ev_score")].copy()
    module_b_scores = paired_scores[paired_scores["score_name"].eq("module_b_ev_score")].copy()
    core_lps = core_scores[core_scores["contrast"].eq("lps_vs_media")].iloc[0]
    core_mtori = core_scores[core_scores["contrast"].eq("mtori_lps_vs_lps")].iloc[0]
    core_combined = core_scores[core_scores["contrast"].eq("mtori_lps_vs_media")].iloc[0]
    module_b_mtori = module_b_scores[module_b_scores["contrast"].eq("mtori_lps_vs_lps")].iloc[0]

    lines = [
        "# GSE187403 Orthogonal Validation Summary",
        "",
        "## Dataset",
        "- GEO series: GSE187403",
        "- Title: Transcriptional effects of catalytic mTOR inhibition among LPS-stimulated primary human monocytes",
        "- Design: six donors; media, LPS, and AZD2014+LPS conditions",
        "",
        "## Key findings",
        (
            f"- LPS alone reproduced the EV down arm better than the up arm "
            f"({int(lps_down['n_concordant'])}/{int(lps_down['n_available'])}, "
            f"{lps_down['concordance_fraction']:.3f} vs "
            f"{int(lps_up['n_concordant'])}/{int(lps_up['n_available'])}, "
            f"{lps_up['concordance_fraction']:.3f})."
        ),
        (
            f"- AZD2014 within LPS did not preferentially reconstruct the EV up arm "
            f"({int(mtori_up['n_concordant'])}/{int(mtori_up['n_available'])}, "
            f"{mtori_up['concordance_fraction']:.3f}) and remained less aligned than the down arm "
            f"({int(mtori_down['n_concordant'])}/{int(mtori_down['n_available'])}, "
            f"{mtori_down['concordance_fraction']:.3f})."
        ),
        (
            f"- The direct combined perturbation (AZD2014+LPS vs media) achieved "
            f"{int(combined_full['n_concordant'])}/{int(combined_full['n_available'])} "
            f"full-signature concordance ({combined_full['concordance_fraction']:.3f})."
        ),
        (
            f"- In leave-one-gene-out prediction, the best reconstruction model was "
            f"`{additive_best['model']}` with Pearson r = {additive_best['loocv_pearson_r']:.3f}, "
            f"Spearman rho = {additive_best['loocv_spearman_rho']:.3f}, and R² = {additive_best['loocv_r2']:.3f}."
        ),
        (
            f"- For comparison, direct combined contrast alone reached LOOCV Pearson r = "
            f"{direct_model['loocv_pearson_r']:.3f}; LPS alone = {lps_model['loocv_pearson_r']:.3f}; "
            f"AZD2014-within-LPS alone = {mtori_model['loocv_pearson_r']:.3f}."
        ),
        "",
        "## Sample-score behavior",
        (
            f"- Core EV-like score: LPS vs media mean delta = {core_lps['mean_delta']:.3f} "
            f"(sign-flip p = {core_lps['signflip_p']:.3f}); "
            f"AZD2014+LPS vs LPS = {core_mtori['mean_delta']:.3f} "
            f"(p = {core_mtori['signflip_p']:.3f}); "
            f"AZD2014+LPS vs media = {core_combined['mean_delta']:.3f} "
            f"(p = {core_combined['signflip_p']:.3f})."
        ),
        (
            f"- Module B EV-like score decreased under AZD2014 within LPS "
            f"(mean delta = {module_b_mtori['mean_delta']:.3f}, p = {module_b_mtori['signflip_p']:.3f}), "
            "which argues against a simple claim that acute catalytic mTOR inhibition in this monocyte context reproduces the MYC-associated arm."
        ),
        "",
        "## Interpretation",
        (
            "- This dataset is useful because it supplies the missing mTOR-inhibition side of the manuscript's "
            "two-component model in primary human monocytes. The most important question is not whether it "
            "perfectly reproduces the EV state, but whether LPS and mTOR inhibition contribute different arms."
        ),
        (
            "- Here, the strongest reproducible signal is still the LPS-like down arm. The mTOR-inhibitor step adds only modest full-signature support and does not independently recover Module B in this acute monocyte setting."
        ),
    ]
    (OUT_DIR / "gse187403_validation_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()

    signature = load_signature()
    contrasts = [load_contrast(name, CONTRAST_LABELS[name], path) for name, path in CONTRAST_PATHS.items()]
    expression = combine_expression_tables(contrasts)

    all_concordance = []
    gene_level = None
    for contrast in contrasts:
        summary_df, merged = summarize_concordance(signature, contrast)
        all_concordance.append(summary_df)
        selected_cols = [
            "gene",
            "direction",
            "arm",
            "module_group",
            "ev_mean_log2fc",
            f"{contrast.name}_log2fc",
            f"{contrast.name}_fdr",
            "concordant",
        ]
        gene_subset = merged[selected_cols].copy()
        gene_level = gene_subset if gene_level is None else gene_level.merge(
            gene_subset, on=["gene", "direction", "arm", "module_group", "ev_mean_log2fc"], how="outer"
        )

    concordance_df = pd.concat(all_concordance, ignore_index=True)
    sample_scores, paired_scores = compute_sample_scores(signature, expression)
    additive_df, additive_predictions = additive_model_summary(signature, contrasts)
    anchor_df = anchor_gene_table(contrasts)

    concordance_df.to_csv(OUT_DIR / "gse187403_concordance_summary.csv", index=False)
    gene_level.sort_values(["direction", "module_group", "gene"]).to_csv(
        OUT_DIR / "gse187403_gene_level_validation.csv", index=False
    )
    sample_scores.sort_values(["score_name", "donor", "condition"]).to_csv(
        OUT_DIR / "gse187403_sample_scores.csv", index=False
    )
    paired_scores.sort_values(["score_name", "contrast"]).to_csv(
        OUT_DIR / "gse187403_paired_score_tests.csv", index=False
    )
    additive_df.to_csv(OUT_DIR / "gse187403_additive_model_summary.csv", index=False)
    additive_predictions.to_csv(OUT_DIR / "gse187403_additive_model_predictions.csv", index=False)
    anchor_df.to_csv(OUT_DIR / "gse187403_anchor_genes.csv", index=False)

    plot_summary(concordance_df, sample_scores, additive_df)
    write_summary(concordance_df, paired_scores, additive_df)


if __name__ == "__main__":
    main()
