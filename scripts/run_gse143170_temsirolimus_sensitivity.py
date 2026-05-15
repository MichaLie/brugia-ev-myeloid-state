#!/usr/bin/env python3

from __future__ import annotations

import gzip
import math
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
COUNT_PATH = RAW_DIR / "GSE143170_count_matrix_RNAseq.txt.gz"


def safe_corr(func, x: pd.Series, y: pd.Series) -> float:
    valid = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(valid) < 3:
        return math.nan
    if valid["x"].nunique() < 2 or valid["y"].nunique() < 2:
        return math.nan
    return float(func(valid["x"], valid["y"])[0])


def exact_paired_signflip_p(diffs: pd.Series) -> float:
    arr = diffs.dropna().astype(float).to_numpy()
    if len(arr) == 0:
        return math.nan
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
    with gzip.open(COUNT_PATH, "rt", encoding="utf-8", errors="replace") as handle:
        return pd.read_csv(handle, sep="\t", index_col=0)


def median_ratio_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    geom_mean = np.exp(np.log(counts.replace(0, np.nan)).mean(axis=1))
    valid = geom_mean.replace(0, np.nan).dropna()
    ratios = counts.loc[valid.index].div(valid, axis=0).replace(0, np.nan)
    size_factors = ratios.median(axis=0)
    return counts.div(size_factors, axis=1)


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1, ddof=0).replace(0, np.nan), axis=0)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    signature = load_signature()
    counts = load_counts()
    norm = median_ratio_normalize(counts)

    ph73 = [col for col in norm.columns if col.startswith("pH73_")]
    tems = [col for col in norm.columns if col.startswith("Tems_")]
    ph65 = [col for col in norm.columns if col.startswith("pH65_")]

    tems_lfc = np.log2(norm[tems].mean(axis=1) + 1) - np.log2(norm[ph73].mean(axis=1) + 1)
    ph65_lfc = np.log2(norm[ph65].mean(axis=1) + 1) - np.log2(norm[ph73].mean(axis=1) + 1)

    merged = signature.merge(tems_lfc.rename("tems_vs_ph73_log2fc"), left_on="gene", right_index=True, how="left")
    merged = merged.merge(ph65_lfc.rename("ph65_vs_ph73_log2fc"), left_on="gene", right_index=True, how="left")
    merged["tems_concordant"] = np.sign(merged["tems_vs_ph73_log2fc"]) == np.sign(merged["ev_mean_log2fc"])
    merged["ph65_concordant"] = np.sign(merged["ph65_vs_ph73_log2fc"]) == np.sign(merged["ev_mean_log2fc"])

    summary_rows = []
    for contrast_name, lfc_col, conc_col in [
        ("tems_vs_ph73", "tems_vs_ph73_log2fc", "tems_concordant"),
        ("ph65_vs_ph73", "ph65_vs_ph73_log2fc", "ph65_concordant"),
    ]:
        for group_name, mask in {
            "full_signature": merged[lfc_col].notna(),
            "down_arm": merged[lfc_col].notna() & merged["arm"].eq("down_arm"),
            "up_arm": merged[lfc_col].notna() & merged["arm"].eq("up_arm"),
            "module_a": merged[lfc_col].notna() & merged["module_group"].eq("module_a"),
            "module_b": merged[lfc_col].notna() & merged["module_group"].eq("module_b"),
            "module_c": merged[lfc_col].notna() & merged["module_group"].eq("module_c"),
        }.items():
            subset = merged.loc[mask].copy()
            n = len(subset)
            concordant = int(subset[conc_col].sum()) if n else 0
            summary_rows.append(
                {
                    "contrast": contrast_name,
                    "group": group_name,
                    "n_available": n,
                    "n_concordant": concordant,
                    "concordance_fraction": float(concordant / n) if n else math.nan,
                    "pearson_r": safe_corr(pearsonr, subset["ev_mean_log2fc"], subset[lfc_col]),
                    "spearman_rho": safe_corr(spearmanr, subset["ev_mean_log2fc"], subset[lfc_col]),
                }
            )

    summary_df = pd.DataFrame(summary_rows)

    common_genes = [gene for gene in signature["gene"] if gene in norm.index]
    zexpr = zscore_rows(norm.loc[common_genes])
    up_genes = [gene for gene in signature.loc[signature["arm"].eq("up_arm"), "gene"] if gene in zexpr.index]
    down_genes = [gene for gene in signature.loc[signature["arm"].eq("down_arm"), "gene"] if gene in zexpr.index]
    module_b_genes = [gene for gene in signature.loc[signature["module_group"].eq("module_b"), "gene"] if gene in zexpr.index]

    score_df = pd.DataFrame(
        {
            "core_ev_score": zexpr.loc[up_genes].mean(axis=0) - zexpr.loc[down_genes].mean(axis=0),
            "up_arm_ev_score": zexpr.loc[up_genes].mean(axis=0),
            "down_arm_ev_score": -zexpr.loc[down_genes].mean(axis=0),
            "module_b_ev_score": zexpr.loc[module_b_genes].mean(axis=0),
        }
    ).reset_index(names="sample")
    score_df["donor"] = score_df["sample"].str.extract(r"_(\d+)$")[0].astype(int)
    score_df["condition"] = np.select(
        [
            score_df["sample"].str.startswith("pH73_"),
            score_df["sample"].str.startswith("Tems_"),
            score_df["sample"].str.startswith("pH65_"),
        ],
        ["pH73", "Tems", "pH65"],
        default="other",
    )

    paired_rows = []
    for score_name in ["core_ev_score", "up_arm_ev_score", "down_arm_ev_score", "module_b_ev_score"]:
        wide = score_df.pivot(index="donor", columns="condition", values=score_name)
        for cond_a, cond_b, contrast_name in [("pH73", "Tems", "tems_vs_ph73"), ("pH73", "pH65", "ph65_vs_ph73")]:
            diffs = wide[cond_b] - wide[cond_a]
            paired_rows.append(
                {
                    "score_name": score_name,
                    "contrast": contrast_name,
                    "n_donors": int(diffs.notna().sum()),
                    "mean_delta": float(diffs.mean()),
                    "median_delta": float(diffs.median()),
                    "signflip_p": exact_paired_signflip_p(diffs),
                }
            )
    paired_df = pd.DataFrame(paired_rows)

    merged.sort_values(["direction", "module_group", "gene"]).to_csv(
        OUT_DIR / "gse143170_gene_level_sensitivity.csv", index=False
    )
    summary_df.to_csv(OUT_DIR / "gse143170_concordance_summary.csv", index=False)
    paired_df.to_csv(OUT_DIR / "gse143170_paired_score_tests.csv", index=False)

    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    plot_df = summary_df[summary_df["contrast"].eq("tems_vs_ph73") & summary_df["group"].isin(["full_signature", "down_arm", "up_arm", "module_b"])].copy()
    plot_df["group"] = plot_df["group"].map(
        {
            "full_signature": "Full signature",
            "down_arm": "Down arm",
            "up_arm": "Up arm",
            "module_b": "Module B",
        }
    )
    sns.barplot(data=plot_df, x="group", y="concordance_fraction", color="#7FA8D1", ax=axes[0])
    axes[0].set_ylim(0, 1)
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Concordance fraction")
    axes[0].set_title("Temsirolimus vs pH 7.3")

    paired_plot = score_df.melt(id_vars=["sample", "donor", "condition"], value_vars=["core_ev_score", "module_b_ev_score"], var_name="score_name", value_name="score_value")
    paired_plot = paired_plot[paired_plot["condition"].isin(["pH73", "Tems"])].copy()
    paired_plot["condition"] = pd.Categorical(paired_plot["condition"], categories=["pH73", "Tems"], ordered=True)
    sns.pointplot(
        data=paired_plot,
        x="condition",
        y="score_value",
        hue="score_name",
        errorbar=None,
        dodge=0.25,
        markers="o",
        linestyles="-",
        palette=["#D6866A", "#5C8A58"],
        ax=axes[1],
    )
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Higher = more EV-like")
    axes[1].set_title("EV-like score shift under temsirolimus")
    axes[1].legend(frameon=False, title="")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "gse143170_temsirolimus_sensitivity.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    tems_full = summary_df[(summary_df["contrast"].eq("tems_vs_ph73")) & (summary_df["group"].eq("full_signature"))].iloc[0]
    tems_down = summary_df[(summary_df["contrast"].eq("tems_vs_ph73")) & (summary_df["group"].eq("down_arm"))].iloc[0]
    tems_up = summary_df[(summary_df["contrast"].eq("tems_vs_ph73")) & (summary_df["group"].eq("up_arm"))].iloc[0]
    tems_core = paired_df[(paired_df["contrast"].eq("tems_vs_ph73")) & (paired_df["score_name"].eq("core_ev_score"))].iloc[0]
    tems_module_b = paired_df[(paired_df["contrast"].eq("tems_vs_ph73")) & (paired_df["score_name"].eq("module_b_ev_score"))].iloc[0]

    lines = [
        "# GSE143170 Temsirolimus Sensitivity Summary",
        "",
        "## Dataset",
        "- GEO series: GSE143170",
        "- System: human monocyte-derived dendritic cells differentiated for 5 days in vitro",
        "- Conditions used here: temsirolimus vs pH 7.3 control (paired by donor, n = 6)",
        "- Note: this is a sensitivity analysis from normalized count data, not a formal DESeq2 rerun",
        "",
        "## Key findings",
        (
            f"- Temsirolimus showed {int(tems_full['n_concordant'])}/{int(tems_full['n_available'])} "
            f"full-signature concordance ({tems_full['concordance_fraction']:.3f})."
        ),
        (
            f"- Concordance was strong for both the down arm "
            f"({int(tems_down['n_concordant'])}/{int(tems_down['n_available'])}, {tems_down['concordance_fraction']:.3f}) "
            f"and the up arm ({int(tems_up['n_concordant'])}/{int(tems_up['n_available'])}, {tems_up['concordance_fraction']:.3f})."
        ),
        (
            f"- Temsirolimus increased the core EV-like score by {tems_core['mean_delta']:.3f} "
            f"(sign-flip p = {tems_core['signflip_p']:.3f}) and the Module B EV-like score by "
            f"{tems_module_b['mean_delta']:.3f} (p = {tems_module_b['signflip_p']:.3f})."
        ),
        "",
        "## Interpretation",
        "- Unlike the acute AZD2014-within-LPS monocyte dataset, this moDC differentiation context shows a strong EV-like response under mTORC1 inhibition.",
        "- This supports the interpretive caution that the MYC-associated arm is compatible with mTOR-linked control, but context dependent rather than universally reproduced by every mTOR perturbation.",
    ]
    (OUT_DIR / "gse143170_temsirolimus_sensitivity_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
