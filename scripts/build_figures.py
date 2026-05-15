#!/usr/bin/env python3
from __future__ import annotations

import math
import re
from pathlib import Path
from string import ascii_uppercase

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
FIGURE_DIR = RESULTS_DIR / "figures"
MAIN_DIR = FIGURE_DIR / "main"
SUPP_DIR = FIGURE_DIR / "supplementary"
MAIN_PANEL_DIR = FIGURE_DIR / "panels" / "main"
SUPP_PANEL_DIR = FIGURE_DIR / "panels" / "supplementary"

try:
    from figure_export import save_publication_figure  # type: ignore  # noqa: E402
    from style_presets import apply_publication_style  # type: ignore  # noqa: E402
except ModuleNotFoundError:  # pragma: no cover
    def apply_publication_style(_: str) -> None:
        return None

    def save_publication_figure(
        fig: plt.Figure,
        outbase: Path,
        formats: list[str] | None = None,
        dpi: int = 600,
        pad_inches: float = 0.05,
    ) -> None:
        formats = formats or ["pdf", "png"]
        for ext in formats:
            save_kwargs = {
                "bbox_inches": "tight",
                "pad_inches": pad_inches,
                "facecolor": "white",
                "edgecolor": "none",
            }
            if ext != "pdf":
                save_kwargs["dpi"] = dpi
            fig.savefig(outbase.with_suffix(f".{ext}"), **save_kwargs)


BLUE = "#0072B2"
LIGHT_BLUE = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
GOLD = "#F0E442"
BLACK = "#000000"
GREY = "#8A8A8A"
LIGHT_GREY = "#D9D9D9"
DARK_GREY = "#4D4D4D"

MODULE_COLORS = {
    "module_a": BLUE,
    "module_b": VERMILLION,
    "module_c": GREEN,
    "down_arm": BLUE,
    "up_arm": VERMILLION,
}

CONTEXT_COLORS = {
    "tlr_context": BLUE,
    "orthogonal_validation": VERMILLION,
    "orthogonal_perturbation": VERMILLION,
    "direct_human_brugia_exposure": GREEN,
    "natural_human_filarial_infection": ORANGE,
    "same_parasite_l3_langerhans_exposure": PURPLE,
    "cross_helminth_modc_resting": LIGHT_BLUE,
    "cross_helminth_modc_lps_background": DARK_GREY,
}

SHORT_EXTERNAL_LABELS = {
    "DC_BrugiaCombined_vs_Baseline": "GSE360 DC combined",
    "Macrophage_BrugiaCombined_vs_Baseline": "GSE360 Mac combined",
    "InfectedPre_vs_Normal": "GSE2135 infected vs normal",
    "LC_L3_vs_LC_Control": "GSE42694 LC vs control",
    "LPS_Al_CPI_0.1uM_vs_LPS": "GSE250463 LPS + Al-CPI 0.1",
    "LPS_Al_CPI_1uM_vs_LPS": "GSE250463 LPS + Al-CPI 1",
    "LPS_ABF_vs_LPS": "GSE250463 LPS + ABF",
    "Al_CPI_0.1uM_vs_Control": "GSE250463 Al-CPI 0.1",
    "Al_CPI_1uM_vs_Control": "GSE250463 Al-CPI 1",
    "DC_Brugia50_vs_Baseline": "GSE360 DC 50",
    "DC_Brugia5_vs_Baseline": "GSE360 DC 5",
    "Macrophage_Brugia50_vs_Baseline": "GSE360 Mac 50",
    "Macrophage_Brugia5_vs_Baseline": "GSE360 Mac 5",
    "InfectedPre_vs_EndemicNegativePre": "GSE2135 infected vs endemic-",
    "InfectedPost_vs_InfectedPre": "GSE2135 post vs pre",
}


def read_csv(rel: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / rel)


def clean_term(term: str) -> str:
    term = term.replace("HALLMARK_", "")
    return term.replace("_", " ").title()


def clean_context(context: str) -> str:
    mapping = {
        "TLR7/8 monocytes vs unstimulated": "TLR7/8 mono",
        "TLR4 monocytes vs unstimulated": "TLR4 mono",
        "LPS vs media": "LPS",
        "AZD2014+LPS vs LPS": "AZD+LPS vs LPS",
        "AZD2014+LPS vs media": "AZD+LPS",
        "Temsirolimus vs pH7.3": "Tems",
        "Acidic pH vs pH7.3": "pH6.5",
        "Direct human Brugia exposure": "Direct Brugia",
        "Natural human filarial infection": "Natural infection",
        "Same-parasite L3 Langerhans exposure": "L3 LC",
        "Cross-helminth moDC resting": "Cross-helminth resting",
        "Cross-helminth moDC LPS background": "Cross-helminth LPS",
    }
    return mapping.get(context, context)


def module_short(name: str) -> str:
    mapping = {
        "A_Myeloid_lineage_state": "Module A",
        "B_Proliferative_metabolic": "Module B",
        "C_Context_sensitive_residual": "Module C",
        "module_a": "Module A",
        "module_b": "Module B",
        "module_c": "Module C",
        "down_arm": "Down arm",
        "up_arm": "Up arm",
    }
    return mapping.get(name, name)


def external_short_label(name: str) -> str:
    return SHORT_EXTERNAL_LABELS.get(name, name.replace("_", " "))


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
        ha="left",
    )


def finish_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3, width=0.5)


def save_panel_exports(
    fig: plt.Figure,
    outbase: Path,
    panel_axes: dict[str, plt.Axes],
) -> None:
    if outbase.parent == MAIN_DIR:
        panel_dir = MAIN_PANEL_DIR
    elif outbase.parent == SUPP_DIR:
        panel_dir = SUPP_PANEL_DIR
    else:
        return

    panel_dir.mkdir(parents=True, exist_ok=True)
    all_axes = list(fig.axes)

    for label, ax in panel_axes.items():
        visibility_state = {other: other.get_visible() for other in all_axes}
        label_text = None
        original_label_state = None
        for text in ax.texts:
            if text.get_text() == label:
                label_text = text
                original_label_state = (
                    text.get_position(),
                    text.get_ha(),
                    text.get_va(),
                )
                break
        if label_text is not None and label in {"B", "C", "D", "E", "F"}:
            label_text.set_position((-0.15, 1.08))
            label_text.set_ha("left")
            label_text.set_va("top")
        for other in all_axes:
            other.set_visible(other is ax)
        fig.canvas.draw()
        panel_base = panel_dir / f"{outbase.stem}_panel_{label}"
        for ext in ["png", "pdf"]:
            fig.savefig(
                panel_base.with_suffix(f".{ext}"),
                dpi=600,
                bbox_inches="tight",
                pad_inches=0.04,
                facecolor="white",
                edgecolor="none",
            )
        for other, visible in visibility_state.items():
            other.set_visible(visible)
        if label_text is not None and original_label_state is not None:
            pos, ha, va = original_label_state
            label_text.set_position(pos)
            label_text.set_ha(ha)
            label_text.set_va(va)
    fig.canvas.draw()


def save_figure(
    fig: plt.Figure,
    outbase: Path,
    panel_axes: dict[str, plt.Axes] | None = None,
) -> None:
    save_publication_figure(fig, outbase, formats=["pdf", "png"], dpi=600, pad_inches=0.05)
    if panel_axes:
        save_panel_exports(fig, outbase, panel_axes)
    plt.close(fig)


def parse_fraction(text: str) -> float:
    match = re.search(r"\((0?\.\d+)\)", str(text))
    if match:
        return float(match.group(1))
    if "/" in str(text):
        num, den = str(text).split(" ")[0].split("/")
        return float(num) / float(den)
    return float(text)


def format_p(p: float | int | None) -> str:
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return "NA"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3f}"


def color_for_context(label: str) -> str:
    if label.startswith("TLR7/8") or label.startswith("TLR4"):
        return BLUE
    if label in {"LPS", "AZD+LPS vs LPS", "AZD+LPS", "Tems", "pH6.5"}:
        return VERMILLION if label in {"AZD+LPS vs LPS", "AZD+LPS", "Tems"} else ORANGE
    if label == "Direct Brugia":
        return GREEN
    if label == "Natural infection":
        return ORANGE
    if label == "L3 LC":
        return PURPLE
    if label == "Cross-helminth resting":
        return LIGHT_BLUE
    if label == "Cross-helminth LPS":
        return DARK_GREY
    return GREY


def draw_heatmap(
    ax: plt.Axes,
    matrix: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    annotate: bool = False,
    fmt: str = ".2f",
) -> None:
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=30, ha="right")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels)
    if annotate:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = matrix[i, j]
                if math.isnan(value):
                    continue
                ax.text(
                    j,
                    i,
                    format(value, fmt),
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="black",
                )
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.outline.set_linewidth(0.4)
    cbar.ax.tick_params(labelsize=6, length=2)
    finish_axis(ax)


def figure1() -> None:
    overlap = read_csv("data/evidence/core/frozen_signature_cross_dataset_overlap.csv")
    signature = read_csv("data/evidence/core/frozen_76_gene_signature.csv")
    pathways = read_csv("data/evidence/core/frozen_signature_hallmark_enrichment.csv")

    fig = plt.figure(figsize=(7.2, 5.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1.0], width_ratios=[1.3, 1.0])

    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])

    add_panel_label(ax1, "A")
    add_panel_label(ax2, "B")
    add_panel_label(ax3, "C")

    ax1.scatter(
        overlap["log2FoldChange_monocytes"],
        overlap["log2FoldChange_moDCs"],
        s=9,
        color=LIGHT_GREY,
        alpha=0.35,
        linewidth=0,
        rasterized=True,
    )
    down = signature[signature["direction"] == "down"]
    up = signature[signature["direction"] == "up"]
    ax1.scatter(
        down["log2FoldChange_monocytes"],
        down["log2FoldChange_moDCs"],
        s=22,
        color=BLUE,
        alpha=0.9,
        label=f"Frozen down arm (n={len(down)})",
    )
    ax1.scatter(
        up["log2FoldChange_monocytes"],
        up["log2FoldChange_moDCs"],
        s=26,
        color=VERMILLION,
        alpha=0.9,
        label=f"Frozen up arm (n={len(up)})",
    )
    lim = 4.5
    ax1.axhline(0, color=LIGHT_GREY, linewidth=0.8)
    ax1.axvline(0, color=LIGHT_GREY, linewidth=0.8)
    ax1.set_xlim(-lim, lim)
    ax1.set_ylim(-lim, lim)
    ax1.set_xlabel("Monocyte log2 fold change (EV vs media)")
    ax1.set_ylabel("moDC log2 fold change (EV vs media)")
    ax1.set_title("Cross-dataset reconstruction identifies a stable shared signature", loc="left", fontsize=9)
    ax1.legend(frameon=False, fontsize=6, loc="upper left")
    finish_axis(ax1)

    counts = signature["direction"].value_counts().reindex(["down", "up"])
    ax2.barh(
        ["Frozen 76-gene signature"],
        [counts["down"]],
        color=BLUE,
        height=0.55,
    )
    ax2.barh(
        ["Frozen 76-gene signature"],
        [counts["up"]],
        left=[counts["down"]],
        color=VERMILLION,
        height=0.55,
    )
    ax2.text(counts["down"] / 2, 0, f"{counts['down']}", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    ax2.text(counts["down"] + counts["up"] / 2, 0, f"{counts['up']}", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    ax2.set_xlim(0, len(signature) + 5)
    ax2.set_xlabel("Genes")
    ax2.set_title("Direction split", loc="left", fontsize=9)
    handles = [
        patches.Patch(facecolor=BLUE, edgecolor="none", label="Down-regulated"),
        patches.Patch(facecolor=VERMILLION, edgecolor="none", label="Up-regulated"),
    ]
    legend = ax2.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.03, 0.98),
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        borderpad=0.28,
        handlelength=0.7,
        handleheight=0.7,
        handletextpad=0.55,
        labelspacing=0.35,
        fontsize=5.8,
    )
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor(LIGHT_GREY)
    legend.get_frame().set_linewidth(0.6)
    finish_axis(ax2)

    preferred_terms = [
        "HALLMARK_G2M_CHECKPOINT",
        "HALLMARK_E2F_TARGETS",
        "HALLMARK_MYC_TARGETS_V1",
        "HALLMARK_MYOGENESIS",
        "HALLMARK_APICAL_JUNCTION",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "HALLMARK_ANGIOGENESIS",
        "HALLMARK_P53_PATHWAY",
        "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
    ]
    psub = pathways[pathways["term"].isin(preferred_terms)].copy()
    psub["order"] = psub["term"].map({term: i for i, term in enumerate(preferred_terms)})
    psub = psub.sort_values("order")
    matrix = psub[["NES_monocytes", "NES_moDCs"]].to_numpy()
    row_labels = [clean_term(t) for t in psub["term"]]
    draw_heatmap(ax3, matrix, row_labels, ["Monocytes", "moDCs"], vmin=-2.2, vmax=2.2, annotate=True)
    ax3.set_title("Consensus Hallmark structure", loc="left", fontsize=9)

    save_figure(
        fig,
        MAIN_DIR / "Figure_1_reproducible_frozen_signature",
        panel_axes={"A": ax1, "B": ax2, "C": ax3},
    )


def figure2() -> None:
    tf_primary = read_csv("data/evidence/core/preregistered_tf_enrichment_primary.csv")
    coverage = read_csv("data/evidence/regulon/regulon_coverage_summary.csv")
    modules = read_csv("data/evidence/core/module_membership.csv")
    corr = read_csv("data/evidence/core/baseline_myeloid_state_correlations.csv")

    fig = plt.figure(figsize=(7.95, 5.95))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.31, 1.23],
        height_ratios=[1.0, 0.90],
        wspace=0.30,
        hspace=0.30,
    )

    ax1 = fig.add_subplot(gs[:, 0])
    right = gs[:, 1].subgridspec(3, 1, height_ratios=[1.0, 0.34, 0.34], hspace=0.48)
    ax2 = fig.add_subplot(right[0, 0])
    subgs = right[1:, 0].subgridspec(2, 1, hspace=0.74)
    ax3 = fig.add_subplot(subgs[0, 0])
    ax4 = fig.add_subplot(subgs[1, 0])

    panel_offsets = {
        "A": (-0.14, 1.03),
        "B": (-0.16, 1.06),
        "C": (-0.16, 1.08),
        "D": (-0.16, 1.08),
    }
    for label, ax in zip(ascii_uppercase[:4], [ax1, ax2, ax3, ax4]):
        xoff, yoff = panel_offsets[label]
        ax.text(
            xoff,
            yoff,
            label,
            transform=ax.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
            ha="left",
        )

    selected_tfs = ["SPI1", "RUNX1", "MECOM", "CEBPB", "LXR", "MYC", "MAX"]
    plot_df = tf_primary[
        tf_primary["TF"].isin(selected_tfs) & tf_primary["Arm"].isin(["Down", "Up", "Combined"])
    ].copy()
    tf_order = [tf for tf in selected_tfs if tf in plot_df["TF"].unique()]
    y_map = {tf: i for i, tf in enumerate(reversed(tf_order))}
    arm_offsets = {"Down": 0.24, "Combined": 0.0, "Up": -0.24}
    arm_colors = {"Down": BLUE, "Combined": DARK_GREY, "Up": VERMILLION}
    threshold = -math.log10(0.05)
    ax1.axvline(threshold, color=LIGHT_GREY, linestyle="--", linewidth=0.8)
    for _, row in plot_df.iterrows():
        x = -math.log10(max(row["Best_FDR"], 1e-8))
        y = y_map[row["TF"]] + arm_offsets[row["Arm"]]
        face = arm_colors[row["Arm"]] if row["Best_FDR"] < 0.05 else "white"
        ax1.scatter(
            x,
            y,
            s=40 + row["Libraries_Significant"] * 70,
            facecolor=face,
            edgecolor=arm_colors[row["Arm"]],
            linewidth=1.0,
            zorder=3,
        )
    ax1.set_yticks(list(y_map.values()))
    ax1.set_yticklabels(list(reversed(tf_order)))
    ax1.set_xlabel("-log10 best FDR")
    ax1.set_title("Arm-specific TF enrichment reveals regulatory heterogeneity", loc="left", fontsize=8.4, pad=2)
    ax1.set_xlim(0, max(5.5, ax1.get_xlim()[1]))
    finish_axis(ax1)
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=arm_colors[arm], markeredgecolor=arm_colors[arm], markersize=5.6, label=arm)
        for arm in ["Down", "Combined", "Up"]
    ]
    legend = ax1.legend(
        handles=handles,
        frameon=True,
        fancybox=False,
        framealpha=1.0,
        fontsize=5.8,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.995),
        title="Arm",
        title_fontsize=6,
        borderpad=0.25,
        handletextpad=0.45,
        labelspacing=0.32,
    )
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor(LIGHT_GREY)
    legend.get_frame().set_linewidth(0.6)

    cov = coverage[coverage["TF"].isin(["SPI1", "RUNX1"])].copy()
    contexts = ["Down", "ModA", "ModC", "Up"]
    context_cols = ["Down_fraction", "ModA_fraction", "ModC_fraction", "Up_fraction"]
    x = np.arange(len(contexts))
    width = 0.32
    for i, tf in enumerate(["SPI1", "RUNX1"]):
        row = cov[cov["TF"] == tf].iloc[0]
        vals = [row[col] for col in context_cols]
        ax2.bar(x + (i - 0.5) * width, vals, width=width, color=[LIGHT_BLUE, GOLD][i], edgecolor=BLACK, linewidth=0.4, label=tf)
    ax2.set_xticks(x)
    ax2.set_xticklabels(["Down arm", "Module A", "Module C", "Up arm"], rotation=20, ha="right")
    ax2.set_ylim(0, 1.0)
    ax2.set_ylabel("Union target coverage fraction")
    ax2.set_title("SPI1/RUNX1 coverage is strongest in the down arm", loc="left", fontsize=8.4, pad=2)
    ax2.legend(frameon=False, fontsize=6, loc="upper right")
    ax2.tick_params(axis="x", labelsize=6.4, pad=1)
    finish_axis(ax2)

    module_counts = modules["module"].value_counts().reindex(
        ["A_Myeloid_lineage_state", "B_Proliferative_metabolic", "C_Context_sensitive_residual"]
    )
    ax3.barh(
        ["Module A", "Module B", "Module C"],
        module_counts.values,
        color=[BLUE, VERMILLION, GREEN],
        edgecolor=BLACK,
        linewidth=0.4,
    )
    for i, v in enumerate(module_counts.values):
        ax3.text(v + 0.6, i, str(v), va="center", fontsize=7)
    ax3.set_xlim(0, max(module_counts.values) + 8)
    ax3.set_xlabel("Genes")
    ax3.set_title("Deterministic module sizes", loc="left", fontsize=8.4, pad=4)
    ax3.tick_params(axis="y", labelsize=6.5, pad=2)
    finish_axis(ax3)

    corr_labels = [
        "Core EV score",
        "Up arm",
        "Down arm",
    ]
    corr_colors = [DARK_GREY, VERMILLION, BLUE]
    ax4.barh(corr_labels, corr["Spearman_r"], color=corr_colors, edgecolor=BLACK, linewidth=0.4)
    for i, v in enumerate(corr["Spearman_r"]):
        if v >= 0:
            xpos = v + 0.03
            ha = "left"
            color = BLACK
        else:
            xpos = -0.05
            ha = "right"
            color = "white"
        ax4.text(xpos, i, f"r={v:.2f}", va="center", ha=ha, fontsize=7, color=color)
    ax4.set_xlim(-1, 1)
    ax4.axvline(0, color=LIGHT_GREY, linewidth=0.8)
    ax4.set_xlabel("Spearman r")
    ax4.set_title("Baseline myeloid geometry supports a host-state model", loc="left", fontsize=8.4, pad=4)
    ax4.tick_params(axis="y", labelsize=6.5, pad=2)
    finish_axis(ax4)

    save_figure(
        fig,
        MAIN_DIR / "Figure_2_two_arm_regulatory_architecture",
        panel_axes={"A": ax1, "B": ax2, "C": ax3, "D": ax4},
    )


def figure3() -> None:
    panel = read_csv("data/evidence/core/regulator_panel_target_rerun_results.csv")
    mtor_control = read_csv("data/evidence/core/mir100_mtor_positive_control.csv").iloc[0]
    motif6 = read_csv("data/evidence/core/frozen_signature_6mer_enrichment.csv")
    motif7 = read_csv("data/evidence/core/frozen_signature_7mer_enrichment.csv")

    fig = plt.figure(figsize=(7.2, 6.25), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1.0], height_ratios=[1.0, 1.0])

    ax1 = fig.add_subplot(gs[:, 0])
    right = gs[:, 1].subgridspec(3, 1, height_ratios=[1.0, 1.0, 0.52], hspace=0.20)
    ax2 = fig.add_subplot(right[0, 0])
    ax3 = fig.add_subplot(right[1, 0])
    ax4 = fig.add_subplot(right[2, 0])

    panel_offsets = {
        "A": (-0.12, 1.05),
        "B": (-0.18, 1.08),
        "C": (-0.18, 1.08),
    }
    for label, ax in zip(ascii_uppercase[:3], [ax1, ax2, ax3]):
        xoff, yoff = panel_offsets[label]
        ax.text(
            xoff,
            yoff,
            label,
            transform=ax.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
            ha="left",
        )

    panel = panel.sort_values(["Category", "Gene"]).reset_index(drop=True)
    hit_matrix = panel[["Targeted_by_TargetScan_rerun", "Targeted_by_miRDB_rerun"]].astype(int).to_numpy()
    im = ax1.imshow(hit_matrix, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(["TargetScan", "miRDB"])
    ax1.set_yticks(range(len(panel)))
    ax1.set_yticklabels(panel["Gene"], fontsize=6)
    ax1.set_title("Preserved regulator panel shows no direct hsa-miR-100-5p hits", loc="left", fontsize=8.7)
    for i in range(hit_matrix.shape[0] + 1):
        ax1.axhline(i - 0.5, color="white", linewidth=0.6)
    for j in range(hit_matrix.shape[1] + 1):
        ax1.axvline(j - 0.5, color="white", linewidth=0.6)
    cat_positions = {}
    for cat, sub in panel.groupby("Category", sort=False):
        cat_positions[cat] = (sub.index.min() + sub.index.max()) / 2.0
    for cat, pos in cat_positions.items():
        ax1.text(-0.8, pos, cat.replace("Master TF", "Master TF"), va="center", ha="right", fontsize=6, color=DARK_GREY)
    ax1.text(
        0.60,
        0.985,
        "0/17 genes hit in either rerun resource",
        transform=ax1.transAxes,
        ha="center",
        va="top",
        fontsize=6.5,
        color=DARK_GREY,
        bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor=LIGHT_GREY, linewidth=0.6),
    )
    finish_axis(ax1)
    plt.colorbar(im, ax=ax1, fraction=0.03, pad=0.02).ax.tick_params(labelsize=6, length=2)

    def motif_panel(ax: plt.Axes, df: pd.DataFrame, title: str, color: str) -> None:
        plot_df = df.copy().sort_values("fold_change", ascending=True)
        labels = [f"{m} ({s})" for m, s in zip(plot_df["motif_name"], plot_df["motif_seq"])]
        ax.barh(labels, plot_df["fold_change"], color=color, edgecolor=BLACK, linewidth=0.4)
        ax.axvline(1.0, color=LIGHT_GREY, linestyle="--", linewidth=0.8)
        for i, row in plot_df.reset_index(drop=True).iterrows():
            ax.text(
                row["fold_change"] + 0.03,
                i,
                f"p={row['pval_depletion']:.2e}",
                va="center",
                fontsize=6,
            )
        ax.set_xlim(0, 1.35)
        ax.set_xlabel("Observed / expected motif count")
        ax.set_title(title, loc="left", fontsize=9)
        finish_axis(ax)

    motif_panel(ax2, motif6, "6-mer seed motifs are depleted in the down arm", LIGHT_BLUE)
    motif_panel(ax3, motif7, "7-mer motifs show no compensatory enrichment", GOLD)

    ax4.set_axis_off()
    ax4.text(
        0.0,
        0.92,
        "Positive control retained",
        transform=ax4.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        fontweight="bold",
        color=DARK_GREY,
    )
    ax4.text(
        0.0,
        0.62,
        f"MTOR remains a predicted {mtor_control['hsa_miRNA']} target\nin both preserved rerun resources.",
        transform=ax4.transAxes,
        ha="left",
        va="top",
        fontsize=6.7,
        color=DARK_GREY,
    )
    ax4.text(
        0.0,
        0.20,
        (
            f"TargetScan percentile {int(mtor_control['TargetScan_context_percentile'])}   |   "
            f"TargetScan score {float(mtor_control['TargetScan_weighted_context_score']):.3f}   |   "
            f"miRDB score {float(mtor_control['miRDB_score']):.2f}"
        ),
        transform=ax4.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.3,
        color=DARK_GREY,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=LIGHT_GREY, linewidth=0.8),
    )

    save_figure(
        fig,
        MAIN_DIR / "Figure_3_anti_canonical_miRNA_evidence",
        panel_axes={"A": ax1, "B": ax2, "C": ax3},
    )


def figure4() -> None:
    orth = read_csv("results/tables/orthogonal_validation_summary.csv")
    ext_context = read_csv("results/tables/external_transport_context_summary.csv")
    gse143 = read_csv("results/tables/gse143170_paired_deseq2_concordance_summary.csv")
    score_187 = read_csv("results/validation/gse187403_paired_score_tests.csv")
    score_143 = read_csv("results/validation/gse143170_paired_score_tests.csv")

    fig = plt.figure(figsize=(7.4, 6.3), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.0, 1.0])
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])

    panel_offsets = {
        "A": (-0.08, 1.12),
        "B": (-0.12, 1.05),
        "C": (-0.12, 1.05),
    }
    for label, ax in zip(ascii_uppercase[:3], [ax1, ax2, ax3]):
        xoff, yoff = panel_offsets[label]
        ax.text(
            xoff,
            yoff,
            label,
            transform=ax.transAxes,
            fontsize=11,
            fontweight="bold",
            va="top",
            ha="left",
        )

    pH65_row = gse143[gse143["contrast"] == "ph65_vs_ph73"]
    contexts = [
        ("TLR7/8 mono", 0.903, 0.176),
        ("TLR4 mono", 0.815, 0.158),
        ("LPS", 0.786, 0.400),
        ("AZD+LPS vs LPS", 0.661, 0.500),
        ("Tems", 0.833, 0.895),
        (
            "pH6.5",
            float(pH65_row[pH65_row["group"] == "down_arm"]["concordance_fraction"].iloc[0]),
            float(pH65_row[pH65_row["group"] == "up_arm"]["concordance_fraction"].iloc[0]),
        ),
    ]
    x = np.arange(len(contexts))
    width = 0.36
    ax1.bar(
        x - width / 2,
        [c[1] for c in contexts],
        width=width,
        color=BLUE,
        edgecolor=BLACK,
        linewidth=0.4,
        label="Down arm or Module A",
    )
    ax1.bar(
        x + width / 2,
        [c[2] for c in contexts],
        width=width,
        color=VERMILLION,
        edgecolor=BLACK,
        linewidth=0.4,
        label="Up arm or Module B",
    )
    ax1.axhline(0.5, color=LIGHT_GREY, linestyle="--", linewidth=0.8)
    ax1.set_ylim(0, 1.0)
    ax1.set_xticks(x)
    ax1.set_xticklabels([c[0] for c in contexts], rotation=20, ha="right")
    ax1.set_ylabel("Directional concordance fraction")
    ax1.set_title("The two axes separate under TLR stimulation and orthogonal perturbation", loc="left", fontsize=9)
    ax1.legend(frameon=False, fontsize=6, loc="upper left")
    finish_axis(ax1)

    module_b_scores = pd.concat(
        [
            score_187[score_187["score_name"] == "module_b_ev_score"].assign(dataset="GSE187403"),
            score_143[score_143["score_name"] == "module_b_ev_score"].assign(dataset="GSE143170"),
        ],
        ignore_index=True,
    )
    label_map = {
        "lps_vs_media": "LPS",
        "mtori_lps_vs_lps": "AZD+LPS vs LPS",
        "mtori_lps_vs_media": "AZD+LPS",
        "tems_vs_ph73": "Tems",
        "ph65_vs_ph73": "pH6.5",
    }
    module_b_scores["label"] = module_b_scores["contrast"].map(label_map)
    module_b_scores = module_b_scores[module_b_scores["label"].notna()].copy()
    order = ["LPS", "AZD+LPS vs LPS", "AZD+LPS", "Tems", "pH6.5"]
    module_b_scores["label"] = pd.Categorical(module_b_scores["label"], order)
    module_b_scores = module_b_scores.sort_values("label")
    ax2.bar(
        module_b_scores["label"].astype(str),
        module_b_scores["mean_delta"],
        color=[color_for_context(x) for x in module_b_scores["label"].astype(str)],
        edgecolor=BLACK,
        linewidth=0.4,
    )
    ax2.axhline(0, color=LIGHT_GREY, linewidth=0.8)
    for i, row in module_b_scores.reset_index(drop=True).iterrows():
        ax2.text(i, row["mean_delta"] + (0.04 if row["mean_delta"] >= 0 else -0.08), f"p={format_p(row['signflip_p'])}", ha="center", va="bottom" if row["mean_delta"] >= 0 else "top", fontsize=6)
    ax2.set_ylabel("Mean Module B EV-like score delta")
    ax2.set_title("Module B responds in a context-dependent direction", loc="left", fontsize=9)
    ax2.tick_params(axis="x", rotation=20)
    finish_axis(ax2)

    state_points = [
        ("TLR7/8 mono", 0.903, 0.176),
        ("TLR4 mono", 0.815, 0.158),
        ("LPS", 0.786, 0.400),
        ("AZD+LPS vs LPS", 0.661, 0.500),
        ("Tems", 0.833, 0.895),
        ("pH6.5", float(pH65_row[pH65_row["group"] == "down_arm"]["concordance_fraction"].iloc[0]), float(pH65_row[pH65_row["group"] == "up_arm"]["concordance_fraction"].iloc[0])),
        ("Direct Brugia", float(ext_context[ext_context["context"] == "Direct human Brugia exposure"]["mean_down_fraction"].iloc[0]), float(ext_context[ext_context["context"] == "Direct human Brugia exposure"]["mean_up_fraction"].iloc[0])),
        ("Natural infection", float(ext_context[ext_context["context"] == "Natural human filarial infection"]["mean_down_fraction"].iloc[0]), float(ext_context[ext_context["context"] == "Natural human filarial infection"]["mean_up_fraction"].iloc[0])),
        ("L3 LC", float(ext_context[ext_context["context"] == "Same-parasite L3 Langerhans exposure"]["mean_down_fraction"].iloc[0]), float(ext_context[ext_context["context"] == "Same-parasite L3 Langerhans exposure"]["mean_up_fraction"].iloc[0])),
        ("Cross-helminth resting", float(ext_context[ext_context["context"] == "Cross-helminth moDC resting"]["mean_down_fraction"].iloc[0]), float(ext_context[ext_context["context"] == "Cross-helminth moDC resting"]["mean_up_fraction"].iloc[0])),
        ("Cross-helminth LPS", float(ext_context[ext_context["context"] == "Cross-helminth moDC LPS background"]["mean_down_fraction"].iloc[0]), float(ext_context[ext_context["context"] == "Cross-helminth moDC LPS background"]["mean_up_fraction"].iloc[0])),
    ]
    label_offsets = {
        "TLR7/8 mono": (0.012, 0.012),
        "TLR4 mono": (0.012, -0.018),
        "LPS": (0.012, 0.016),
        "AZD+LPS vs LPS": (0.012, -0.018),
        "Tems": (0.012, 0.016),
        "pH6.5": (0.012, 0.012),
        "Direct Brugia": (0.012, 0.022),
        "Natural infection": (0.012, -0.022),
        "L3 LC": (0.012, 0.012),
        "Cross-helminth resting": (0.012, 0.000),
        "Cross-helminth LPS": (0.012, 0.020),
    }
    for label, xval, yval in state_points:
        color = color_for_context(label)
        ax3.scatter(xval, yval, s=52, color=color, edgecolor="white", linewidth=0.5, zorder=3)
        dx, dy = label_offsets.get(label, (0.012, 0.012))
        ax3.text(xval + dx, yval + dy, label, fontsize=5.8, va="center")
    ax3.axvline(0.5, color=LIGHT_GREY, linestyle="--", linewidth=0.8)
    ax3.axhline(0.5, color=LIGHT_GREY, linestyle="--", linewidth=0.8)
    ax3.set_xlim(0, 1.0)
    ax3.set_ylim(0, 1.0)
    ax3.set_xlabel("Down-arm or Module A concordance")
    ax3.set_ylabel("Up-arm or Module B concordance")
    ax3.set_title("A two-axis state-space resolves the central pattern", loc="left", fontsize=9)
    finish_axis(ax3)

    save_figure(
        fig,
        MAIN_DIR / "Figure_4_tlr_and_orthogonal_validation",
        panel_axes={"A": ax1, "B": ax2, "C": ax3},
    )


def figure5() -> None:
    ext_contrast = read_csv("data/evidence/external/external_arm_transport_contrast_summary.csv")
    ext_context = read_csv("results/tables/external_transport_context_summary.csv")

    fig = plt.figure(figsize=(7.4, 6.35), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1.0], height_ratios=[1.0, 1.0])
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])

    for label, ax in zip(ascii_uppercase[:3], [ax1, ax2, ax3]):
        add_panel_label(ax, label)

    plot_df = ext_contrast.sort_values("concordance_fraction").copy()
    plot_df["short_label"] = plot_df["contrast"].map(external_short_label)
    y = np.arange(len(plot_df))
    colors = []
    for cgroup in plot_df["context_group"]:
        colors.append(CONTEXT_COLORS.get(cgroup, GREY))
    ax1.hlines(y, 0, plot_df["concordance_fraction"], color=colors, linewidth=2.0)
    ax1.scatter(plot_df["concordance_fraction"], y, color=colors, s=35, zorder=3)
    ax1.axvline(0.5, color=LIGHT_GREY, linestyle="--", linewidth=0.8)
    ax1.set_yticks(y)
    ax1.set_yticklabels(plot_df["short_label"], fontsize=6)
    ax1.set_xlim(0, 1.0)
    ax1.set_xlabel("Total concordance fraction")
    ax1.set_title("External projection defines context limits", loc="left", fontsize=8.6)
    ax1.text(
        0.98,
        0.03,
        "Supportive contrasts: 0 / 15",
        transform=ax1.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=LIGHT_GREY, linewidth=0.8),
    )
    finish_axis(ax1)

    order = [
        "Direct human Brugia exposure",
        "Natural human filarial infection",
        "Same-parasite L3 Langerhans exposure",
        "Cross-helminth moDC resting",
        "Cross-helminth moDC LPS background",
    ]
    cdf = ext_context.copy()
    cdf["context"] = pd.Categorical(cdf["context"], order)
    cdf = cdf.sort_values("context")
    x = np.arange(len(cdf))
    width = 0.35
    ax2.bar(x - width / 2, cdf["mean_down_fraction"], width=width, color=BLUE, edgecolor=BLACK, linewidth=0.4, label="Down arm")
    ax2.bar(x + width / 2, cdf["mean_up_fraction"], width=width, color=VERMILLION, edgecolor=BLACK, linewidth=0.4, label="Up arm")
    ax2.axhline(0.5, color=LIGHT_GREY, linestyle="--", linewidth=0.8)
    ax2.set_ylim(0, 1.0)
    ax2.set_xticks(x)
    ax2.set_xticklabels([clean_context(v) for v in cdf["context"]], rotation=32, ha="right")
    ax2.set_ylabel("Mean directional concordance")
    ax2.tick_params(axis="x", labelsize=6.3)
    ax2.set_title("Arm projection across external contexts", loc="left", fontsize=8.6)
    ax2.legend(frameon=False, fontsize=6, loc="upper right")
    finish_axis(ax2)

    ax3.set_axis_off()
    ax3.set_title("Bounded composite-state model", loc="left", fontsize=9, pad=6)

    left = patches.FancyBboxPatch((0.05, 0.55), 0.21, 0.24, boxstyle="round,pad=0.02", facecolor="white", edgecolor=BLACK, linewidth=0.8)
    down_box = patches.FancyBboxPatch((0.38, 0.64), 0.57, 0.20, boxstyle="round,pad=0.02", facecolor="#EAF3FA", edgecolor=BLUE, linewidth=1.0)
    up_box = patches.FancyBboxPatch((0.38, 0.34), 0.57, 0.20, boxstyle="round,pad=0.02", facecolor="#FCEFEA", edgecolor=VERMILLION, linewidth=1.0)
    boundary = patches.FancyBboxPatch((0.04, 0.06), 0.92, 0.20, boxstyle="round,pad=0.02", facecolor="#F7F7F7", edgecolor=GREY, linewidth=0.8)
    for patch in [left, down_box, up_box, boundary]:
        ax3.add_patch(patch)

    ax3.text(0.155, 0.67, "Brugia EV\nexposure", ha="center", va="center", fontsize=8, fontweight="bold")
    ax3.text(0.665, 0.775, "Innate-like down arm", ha="center", va="center", fontsize=7.8, color=BLUE, fontweight="bold")
    ax3.text(
        0.665,
        0.712,
        "SPI1/RUNX1-linked myeloid contraction\nInnate-concordant; recovered by LPS",
        ha="center",
        va="center",
        fontsize=6.2,
        color=DARK_GREY,
    )
    ax3.text(0.665, 0.465, "MYC-associated up arm", ha="center", va="center", fontsize=7.8, color=VERMILLION, fontweight="bold")
    ax3.text(
        0.665,
        0.395,
        "Context-dependent biosynthetic program\nRecovered by temsirolimus during\nmoDC differentiation",
        ha="center",
        va="center",
        fontsize=6.0,
        color=DARK_GREY,
    )
    ax3.text(
        0.50,
        0.16,
        "Canonical seed-based miRNA targeting\n"
        "does not explain the full state.\n"
        "External datasets define context limits\n"
        "rather than uniform confirmation.",
        ha="center",
        va="center",
        fontsize=5.9,
    )
    ax3.annotate("", xy=(0.38, 0.74), xytext=(0.26, 0.70), arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.2))
    ax3.annotate("", xy=(0.38, 0.44), xytext=(0.26, 0.64), arrowprops=dict(arrowstyle="->", color=VERMILLION, lw=1.2))

    save_figure(
        fig,
        MAIN_DIR / "Figure_5_external_boundary_and_composite_state",
        panel_axes={"A": ax1, "B": ax2, "C": ax3},
    )


def supplement_s1() -> None:
    tlr_all = read_csv("data/evidence/tlr_projection/tlr78_projection_all_genes.csv")
    tlr_summary = read_csv("data/evidence/tlr_projection/tlr78_concordance_summary.csv")

    fig = plt.figure(figsize=(7.2, 6.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1.0], height_ratios=[1.1, 1.0])
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])
    for label, ax in zip(ascii_uppercase[:3], [ax1, ax2, ax3]):
        add_panel_label(ax, label)

    mono = tlr_all[tlr_all["cell_type"] == "mono"].copy()
    module_order = {"A": 0, "C": 1, "B": 2}
    mono["module_order"] = mono["module"].map(module_order).fillna(9)
    mono["ev_mean_lfc"] = mono[["ev_log2FC_monocytes", "ev_log2FC_moDCs"]].mean(axis=1)
    mono = mono.sort_values(["module_order", "ev_direction", "ev_mean_lfc"])
    heat_cols = ["log2FC_TLR7_8", "log2FC_TLR4"]
    heat = mono[heat_cols].to_numpy()
    row_labels = mono["gene"].tolist()
    draw_heatmap(ax1, heat, row_labels, ["Monocyte TLR7/8", "Monocyte TLR4"], vmin=-3.0, vmax=3.0, annotate=False)
    ax1.set_title("Gene-level TLR projection in monocytes", loc="left", fontsize=9)
    if len(row_labels) > 30:
        step = 4
        ax1.set_yticks(range(0, len(row_labels), step))
        ax1.set_yticklabels([row_labels[i] for i in range(0, len(row_labels), step)], fontsize=5)

    mono["ev_mean_lfc"] = mono[["ev_log2FC_monocytes", "ev_log2FC_moDCs"]].mean(axis=1)
    color_map = {"A": BLUE, "B": VERMILLION, "C": GREEN}
    for module, sub in mono.groupby("module"):
        ax2.scatter(sub["ev_mean_lfc"], sub["log2FC_TLR7_8"], s=18, color=color_map.get(module, GREY), alpha=0.8, label=f"Module {module}")
    ax2.axhline(0, color=LIGHT_GREY, linewidth=0.8)
    ax2.axvline(0, color=LIGHT_GREY, linewidth=0.8)
    ax2.set_xlabel("Mean EV log2 fold change")
    ax2.set_ylabel("Monocyte TLR7/8 log2 fold change")
    ax2.set_title("Module B occupies the discordant quadrant under acute TLR7/8", loc="left", fontsize=9)
    ax2.legend(frameon=False, fontsize=6, loc="lower right")
    finish_axis(ax2)

    summary = tlr_summary.copy()
    summary["context"] = summary["stimulus"] + " " + summary["cell_type"]
    pivot = summary.pivot(index="module", columns="context", values="concordance_fraction").reindex(["A", "B", "C"])
    draw_heatmap(ax3, pivot.to_numpy(), [f"Module {m}" for m in pivot.index], list(pivot.columns), cmap="viridis", vmin=0, vmax=1, annotate=True)
    ax3.set_title("Concordance fractions across stimuli and cell types", loc="left", fontsize=9)

    save_figure(
        fig,
        SUPP_DIR / "Figure_S1_tlr_projection_detail",
        panel_axes={"A": ax1, "B": ax2, "C": ax3},
    )


def supplement_s2() -> None:
    tf_primary = read_csv("data/evidence/core/preregistered_tf_enrichment_primary.csv")
    module_eval = read_csv("data/evidence/core/module_gene_level_concordance.csv")

    fig = plt.figure(figsize=(7.2, 5.8), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    add_panel_label(ax1, "A")
    add_panel_label(ax2, "B")

    tf_primary = tf_primary.copy()
    tf_primary = tf_primary[tf_primary["Best_FDR"] < 0.1].copy()
    tf_primary["minus_log10_fdr"] = -np.log10(tf_primary["Best_FDR"].clip(lower=1e-8))
    tf_primary["label"] = tf_primary["TF"] + " (" + tf_primary["Arm"] + ")"
    tf_primary = tf_primary.sort_values("minus_log10_fdr", ascending=True).tail(18)
    colors = tf_primary["Arm"].map({"Down": BLUE, "Up": VERMILLION, "Combined": DARK_GREY})
    ax1.barh(tf_primary["label"], tf_primary["minus_log10_fdr"], color=colors, edgecolor=BLACK, linewidth=0.4)
    ax1.axvline(-math.log10(0.05), color=LIGHT_GREY, linestyle="--", linewidth=0.8)
    ax1.set_xlabel("-log10 best FDR")
    ax1.set_title("Expanded TF signal remains arm-specific", loc="left", fontsize=9)
    finish_axis(ax1)

    mod = module_eval.copy()
    mod["module_short"] = mod["module"].map({
        "A_Myeloid_lineage_state": "Module A",
        "B_Proliferative_metabolic": "Module B",
        "C_Context_sensitive_residual": "Module C",
    })
    xpos = {"Module A": 0, "Module B": 1, "Module C": 2}
    colors = {"resting": BLUE, "lps_background": VERMILLION}
    for ctx, sub in mod.groupby("context_type"):
        jitter = np.random.default_rng(4).uniform(-0.08, 0.08, size=len(sub))
        x = [xpos[m] for m in sub["module_short"]] + jitter
        ax2.scatter(x, sub["concordance_fraction"], s=26, alpha=0.8, color=colors.get(ctx, GREY), label=ctx.replace("_", " "))
    mean_df = mod.groupby("module_short", as_index=False)["concordance_fraction"].mean()
    ax2.plot([xpos[m] for m in mean_df["module_short"]], mean_df["concordance_fraction"], color=BLACK, linewidth=1.0, marker="o", markersize=3)
    ax2.set_xticks([0, 1, 2])
    ax2.set_xticklabels(["Module A", "Module B", "Module C"])
    ax2.set_ylim(0, 1.0)
    ax2.set_ylabel("Held-out gene-level concordance")
    ax2.set_title("Module B is the most independently robust component", loc="left", fontsize=9)
    ax2.legend(frameon=False, fontsize=6, loc="upper right")
    finish_axis(ax2)

    save_figure(
        fig,
        SUPP_DIR / "Figure_S2_regulatory_detail_and_module_robustness",
        panel_axes={"A": ax1, "B": ax2},
    )


def supplement_s3() -> None:
    val187 = read_csv("results/validation/gse187403_gene_level_validation.csv")
    val143 = read_csv("results/tables/gse143170_paired_deseq2_gene_level.csv")

    fig = plt.figure(figsize=(7.2, 6.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    axes = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]
    titles = [
        ("lps_vs_media_log2fc", "GSE187403: LPS vs media"),
        ("mtori_lps_vs_lps_log2fc", "GSE187403: AZD2014+LPS vs LPS"),
        ("tems_vs_ph73_log2fc", "GSE143170: temsirolimus vs pH7.3"),
        ("ph65_vs_ph73_log2fc", "GSE143170: pH6.5 vs pH7.3"),
    ]
    for label, ax in zip(ascii_uppercase[:4], axes):
        add_panel_label(ax, label)

    color_by_arm = {"down_arm": BLUE, "up_arm": VERMILLION}
    for ax, (col, title) in zip(axes, titles):
        if "GSE187403" in title:
            df = val187.copy()
        else:
            df = val143.copy()
        for arm, sub in df.groupby("arm"):
            ax.scatter(sub["ev_mean_log2fc"], sub[col], s=18, alpha=0.8, color=color_by_arm.get(arm, GREY), label=arm.replace("_", " "))
        ax.axhline(0, color=LIGHT_GREY, linewidth=0.8)
        ax.axvline(0, color=LIGHT_GREY, linewidth=0.8)
        ax.set_xlabel("Mean EV log2 fold change")
        ax.set_ylabel("Orthogonal perturbation log2 fold change")
        ax.set_title(title, loc="left", fontsize=8.5)
        finish_axis(ax)
    axes[0].legend(frameon=False, fontsize=6, loc="lower right")

    save_figure(
        fig,
        SUPP_DIR / "Figure_S3_orthogonal_validation_detail",
        panel_axes={"A": axes[0], "B": axes[1], "C": axes[2], "D": axes[3]},
    )


def supplement_s4() -> None:
    ext_contrast = read_csv("data/evidence/external/external_arm_transport_contrast_summary.csv")
    ext_scores = read_csv("data/evidence/external/external_standardized_module_contrast_summary.csv")

    fig = plt.figure(figsize=(7.2, 6.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    add_panel_label(ax1, "A")
    add_panel_label(ax2, "B")

    plot_df = ext_contrast.sort_values(["context_group", "concordance_fraction"]).copy()
    plot_df["label"] = plot_df["contrast"].map(external_short_label)
    y = np.arange(len(plot_df))
    colors = [CONTEXT_COLORS.get(c, GREY) for c in plot_df["context_group"]]
    ax1.barh(y, plot_df["concordance_fraction"], color=colors, edgecolor=BLACK, linewidth=0.3)
    ax1.axvline(0.5, color=LIGHT_GREY, linestyle="--", linewidth=0.8)
    ax1.set_yticks(y)
    ax1.set_yticklabels(plot_df["label"], fontsize=6)
    ax1.set_xlim(0, 1.0)
    ax1.set_xlabel("Total concordance fraction")
    ax1.set_title("All external contrasts remain below the formal supportive threshold", loc="left", fontsize=9)
    finish_axis(ax1)

    wanted_scores = [
        "core_ev_signed_score_z",
        "module_a_signed_score_z",
        "module_b_signed_score_z",
        "module_c_signed_score_z",
    ]
    hdf = ext_scores[ext_scores["score"].isin(wanted_scores)].copy()
    score_map = {
        "core_ev_signed_score_z": "Core",
        "module_a_signed_score_z": "Module A",
        "module_b_signed_score_z": "Module B",
        "module_c_signed_score_z": "Module C",
    }
    hdf["score_label"] = hdf["score"].map(score_map)
    pivot = hdf.pivot(index="score_label", columns="contrast", values="delta_case_minus_control")
    short_cols = [external_short_label(col) for col in pivot.columns]
    draw_heatmap(ax2, pivot.to_numpy(), list(pivot.index), short_cols, vmin=-2.0, vmax=2.0, annotate=False)
    ax2.set_title("Standardized score deltas vary strongly by context", loc="left", fontsize=9)
    ax2.tick_params(axis="x", labelsize=5.5)

    save_figure(
        fig,
        SUPP_DIR / "Figure_S4_external_transport_detail",
        panel_axes={"A": ax1, "B": ax2},
    )


def supplement_s5() -> None:
    tol_overlap = read_csv("data/evidence/tolerance/tolerance_signature_overlap.csv")
    tol_external = read_csv("data/evidence/tolerance/tolerance_markers_in_external.csv")
    pdac = read_csv("data/evidence/comparators/pdac_cross_cohort_state_geometry.csv")
    macrophage = read_csv("data/evidence/comparators/macrophage_core_by_polarization.csv")

    fig = plt.figure(figsize=(7.2, 6.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])
    for label, ax in zip(ascii_uppercase[:3], [ax1, ax2, ax3]):
        add_panel_label(ax, label)

    tol_mat = tol_overlap[["in_frozen_signature"]].astype(int).to_numpy()
    im = ax1.imshow(tol_mat, aspect="auto", cmap="Greys", vmin=0, vmax=1)
    ax1.set_xticks([0])
    ax1.set_xticklabels(["Present in\nfixed signature"])
    ax1.set_yticks(range(len(tol_overlap)))
    ax1.set_yticklabels(tol_overlap["gene"], fontsize=6)
    ax1.set_xlim(-0.5, 1.55)
    ax1.set_title("Only IRAK1 overlaps the curated tolerance marker set", loc="left", fontsize=9)
    for i, gene in enumerate(tol_overlap["gene"]):
        row = tol_overlap.iloc[i]
        if bool(row["in_frozen_signature"]):
            module_label = module_short(str(row["frozen_module"]))
            direction_label = "down arm" if str(row["frozen_direction"]) == "down" else str(row["frozen_direction"])
            ax1.text(0.62, i, f"{module_label} / {direction_label}", va="center", ha="left", fontsize=6, color=BLACK)
    finish_axis(ax1)
    plt.colorbar(im, ax=ax1, fraction=0.03, pad=0.02).ax.tick_params(labelsize=6, length=2)

    tol_external = tol_external.copy()
    tol_external["consistent"] = (tol_external["consistent_with_tolerance"] == "YES").astype(int)
    agg = tol_external.groupby("dataset", as_index=False)["consistent"].mean()
    ax2.bar(agg["dataset"], agg["consistent"], color=[GREEN, ORANGE, PURPLE], edgecolor=BLACK, linewidth=0.4)
    ax2.axhline(0.5, color=LIGHT_GREY, linestyle="--", linewidth=0.8)
    ax2.set_ylim(0, 1.0)
    ax2.set_ylabel("Mean tolerance-marker consistency")
    ax2.set_title("Tolerance consistency does not track the EV signature coherently", loc="left", fontsize=9)
    finish_axis(ax2)

    ax3.bar(pdac["coarse_group"], pdac["mean"], color=LIGHT_BLUE, edgecolor=BLACK, linewidth=0.4, label="PDAC myeloid cohorts")
    ax3.plot(range(len(macrophage)), macrophage["mean"], color=VERMILLION, marker="o", linewidth=1.4, label="Macrophage polarization")
    ax3.set_xticks(range(len(macrophage)))
    ax3.set_xticklabels(["Dendritic/\nM2-like", "Macrophage/\nM0-like", "Monocyte/\nM1-like"])
    ax3.set_ylabel("Mean Core EV score")
    ax3.set_title("Comparator systems remain contextual rather than mechanistic", loc="left", fontsize=9)
    ax3.legend(frameon=False, fontsize=6, loc="upper right")
    finish_axis(ax3)

    save_figure(
        fig,
        SUPP_DIR / "Figure_S5_tolerance_and_comparator_context",
        panel_axes={"A": ax1, "B": ax2, "C": ax3},
    )


def main() -> None:
    apply_publication_style("nature")
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    supplement_s1()
    supplement_s2()
    supplement_s3()
    supplement_s4()
    supplement_s5()


if __name__ == "__main__":
    main()
