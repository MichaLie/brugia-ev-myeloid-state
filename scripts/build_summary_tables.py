#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
OUT_DIR = RESULTS_DIR / "tables"
VALIDATION_DIR = RESULTS_DIR / "validation"
DOCS_DIR = ROOT / "docs"


def build_orthogonal_validation_table() -> pd.DataFrame:
    tlr = pd.read_csv(DATA_DIR / "evidence" / "tlr_projection" / "tlr78_concordance_summary.csv")
    gse187403 = pd.read_csv(VALIDATION_DIR / "gse187403_concordance_summary.csv")
    gse187403_scores = pd.read_csv(VALIDATION_DIR / "gse187403_paired_score_tests.csv")
    gse143170 = pd.read_csv(OUT_DIR / "gse143170_paired_deseq2_concordance_summary.csv")
    gse143170_scores = pd.read_csv(VALIDATION_DIR / "gse143170_paired_score_tests.csv")

    rows = []

    for stimulus in ("TLR7/8", "TLR4"):
        a = tlr[(tlr["stimulus"].eq(stimulus)) & (tlr["cell_type"].eq("mono")) & (tlr["module"].eq("A"))].iloc[0]
        b = tlr[(tlr["stimulus"].eq(stimulus)) & (tlr["cell_type"].eq("mono")) & (tlr["module"].eq("B"))].iloc[0]
        rows.append(
            {
                "analysis_layer": "tlr_context",
                "dataset": "GSE147310",
                "contrast": f"{stimulus} monocytes vs unstimulated",
                "full_signature_concordance": "",
                "down_or_module_a_concordance": f"{int(a['concordant'])}/{int(a['evaluable'])} ({a['concordance_fraction']:.3f})",
                "up_or_module_b_concordance": f"{int(b['concordant'])}/{int(b['evaluable'])} ({b['concordance_fraction']:.3f})",
                "module_b_score_delta": "",
                "module_b_score_p": "",
                "note": "Module A strongly concordant while Module B is discordant under acute innate stimulation.",
            }
        )

    for contrast_name, label in [
        ("lps_vs_media", "LPS vs media"),
        ("mtori_lps_vs_lps", "AZD2014+LPS vs LPS"),
        ("mtori_lps_vs_media", "AZD2014+LPS vs media"),
    ]:
        full = gse187403[(gse187403["contrast"].eq(contrast_name)) & (gse187403["group"].eq("full_signature"))].iloc[0]
        down = gse187403[(gse187403["contrast"].eq(contrast_name)) & (gse187403["group"].eq("down_arm"))].iloc[0]
        modb = gse187403[(gse187403["contrast"].eq(contrast_name)) & (gse187403["group"].eq("module_b"))].iloc[0]
        score = gse187403_scores[
            (gse187403_scores["contrast"].eq(contrast_name))
            & (gse187403_scores["score_name"].eq("module_b_ev_score"))
        ].iloc[0]
        rows.append(
            {
                "analysis_layer": "orthogonal_perturbation",
                "dataset": "GSE187403",
                "contrast": label,
                "full_signature_concordance": f"{int(full['n_concordant'])}/{int(full['n_available'])} ({full['concordance_fraction']:.3f})",
                "down_or_module_a_concordance": f"{int(down['n_concordant'])}/{int(down['n_available'])} ({down['concordance_fraction']:.3f})",
                "up_or_module_b_concordance": f"{int(modb['n_concordant'])}/{int(modb['n_available'])} ({modb['concordance_fraction']:.3f})",
                "module_b_score_delta": f"{score['mean_delta']:.3f}",
                "module_b_score_p": f"{score['signflip_p']:.5f}",
                "note": "Primary monocytes; acute catalytic mTOR inhibition is tested only within the LPS setting.",
            }
        )

    for contrast_name, label in [
        ("tems_vs_ph73", "Temsirolimus vs pH7.3"),
        ("ph65_vs_ph73", "Acidic pH vs pH7.3"),
    ]:
        full = gse143170[(gse143170["contrast"].eq(contrast_name)) & (gse143170["group"].eq("full_signature"))].iloc[0]
        down = gse143170[(gse143170["contrast"].eq(contrast_name)) & (gse143170["group"].eq("down_arm"))].iloc[0]
        modb = gse143170[(gse143170["contrast"].eq(contrast_name)) & (gse143170["group"].eq("module_b"))].iloc[0]
        score = gse143170_scores[
            (gse143170_scores["contrast"].eq(contrast_name))
            & (gse143170_scores["score_name"].eq("module_b_ev_score"))
        ].iloc[0]
        rows.append(
            {
                "analysis_layer": "orthogonal_perturbation",
                "dataset": "GSE143170",
                "contrast": label,
                "full_signature_concordance": f"{int(full['n_concordant'])}/{int(full['n_available'])} ({full['concordance_fraction']:.3f})",
                "down_or_module_a_concordance": f"{int(down['n_concordant'])}/{int(down['n_available'])} ({down['concordance_fraction']:.3f})",
                "up_or_module_b_concordance": f"{int(modb['n_concordant'])}/{int(modb['n_available'])} ({modb['concordance_fraction']:.3f})",
                "module_b_score_delta": f"{score['mean_delta']:.3f}",
                "module_b_score_p": f"{score['signflip_p']:.5f}",
                "note": "Paired DESeq2 rerun in differentiating moDCs; acidic pH serves as a negative control.",
            }
        )

    return pd.DataFrame(rows)


def build_external_transport_contrast_table() -> pd.DataFrame:
    overview = pd.read_csv(DATA_DIR / "evidence" / "external" / "external_contrast_overview.csv")
    gse250463 = pd.read_csv(DATA_DIR / "evidence" / "external" / "gse250463_contrast_results.csv")
    gse42694 = pd.read_csv(DATA_DIR / "evidence" / "external" / "gse42694_contrast_results.csv")

    base = overview.copy()
    base["context"] = base["dataset"].map(
        {
            "GSE360": "Direct human Brugia exposure",
            "GSE2135": "Natural human filarial infection",
        }
    )
    base["n_samples_case"] = base["n_case"]
    base["n_samples_control"] = base["n_control"]

    extra_250463 = gse250463.copy()
    extra_250463["context"] = extra_250463["contrast"].map(
        lambda x: "Cross-helminth moDC LPS background" if "LPS" in x else "Cross-helminth moDC resting"
    )
    extra_250463["n_samples_case"] = extra_250463["n_pairs"]
    extra_250463["n_samples_control"] = extra_250463["n_pairs"]

    extra_42694 = gse42694.copy()
    extra_42694["context"] = "Same-parasite L3 Langerhans exposure"
    extra_42694["n_samples_case"] = extra_42694["n_pairs"]
    extra_42694["n_samples_control"] = extra_42694["n_pairs"]

    combined = pd.concat(
        [
            base[
                [
                    "dataset",
                    "contrast",
                    "context",
                    "n_samples_case",
                    "n_samples_control",
                    "core_mean_difference",
                    "core_score_p_value",
                    "available_up_genes",
                    "available_down_genes",
                    "up_concordant",
                    "down_concordant",
                    "total_available_genes",
                    "total_concordant_genes",
                    "concordance_fraction",
                    "direction_empirical_p_value",
                ]
            ],
            extra_250463[
                [
                    "dataset",
                    "contrast",
                    "context",
                    "n_samples_case",
                    "n_samples_control",
                    "core_mean_difference",
                    "core_score_p_value",
                    "available_up_genes",
                    "available_down_genes",
                    "up_concordant",
                    "down_concordant",
                    "total_available_genes",
                    "total_concordant_genes",
                    "concordance_fraction",
                    "direction_empirical_p_value",
                ]
            ],
            extra_42694[
                [
                    "dataset",
                    "contrast",
                    "context",
                    "n_samples_case",
                    "n_samples_control",
                    "core_mean_difference",
                    "core_score_p_value",
                    "available_up_genes",
                    "available_down_genes",
                    "up_concordant",
                    "down_concordant",
                    "total_available_genes",
                    "total_concordant_genes",
                    "concordance_fraction",
                    "direction_empirical_p_value",
                ]
            ],
        ],
        ignore_index=True,
    )

    combined["up_fraction"] = combined["up_concordant"] / combined["available_up_genes"]
    combined["down_fraction"] = combined["down_concordant"] / combined["available_down_genes"]
    combined["arm_gap_down_minus_up"] = combined["down_fraction"] - combined["up_fraction"]
    return combined.sort_values(["dataset", "contrast"]).reset_index(drop=True)


def build_external_transport_context_table() -> pd.DataFrame:
    return pd.read_csv(DATA_DIR / "evidence" / "external" / "external_arm_transport_context_summary.csv")


def build_comparator_summary_table() -> pd.DataFrame:
    pdac_h = pd.read_csv(DATA_DIR / "evidence" / "comparators" / "pdac_cross_cohort_hallmark_concordance.csv")
    pdac_m = pd.read_csv(DATA_DIR / "evidence" / "comparators" / "pdac_cross_cohort_marker_concordance.csv")
    pdac_geom = pd.read_csv(DATA_DIR / "evidence" / "comparators" / "pdac_cross_cohort_state_geometry.csv")
    mac_mark = pd.read_csv(DATA_DIR / "evidence" / "comparators" / "macrophage_marker_correlations.csv")
    mac_hall = pd.read_csv(DATA_DIR / "evidence" / "comparators" / "macrophage_hallmark_correlations.csv")
    mac_pair = pd.read_csv(DATA_DIR / "evidence" / "comparators" / "macrophage_paired_tests.csv")
    mac_core = pd.read_csv(DATA_DIR / "evidence" / "comparators" / "macrophage_core_by_polarization.csv")

    pdac_h_rho = spearmanr(pdac_h["spearman_rho_task36"], pdac_h["spearman_rho_gse111672"]).statistic
    pdac_m_rho = spearmanr(pdac_m["spearman_rho_task14"], pdac_m["spearman_rho_gse111672"]).statistic

    geom_task36 = pdac_geom[pdac_geom["cohort"].eq("GSE155698_task36")].sort_values("mean", ascending=False)
    geom_111672 = pdac_geom[pdac_geom["cohort"].eq("GSE111672")].sort_values("mean", ascending=False)
    mac_order = mac_core.sort_values("mean", ascending=False)
    myc_row = mac_mark[mac_mark["gene"].eq("MYC")].iloc[0]
    stat3_row = mac_mark[mac_mark["gene"].eq("STAT3")].iloc[0]
    il6_row = mac_hall[mac_hall["pathway"].eq("HALLMARK_IL6_JAK_STAT3_SIGNALING")].iloc[0]
    mtorc1_row = mac_hall[mac_hall["pathway"].eq("HALLMARK_MTORC1_SIGNALING")].iloc[0]
    m2m1_row = mac_pair[mac_pair["contrast"].eq("M2_minus_M1")].iloc[0]

    rows = [
        {
            "comparator": "PDAC",
            "summary_metric": "hallmark_vector_rho",
            "value": f"{pdac_h_rho:.4f}",
            "supportive_direction": "positive",
            "interpretation": "Strong cross-cohort pathway-level concordance.",
        },
        {
            "comparator": "PDAC",
            "summary_metric": "marker_vector_rho",
            "value": f"{pdac_m_rho:.4f}",
            "supportive_direction": "weak",
            "interpretation": "Marker-level concordance is weak and should not be mechanistically over-read.",
        },
        {
            "comparator": "PDAC",
            "summary_metric": "state_geometry_order",
            "value": (
                f"GSE155698_task36: {' > '.join(geom_task36['coarse_group'])}; "
                f"GSE111672: {' > '.join(geom_111672['coarse_group'])}"
            ),
            "supportive_direction": "dendritic_high",
            "interpretation": "Broad myeloid state geometry is consistent across cohorts.",
        },
        {
            "comparator": "Macrophage",
            "summary_metric": "core_score_order",
            "value": " > ".join(mac_order["polarization_state"]) + " by mean score",
            "supportive_direction": "M2_adjacent",
            "interpretation": "EV state is non-M1 and M2-adjacent, but not reducible to a canonical polarization label.",
        },
        {
            "comparator": "Macrophage",
            "summary_metric": "MYC_marker_rho",
            "value": f"rho={myc_row['spearman_rho']:.4f}; p={myc_row['p_value']:.4g}; fdr={myc_row['fdr']:.4g}",
            "supportive_direction": "positive",
            "interpretation": "MYC is the strongest positive marker-level correlate.",
        },
        {
            "comparator": "Macrophage",
            "summary_metric": "MTORC1_hallmark_rho",
            "value": f"rho={mtorc1_row['spearman_rho']:.4f}; fdr={mtorc1_row['fdr']:.4g}",
            "supportive_direction": "positive",
            "interpretation": "Metabolic and growth-associated hallmark structure is the strongest shared layer.",
        },
        {
            "comparator": "Macrophage",
            "summary_metric": "IL6_JAK_STAT3_hallmark_rho",
            "value": f"rho={il6_row['spearman_rho']:.4f}; p={il6_row['p_value']:.4g}",
            "supportive_direction": "null",
            "interpretation": "STAT3-linked inflammatory structure is not a stable organizing feature.",
        },
        {
            "comparator": "Macrophage",
            "summary_metric": "STAT3_marker_rho",
            "value": f"rho={stat3_row['spearman_rho']:.4f}; p={stat3_row['p_value']:.4g}",
            "supportive_direction": "null",
            "interpretation": "STAT3 expression does not provide a stable marker-level bridge.",
        },
        {
            "comparator": "Macrophage",
            "summary_metric": "M2_minus_M1_paired_test",
            "value": f"mean_delta={m2m1_row['mean_delta']:.4f}; p={m2m1_row['p_value']:.4g}; n={int(m2m1_row['n_donors'])}",
            "supportive_direction": "positive",
            "interpretation": "Ordering is directionally consistent across donors but limited by n=4.",
        },
    ]

    return pd.DataFrame(rows)


def build_manifest() -> str:
    return "\n".join(
        [
            "# Derived Table Manifest",
            "",
            "## Included tables",
            "- `claim_support_stats.csv`",
            "  Suggested caption: Claim-facing statistical anchors for the main text and extended materials.",
            "- `orthogonal_validation_summary.csv`",
            "  Suggested caption: Summary of TLR-context projection and orthogonal perturbation testing across GSE147310, GSE187403, and GSE143170.",
            "- `gse143170_paired_deseq2_concordance_summary.csv`",
            "  Suggested caption: Paired DESeq2 concordance summary for GSE143170 temsirolimus and acidic-pH contrasts.",
            "- `gse143170_paired_deseq2_crosscheck.csv`",
            "  Suggested caption: Crosscheck between the paired DESeq2 rerun and the normalized-count sensitivity rerun for GSE143170.",
            "- `gse143170_paired_deseq2_gene_level.csv`",
            "  Suggested caption: Gene-level paired DESeq2 output for fixed-signature genes in GSE143170.",
            "- `external_transport_contrast_summary.csv`",
            "  Suggested caption: Contrast-level external projection summary across direct Brugia exposure, natural infection, same-parasite stress-test, and cross-helminth datasets.",
            "- `external_transport_context_summary.csv`",
            "  Suggested caption: Context-level external-projection summary showing differential down-arm and up-arm behavior across external settings.",
            "- `comparator_state_summary.csv`",
            "  Suggested caption: Comparator-state summary for PDAC myeloid cohorts and canonical macrophage polarization, emphasizing pathway-level neighborhood rather than shared mechanism.",
            "",
            "## External support rule",
            "- The preserved `overall_supportive` flag for external contrasts requires total directional concordance >= 0.60 and empirical direction-test P < 0.05.",
            "- The context-level `n_supportive_contrasts` field counts contrasts satisfying that fixed rule.",
            "",
            "## Associated summaries",
            "- `gse143170_paired_deseq2_summary.md`",
            "  Suggested use: human-readable summary of the paired DESeq2 rerun.",
            "",
            "## Repository use",
            "- These tables are derived summaries built from the repository evidence and perturbation-output layers.",
            "- They are suitable for supplements, appendices, or public repository browsing without publication-package files.",
        ]
    ) + "\n"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    build_orthogonal_validation_table().to_csv(OUT_DIR / "orthogonal_validation_summary.csv", index=False)
    build_external_transport_contrast_table().to_csv(OUT_DIR / "external_transport_contrast_summary.csv", index=False)
    build_external_transport_context_table().to_csv(OUT_DIR / "external_transport_context_summary.csv", index=False)
    build_comparator_summary_table().to_csv(OUT_DIR / "comparator_state_summary.csv", index=False)
    (DOCS_DIR / "derived_table_manifest.md").write_text(build_manifest(), encoding="utf-8")


if __name__ == "__main__":
    main()
