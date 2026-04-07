#!/usr/bin/env python3

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
from scipy.stats import binomtest, fisher_exact


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
OUT_DIR = RESULTS_DIR / "tables"
VALIDATION_DIR = RESULTS_DIR / "validation"
DOCS_DIR = ROOT / "docs"


def fisher_greater(concordant_a: int, discordant_a: int, concordant_b: int, discordant_b: int) -> tuple[float, float]:
    odds_ratio, p_value = fisher_exact(
        [[concordant_a, discordant_a], [concordant_b, discordant_b]],
        alternative="greater",
    )
    return float(odds_ratio), float(p_value)


def format_fraction(concordant: int, total: int) -> str:
    return f"{concordant}/{total} ({concordant / total:.3f})"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DOCS_DIR.mkdir(parents=True, exist_ok=True)

    tlr = pd.read_csv(DATA_DIR / "evidence" / "tlr_projection" / "tlr78_concordance_summary.csv")
    gse187403 = pd.read_csv(VALIDATION_DIR / "gse187403_concordance_summary.csv")
    gse187403_scores = pd.read_csv(VALIDATION_DIR / "gse187403_paired_score_tests.csv")
    gse143170 = pd.read_csv(OUT_DIR / "gse143170_paired_deseq2_concordance_summary.csv")
    gse143170_scores = pd.read_csv(VALIDATION_DIR / "gse143170_paired_score_tests.csv")
    external_overview = pd.read_csv(DATA_DIR / "evidence" / "external" / "external_contrast_overview.csv")
    external_context = pd.read_csv(DATA_DIR / "evidence" / "external" / "external_arm_transport_context_summary.csv")
    gse42694_results = pd.read_csv(DATA_DIR / "evidence" / "external" / "gse42694_contrast_results.csv")
    gse250463_results = pd.read_csv(DATA_DIR / "evidence" / "external" / "gse250463_contrast_results.csv")

    rows = []

    mono_a = tlr[(tlr["stimulus"].eq("TLR7/8")) & (tlr["cell_type"].eq("mono")) & (tlr["module"].eq("A"))].iloc[0]
    mono_b = tlr[(tlr["stimulus"].eq("TLR7/8")) & (tlr["cell_type"].eq("mono")) & (tlr["module"].eq("B"))].iloc[0]
    mono_or, mono_p = fisher_greater(
        int(mono_a["concordant"]),
        int(mono_a["discordant"]),
        int(mono_b["concordant"]),
        int(mono_b["discordant"]),
    )
    rows.append(
        {
            "claim_code": "tlr78_mono_module_a_vs_b",
            "claim_area": "mechanistic_context",
            "dataset": "GSE147310",
            "summary": (
                f"TLR7/8 monocytes: Module A {format_fraction(int(mono_a['concordant']), int(mono_a['evaluable']))} "
                f"vs Module B {format_fraction(int(mono_b['concordant']), int(mono_b['evaluable']))}"
            ),
            "odds_ratio": mono_or,
            "p_value": mono_p,
            "destination": "main_text",
        }
    )

    mono_tlr4_a = tlr[(tlr["stimulus"].eq("TLR4")) & (tlr["cell_type"].eq("mono")) & (tlr["module"].eq("A"))].iloc[0]
    mono_tlr4_b = tlr[(tlr["stimulus"].eq("TLR4")) & (tlr["cell_type"].eq("mono")) & (tlr["module"].eq("B"))].iloc[0]
    tlr4_or, tlr4_p = fisher_greater(
        int(mono_tlr4_a["concordant"]),
        int(mono_tlr4_a["discordant"]),
        int(mono_tlr4_b["concordant"]),
        int(mono_tlr4_b["discordant"]),
    )
    rows.append(
        {
            "claim_code": "tlr4_mono_module_a_vs_b",
            "claim_area": "mechanistic_context",
            "dataset": "GSE147310",
            "summary": (
                f"TLR4 monocytes: Module A {format_fraction(int(mono_tlr4_a['concordant']), int(mono_tlr4_a['evaluable']))} "
                f"vs Module B {format_fraction(int(mono_tlr4_b['concordant']), int(mono_tlr4_b['evaluable']))}"
            ),
            "odds_ratio": tlr4_or,
            "p_value": tlr4_p,
            "destination": "supplement",
        }
    )

    lps_a = gse187403[(gse187403["contrast"].eq("lps_vs_media")) & (gse187403["group"].eq("module_a"))].iloc[0]
    lps_b = gse187403[(gse187403["contrast"].eq("lps_vs_media")) & (gse187403["group"].eq("module_b"))].iloc[0]
    lps_or, lps_p = fisher_greater(
        int(lps_a["n_concordant"]),
        int(lps_a["n_available"] - lps_a["n_concordant"]),
        int(lps_b["n_concordant"]),
        int(lps_b["n_available"] - lps_b["n_concordant"]),
    )
    rows.append(
        {
            "claim_code": "gse187403_lps_module_a_vs_b",
            "claim_area": "orthogonal_validation",
            "dataset": "GSE187403",
            "summary": (
                f"LPS vs media: Module A {format_fraction(int(lps_a['n_concordant']), int(lps_a['n_available']))} "
                f"vs Module B {format_fraction(int(lps_b['n_concordant']), int(lps_b['n_available']))}"
            ),
            "odds_ratio": lps_or,
            "p_value": lps_p,
            "destination": "main_text",
        }
    )

    azd_b = gse187403[(gse187403["contrast"].eq("mtori_lps_vs_lps")) & (gse187403["group"].eq("module_b"))].iloc[0]
    rows.append(
        {
            "claim_code": "gse187403_mtori_module_b_binomial",
            "claim_area": "orthogonal_validation",
            "dataset": "GSE187403",
            "summary": (
                f"AZD2014+LPS vs LPS: Module B {format_fraction(int(azd_b['n_concordant']), int(azd_b['n_available']))}"
            ),
            "odds_ratio": math.nan,
            "p_value": float(binomtest(int(azd_b["n_concordant"]), int(azd_b["n_available"]), 0.5, alternative="greater").pvalue),
            "destination": "main_text",
        }
    )

    azd_modb_score = gse187403_scores[
        (gse187403_scores["contrast"].eq("mtori_lps_vs_lps"))
        & (gse187403_scores["score_name"].eq("module_b_ev_score"))
    ].iloc[0]
    rows.append(
        {
            "claim_code": "gse187403_mtori_module_b_score_shift",
            "claim_area": "orthogonal_validation",
            "dataset": "GSE187403",
            "summary": f"Module B EV-like score delta = {azd_modb_score['mean_delta']:.3f}",
            "odds_ratio": math.nan,
            "p_value": float(azd_modb_score["signflip_p"]),
            "destination": "main_text",
        }
    )

    tems_b = gse143170[(gse143170["contrast"].eq("tems_vs_ph73")) & (gse143170["group"].eq("module_b"))].iloc[0]
    rows.append(
        {
            "claim_code": "gse143170_tems_module_b_binomial",
            "claim_area": "orthogonal_validation",
            "dataset": "GSE143170",
            "summary": f"Temsirolimus vs pH7.3: Module B {format_fraction(int(tems_b['n_concordant']), int(tems_b['n_available']))}",
            "odds_ratio": math.nan,
            "p_value": float(tems_b["binom_p_greater"]),
            "destination": "main_text",
        }
    )

    ph65_b = gse143170[(gse143170["contrast"].eq("ph65_vs_ph73")) & (gse143170["group"].eq("module_b"))].iloc[0]
    rows.append(
        {
            "claim_code": "gse143170_ph65_module_b_binomial",
            "claim_area": "orthogonal_validation",
            "dataset": "GSE143170",
            "summary": f"Acidic pH vs pH7.3: Module B {format_fraction(int(ph65_b['n_concordant']), int(ph65_b['n_available']))}",
            "odds_ratio": math.nan,
            "p_value": float(ph65_b["binom_p_greater"]),
            "destination": "supplement",
        }
    )

    tems_modb_score = gse143170_scores[
        (gse143170_scores["contrast"].eq("tems_vs_ph73"))
        & (gse143170_scores["score_name"].eq("module_b_ev_score"))
    ].iloc[0]
    rows.append(
        {
            "claim_code": "gse143170_tems_module_b_score_shift",
            "claim_area": "orthogonal_validation",
            "dataset": "GSE143170",
            "summary": f"Module B EV-like score delta = {tems_modb_score['mean_delta']:.3f}",
            "odds_ratio": math.nan,
            "p_value": float(tems_modb_score["signflip_p"]),
            "destination": "main_text",
        }
    )

    total_contrasts = int(external_context["n_contrasts"].sum())
    total_supportive = int(external_context["n_supportive_contrasts"].sum())
    rows.append(
        {
            "claim_code": "external_no_uniform_support",
            "claim_area": "external_transport",
            "dataset": "external_multi_cohort",
            "summary": f"Supportive external contrasts: {total_supportive}/{total_contrasts}",
            "odds_ratio": math.nan,
            "p_value": math.nan,
            "destination": "main_text",
        }
    )

    best_gse360 = external_overview[external_overview["dataset"].eq("GSE360")].sort_values(
        "concordance_fraction", ascending=False
    ).iloc[0]
    rows.append(
        {
            "claim_code": "external_gse360_best_directional",
            "claim_area": "external_transport",
            "dataset": "GSE360",
            "summary": (
                f"{best_gse360['contrast']}: {int(best_gse360['total_concordant_genes'])}/"
                f"{int(best_gse360['total_available_genes'])} ({best_gse360['concordance_fraction']:.3f})"
            ),
            "odds_ratio": math.nan,
            "p_value": float(best_gse360["direction_empirical_p_value"]),
            "destination": "supplement",
        }
    )

    best_gse250463 = gse250463_results.sort_values("concordance_fraction", ascending=False).iloc[0]
    rows.append(
        {
            "claim_code": "external_gse250463_best_directional",
            "claim_area": "external_transport",
            "dataset": "GSE250463",
            "summary": (
                f"{best_gse250463['contrast']}: {int(best_gse250463['total_concordant_genes'])}/"
                f"{int(best_gse250463['total_available_genes'])} ({best_gse250463['concordance_fraction']:.3f})"
            ),
            "odds_ratio": math.nan,
            "p_value": float(best_gse250463["direction_empirical_p_value"]),
            "destination": "supplement",
        }
    )

    gse42694 = gse42694_results.iloc[0]
    rows.append(
        {
            "claim_code": "external_gse42694_boundary",
            "claim_area": "external_transport",
            "dataset": "GSE42694",
            "summary": (
                f"{gse42694['contrast']}: {int(gse42694['total_concordant_genes'])}/"
                f"{int(gse42694['total_available_genes'])} ({gse42694['concordance_fraction']:.3f})"
            ),
            "odds_ratio": math.nan,
            "p_value": float(gse42694["direction_empirical_p_value"]),
            "destination": "supplement",
        }
    )

    direct_ctx = external_context[external_context["context"].eq("Direct human Brugia exposure")].iloc[0]
    lps_ctx = external_context[external_context["context"].eq("Cross-helminth moDC LPS background")].iloc[0]
    rows.append(
        {
            "claim_code": "external_context_structure",
            "claim_area": "external_transport",
            "dataset": "external_multi_cohort",
            "summary": (
                f"Direct exposure down-up gap = {direct_ctx['mean_arm_gap_down_minus_up']:.3f}; "
                f"LPS-background gap = {lps_ctx['mean_arm_gap_down_minus_up']:.3f}"
            ),
            "odds_ratio": math.nan,
            "p_value": math.nan,
            "destination": "supplement",
        }
    )

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(OUT_DIR / "claim_support_stats.csv", index=False)

    md_lines = [
        "# Claim Support Summary",
        "",
        "## Primary Narrative Anchors",
    ]
    for _, row in summary_df[summary_df["destination"].eq("main_text")].iterrows():
        pval = "NA" if pd.isna(row["p_value"]) else f"{row['p_value']:.4g}"
        md_lines.append(f"- `{row['claim_code']}`: {row['summary']} (p = {pval})")

    md_lines.extend(
        [
            "",
            "## Extended-Materials Anchors",
        ]
    )
    for _, row in summary_df[summary_df["destination"].eq("supplement")].iterrows():
        pval = "NA" if pd.isna(row["p_value"]) else f"{row['p_value']:.4g}"
        md_lines.append(f"- `{row['claim_code']}`: {row['summary']} (p = {pval})")

    md_lines.extend(
        [
        ]
    )

    (DOCS_DIR / "claim_support_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
