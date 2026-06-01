#!/usr/bin/env python3
"""Regenerate reference-style robustness supplementary figures for v5 selected death-mode.

Outputs:
- Supplementary_Figure_S2_hub_rank_stability_across_pair_count_thresholds_selected_death_modes_2000_2025
- Supplementary_Figure_S3_observed_vs_null_selected_death_modes_2000_2025

The visual layouts and metric content follow the Cell_Death_Select S2/S3
reference figures, while all values are recomputed from the v5 selected death-mode
2000-2025 pair-metric table.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(os.environ.get("CELLDEATH_ATLAS_ROOT", Path(__file__).resolve().parents[2]))
V5 = Path(os.environ.get("CELLDEATH_ATLAS_PACKAGE_DIR", ROOT / os.environ.get("CELLDEATH_ATLAS_PACKAGE", "selected_death_modes_2000_2025")))
SRC = V5 / "04_supplementary_figures_regenerated" / "source_data"
OUT_ROOT = V5 / "04_supplementary_figures_regenerated"
SCRIPT_NAME = Path(__file__).name
THRESHOLDS = [1, 3, 5, 10]
TOP_NS = [10, 25, 50]
PERMUTATION_ITERATIONS = 200
SEED = 20260513


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.frameon": False,
    }
)


def ensure_dirs() -> None:
    for subdir in ["pdf", "png", "tiff", "source_data"]:
        (OUT_ROOT / subdir).mkdir(parents=True, exist_ok=True)


def save_all(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT_ROOT / "pdf" / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT_ROOT / "png" / f"{stem}.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT_ROOT / "tiff" / f"{stem}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def require_csv(name: str) -> pd.DataFrame:
    path = SRC / name
    if not path.exists():
        raise SystemExit(f"Missing required source data: {path}")
    return pd.read_csv(path)


def load_pair_data() -> pd.DataFrame:
    path = V5 / "01_source_data_filtered" / "pair_metrics_selected_death_modes_2000_2025.csv"
    if not path.exists():
        path = SRC / "pair_metrics_selected_death_modes_2000_2025.csv"
    if not path.exists():
        raise SystemExit(f"Missing selected death-mode pair metrics source: {path}")
    pair = pd.read_csv(path)
    required = {
        "Death_Mode",
        "Disease",
        "Pair_Count",
        "Disease_Category",
        "Mechanistic_Convergence_Score",
        "Disease_Specificity_Score",
        "Gap_Score",
    }
    missing = required - set(pair.columns)
    if missing:
        raise SystemExit(f"pair metrics missing columns: {sorted(missing)}")
    if "Pair_ID" not in pair.columns:
        pair["Pair_ID"] = pair["Death_Mode"].astype(str) + "|||" + pair["Disease"].astype(str)
    numeric_cols = [
        "Pair_Count",
        "Mechanistic_Convergence_Score",
        "Disease_Specificity_Score",
        "Gap_Score",
    ]
    for col in numeric_cols:
        pair[col] = pd.to_numeric(pair[col], errors="coerce").fillna(0.0)
    pair["Disease_Category"] = pair["Disease_Category"].fillna("Other or mixed")
    return pair


def zscore(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce").fillna(0.0)
    sd = values.std(ddof=0)
    if sd == 0 or pd.isna(sd):
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - values.mean()) / sd


def rank_correlation(ref: pd.Series, comp: pd.Series, method: str) -> float:
    aligned = pd.concat([ref, comp], axis=1, join="inner").dropna()
    if len(aligned) < 3:
        return np.nan
    if aligned.iloc[:, 0].nunique() < 2 or aligned.iloc[:, 1].nunique() < 2:
        return np.nan
    return float(aligned.iloc[:, 0].corr(aligned.iloc[:, 1], method=method))


def top_jaccard(ref_order: pd.Series, comp_order: pd.Series, n: int) -> float:
    ref_top = set(ref_order.sort_values(ascending=False).head(n).index)
    comp_top = set(comp_order.sort_values(ascending=False).head(n).index)
    union = ref_top | comp_top
    return len(ref_top & comp_top) / len(union) if union else np.nan


def subset_by_threshold(pair: pd.DataFrame, threshold: int) -> pd.DataFrame:
    return pair[pair["Pair_Count"] >= threshold].copy()


def compute_disease_hubs(sub: pd.DataFrame) -> pd.DataFrame:
    if sub.empty:
        return pd.DataFrame(
            columns=[
                "Disease",
                "Disease_Category",
                "death_mode_n",
                "pair_n",
                "pair_count_sum",
                "mean_cleaned_mcs",
                "mean_dss",
                "mean_frontier",
                "literature_hub_score",
                "literature_hub_rank",
            ]
        )
    grouped = (
        sub.groupby(["Disease", "Disease_Category"], dropna=False)
        .agg(
            death_mode_n=("Death_Mode", "nunique"),
            pair_n=("Pair_ID", "nunique"),
            pair_count_sum=("Pair_Count", "sum"),
            mean_cleaned_mcs=("Mechanistic_Convergence_Score", "mean"),
            mean_dss=("Disease_Specificity_Score", "mean"),
            mean_frontier=("Gap_Score", "mean"),
        )
        .reset_index()
    )
    grouped["literature_hub_score"] = (
        zscore(np.log1p(grouped["pair_count_sum"]))
        + zscore(grouped["death_mode_n"])
        + zscore(grouped["mean_cleaned_mcs"])
        + zscore(grouped["mean_dss"])
    )
    grouped = grouped.sort_values("literature_hub_score", ascending=False).reset_index(drop=True)
    grouped["literature_hub_rank"] = np.arange(1, len(grouped) + 1)
    return grouped


def compute_system_composition(sub: pd.DataFrame) -> pd.DataFrame:
    if sub.empty:
        return pd.DataFrame(columns=["Disease_Category", "pair_n", "disease_n", "pair_count_sum", "volume_share", "system_rank_score"])
    grouped = (
        sub.groupby("Disease_Category", dropna=False)
        .agg(
            pair_n=("Pair_ID", "nunique"),
            disease_n=("Disease", "nunique"),
            pair_count_sum=("Pair_Count", "sum"),
        )
        .reset_index()
    )
    grouped["volume_share"] = grouped["pair_count_sum"] / max(grouped["pair_count_sum"].sum(), 1)
    grouped["system_rank_score"] = grouped["volume_share"]
    return grouped.sort_values("system_rank_score", ascending=False).reset_index(drop=True)


def series_for_metric(details: dict[int, dict[str, pd.DataFrame]], threshold: int, metric: str) -> pd.Series:
    if metric == "hub_score":
        hubs = details[threshold]["hubs"]
        return hubs.set_index("Disease")["literature_hub_score"] if not hubs.empty else pd.Series(dtype=float)
    if metric == "cleaned_mcs":
        sub = details[threshold]["pairs"]
        return sub.set_index("Pair_ID")["Mechanistic_Convergence_Score"]
    if metric == "dss":
        sub = details[threshold]["pairs"]
        return sub.set_index("Pair_ID")["Disease_Specificity_Score"]
    if metric == "frontier":
        sub = details[threshold]["pairs"]
        return sub.set_index("Pair_ID")["Gap_Score"]
    if metric == "disease_system_enrichment":
        systems = details[threshold]["systems"]
        return systems.set_index("Disease_Category")["system_rank_score"] if not systems.empty else pd.Series(dtype=float)
    raise ValueError(metric)


def compute_rank_stability_tables(pair: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    details: dict[int, dict[str, pd.DataFrame]] = {}
    for threshold in THRESHOLDS:
        sub = subset_by_threshold(pair, threshold)
        details[threshold] = {
            "pairs": sub,
            "hubs": compute_disease_hubs(sub),
            "systems": compute_system_composition(sub),
        }

    hub_rows = []
    metric_rows = []
    for metric in ["hub_score", "cleaned_mcs", "dss", "frontier", "disease_system_enrichment"]:
        ref = series_for_metric(details, 1, metric)
        for threshold in THRESHOLDS:
            comp = series_for_metric(details, threshold, metric)
            row = {
                "metric_name": metric,
                "reference_threshold": 1,
                "comparison_threshold": threshold,
                "intersection_entity_n": int(len(ref.index.intersection(comp.index))),
                "spearman_correlation": rank_correlation(ref, comp, "spearman"),
                "kendall_tau": rank_correlation(ref, comp, "kendall"),
            }
            for n in TOP_NS:
                row[f"top_{n}_jaccard"] = top_jaccard(ref, comp, n)
            if metric == "hub_score":
                hub_rows.append(row)
            else:
                metric_rows.append(row)
    return pd.DataFrame(hub_rows), pd.DataFrame(metric_rows)


def scalar_structure_metrics(sub: pd.DataFrame) -> dict[str, float]:
    hubs = compute_disease_hubs(sub)
    systems = compute_system_composition(sub)
    system_var = float(systems["volume_share"].var(ddof=0)) if not systems.empty else 0.0
    return {
        "mean_cleaned_MCS": float(sub["Mechanistic_Convergence_Score"].mean()) if not sub.empty else np.nan,
        "max_literature_hub_score": float(hubs["literature_hub_score"].max()) if not hubs.empty else np.nan,
        "mean_literature_hub_score_top25": float(hubs.head(25)["literature_hub_score"].mean()) if not hubs.empty else np.nan,
        "disease_system_enrichment_variance": system_var,
        "mean_frontier_score": float(sub["Gap_Score"].mean()) if not sub.empty else np.nan,
        "max_frontier_score": float(sub["Gap_Score"].max()) if not sub.empty else np.nan,
    }


def compute_permutation_null(pair: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    base = subset_by_threshold(pair, 1).reset_index(drop=True)
    observed = scalar_structure_metrics(base)
    records = []
    for null_model in ["shuffle_disease_labels", "shuffle_death_mode_labels"]:
        for iteration in range(PERMUTATION_ITERATIONS):
            perm = base.copy()
            if null_model == "shuffle_disease_labels":
                idx = rng.permutation(len(perm))
                perm["Disease"] = base["Disease"].iloc[idx].to_numpy()
                perm["Disease_Category"] = base["Disease_Category"].iloc[idx].to_numpy()
            else:
                perm["Death_Mode"] = base["Death_Mode"].iloc[rng.permutation(len(perm))].to_numpy()
            metrics = scalar_structure_metrics(perm)
            for metric_name, value in metrics.items():
                records.append(
                    {
                        "null_model": null_model,
                        "iteration": iteration + 1,
                        "metric_name": metric_name,
                        "null_value": value,
                    }
                )
    null_long = pd.DataFrame(records)
    summary_rows = []
    for (null_model, metric_name), sub in null_long.groupby(["null_model", "metric_name"], dropna=False):
        vals = pd.to_numeric(sub["null_value"], errors="coerce")
        obs = observed.get(metric_name, np.nan)
        summary_rows.append(
            {
                "null_model": null_model,
                "metric_name": metric_name,
                "observed_value": obs,
                "null_mean": float(vals.mean()),
                "null_sd": float(vals.std(ddof=1)) if len(vals) > 1 else 0.0,
                "null_ci_lower": float(vals.quantile(0.025)),
                "null_ci_upper": float(vals.quantile(0.975)),
                "observed_minus_null_mean": float(obs - vals.mean()) if pd.notna(obs) else np.nan,
                "permutation_iterations": PERMUTATION_ITERATIONS,
                "interpretation_note": "Exploratory null control only; do not claim formal statistical significance.",
            }
        )
    return pd.DataFrame(summary_rows), null_long


def pretty_number(value: float) -> str:
    if pd.isna(value):
        return ""
    if abs(value - 1.0) < 1e-9:
        return "1"
    if abs(value) < 0.005 and value != 0:
        return f"{value:.2g}"
    return f"{value:.2f}".rstrip("0").rstrip(".")


def heatmap_with_text(
    ax: plt.Axes,
    data: pd.DataFrame,
    title: str,
    cbar_label: str,
    cmap: str,
    tick_fontsize: float = 8,
) -> None:
    sns.heatmap(
        data,
        ax=ax,
        cmap=cmap,
        vmin=0,
        vmax=1,
        linewidths=0.7,
        linecolor="white",
        cbar=True,
        cbar_kws={"label": cbar_label, "shrink": 0.62, "pad": 0.05},
        annot=False,
        square=False,
    )
    for row_i in range(data.shape[0]):
        for col_i in range(data.shape[1]):
            value = data.iloc[row_i, col_i]
            if pd.notna(value):
                ax.text(
                    col_i + 0.5,
                    row_i + 0.5,
                    pretty_number(float(value)),
                    ha="center",
                    va="center",
                    fontsize=8.2,
                    color="black",
                )
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold", pad=14)
    ax.set_xlabel("Comparison threshold", fontsize=9, fontweight="bold")
    ax.set_ylabel("")
    ax.tick_params(axis="both", length=0, labelsize=tick_fontsize)
    ax.set_xticklabels([str(x.get_text()) for x in ax.get_xticklabels()], rotation=0)
    ax.set_yticklabels([str(y.get_text()) for y in ax.get_yticklabels()], rotation=0)


def regenerate_s2() -> None:
    pair = load_pair_data()
    hub, metric = compute_rank_stability_tables(pair)
    hub.to_csv(SRC / "hub_rank_stability_selected_death_modes_2000_2025.csv", index=False)
    metric.to_csv(SRC / "metric_rank_stability_selected_death_modes_2000_2025.csv", index=False)

    hub_matrix = (
        hub.pivot_table(
            index="metric_name",
            columns="comparison_threshold",
            values="spearman_correlation",
            aggfunc="mean",
        )
        .reindex(index=["hub_score"])
        .reindex(columns=THRESHOLDS)
    )

    metric_order = [
        m
        for m in ["frontier", "dss", "cleaned_mcs", "disease_system_enrichment"]
        if m in set(metric["metric_name"])
    ]
    metric_matrix = (
        metric.pivot_table(
            index="metric_name",
            columns="comparison_threshold",
            values="top_25_jaccard",
            aggfunc="mean",
        )
        .reindex(index=metric_order)
        .reindex(columns=THRESHOLDS)
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(9.0, 6.5),
        gridspec_kw={"width_ratios": [1.0, 1.35], "wspace": 0.72},
    )
    heatmap_with_text(
        axes[0],
        hub_matrix,
        "A. Hub-rank stability",
        "Spearman",
        "viridis",
        tick_fontsize=8.5,
    )
    heatmap_with_text(
        axes[1],
        metric_matrix,
        "B. Metric top-rank overlap",
        "Top-25 Jaccard",
        "magma",
        tick_fontsize=8.5,
    )
    fig.text(
        0.02,
        0.025,
        "Reference-style selected death-mode robustness diagnostic. Thresholds are Pair_Count >= 1, 3, 5, and 10.",
        fontsize=6.3,
        color="#52616B",
    )
    stem = "Supplementary_Figure_S2_hub_rank_stability_across_pair_count_thresholds_selected_death_modes_2000_2025"
    save_all(fig, stem)

    qc = pd.DataFrame(
        [
            {
                "figure": stem,
                "source_hub_rows": len(hub),
                "source_metric_rows": len(metric),
                "thresholds": "|".join(map(str, THRESHOLDS)),
                "metric_rows_expected_reference_style": 16,
                "script": SCRIPT_NAME,
            }
        ]
    )
    qc.to_csv(SRC / f"{stem}_QC.csv", index=False)


def regenerate_s3() -> None:
    pair = load_pair_data()
    nulls, null_long = compute_permutation_null(pair)
    nulls.to_csv(SRC / "permutation_null_summary_selected_death_modes_2000_2025.csv", index=False)
    null_long.to_csv(SRC / "permutation_null_long_selected_death_modes_2000_2025.csv", index=False)
    required = {
        "null_model",
        "metric_name",
        "observed_value",
        "null_mean",
        "null_ci_lower",
        "null_ci_upper",
        "permutation_iterations",
    }
    missing = required - set(nulls.columns)
    if missing:
        raise SystemExit(f"permutation_null_summary missing columns: {sorted(missing)}")

    nulls["entity"] = nulls["null_model"].astype(str) + " / " + nulls["metric_name"].astype(str)
    nulls = nulls.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(nulls))

    fig, ax = plt.subplots(figsize=(9.0, 6.5))
    fig.subplots_adjust(left=0.34, right=0.96, top=0.78, bottom=0.14)
    ax.hlines(
        y,
        nulls["null_ci_lower"],
        nulls["null_ci_upper"],
        color="#9AA7B1",
        lw=2.2,
        label="Null 95% interval",
        zorder=1,
    )
    ax.scatter(nulls["null_mean"], y, color="#5B6770", s=38, label="Null mean", zorder=3)
    ax.scatter(nulls["observed_value"], y, color="#B24745", s=42, label="Observed", zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(nulls["entity"], fontsize=7.8)
    ax.set_xlabel("Metric value", fontsize=9, fontweight="bold")
    fig.suptitle(
        "Observed versus shuffled null structure",
        x=0.34,
        y=0.96,
        ha="left",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.34,
        0.90,
        "Diagnostics only; no formal statistical significance is claimed.",
        ha="left",
        va="bottom",
        fontsize=8.4,
        color="#52616B",
    )
    ax.grid(axis="x", color="#E1E5EA", lw=0.8)
    ax.grid(axis="y", color="#E1E5EA", lw=0.8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=7.5)
    ax.margins(x=0.04, y=0.12)

    stem = "Supplementary_Figure_S3_observed_vs_null_selected_death_modes_2000_2025"
    save_all(fig, stem)

    qc = pd.DataFrame(
        [
            {
                "figure": stem,
                "source_null_rows": len(nulls),
                "source_null_long_rows": len(null_long),
                "permutation_iterations_min": int(nulls["permutation_iterations"].min()),
                "permutation_iterations_max": int(nulls["permutation_iterations"].max()),
                "expected_reference_style_rows": 12,
                "script": SCRIPT_NAME,
            }
        ]
    )
    qc.to_csv(SRC / f"{stem}_QC.csv", index=False)


def main() -> int:
    ensure_dirs()
    regenerate_s2()
    regenerate_s3()
    print("Regenerated reference-style selected death-mode S2/S3 robustness figures.")
    print(OUT_ROOT / "pdf" / "Supplementary_Figure_S2_hub_rank_stability_across_pair_count_thresholds_selected_death_modes_2000_2025.pdf")
    print(OUT_ROOT / "pdf" / "Supplementary_Figure_S3_observed_vs_null_selected_death_modes_2000_2025.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
