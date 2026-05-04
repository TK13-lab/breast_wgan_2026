#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]

PATIENT_PAIR_FOLDS = (
    ROOT
    / "outputs"
    / "revision"
    / "patient_pair_equal_budget"
    / "tables"
    / "patient_pair_equal_budget_fold_metrics.csv"
)
SCVI_FOLDS = (
    ROOT
    / "outputs"
    / "revision"
    / "latent_patient_pair"
    / "tables"
    / "scvi_latent_patient_pair_fold_metrics.csv"
)
COX_GRID = ROOT / "outputs" / "figures" / "fig05_bulk_projection_survival" / "tables" / "cox_grid_robust.csv"
CINDEX_GRID = ROOT / "outputs" / "figures" / "fig05_bulk_projection_survival" / "tables" / "cindex_grid.csv"

OUT_ROOT = ROOT / "outputs" / "revision" / "statistical_reporting"
OUT_TABLES = OUT_ROOT / "tables"
OUT_MARKDOWN = OUT_ROOT / "REVISION_STATISTICAL_REPORTING.md"

METRICS = ("AUROC", "AUPRC", "F1", "ACC")
METHOD_LABELS = {
    "scvi_latent_lr": "scVI latent + LR",
}


def mean_ci(series: pd.Series) -> tuple[float, float, float, int]:
    values = series.dropna().astype(float)
    n = int(values.shape[0])
    if n == 0:
        return float("nan"), float("nan"), float("nan"), 0
    mean = float(values.mean())
    if n == 1:
        return mean, mean, mean, 1
    se = float(values.std(ddof=1) / math.sqrt(n))
    delta = 1.96 * se
    return mean, mean - delta, mean + delta, n


def summarize_fold_metrics(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for keys, group in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        record = {col: value for col, value in zip(group_cols, keys)}
        n_folds = 0
        for metric in METRICS:
            mean, ci_low, ci_high, n_obs = mean_ci(group[metric])
            record[f"mean_{metric}"] = mean
            record[f"ci95_low_{metric}"] = ci_low
            record[f"ci95_high_{metric}"] = ci_high
            n_folds = max(n_folds, n_obs)
        record["n_folds"] = n_folds
        records.append(record)
    return pd.DataFrame.from_records(records)


def ensure_method_label(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "method_label" not in out.columns and "method" in out.columns:
        out["method_label"] = out["method"].map(METHOD_LABELS).fillna(out["method"])
    return out


def build_sample_size_summary(patient_pair_df: pd.DataFrame, cindex_df: pd.DataFrame) -> pd.DataFrame:
    normals = sorted(patient_pair_df["normal_holdout"].dropna().unique().tolist())
    tumors = sorted(patient_pair_df["tumor_holdout"].dropna().unique().tolist())
    rows: list[dict[str, object]] = [
        {
            "scope": "scRNA_discovery",
            "group": "overall",
            "n": len(normals) + len(tumors),
            "events": pd.NA,
            "detail": f"{len(normals)} normal samples ({', '.join(normals)}) and {len(tumors)} tumor samples ({', '.join(tumors)})",
        }
    ]

    cindex_clean = cindex_df.drop_duplicates(subset=["subtype", "n", "events"]).copy()
    for row in cindex_clean.itertuples(index=False):
        rows.append(
            {
                "scope": "TCGA_BRCA_survival",
                "group": row.subtype,
                "n": int(row.n),
                "events": int(row.events),
                "detail": f"{int(row.events)} events in {int(row.n)} patients",
            }
        )

    return pd.DataFrame(rows)


def build_survival_hr_summary(cox_df: pd.DataFrame, cindex_df: pd.DataFrame) -> pd.DataFrame:
    merged = cox_df.merge(
        cindex_df[["subtype", "lineage", "cindex"]],
        on=["subtype", "lineage"],
        how="left",
    ).copy()
    merged["hr"] = merged["beta"].apply(lambda value: math.exp(value) if pd.notna(value) else pd.NA)
    merged["hr_ci95_low"] = merged.apply(
        lambda row: math.exp(row.beta - 1.96 * row.se) if pd.notna(row.beta) and pd.notna(row.se) else pd.NA,
        axis=1,
    )
    merged["hr_ci95_high"] = merged.apply(
        lambda row: math.exp(row.beta + 1.96 * row.se) if pd.notna(row.beta) and pd.notna(row.se) else pd.NA,
        axis=1,
    )
    return merged[
        [
            "subtype",
            "lineage",
            "n",
            "events",
            "beta",
            "se",
            "hr",
            "hr_ci95_low",
            "hr_ci95_high",
            "cindex",
            "p",
            "FDR",
            "method",
            "sig",
        ]
    ].sort_values(["subtype", "lineage"])


def format_metric_line(df: pd.DataFrame, method_label: str, metric: str) -> str:
    row = df.loc[df["method_label"] == method_label].iloc[0]
    return (
        f"- {method_label} {metric}: {row[f'mean_{metric}']:.3f} "
        f"(95% CI {row[f'ci95_low_{metric}']:.3f}-{row[f'ci95_high_{metric}']:.3f}; "
        f"n = {int(row['n_folds'])} held-out folds)"
    )


def write_markdown(
    overall_df: pd.DataFrame,
    scvi_df: pd.DataFrame | None,
    sample_size_df: pd.DataFrame,
) -> None:
    lines: list[str] = []
    lines.append("# Revision Statistical Reporting Summary")
    lines.append("")
    lines.append("This note consolidates the quantitative additions used to strengthen the reviewer response on statistical reporting.")
    lines.append("")
    lines.append("## Confidence Intervals")
    lines.append("")
    lines.append("Overall patient-pair benchmark means with descriptive 95% Wald intervals across held-out folds:")
    lines.append("")
    lines.append(format_metric_line(overall_df, "GAN fixed top-k + LR", "AUROC"))
    lines.append(format_metric_line(overall_df, "GAN fixed top-k + LR", "AUPRC"))
    lines.append(format_metric_line(overall_df, "DEG fixed top-k + LR", "AUROC"))
    lines.append(format_metric_line(overall_df, "DEG fixed top-k + LR", "AUPRC"))
    if scvi_df is not None and not scvi_df.empty:
        lines.append(format_metric_line(scvi_df, "scVI latent + LR", "AUROC"))
        lines.append(format_metric_line(scvi_df, "scVI latent + LR", "AUPRC"))
    lines.append("")
    lines.append("Subtype-stratified Cox results are paired with 95% Wald confidence intervals for hazard ratios in `tables/survival_hr_95ci.csv`.")
    lines.append("")
    lines.append("## Sample Size Limitations")
    lines.append("")
    for row in sample_size_df.itertuples(index=False):
        lines.append(f"- {row.scope}: {row.group} -> {row.detail}")
    lines.append("")
    lines.append("## Class Imbalance Handling")
    lines.append("")
    lines.append("- Top-k sensitivity analysis uses stratified tumor/normal holdout splits.")
    lines.append("- The revision benchmark holds out one normal sample and one tumor sample per fold.")
    lines.append('- Downstream logistic regression baselines use `class_weight=\"balanced\"`.')
    lines.append("")
    OUT_MARKDOWN.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main() -> None:
    OUT_TABLES.mkdir(parents=True, exist_ok=True)

    patient_pair_df = ensure_method_label(pd.read_csv(PATIENT_PAIR_FOLDS))
    overall_df = summarize_fold_metrics(patient_pair_df, ["method", "method_label"]).sort_values("mean_AUROC", ascending=False)
    by_celltype_df = summarize_fold_metrics(
        patient_pair_df,
        ["cell_type_key", "cell_type", "method", "method_label"],
    ).sort_values(["cell_type", "mean_AUROC"], ascending=[True, False])

    scvi_df: pd.DataFrame | None = None
    scvi_by_celltype_df: pd.DataFrame | None = None
    if SCVI_FOLDS.exists():
        scvi_folds_df = ensure_method_label(pd.read_csv(SCVI_FOLDS))
        scvi_df = summarize_fold_metrics(scvi_folds_df, ["method", "method_label"]).sort_values("mean_AUROC", ascending=False)
        scvi_by_celltype_df = summarize_fold_metrics(
            scvi_folds_df,
            ["cell_type_key", "cell_type", "method", "method_label"],
        ).sort_values(["cell_type", "mean_AUROC"], ascending=[True, False])

    cox_df = pd.read_csv(COX_GRID)
    cindex_df = pd.read_csv(CINDEX_GRID)

    sample_size_df = build_sample_size_summary(patient_pair_df, cindex_df)
    survival_df = build_survival_hr_summary(cox_df, cindex_df)

    overall_df.to_csv(OUT_TABLES / "patient_pair_metric_95ci_overall.csv", index=False)
    by_celltype_df.to_csv(OUT_TABLES / "patient_pair_metric_95ci_by_celltype.csv", index=False)
    sample_size_df.to_csv(OUT_TABLES / "sample_size_summary.csv", index=False)
    survival_df.to_csv(OUT_TABLES / "survival_hr_95ci.csv", index=False)

    if scvi_df is not None and scvi_by_celltype_df is not None:
        scvi_df.to_csv(OUT_TABLES / "scvi_patient_pair_metric_95ci_overall.csv", index=False)
        scvi_by_celltype_df.to_csv(OUT_TABLES / "scvi_patient_pair_metric_95ci_by_celltype.csv", index=False)

    write_markdown(overall_df, scvi_df, sample_size_df)
    print(f"Wrote tables to: {OUT_TABLES}")
    print(f"Wrote summary note: {OUT_MARKDOWN}")


if __name__ == "__main__":
    main()
