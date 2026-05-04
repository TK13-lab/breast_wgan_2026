#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import os
from pathlib import Path
import re
import sys
import warnings

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, roc_auc_score


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

import project_paths as pp


MPL_CACHE_DIR = pp.PROJECT_ROOT / "cache" / "matplotlib"
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE_DIR))

import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn.linear_model._logistic")


LABEL_PATTERN = re.compile(r"^(BN1|BN2|DN|EN|DT|ET|FT|GT|BT)")
LABEL_MAP = {
    "BN1": 0,
    "BN2": 0,
    "DN": 0,
    "EN": 0,
    "DT": 1,
    "ET": 1,
    "FT": 1,
    "GT": 1,
    "BT": 1,
}
NORMAL_SAMPLES = ("BN1", "BN2", "DN", "EN")
TUMOR_SAMPLES = ("BT", "DT", "ET", "FT", "GT")

METHOD_LABELS = {
    "gan_fixed_topk_lr": "GAN fixed top-k + LR",
    "deg_fixed_topk_lr": "DEG fixed top-k + LR",
    "variance_fixed_topk_lr": "Variance fixed top-k + LR",
    "pca_latent_lr": "PCA latent + LR",
    "all_genes_lr": "All candidate genes + LR",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Patient-pair holdout benchmark for the Paper3 revision. "
            "Each fold holds out one normal sample and one tumor sample, then evaluates "
            "fixed equal-budget signatures discovered on the full atlas."
        )
    )
    parser.add_argument("--cell-types", nargs="+", default=["all"], help="Cell type keys or 'all'.")
    parser.add_argument("--count-dir", default=None, help="Directory containing *_sigGene_expression_counts.csv files.")
    parser.add_argument("--gan-raw-dir", default=None, help="Directory containing ranked_genes_*.csv files.")
    parser.add_argument("--best-k-table", default=None, help="Optional best_k_by_celltype.csv from the Fig. 3 workflow.")
    parser.add_argument("--out-dir", default=None, help="Output directory for patient-pair benchmark tables and plots.")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--pca-components", type=int, default=20)
    return parser.parse_args()


def first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def resolve_count_dir(arg_value: str | None) -> Path:
    candidates: list[Path] = []
    if arg_value:
        candidates.append(Path(arg_value).expanduser())
    if "PAPER3_COUNT_DIR" in os.environ:
        candidates.append(Path(os.environ["PAPER3_COUNT_DIR"]).expanduser())
    candidates.extend(
        [
            pp.LEGACY_ROOT,
            pp.LEGACY_ROOT / "objects" / "csv",
        ]
    )
    found = first_existing(candidates)
    if found is None:
        tried = "\n".join(str(path) for path in candidates)
        raise FileNotFoundError(f"Could not resolve count directory. Tried:\n{tried}")
    return found.resolve()


def resolve_gan_raw_dir(arg_value: str | None) -> Path:
    candidates: list[Path] = []
    if arg_value:
        candidates.append(Path(arg_value).expanduser())
    if "PAPER3_FIG03_RAW_DIR" in os.environ:
        candidates.append(Path(os.environ["PAPER3_FIG03_RAW_DIR"]).expanduser())
    candidates.extend(
        [
            pp.FIG03_RAW,
            pp.LEGACY_ROOT / "output-GAN2" / "output2",
        ]
    )
    found = first_existing(candidates)
    if found is None:
        tried = "\n".join(str(path) for path in candidates)
        raise FileNotFoundError(f"Could not resolve GAN raw directory. Tried:\n{tried}")
    return found.resolve()


def resolve_best_k_table(arg_value: str | None) -> Path:
    candidates: list[Path] = []
    if arg_value:
        candidates.append(Path(arg_value).expanduser())
    candidates.extend(
        [
            pp.FIG03_TABLES / "best_k_by_celltype.csv",
        ]
    )
    found = first_existing(candidates)
    if found is None:
        tried = "\n".join(str(path) for path in candidates)
        raise FileNotFoundError(f"Could not resolve best-k table. Tried:\n{tried}")
    return found.resolve()


def count_path(count_dir: Path, spec: pp.CellTypeSpec) -> Path:
    direct = count_dir / spec.input_filename
    nested = count_dir / "objects" / "csv" / spec.input_filename
    found = first_existing([direct, nested])
    if found is None:
        raise FileNotFoundError(f"Missing count matrix for {spec.key}: tried {direct} and {nested}")
    return found


def read_prepared_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0, low_memory=False).T.reset_index().rename(columns={"index": "sample_id"})
    df["sample_group"] = df["sample_id"].str.extract(LABEL_PATTERN.pattern, expand=False)
    df["label"] = df["sample_group"].map(LABEL_MAP)
    df = df.dropna(subset=["sample_group", "label"]).copy()
    df["label"] = df["label"].astype(int)
    gene_cols = [
        col for col in df.columns
        if col not in {"sample_id", "sample_group", "label"} and not col.startswith(("RPS", "RPL"))
    ]
    out = df[["sample_id", "sample_group", "label"] + gene_cols].copy()
    return out


def fixed_effect_rank(df: pd.DataFrame, gene_cols: list[str]) -> list[str]:
    tumor_mean = df.loc[df["label"] == 1, gene_cols].mean(axis=0)
    normal_mean = df.loc[df["label"] == 0, gene_cols].mean(axis=0)
    return list((tumor_mean - normal_mean).abs().sort_values(ascending=False).index)


def fixed_variance_rank(df: pd.DataFrame, gene_cols: list[str]) -> list[str]:
    return list(df[gene_cols].var(axis=0).sort_values(ascending=False).index)


def fixed_gan_rank(gan_raw_dir: Path, spec: pp.CellTypeSpec) -> list[str]:
    path = gan_raw_dir / f"ranked_genes_{spec.output_stub}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing GAN ranked genes for {spec.key}: {path}")
    rank_df = pd.read_csv(path)
    return rank_df["gene"].astype(str).tolist()


def evaluate_lr(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    y_test: np.ndarray,
    seed: int,
) -> dict[str, float]:
    clf = LogisticRegression(
        solver="liblinear",
        penalty="l2",
        max_iter=2000,
        class_weight="balanced",
        random_state=seed,
    )
    clf.fit(x_train, y_train)
    probs = clf.predict_proba(x_test)[:, 1]
    preds = (probs >= 0.5).astype(int)
    return {
        "AUROC": float(roc_auc_score(y_test, probs)) if len(np.unique(y_test)) > 1 else float("nan"),
        "AUPRC": float(average_precision_score(y_test, probs)) if len(np.unique(y_test)) > 1 else float("nan"),
        "F1": float(f1_score(y_test, preds)),
        "ACC": float(accuracy_score(y_test, preds)),
    }


def evaluate_pca(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    gene_cols: list[str],
    n_components: int,
    seed: int,
) -> tuple[dict[str, float], int]:
    max_components = min(n_components, train_df.shape[0] - 1, len(gene_cols))
    if max_components < 1:
        return {"AUROC": float("nan"), "AUPRC": float("nan"), "F1": float("nan"), "ACC": float("nan")}, 0
    train_x = train_df[gene_cols].values
    test_x = test_df[gene_cols].values
    train_mean = train_x.mean(axis=0)
    train_sd = train_x.std(axis=0)
    train_sd[train_sd == 0] = 1.0
    train_scaled = (train_x - train_mean) / train_sd
    test_scaled = (test_x - train_mean) / train_sd
    pca = PCA(n_components=max_components, random_state=seed)
    train_latent = pca.fit_transform(train_scaled)
    test_latent = pca.transform(test_scaled)
    return evaluate_lr(train_latent, train_df["label"].values, test_latent, test_df["label"].values, seed), max_components


def fold_pairs(df: pd.DataFrame) -> list[tuple[str, str]]:
    normals = [sample for sample in NORMAL_SAMPLES if (df["sample_group"] == sample).any()]
    tumors = [sample for sample in TUMOR_SAMPLES if (df["sample_group"] == sample).any()]
    return list(itertools.product(normals, tumors))


def save_metric_plot(summary_df: pd.DataFrame, metric: str, plot_dir: Path) -> None:
    methods = [
        "gan_fixed_topk_lr",
        "deg_fixed_topk_lr",
        "variance_fixed_topk_lr",
        "pca_latent_lr",
        "all_genes_lr",
    ]
    plot_df = summary_df[summary_df["method"].isin(methods)].copy()
    x = np.arange(plot_df.shape[0])
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.bar(x, plot_df[f"mean_{metric}"], color="#4c78a8")
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df["method_label"], rotation=25, ha="right")
    ax.set_ylabel(f"Mean {metric}")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Patient-pair holdout benchmark: {metric}")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(plot_dir / f"patient_pair_equal_budget_{metric.lower()}.png", dpi=300)
    fig.savefig(plot_dir / f"patient_pair_equal_budget_{metric.lower()}.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    count_dir = resolve_count_dir(args.count_dir)
    gan_raw_dir = resolve_gan_raw_dir(args.gan_raw_dir)
    best_k_df = pd.read_csv(resolve_best_k_table(args.best_k_table))

    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else pp.OUTPUTS_ROOT / "revision" / "patient_pair_equal_budget"
    table_dir = out_dir / "tables"
    plot_dir = out_dir / "plots"
    table_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    fold_rows: list[dict[str, object]] = []

    for spec in pp.resolve_celltypes(args.cell_types):
        print(f"[patient-pair] {spec.display_name}", flush=True)
        full_df = read_prepared_counts(count_path(count_dir, spec))
        gene_cols = [col for col in full_df.columns if col not in {"sample_id", "sample_group", "label"}]
        best_k_match = best_k_df[best_k_df["cell_type"].astype(str) == spec.display_name]
        if best_k_match.empty:
            raise ValueError(f"Missing best_k for {spec.display_name}")
        selected_k = int(best_k_match.iloc[0]["best_k"])

        gan_genes = [g for g in fixed_gan_rank(gan_raw_dir, spec) if g in gene_cols][:selected_k]
        deg_genes = [g for g in fixed_effect_rank(full_df, gene_cols) if g in gene_cols][:selected_k]
        var_genes = [g for g in fixed_variance_rank(full_df, gene_cols) if g in gene_cols][:selected_k]

        for normal_sample, tumor_sample in fold_pairs(full_df):
            test_mask = full_df["sample_group"].isin([normal_sample, tumor_sample])
            train_df = full_df.loc[~test_mask].copy()
            test_df = full_df.loc[test_mask].copy()
            if train_df["label"].nunique() < 2 or test_df["label"].nunique() < 2:
                continue

            for method, genes in [
                ("gan_fixed_topk_lr", gan_genes),
                ("deg_fixed_topk_lr", deg_genes),
                ("variance_fixed_topk_lr", var_genes),
            ]:
                metrics = evaluate_lr(
                    train_df[genes].values,
                    train_df["label"].values,
                    test_df[genes].values,
                    test_df["label"].values,
                    seed=args.seed,
                )
                fold_rows.append(
                    {
                        "cell_type_key": spec.key,
                        "cell_type": spec.display_name,
                        "method": method,
                        "method_label": METHOD_LABELS[method],
                        "selected_k": selected_k,
                        "normal_holdout": normal_sample,
                        "tumor_holdout": tumor_sample,
                        "n_train": int(train_df.shape[0]),
                        "n_test": int(test_df.shape[0]),
                        **metrics,
                    }
                )

            for method, metrics, effective_k in [
                (
                    "all_genes_lr",
                    evaluate_lr(
                        train_df[gene_cols].values,
                        train_df["label"].values,
                        test_df[gene_cols].values,
                        test_df["label"].values,
                        seed=args.seed,
                    ),
                    len(gene_cols),
                ),
                (
                    "pca_latent_lr",
                    evaluate_pca(train_df, test_df, gene_cols, args.pca_components, args.seed)[0],
                    min(args.pca_components, train_df.shape[0] - 1, len(gene_cols)),
                ),
            ]:
                fold_rows.append(
                    {
                        "cell_type_key": spec.key,
                        "cell_type": spec.display_name,
                        "method": method,
                        "method_label": METHOD_LABELS[method],
                        "selected_k": int(effective_k),
                        "normal_holdout": normal_sample,
                        "tumor_holdout": tumor_sample,
                        "n_train": int(train_df.shape[0]),
                        "n_test": int(test_df.shape[0]),
                        **metrics,
                    }
                )

    fold_df = pd.DataFrame(fold_rows)
    fold_df.to_csv(table_dir / "patient_pair_equal_budget_fold_metrics.csv", index=False)

    summary_df = (
        fold_df.groupby(["method", "method_label"], as_index=False)
        .agg(
            mean_AUROC=("AUROC", "mean"),
            median_AUROC=("AUROC", "median"),
            mean_AUPRC=("AUPRC", "mean"),
            median_AUPRC=("AUPRC", "median"),
            mean_F1=("F1", "mean"),
            mean_ACC=("ACC", "mean"),
            n_folds=("AUROC", "size"),
            n_cell_types=("cell_type", "nunique"),
        )
        .sort_values(["mean_AUROC", "mean_AUPRC"], ascending=False)
    )
    summary_df.to_csv(table_dir / "patient_pair_equal_budget_summary.csv", index=False)

    by_celltype_df = (
        fold_df.groupby(["cell_type", "method", "method_label"], as_index=False)
        .agg(
            mean_AUROC=("AUROC", "mean"),
            mean_AUPRC=("AUPRC", "mean"),
            mean_F1=("F1", "mean"),
            mean_ACC=("ACC", "mean"),
            n_folds=("AUROC", "size"),
        )
        .sort_values(["cell_type", "mean_AUROC"], ascending=[True, False])
    )
    by_celltype_df.to_csv(table_dir / "patient_pair_equal_budget_by_celltype.csv", index=False)

    save_metric_plot(summary_df, "AUROC", plot_dir)
    save_metric_plot(summary_df, "AUPRC", plot_dir)

    print(f"Count matrices: {count_dir}")
    print(f"GAN ranked genes: {gan_raw_dir}")
    print(f"Wrote fold metrics to: {table_dir / 'patient_pair_equal_budget_fold_metrics.csv'}")
    print(f"Wrote summary to: {table_dir / 'patient_pair_equal_budget_summary.csv'}")


if __name__ == "__main__":
    main()
