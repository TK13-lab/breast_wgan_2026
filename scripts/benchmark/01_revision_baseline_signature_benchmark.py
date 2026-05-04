#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path
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

import gan_signature_utils as gsu
import project_paths as pp


MPL_CACHE_DIR = pp.PROJECT_ROOT / "cache" / "matplotlib"
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE_DIR))

import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn.linear_model._logistic")


METHOD_LABELS = {
    "gan_rank_lr": "GAN rank + LR",
    "gan_rank_knn": "GAN rank + kNN",
    "deg_effect_rank_lr": "DEG effect rank + LR",
    "deg_effect_rank_knn": "DEG effect rank + kNN",
    "variance_rank_lr": "Variance rank + LR",
    "variance_rank_knn": "Variance rank + kNN",
    "pca_latent_lr": "PCA latent + LR",
    "all_genes_lr": "All candidate genes + LR",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark Paper3 GAN-derived gene rankings against simpler feature-ranking "
            "and latent-representation baselines using the same cell-type input matrices."
        )
    )
    parser.add_argument("--cell-types", nargs="+", default=["all"], help="Cell type keys or 'all'.")
    parser.add_argument("--count-dir", default=None, help="Directory containing *_sigGene_expression_counts.csv files.")
    parser.add_argument("--gan-raw-dir", default=None, help="Directory containing Fig. 3 GAN k-curve/ranking CSV files.")
    parser.add_argument("--best-k-table", default=None, help="Optional best_k_by_celltype.csv from the Fig. 3 workflow.")
    parser.add_argument("--out-dir", default=None, help="Output directory for revision benchmark tables and plots.")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--test-size", type=float, default=0.2)
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
            pp.LEGACY_ROOT / "results" / "gan_output_v2",
        ]
    )
    found = first_existing(candidates)
    if found is None:
        tried = "\n".join(str(path) for path in candidates)
        raise FileNotFoundError(f"Could not resolve GAN raw result directory. Tried:\n{tried}")
    return found.resolve()


def resolve_best_k_table(arg_value: str | None) -> Path | None:
    candidates: list[Path] = []
    if arg_value:
        candidates.append(Path(arg_value).expanduser())
    candidates.extend(
        [
            pp.FIG03_TABLES / "best_k_by_celltype.csv",
        ]
    )
    return first_existing(candidates)


def count_path(count_dir: Path, spec: pp.CellTypeSpec) -> Path:
    direct = count_dir / spec.input_filename
    nested = count_dir / "objects" / "csv" / spec.input_filename
    found = first_existing([direct, nested])
    if found is None:
        raise FileNotFoundError(f"Missing count matrix for {spec.key}: tried {direct} and {nested}")
    return found


def safe_metric(y_true: np.ndarray, y_score: np.ndarray, metric: str) -> float:
    if metric == "AUROC":
        return float(roc_auc_score(y_true, y_score)) if len(np.unique(y_true)) > 1 else float("nan")
    if metric == "AUPRC":
        return float(average_precision_score(y_true, y_score)) if len(np.unique(y_true)) > 1 else float("nan")
    raise ValueError(metric)


def evaluate_matrix_classifier(
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
        "AUROC": safe_metric(y_test, probs, "AUROC"),
        "AUPRC": safe_metric(y_test, probs, "AUPRC"),
        "F1": float(f1_score(y_test, preds)),
        "ACC": float(accuracy_score(y_test, preds)),
    }


def rank_by_effect_size(x_train_df: pd.DataFrame, y_train: np.ndarray) -> pd.DataFrame:
    tumor_mean = x_train_df.loc[y_train == 1].mean(axis=0)
    normal_mean = x_train_df.loc[y_train == 0].mean(axis=0)
    scores = (tumor_mean - normal_mean).abs().sort_values(ascending=False)
    return pd.DataFrame({"gene": scores.index, "score": scores.values, "rank": np.arange(1, len(scores) + 1)})


def rank_by_variance(x_train_df: pd.DataFrame) -> pd.DataFrame:
    scores = x_train_df.var(axis=0).sort_values(ascending=False)
    return pd.DataFrame({"gene": scores.index, "score": scores.values, "rank": np.arange(1, len(scores) + 1)})


def read_gan_curve(raw_dir: Path, spec: pp.CellTypeSpec, kind: str) -> pd.DataFrame:
    path = raw_dir / f"k_curve_{kind}_{spec.output_stub}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing GAN {kind} k-curve for {spec.key}: {path}")
    df = pd.read_csv(path)
    if "k" in df.columns:
        df = df.rename(columns={"k": "top_k"})
    df["top_k"] = df["top_k"].astype(int)
    return df[["top_k", "AUROC", "AUPRC", "F1", "ACC"]].drop_duplicates("top_k")


def baseline_curve(
    rank_df: pd.DataFrame,
    x_train_df: pd.DataFrame,
    y_train: np.ndarray,
    x_test_df: pd.DataFrame,
    y_test: np.ndarray,
    clf_kind: str,
    seed: int,
) -> pd.DataFrame:
    curve = gsu.evaluate_k_curve(
        rank_df=rank_df,
        x_train_df=x_train_df,
        y_train=y_train,
        x_test_df=x_test_df,
        y_test=y_test,
        ks=pp.TOP_K_GRID,
        clf_kind=clf_kind,
        seed=seed,
    )
    return curve.rename(columns={"k": "top_k"})


def pca_latent_metrics(
    x_train_df: pd.DataFrame,
    y_train: np.ndarray,
    x_test_df: pd.DataFrame,
    y_test: np.ndarray,
    n_components: int,
    seed: int,
) -> tuple[dict[str, float], int]:
    max_components = min(n_components, x_train_df.shape[0] - 1, x_train_df.shape[1])
    if max_components < 1:
        return {"AUROC": float("nan"), "AUPRC": float("nan"), "F1": float("nan"), "ACC": float("nan")}, 0
    pca = PCA(n_components=max_components, random_state=seed)
    x_train_latent = pca.fit_transform(x_train_df.values)
    x_test_latent = pca.transform(x_test_df.values)
    return evaluate_matrix_classifier(x_train_latent, y_train, x_test_latent, y_test, seed=seed), max_components


def add_curve_rows(rows: list[dict[str, object]], curve: pd.DataFrame, spec: pp.CellTypeSpec, method: str) -> None:
    for rec in curve.to_dict(orient="records"):
        rows.append(
            {
                "cell_type_key": spec.key,
                "cell_type": spec.display_name,
                "method": method,
                "method_label": METHOD_LABELS[method],
                "top_k": int(rec["top_k"]),
                "AUROC": float(rec["AUROC"]),
                "AUPRC": float(rec["AUPRC"]),
                "F1": float(rec["F1"]),
                "ACC": float(rec["ACC"]),
            }
        )


def selected_k_for_cell(best_k_df: pd.DataFrame | None, spec: pp.CellTypeSpec, fallback_curve: pd.DataFrame) -> int:
    if best_k_df is not None:
        candidates = best_k_df[
            best_k_df["cell_type"].astype(str).isin([spec.display_name, spec.key, spec.output_stub])
        ]
        if not candidates.empty:
            return int(candidates.iloc[0]["best_k"])
    best_idx = fallback_curve["AUROC"].astype(float).idxmax()
    return int(fallback_curve.loc[best_idx, "top_k"])


def write_summary_tables(long_df: pd.DataFrame, selected_df: pd.DataFrame, table_dir: Path) -> None:
    long_df.to_csv(table_dir / "celltype_benchmark_long.csv", index=False)
    selected_df.to_csv(table_dir / "celltype_selected_k_comparison.csv", index=False)

    summary = (
        selected_df.groupby(["method", "method_label"], as_index=False)
        .agg(
            mean_AUROC=("AUROC", "mean"),
            median_AUROC=("AUROC", "median"),
            mean_AUPRC=("AUPRC", "mean"),
            median_AUPRC=("AUPRC", "median"),
            mean_F1=("F1", "mean"),
            mean_ACC=("ACC", "mean"),
            n_cell_types=("cell_type", "nunique"),
        )
        .sort_values(["mean_AUROC", "mean_AUPRC"], ascending=False)
    )
    summary.to_csv(table_dir / "celltype_benchmark_summary.csv", index=False)


def save_metric_plot(selected_df: pd.DataFrame, metric: str, plot_dir: Path) -> None:
    methods = [
        "gan_rank_lr",
        "deg_effect_rank_lr",
        "variance_rank_lr",
        "pca_latent_lr",
        "all_genes_lr",
    ]
    plot_df = selected_df[selected_df["method"].isin(methods)].copy()
    cell_types = list(dict.fromkeys(plot_df["cell_type"].tolist()))

    x = np.arange(len(cell_types))
    width = 0.15
    fig, ax = plt.subplots(figsize=(12, 5.5))
    for idx, method in enumerate(methods):
        vals = []
        for cell_type in cell_types:
            sub = plot_df[(plot_df["cell_type"] == cell_type) & (plot_df["method"] == method)]
            vals.append(float(sub[metric].iloc[0]) if not sub.empty else np.nan)
        ax.bar(x + (idx - 2) * width, vals, width=width, label=METHOD_LABELS[method])
    ax.set_xticks(x)
    ax.set_xticklabels(cell_types, rotation=35, ha="right")
    ax.set_ylabel(metric)
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, ncol=2)
    ax.set_title(f"Revision benchmark: {metric} at selected k")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(plot_dir / f"revision_benchmark_{metric.lower()}.png", dpi=300)
    fig.savefig(plot_dir / f"revision_benchmark_{metric.lower()}.pdf")
    plt.close(fig)


def save_delta_plot(selected_df: pd.DataFrame, plot_dir: Path) -> None:
    rows = []
    for cell_type, group in selected_df.groupby("cell_type"):
        gan = group[group["method"] == "gan_rank_lr"]
        deg = group[group["method"] == "deg_effect_rank_lr"]
        if gan.empty or deg.empty:
            continue
        rows.append(
            {
                "cell_type": cell_type,
                "delta_AUROC": float(gan["AUROC"].iloc[0] - deg["AUROC"].iloc[0]),
                "delta_AUPRC": float(gan["AUPRC"].iloc[0] - deg["AUPRC"].iloc[0]),
            }
        )
    delta_df = pd.DataFrame(rows)
    if delta_df.empty:
        return

    x = np.arange(delta_df.shape[0])
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.axhline(0, color="black", linewidth=0.8)
    ax.bar(x - 0.18, delta_df["delta_AUROC"], width=0.36, label="Delta AUROC")
    ax.bar(x + 0.18, delta_df["delta_AUPRC"], width=0.36, label="Delta AUPRC")
    ax.set_xticks(x)
    ax.set_xticklabels(delta_df["cell_type"], rotation=35, ha="right")
    ax.set_ylabel("GAN rank + LR minus DEG effect rank + LR")
    ax.legend(frameon=False)
    ax.set_title("Improvement over DEG effect-size ranking")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(plot_dir / "revision_benchmark_delta_vs_deg.png", dpi=300)
    fig.savefig(plot_dir / "revision_benchmark_delta_vs_deg.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    count_dir = resolve_count_dir(args.count_dir)
    gan_raw_dir = resolve_gan_raw_dir(args.gan_raw_dir)
    best_k_path = resolve_best_k_table(args.best_k_table)
    best_k_df = pd.read_csv(best_k_path) if best_k_path else None

    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else pp.OUTPUTS_ROOT / "revision" / "benchmark"
    table_dir = out_dir / "tables"
    plot_dir = out_dir / "plots"
    table_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    selected_rows: list[dict[str, object]] = []

    for spec in pp.resolve_celltypes(args.cell_types):
        prepared = gsu.read_and_prepare_inputs(
            raw_counts_csv=count_path(count_dir, spec),
            seed=args.seed,
            test_size=args.test_size,
            exclude_prefixes=pp.EXCLUDED_GENE_PREFIXES,
        )
        x_train_df = prepared["x_train_df"]
        x_test_df = prepared["x_test_df"]
        y_train = prepared["y_train"]
        y_test = prepared["y_test"]

        gan_lr = read_gan_curve(gan_raw_dir, spec, "lr")
        gan_knn = read_gan_curve(gan_raw_dir, spec, "knn")
        deg_rank = rank_by_effect_size(x_train_df, y_train)
        var_rank = rank_by_variance(x_train_df)

        curves = {
            "gan_rank_lr": gan_lr,
            "gan_rank_knn": gan_knn,
            "deg_effect_rank_lr": baseline_curve(deg_rank, x_train_df, y_train, x_test_df, y_test, "lr", args.seed),
            "deg_effect_rank_knn": baseline_curve(deg_rank, x_train_df, y_train, x_test_df, y_test, "knn", args.seed),
            "variance_rank_lr": baseline_curve(var_rank, x_train_df, y_train, x_test_df, y_test, "lr", args.seed),
            "variance_rank_knn": baseline_curve(var_rank, x_train_df, y_train, x_test_df, y_test, "knn", args.seed),
        }

        for method, curve in curves.items():
            add_curve_rows(rows, curve, spec, method)

        selected_k = selected_k_for_cell(best_k_df, spec, gan_lr)
        for method, curve in curves.items():
            selected = curve[curve["top_k"] == selected_k]
            if selected.empty:
                selected = curve.iloc[[curve["top_k"].sub(selected_k).abs().idxmin()]]
            rec = selected.iloc[0].to_dict()
            selected_rows.append(
                {
                    "cell_type_key": spec.key,
                    "cell_type": spec.display_name,
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "selected_k": selected_k,
                    "evaluated_top_k": int(rec["top_k"]),
                    "AUROC": float(rec["AUROC"]),
                    "AUPRC": float(rec["AUPRC"]),
                    "F1": float(rec["F1"]),
                    "ACC": float(rec["ACC"]),
                }
            )

        all_gene_metrics = evaluate_matrix_classifier(
            x_train_df.values, y_train, x_test_df.values, y_test, seed=args.seed
        )
        pca_metrics, pca_components = pca_latent_metrics(
            x_train_df, y_train, x_test_df, y_test, args.pca_components, seed=args.seed
        )
        for method, metrics, effective_k in [
            ("all_genes_lr", all_gene_metrics, x_train_df.shape[1]),
            ("pca_latent_lr", pca_metrics, pca_components),
        ]:
            selected_rows.append(
                {
                    "cell_type_key": spec.key,
                    "cell_type": spec.display_name,
                    "method": method,
                    "method_label": METHOD_LABELS[method],
                    "selected_k": selected_k,
                    "evaluated_top_k": int(effective_k),
                    **metrics,
                }
            )

    long_df = pd.DataFrame(rows)
    selected_df = pd.DataFrame(selected_rows)
    write_summary_tables(long_df, selected_df, table_dir)
    save_metric_plot(selected_df, "AUROC", plot_dir)
    save_metric_plot(selected_df, "AUPRC", plot_dir)
    save_delta_plot(selected_df, plot_dir)

    print(f"Count matrices: {count_dir}")
    print(f"GAN raw curves: {gan_raw_dir}")
    print(f"Best-k table: {best_k_path if best_k_path else 'not found; used GAN LR max AUROC fallback'}")
    print(f"Wrote benchmark tables to: {table_dir}")
    print(f"Wrote benchmark plots to: {plot_dir}")


if __name__ == "__main__":
    main()
