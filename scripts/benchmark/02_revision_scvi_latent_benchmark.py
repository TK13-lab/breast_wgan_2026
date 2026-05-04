#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, roc_auc_score
from sklearn.model_selection import train_test_split


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

import project_paths as pp


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Optional scVI latent-representation benchmark for the Paper3 revision. "
            "This compares scVI latent embeddings with the proposed compact gene-ranking workflow "
            "on the same cell-type count matrices."
        )
    )
    parser.add_argument("--cell-types", nargs="+", default=["all"], help="Cell type keys or 'all'.")
    parser.add_argument("--count-dir", default=None, help="Directory containing *_sigGene_expression_counts.csv files.")
    parser.add_argument("--out-dir", default=None, help="Output directory for scVI benchmark tables.")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--test-size", type=float, default=0.2)
    parser.add_argument("--latent-dim", type=int, default=20)
    parser.add_argument("--max-epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=256)
    return parser.parse_args()


def require_scvi():
    try:
        import anndata as ad
        import scvi
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "Missing optional scVI dependencies. Install scvi-tools and anndata to run this benchmark:\n"
            "  pip install scvi-tools anndata\n"
            f"Original import error: {exc}"
        ) from exc
    return ad, scvi


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


def count_path(count_dir: Path, spec: pp.CellTypeSpec) -> Path:
    direct = count_dir / spec.input_filename
    nested = count_dir / "objects" / "csv" / spec.input_filename
    found = first_existing([direct, nested])
    if found is None:
        raise FileNotFoundError(f"Missing count matrix for {spec.key}: tried {direct} and {nested}")
    return found


def read_counts(path: Path) -> tuple[pd.DataFrame, pd.Series]:
    df = pd.read_csv(path, index_col=0, low_memory=False)
    df_t = df.T.reset_index().rename(columns={"index": "sample_id"})
    labels = df_t["sample_id"].str.extract(LABEL_PATTERN.pattern, expand=False).map(LABEL_MAP)
    keep = labels.notna()
    labels = labels.loc[keep].astype(int)
    x_df = df_t.loc[keep].drop(columns=["sample_id"]).copy()
    x_df = x_df.loc[:, [col for col in x_df.columns if not col.startswith(("RPS", "RPL"))]]
    return x_df.astype(np.float32), labels


def evaluate_latent(x_latent: np.ndarray, labels: np.ndarray, seed: int, test_size: float) -> dict[str, float]:
    train_idx, test_idx = train_test_split(
        np.arange(labels.shape[0]),
        test_size=test_size,
        stratify=labels,
        random_state=seed,
    )
    clf = LogisticRegression(
        solver="liblinear",
        penalty="l2",
        max_iter=2000,
        class_weight="balanced",
        random_state=seed,
    )
    clf.fit(x_latent[train_idx], labels[train_idx])
    probs = clf.predict_proba(x_latent[test_idx])[:, 1]
    preds = (probs >= 0.5).astype(int)
    return {
        "n_train": int(train_idx.shape[0]),
        "n_test": int(test_idx.shape[0]),
        "AUROC": float(roc_auc_score(labels[test_idx], probs)) if len(np.unique(labels[test_idx])) > 1 else float("nan"),
        "AUPRC": float(average_precision_score(labels[test_idx], probs))
        if len(np.unique(labels[test_idx])) > 1
        else float("nan"),
        "F1": float(f1_score(labels[test_idx], preds)),
        "ACC": float(accuracy_score(labels[test_idx], preds)),
    }


def run_scvi_one(ad, scvi, spec: pp.CellTypeSpec, count_dir: Path, args: argparse.Namespace) -> dict[str, object]:
    x_df, labels = read_counts(count_path(count_dir, spec))
    adata = ad.AnnData(X=x_df.values)
    adata.obs["tumor_status"] = pd.Categorical(np.where(labels.values == 1, "Tumor", "Normal"))
    adata.obs_names = x_df.index.astype(str)
    adata.var_names = x_df.columns.astype(str)

    scvi.settings.seed = args.seed
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata, n_latent=args.latent_dim)
    model.train(max_epochs=args.max_epochs, batch_size=args.batch_size)
    latent = model.get_latent_representation()
    metrics = evaluate_latent(latent, labels.values, seed=args.seed, test_size=args.test_size)
    return {
        "cell_type_key": spec.key,
        "cell_type": spec.display_name,
        "method": "scvi_latent_lr",
        "latent_dim": int(args.latent_dim),
        "max_epochs": int(args.max_epochs),
        "n_cells": int(x_df.shape[0]),
        "n_genes": int(x_df.shape[1]),
        **metrics,
    }


def main() -> None:
    args = parse_args()
    ad, scvi = require_scvi()
    count_dir = resolve_count_dir(args.count_dir)
    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else pp.OUTPUTS_ROOT / "revision" / "benchmark"
    table_dir = out_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for spec in pp.resolve_celltypes(args.cell_types):
        print(f"[scVI] {spec.display_name}", flush=True)
        rows.append(run_scvi_one(ad, scvi, spec, count_dir=count_dir, args=args))
        pd.DataFrame(rows).to_csv(table_dir / "celltype_scvi_latent_benchmark.csv", index=False)

    print(f"Wrote scVI latent benchmark to: {table_dir / 'celltype_scvi_latent_benchmark.csv'}")


if __name__ == "__main__":
    main()
