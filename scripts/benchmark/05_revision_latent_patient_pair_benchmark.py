#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import io
import os
from pathlib import Path
import re
import sys
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, roc_auc_score


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

import project_paths as pp


MPL_CACHE_DIR = pp.PROJECT_ROOT / "cache" / "matplotlib"
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE_DIR))

warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn.linear_model._logistic")
warnings.filterwarnings("ignore", message=".*NVML.*")
warnings.filterwarnings("ignore", message=".*train_dataloader.*does not have many workers.*")
warnings.filterwarnings("ignore", message=".*logging interval.*")


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Patient-pair latent benchmark for scVI and optionally scGen. "
            "Each fold trains on all but one normal and one tumor sample, then evaluates "
            "tumor-normal discrimination using latent coordinates and logistic regression."
        )
    )
    parser.add_argument("--model", choices=["scvi", "scgen"], required=True)
    parser.add_argument("--cell-types", nargs="+", default=["all"], help="Cell type keys or 'all'.")
    parser.add_argument("--count-dir", default=None, help="Directory containing *_sigGene_expression_counts.csv files.")
    parser.add_argument("--out-dir", default=None, help="Output directory for latent patient-pair benchmark tables.")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--latent-dim", type=int, default=20)
    parser.add_argument("--max-epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=256)
    return parser.parse_args()


def require_backend(model_name: str):
    import anndata as ad

    if model_name == "scvi":
        import scvi
        return ad, scvi, None

    try:
        import scgen
    except Exception as exc:
        raise SystemExit(
            "scGen could not be imported in the current environment.\n"
            "This is often caused by a version mismatch between scgen and scvi-tools.\n"
            "Current environment note: the installed scgen package expects APIs that may not exist in scvi-tools 1.4.x.\n"
            "If you want a runnable scGen benchmark, we will likely need a separate pinned environment.\n"
            f"Original import error: {exc}"
        ) from exc
    return ad, None, scgen


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
    return df[["sample_id", "sample_group", "label"] + gene_cols].copy()


def fold_pairs(df: pd.DataFrame) -> list[tuple[str, str]]:
    normals = [sample for sample in NORMAL_SAMPLES if (df["sample_group"] == sample).any()]
    tumors = [sample for sample in TUMOR_SAMPLES if (df["sample_group"] == sample).any()]
    return [(n, t) for n in normals for t in tumors]


def evaluate_lr(x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray, y_test: np.ndarray, seed: int) -> dict[str, float]:
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


def scvi_latent(train_df: pd.DataFrame, test_df: pd.DataFrame, latent_dim: int, max_epochs: int, batch_size: int, seed: int):
    import anndata as ad
    import scvi

    gene_cols = [col for col in train_df.columns if col not in {"sample_id", "sample_group", "label"}]
    adata_train = ad.AnnData(X=train_df[gene_cols].values.astype(np.float32))
    adata_train.obs["sample_group"] = train_df["sample_group"].astype(str).values
    adata_train.obs["tumor_status"] = np.where(train_df["label"].values == 1, "Tumor", "Normal")
    adata_train.var_names = gene_cols

    adata_test = ad.AnnData(X=test_df[gene_cols].values.astype(np.float32))
    adata_test.obs["sample_group"] = test_df["sample_group"].astype(str).values
    adata_test.obs["tumor_status"] = np.where(test_df["label"].values == 1, "Tumor", "Normal")
    adata_test.var_names = gene_cols

    scvi.settings.seed = seed
    scvi.model.SCVI.setup_anndata(adata_train)
    model = scvi.model.SCVI(adata_train, n_latent=min(latent_dim, len(gene_cols)))
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        model.train(max_epochs=max_epochs, batch_size=batch_size, enable_progress_bar=False)
    train_latent = model.get_latent_representation(adata_train)
    test_latent = model.get_latent_representation(adata_test)
    return train_latent, test_latent


def scgen_latent(train_df: pd.DataFrame, test_df: pd.DataFrame, latent_dim: int, max_epochs: int, batch_size: int, seed: int):
    import anndata as ad
    import scgen

    gene_cols = [col for col in train_df.columns if col not in {"sample_id", "sample_group", "label"}]
    adata_train = ad.AnnData(X=train_df[gene_cols].values.astype(np.float32))
    adata_train.obs["condition"] = np.where(train_df["label"].values == 1, "Tumor", "Normal")
    adata_train.obs["celltype"] = "shared_celltype"
    adata_train.var_names = gene_cols

    adata_test = ad.AnnData(X=test_df[gene_cols].values.astype(np.float32))
    adata_test.obs["condition"] = np.where(test_df["label"].values == 1, "Tumor", "Normal")
    adata_test.obs["celltype"] = "shared_celltype"
    adata_test.var_names = gene_cols

    np.random.seed(seed)
    if hasattr(scgen, "SCGEN"):
        model_cls = scgen.SCGEN
    elif hasattr(scgen, "model") and hasattr(scgen.model, "SCGEN"):
        model_cls = scgen.model.SCGEN
    else:
        raise RuntimeError("Could not find SCGEN class in installed scgen package.")

    model_cls.setup_anndata(adata_train, batch_key=None, labels_key="celltype")
    model = model_cls(adata_train, n_latent=min(latent_dim, len(gene_cols)))
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        model.train(max_epochs=max_epochs, batch_size=batch_size)

    if hasattr(model, "get_latent_representation"):
        train_latent = model.get_latent_representation(adata_train)
        test_latent = model.get_latent_representation(adata_test)
    elif hasattr(model, "module") and hasattr(model.module, "encoder"):
        train_latent = model.module.encoder(np.asarray(adata_train.X)).detach().cpu().numpy()
        test_latent = model.module.encoder(np.asarray(adata_test.X)).detach().cpu().numpy()
    else:
        raise RuntimeError("Installed scgen model does not expose a latent representation API.")
    return train_latent, test_latent


def main() -> None:
    args = parse_args()
    require_backend(args.model)
    count_dir = resolve_count_dir(args.count_dir)

    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else pp.OUTPUTS_ROOT / "revision" / "latent_patient_pair"
    table_dir = out_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []

    for spec in pp.resolve_celltypes(args.cell_types):
        print(f"[{args.model}] {spec.display_name}", flush=True)
        full_df = read_prepared_counts(count_path(count_dir, spec))
        for normal_sample, tumor_sample in fold_pairs(full_df):
            test_mask = full_df["sample_group"].isin([normal_sample, tumor_sample])
            train_df = full_df.loc[~test_mask].copy()
            test_df = full_df.loc[test_mask].copy()
            if train_df["label"].nunique() < 2 or test_df["label"].nunique() < 2:
                continue

            if args.model == "scvi":
                train_latent, test_latent = scvi_latent(
                    train_df, test_df, args.latent_dim, args.max_epochs, args.batch_size, args.seed
                )
            else:
                train_latent, test_latent = scgen_latent(
                    train_df, test_df, args.latent_dim, args.max_epochs, args.batch_size, args.seed
                )

            metrics = evaluate_lr(
                train_latent,
                train_df["label"].values,
                test_latent,
                test_df["label"].values,
                args.seed,
            )
            rows.append(
                {
                    "cell_type_key": spec.key,
                    "cell_type": spec.display_name,
                    "method": f"{args.model}_latent_lr",
                    "latent_dim": args.latent_dim,
                    "max_epochs": args.max_epochs,
                    "normal_holdout": normal_sample,
                    "tumor_holdout": tumor_sample,
                    "n_train": int(train_df.shape[0]),
                    "n_test": int(test_df.shape[0]),
                    **metrics,
                }
            )
            pd.DataFrame(rows).to_csv(table_dir / f"{args.model}_latent_patient_pair_fold_metrics.csv", index=False)

    fold_df = pd.DataFrame(rows)
    summary_df = (
        fold_df.groupby(["method"], as_index=False)
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
    )
    summary_df.to_csv(table_dir / f"{args.model}_latent_patient_pair_summary.csv", index=False)

    by_celltype_df = (
        fold_df.groupby(["cell_type", "method"], as_index=False)
        .agg(
            mean_AUROC=("AUROC", "mean"),
            mean_AUPRC=("AUPRC", "mean"),
            mean_F1=("F1", "mean"),
            mean_ACC=("ACC", "mean"),
            n_folds=("AUROC", "size"),
        )
    )
    by_celltype_df.to_csv(table_dir / f"{args.model}_latent_patient_pair_by_celltype.csv", index=False)

    print(f"Wrote latent patient-pair fold metrics to: {table_dir / f'{args.model}_latent_patient_pair_fold_metrics.csv'}")
    print(f"Wrote latent patient-pair summary to: {table_dir / f'{args.model}_latent_patient_pair_summary.csv'}")


if __name__ == "__main__":
    main()
