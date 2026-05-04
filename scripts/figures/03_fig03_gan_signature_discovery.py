#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys

import pandas as pd


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

import gan_signature_utils as gsu
import project_paths as pp


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the Paper3 breast WGAN signature workflow for one or more cell types."
    )
    parser.add_argument("--cell-types", nargs="+", default=["all"], help="Cell type keys or 'all'.")
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--gamma", type=float, default=0.1)
    parser.add_argument("--weight-decay", type=float, default=1e-5)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--test-size", type=float, default=0.2)
    parser.add_argument("--disease-index", type=int, default=1)
    parser.add_argument("--device", default="auto", help="auto, cpu, or cuda")
    parser.add_argument("--raw-dir", default=None, help="Optional override for the raw output directory.")
    return parser.parse_args()


def run_one_celltype(spec: pp.CellTypeSpec, args: argparse.Namespace, raw_dir: Path) -> dict[str, object]:
    input_csv = pp.legacy_count_path(spec)
    if not input_csv.exists():
        raise FileNotFoundError(f"Missing input counts for {spec.key}: {input_csv}")

    print(f"[START] {spec.key}: input={input_csv}", flush=True)
    prepared = gsu.read_and_prepare_inputs(
        raw_counts_csv=input_csv,
        seed=args.seed,
        test_size=args.test_size,
        exclude_prefixes=pp.EXCLUDED_GENE_PREFIXES,
    )

    x_train_df = prepared["x_train_df"]
    x_test_df = prepared["x_test_df"]
    y_train = prepared["y_train"]
    y_test = prepared["y_test"]
    metadata = prepared["metadata"]

    enc_hidden, latent_dim = gsu.propose_architecture(x_train_df.shape[1])
    cls_hidden = [128]
    cfg = gsu.TrainConfig(
        n_epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        weight_decay=args.weight_decay,
        alpha=args.alpha,
        beta=args.beta,
        gamma=args.gamma,
        n_critic=3,
        device=gsu.resolve_device(args.device),
    )

    gsu.set_seed(args.seed)
    best_model, history_df = gsu.train_cv_return_best(
        x=x_train_df.values,
        y=y_train,
        enc_hidden=enc_hidden,
        latent_dim=latent_dim,
        cls_hidden=cls_hidden,
        cfg=cfg,
        n_splits=args.folds,
        seed=args.seed,
    )

    weight_matrix = gsu.compute_total_weight(best_model)
    rank_df = gsu.rank_genes(
        weight_matrix=weight_matrix,
        gene_names=x_train_df.columns.tolist(),
        disease_index=args.disease_index,
        use_abs=True,
    )

    curve_knn = gsu.evaluate_k_curve(
        rank_df=rank_df,
        x_train_df=x_train_df,
        y_train=y_train,
        x_test_df=x_test_df,
        y_test=y_test,
        ks=pp.TOP_K_GRID,
        clf_kind="knn",
        seed=args.seed,
    )
    curve_lr = gsu.evaluate_k_curve(
        rank_df=rank_df,
        x_train_df=x_train_df,
        y_train=y_train,
        x_test_df=x_test_df,
        y_test=y_test,
        ks=pp.TOP_K_GRID,
        clf_kind="lr",
        seed=args.seed,
    )

    history_path = raw_dir / f"training_history_{spec.output_stub}.csv"
    rank_path = raw_dir / f"ranked_genes_{spec.output_stub}.csv"
    curve_knn_path = raw_dir / f"k_curve_knn_{spec.output_stub}.csv"
    curve_lr_path = raw_dir / f"k_curve_lr_{spec.output_stub}.csv"

    history_df.to_csv(history_path, index=False)
    rank_df.to_csv(rank_path, index=False)
    curve_knn.to_csv(curve_knn_path, index=False)
    curve_lr.to_csv(curve_lr_path, index=False)

    gsu.plot_metric_curve(curve_knn, "AUROC", pp.FIG03_PLOTS / f"AUROC_vs_k_knn_{spec.output_stub}.png", spec.display_name)
    gsu.plot_metric_curve(curve_knn, "F1", pp.FIG03_PLOTS / f"F1_vs_k_knn_{spec.output_stub}.png", spec.display_name)
    gsu.plot_metric_curve(curve_lr, "AUROC", pp.FIG03_PLOTS / f"AUROC_vs_k_lr_{spec.output_stub}.png", spec.display_name)
    gsu.plot_metric_curve(curve_lr, "F1", pp.FIG03_PLOTS / f"F1_vs_k_lr_{spec.output_stub}.png", spec.display_name)

    result = {
        "cell_type_key": spec.key,
        "cell_type_display": spec.display_name,
        "output_stub": spec.output_stub,
        "input_csv": str(input_csv),
        "device": cfg.device,
        "n_cells_total": metadata["n_cells_total"],
        "n_train": metadata["n_train"],
        "n_test": metadata["n_test"],
        "n_genes_raw": metadata["n_genes_raw"],
        "n_genes_model": metadata["n_genes_model"],
        "train_label_0": metadata["train_label_0"],
        "train_label_1": metadata["train_label_1"],
        "test_label_0": metadata["test_label_0"],
        "test_label_1": metadata["test_label_1"],
        "max_val_auroc": float(history_df["val_auroc"].max(skipna=True)),
        "max_val_auprc": float(history_df["val_auprc"].max(skipna=True)),
        "history_csv": str(history_path),
        "rank_csv": str(rank_path),
        "k_curve_knn_csv": str(curve_knn_path),
        "k_curve_lr_csv": str(curve_lr_path),
    }
    print(
        f"[DONE] {spec.key}: train={metadata['n_train']} test={metadata['n_test']} "
        f"genes={metadata['n_genes_model']} max_val_auroc={result['max_val_auroc']:.4f}",
        flush=True,
    )
    return result


def main() -> None:
    args = parse_args()
    pp.ensure_project_dirs()
    raw_dir = Path(args.raw_dir).expanduser().resolve() if args.raw_dir else pp.FIG03_RAW
    raw_dir.mkdir(parents=True, exist_ok=True)

    selected_specs = pp.resolve_celltypes(args.cell_types)
    manifest_rows = []
    manifest_path = pp.FIG03_TABLES / "celltype_run_manifest.csv"

    for spec in selected_specs:
        manifest_rows.append(run_one_celltype(spec, args=args, raw_dir=raw_dir))
        pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(manifest_path, index=False)

    config_payload = {
        "selected_cell_types": [spec.key for spec in selected_specs],
        "raw_dir": str(raw_dir),
        "fig03_root": str(pp.FIG03_ROOT),
        "args": vars(args),
        "python": sys.version,
        "platform": platform.platform(),
    }
    config_path = pp.FIG03_LOGS / "fig03_python_run_config.json"
    with config_path.open("w", encoding="utf-8") as handle:
        json.dump(config_payload, handle, indent=2)

    print(manifest_df.to_string(index=False))
    print(f"\nWrote manifest to: {manifest_path}")
    print(f"Wrote config to: {config_path}")


if __name__ == "__main__":
    main()
