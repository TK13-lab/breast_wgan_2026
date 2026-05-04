#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

import project_paths as pp


MPL_CACHE_DIR = pp.PROJECT_ROOT / "cache" / "matplotlib"
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE_DIR))

import matplotlib.pyplot as plt


CAPABILITIES = [
    "Cell-type-aware supervision",
    "Compact gene ranking",
    "Latent representation",
    "Batch-aware modeling",
    "Perturbation prediction",
    "Direct bulk projection",
    "Survival endpoint",
]

METHOD_MATRIX = {
    "Proposed lineage-aware GAN": [1, 1, 1, 0, 0, 1, 1],
    "scVI": [0, 0, 1, 1, 0, 0, 0],
    "scGen": [0, 0, 1, 1, 1, 0, 0],
    "DEG effect ranking": [1, 1, 0, 0, 0, 1, 0],
    "LASSO Cox": [0, 1, 0, 0, 0, 1, 1],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a method-family comparison heatmap for the revision. "
            "This is a conceptual companion to the quantitative benchmark tables."
        )
    )
    parser.add_argument("--out-dir", default=None, help="Output directory for method-family table and plot.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else pp.OUTPUTS_ROOT / "revision" / "benchmark"
    table_dir = out_dir / "tables"
    plot_dir = out_dir / "plots"
    table_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame.from_dict(METHOD_MATRIX, orient="index", columns=CAPABILITIES)
    df.index.name = "method"
    df.to_csv(table_dir / "method_family_comparison.csv")

    fig, ax = plt.subplots(figsize=(11, 4.8))
    matrix = df.values.astype(float)
    ax.imshow(matrix, cmap="Greens", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(CAPABILITIES)))
    ax.set_xticklabels(CAPABILITIES, rotation=35, ha="right")
    ax.set_yticks(np.arange(df.shape[0]))
    ax.set_yticklabels(df.index)

    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            ax.text(col, row, "Yes" if matrix[row, col] == 1 else "-", ha="center", va="center", fontsize=9)

    ax.set_title("Method-family comparison for revision response")
    ax.tick_params(length=0)
    fig.tight_layout()
    fig.savefig(plot_dir / "revision_method_family_comparison.png", dpi=300)
    fig.savefig(plot_dir / "revision_method_family_comparison.pdf")
    plt.close(fig)

    print(f"Wrote method comparison table to: {table_dir / 'method_family_comparison.csv'}")
    print(f"Wrote method comparison plot to: {plot_dir}")


if __name__ == "__main__":
    main()
