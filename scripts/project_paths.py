from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_LINKS = PROJECT_ROOT / "data_links"
LEGACY_ROOT = (DATA_LINKS / "paper3_breast_wgan_win").resolve()
OUTPUTS_ROOT = PROJECT_ROOT / "outputs"
LOGS_ROOT = PROJECT_ROOT / "logs"
CACHE_ROOT = PROJECT_ROOT / "cache"

FIG03_ROOT = OUTPUTS_ROOT / "figures" / "fig03_gan_signature_discovery"
FIG03_RAW = FIG03_ROOT / "raw"
FIG03_PLOTS = FIG03_ROOT / "plots"
FIG03_TABLES = FIG03_ROOT / "tables"
FIG03_RDS = FIG03_ROOT / "rds"
FIG03_LOGS = FIG03_ROOT / "logs"
FIG03_CHECKPOINTS = OUTPUTS_ROOT / "checkpoints" / "fig03_gan_signature_discovery"

TOP_K_GRID = (10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400, 500)
EXCLUDED_GENE_PREFIXES = ("RPS", "RPL")


@dataclass(frozen=True)
class CellTypeSpec:
    key: str
    display_name: str
    input_filename: str
    output_stub: str


CELLTYPE_SPECS = (
    CellTypeSpec("B_cells", "B cells", "B_cells_sigGene_expression_counts.csv", "B_cells"),
    CellTypeSpec("Endothelial", "Endothelial", "Endothelial_sigGene_expression_counts.csv", "Endothelial"),
    CellTypeSpec("Epithelial", "Epithelial", "Epithelial_sigGene_expression_counts.csv", "Epithelial"),
    CellTypeSpec("Fibroblast", "Fibroblast", "Fibroblast_sigGene_expression_counts.csv", "Fibroblast"),
    CellTypeSpec("Mast_cells", "Mast cells", "Mast_cells_sigGene_expression_counts.csv", "Mast_cells"),
    CellTypeSpec("Myeloid", "Myeloid", "Myeloid_sigGene_expression_counts.csv", "Myeloid"),
    CellTypeSpec("Pericytes_SMC", "Pericytes_SMC", "Pericytes_SMC_sigGene_expression_counts.csv", "Pericytes"),
    CellTypeSpec("T_cells", "T cell", "T_cell_sigGene_expression_counts.csv", "T_cells"),
)

PANEL_ORDER = (
    "Endothelial",
    "Fibroblast",
    "Pericytes_SMC",
    "Epithelial",
    "Myeloid",
    "T cell",
    "Mast cells",
    "B cells",
)

CELL_COLORS = {
    "Endothelial": "#A6CEE3",
    "Fibroblast": "#1F78B4",
    "Pericytes_SMC": "#B2DF8A",
    "Epithelial": "#33A02C",
    "Myeloid": "#FDBF6F",
    "T cell": "#FF7F00",
    "Mast cells": "#FB9A99",
    "B cells": "#E31A1C",
}

_CELLTYPE_INDEX = {}
for _spec in CELLTYPE_SPECS:
    for token in (_spec.key, _spec.display_name, _spec.output_stub):
        _CELLTYPE_INDEX[token.lower()] = _spec


def ensure_project_dirs() -> None:
    for path in (
        OUTPUTS_ROOT,
        LOGS_ROOT,
        CACHE_ROOT,
        FIG03_ROOT,
        FIG03_RAW,
        FIG03_PLOTS,
        FIG03_TABLES,
        FIG03_RDS,
        FIG03_LOGS,
        FIG03_CHECKPOINTS,
    ):
        path.mkdir(parents=True, exist_ok=True)


def resolve_celltypes(selected: Iterable[str] | None) -> list[CellTypeSpec]:
    if not selected:
        return list(CELLTYPE_SPECS)

    tokens = [token.strip() for token in selected if token and token.strip()]
    if not tokens or any(token.lower() == "all" for token in tokens):
        return list(CELLTYPE_SPECS)

    resolved: list[CellTypeSpec] = []
    seen: set[str] = set()

    for token in tokens:
        spec = _CELLTYPE_INDEX.get(token.lower())
        if spec is None:
            allowed = ", ".join(spec.key for spec in CELLTYPE_SPECS)
            raise ValueError(f"Unknown cell type '{token}'. Expected one of: {allowed}, all")
        if spec.key not in seen:
            resolved.append(spec)
            seen.add(spec.key)

    return resolved


def legacy_count_path(spec: CellTypeSpec) -> Path:
    return LEGACY_ROOT / spec.input_filename
