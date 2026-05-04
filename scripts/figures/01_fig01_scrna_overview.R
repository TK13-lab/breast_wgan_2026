source(file.path(getwd(), "scripts", "00_project_paths.R"))

paths <- paper3_init()
fig_paths <- paper3_init_fig("fig01_scrna_overview", paths = paths)

paper3_require_packages(c(
  "Seurat", "ggplot2", "dplyr", "tidyr", "ggpubr", "rstatix",
  "RColorBrewer", "scales", "patchwork", "cowplot", "tibble"
))
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(rstatix)
  library(RColorBrewer)
  library(scales)
  library(patchwork)
  library(cowplot)
  library(tibble)
})

save_plot_bundle <- function(plot_obj, stem, width, height, dpi = 300) {
  pdf_path <- file.path(fig_paths$plots, paste0(stem, ".pdf"))
  png_path <- file.path(fig_paths$plots, paste0(stem, ".png"))
  tiff_path <- file.path(fig_paths$plots, paste0(stem, ".tiff"))
  ggsave(pdf_path, plot = plot_obj, width = width, height = height, useDingbats = FALSE)
  ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = dpi)
  ggsave(tiff_path, plot = plot_obj, width = width, height = height, dpi = 600, compression = "lzw")
}

choose_reduction <- function(obj) {
  candidates <- c("UMAP", "umap", "HarmonyUMAP2D", "HarmonyUMAP", "tsne", "TSNE")
  available <- names(obj@reductions)
  hit <- candidates[candidates %in% available]
  if (length(hit) == 0) {
    stop(
      sprintf("No supported dimensional reduction found. Available reductions: %s", paste(available, collapse = ", ")),
      call. = FALSE
    )
  }
  hit[[1]]
}

tag_panel <- function(plot_obj, tag) {
  cowplot::ggdraw(plot_obj) +
    cowplot::draw_label(
      tag,
      x = 0.02,
      y = 0.98,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = 14
    )
}

input_rds <- paper3_legacy_file("breast_268662_scRNA_filter_SCP.rds", paths = paths)
rds_integrated <- readRDS(input_rds)

mc_levels <- c(
  "Endothelial", "Fibroblast", "Pericytes_SMC", "Epithelial",
  "Myeloid", "T cell", "Mast cells", "B cells"
)
cell_cols <- c(
  "Endothelial" = "#A6CEE3",
  "Fibroblast" = "#1F78B4",
  "Pericytes_SMC" = "#B2DF8A",
  "Epithelial" = "#33A02C",
  "Myeloid" = "#FDBF6F",
  "T cell" = "#FF7F00",
  "Mast cells" = "#FB9A99",
  "B cells" = "#E31A1C"
)
group_cols <- c(normal = "#66c2a5", tumor = "#fc8d62")
patient_levels <- c("BN1", "BN2", "DN", "EN", "BT", "DT", "ET", "FT", "GT")

if (!("major_celltype" %in% colnames(rds_integrated@meta.data))) {
  stop("Metadata column 'major_celltype' not found.", call. = FALSE)
}
if (!("group" %in% colnames(rds_integrated@meta.data))) {
  stop("Metadata column 'group' not found.", call. = FALSE)
}
if (!("Patient" %in% colnames(rds_integrated@meta.data))) {
  stop("Metadata column 'Patient' not found.", call. = FALSE)
}

mc_levels_present <- intersect(mc_levels, unique(as.character(rds_integrated$major_celltype)))
if (length(mc_levels_present) < 2) {
  stop("Unexpected major_celltype labels in the integrated object.", call. = FALSE)
}

rds_integrated$major_celltype <- factor(rds_integrated$major_celltype, levels = mc_levels_present)
rds_integrated$group <- factor(rds_integrated$group, levels = c("normal", "tumor"))
reduction_name <- choose_reduction(rds_integrated)

pA <- cowplot::ggdraw() +
  theme_void() +
  cowplot::draw_label(
    "(A) BioRender panel placeholder\nnot generated inside the workspace",
    x = 0.01,
    y = 0.98,
    hjust = 0,
    vjust = 1,
    size = 12
  )

pB_base <- DimPlot(
  rds_integrated,
  reduction = reduction_name,
  group.by = "major_celltype",
  cols = cell_cols[mc_levels_present],
  label = FALSE
) +
  theme_void(base_size = 12) +
  theme(legend.position = "right", plot.margin = margin(3, 3, 3, 3))

pC_base <- DimPlot(
  rds_integrated,
  reduction = reduction_name,
  group.by = "group",
  cols = group_cols,
  label = FALSE
) +
  theme_void(base_size = 12) +
  theme(legend.position = "right", legend.title = element_blank(), plot.margin = margin(3, 3, 3, 3))

meta <- rds_integrated@meta.data %>%
  dplyr::select(Patient, major_celltype, group) %>%
  mutate(
    Patient = factor(Patient, levels = intersect(patient_levels, unique(as.character(Patient)))),
    major_celltype = factor(major_celltype, levels = mc_levels_present),
    group = factor(group, levels = c("normal", "tumor"))
  )

cell_frac <- meta %>%
  count(Patient, group, major_celltype, name = "n_celltype") %>%
  group_by(Patient, group) %>%
  mutate(total = sum(n_celltype), fraction = n_celltype / total) %>%
  ungroup()

cell_counts <- meta %>%
  count(group, major_celltype, name = "n_cells") %>%
  mutate(major_celltype = factor(major_celltype, levels = mc_levels_present))

stat_test <- cell_frac %>%
  filter(group %in% c("normal", "tumor")) %>%
  group_by(major_celltype) %>%
  wilcox_test(fraction ~ group) %>%
  add_significance("p") %>%
  add_xy_position(x = "major_celltype", dodge = 0.8)

paper3_write_csv(cell_frac, file.path(fig_paths$tables, "cell_fraction_by_patient.csv"))
paper3_write_csv(cell_counts, file.path(fig_paths$tables, "cell_count_by_group.csv"))
paper3_write_csv(stat_test, file.path(fig_paths$tables, "cell_fraction_wilcox.csv"))

pD_base <- ggplot(cell_frac, aes(x = Patient, y = fraction, fill = major_celltype)) +
  geom_col(width = 0.95) +
  facet_wrap(~ group, scales = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = cell_cols[mc_levels_present], breaks = mc_levels_present, name = "Major cell type") +
  labs(x = "Patient", y = "Cell fraction") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    plot.margin = margin(3, 3, 3, 3)
  )

pE_base <- ggplot(cell_counts, aes(x = major_celltype, y = n_cells, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = group_cols) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "Number of cells") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(3, 3, 3, 3)
  )

pF_base <- ggplot(cell_frac, aes(x = major_celltype, y = fraction, fill = group)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.6,
    outlier.shape = NA,
    colour = "black"
  ) +
  geom_point(
    aes(color = group),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 1.7,
    alpha = 0.8
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(x = NULL, y = "Fraction") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(3, 3, 3, 3)
  ) +
  stat_pvalue_manual(
    stat_test,
    label = "p.signif",
    bracket.size = 0.35,
    tip.length = 0.01
  )

marker_features <- c(
  "PECAM1", "VWF", "COL1A1", "DCN", "RGS5", "ACTA2", "EPCAM", "KRT8",
  "LST1", "C5AR1", "CD3D", "CD8A", "TPSAB1", "CPA3", "MS4A1", "CD79A"
)
marker_features <- intersect(marker_features, rownames(rds_integrated))
if (length(marker_features) < 4) {
  stop("Too few marker genes are available for the Fig. 1 dot plot.", call. = FALSE)
}

Idents(rds_integrated) <- rds_integrated$major_celltype
pG_base <- DotPlot(rds_integrated, features = marker_features, dot.scale = 5) +
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    plot.margin = margin(3, 3, 3, 3)
  )

pB <- tag_panel(pB_base, "(B)")
pC <- tag_panel(pC_base, "(C)")
pD <- tag_panel(pD_base, "(D)")
pE <- tag_panel(pE_base, "(E)")
pF <- tag_panel(pF_base, "(F)")
pG <- tag_panel(pG_base, "(G)")

save_plot_bundle(pB_base, "fig01_panel_B_umap_celltype", width = 7, height = 6)
save_plot_bundle(pC_base, "fig01_panel_C_umap_group", width = 7, height = 6)
save_plot_bundle(pD_base, "fig01_panel_D_cell_fraction", width = 8, height = 5)
save_plot_bundle(pE_base, "fig01_panel_E_cell_counts", width = 8, height = 5)
save_plot_bundle(pF_base, "fig01_panel_F_boxplot_fraction", width = 8, height = 5)
save_plot_bundle(pG_base, "fig01_panel_G_dotplot_markers", width = 9, height = 6)

panel_inventory <- tibble::tribble(
  ~panel, ~status, ~source,
  "A", "placeholder", "workspace placeholder because original BioRender asset was not included in the legacy folder",
  "B", "rerun", sprintf("DimPlot from %s using reduction %s", basename(input_rds), reduction_name),
  "C", "rerun", sprintf("DimPlot from %s using reduction %s", basename(input_rds), reduction_name),
  "D", "rerun", "cell fractions per patient/group from integrated metadata",
  "E", "rerun", "group-level cell counts from integrated metadata",
  "F", "rerun", "Wilcoxon comparison of cell fractions by group",
  "G", "rerun", "marker gene dot plot generated from the integrated object"
)
paper3_write_csv(panel_inventory, file.path(fig_paths$tables, "panel_inventory.csv"))

paper3_write_lines(
  c(
    "Fig. 1 workspace direct-run summary",
    sprintf("Input object: %s", input_rds),
    sprintf("Reduction used: %s", reduction_name),
    "Panel A is a placeholder because the original BioRender panel was not available in the protected legacy directory.",
    "Panels B-G were regenerated directly in the workspace from the integrated single-cell object.",
    "The public release keeps panel-level exports only and does not write a combined figure canvas."
  ),
  file.path(fig_paths$logs, "provenance_notes.txt")
)
paper3_capture_session_info(file.path(fig_paths$logs, "session_info.txt"))

message("DONE: Fig. 1 regenerated in workspace outputs/figures/fig01_scrna_overview/")
