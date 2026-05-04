source(file.path(getwd(), "scripts", "00_project_paths.R"))

paths <- paper3_init()
fig_paths <- paper3_init_fig("fig04_pathway_lineage_analysis", paths = paths)

paper3_require_packages(c(
  "clusterProfiler", "msigdbr", "dplyr", "purrr", "tidyr", "tibble",
  "ComplexHeatmap", "circlize", "grid", "gridExtra", "ggplot2", "patchwork", "png"
))
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(patchwork)
  library(png)
})

save_plot_bundle <- function(plot_obj, stem, width, height, dpi = 300) {
  pdf_path <- file.path(fig_paths$plots, paste0(stem, ".pdf"))
  png_path <- file.path(fig_paths$plots, paste0(stem, ".png"))
  tiff_path <- file.path(fig_paths$plots, paste0(stem, ".tiff"))
  ggsave(pdf_path, plot = plot_obj, width = width, height = height, useDingbats = FALSE)
  ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = dpi)
  ggsave(tiff_path, plot = plot_obj, width = width, height = height, dpi = 600, compression = "lzw")
}

save_grob_bundle <- function(grob_obj, stem, width, height, dpi = 300) {
  pdf_path <- file.path(fig_paths$plots, paste0(stem, ".pdf"))
  png_path <- file.path(fig_paths$plots, paste0(stem, ".png"))
  tiff_path <- file.path(fig_paths$plots, paste0(stem, ".tiff"))

  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  grid::grid.newpage()
  grid::grid.draw(grob_obj)
  grDevices::dev.off()

  grDevices::png(png_path, width = width, height = height, units = "in", res = dpi)
  grid::grid.newpage()
  grid::grid.draw(grob_obj)
  grDevices::dev.off()

  grDevices::tiff(tiff_path, width = width, height = height, units = "in", res = 600, compression = "lzw")
  grid::grid.newpage()
  grid::grid.draw(grob_obj)
  grDevices::dev.off()
}

render_chord_panel_bundle <- function(stem, width, height, dpi = 300) {
  pdf_path <- file.path(fig_paths$plots, paste0(stem, ".pdf"))
  png_path <- file.path(fig_paths$plots, paste0(stem, ".png"))
  tiff_path <- file.path(fig_paths$plots, paste0(stem, ".tiff"))

  draw_chord_panel <- function() {
    grid::grid.newpage()
    circos.clear()
    circos.par(
      start.degree = 90,
      gap.after = 3,
      track.margin = c(0.01, 0.01),
      cell.padding = c(0, 0, 0, 0)
    )
    chordDiagram(
      x = links_all[, c("from", "to", "value")],
      order = sector_order,
      grid.col = grid_col,
      transparency = 0.7,
      annotationTrack = "grid",
      preAllocateTracks = list(
        list(track.height = 0.04),
        list(track.height = 0.06)
      )
    )
    circos.trackPlotRegion(
      track.index = 1,
      bg.border = NA,
      bg.col = lineage_ring_cols,
      panel.fun = function(x, y) {}
    )
    circos.trackPlotRegion(
      track.index = 2,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector_name <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        if (sector_name %in% pathways) {
          circos.text(
            x = mean(xlim),
            y = ylim[1] + 0.5,
            labels = as.character(path_ids[sector_name]),
            facing = "bending.inside",
            niceFacing = TRUE,
            cex = 1.2
          )
        }
      }
    )
    lgd_cells <- Legend(
      title = "Cell type",
      labels = ct_all,
      legend_gp = gpar(fill = cell_cols[ct_all]),
      nrow = 2,
      by_row = TRUE,
      title_position = "lefttop"
    )
    lgd_lineage <- Legend(
      title = "Lineage",
      labels = names(lineage_cols),
      legend_gp = gpar(fill = lineage_cols),
      nrow = 1,
      by_row = TRUE,
      title_position = "lefttop"
    )
    pushViewport(viewport(x = 0.5, y = 0.995, width = 1, height = 0.08, just = c("center", "top")))
    grid.draw(lgd_cells)
    popViewport()
    pushViewport(viewport(x = 0.5, y = 0.01, width = 1, height = 0.1, just = c("center", "bottom")))
    grid.draw(lgd_lineage)
    popViewport()
    circos.clear()
  }

  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  draw_chord_panel()
  grDevices::dev.off()

  grDevices::png(png_path, width = width, height = height, units = "in", res = dpi)
  draw_chord_panel()
  grDevices::dev.off()

  grDevices::tiff(tiff_path, width = width, height = height, units = "in", res = 600, compression = "lzw")
  draw_chord_panel()
  grDevices::dev.off()

  normalizePath(png_path, winslash = "/", mustWork = TRUE)
}

gene_list_path <- paper3_fig03_gene_list_path(paths = paths, required = TRUE)
gene_lists <- readRDS(gene_list_path)

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
ct_all <- intersect(names(cell_cols), names(gene_lists))
if (length(ct_all) < 4) {
  stop("Fig. 4 needs the canonical Fig. 3 top-k gene list object in the workspace.", call. = FALSE)
}

lineage_map <- c(
  "Epithelial" = "Epithelial",
  "B cells" = "Immune",
  "T cell" = "Immune",
  "Myeloid" = "Immune",
  "Mast cells" = "Immune",
  "Fibroblast" = "Stromal",
  "Pericytes_SMC" = "Vascular",
  "Endothelial" = "Vascular"
)
lineage_cols <- c(
  "Epithelial" = "#1b9e77",
  "Immune" = "#d95f02",
  "Stromal" = "#7570b3",
  "Vascular" = "#e7298a"
)
lineage_order <- c("Epithelial", "Immune", "Stromal", "Vascular")

msig_all <- bind_rows(
  msigdbr(species = "Homo sapiens", collection = "H"),
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
)
term2gene <- msig_all %>% select(gs_name, gene_symbol)
term2name <- msig_all %>% distinct(gs_name, gs_description)

run_enrich_cell <- function(genes) {
  enricher(
    gene = genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5
  )
}

enrich_tbl <- imap(gene_lists[ct_all], function(genes, cell_type) {
  eg <- run_enrich_cell(genes)
  if (is.null(eg)) {
    return(NULL)
  }
  as.data.frame(eg) %>% mutate(cell_type = cell_type)
}) %>%
  bind_rows() %>%
  mutate(score = -log10(qvalue))

if (nrow(enrich_tbl) == 0) {
  stop("Fig. 4 enrichment produced no rows.", call. = FALSE)
}

paper3_write_csv(enrich_tbl, file.path(fig_paths$tables, "enrichment_table_full.csv"))

top2_per_ct <- enrich_tbl %>%
  filter(cell_type %in% ct_all) %>%
  group_by(cell_type) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup()

enrich_sig <- enrich_tbl %>%
  filter(cell_type %in% ct_all, p.adjust < 0.05, qvalue < 0.25)

shared_pathways <- enrich_sig %>%
  group_by(ID) %>%
  summarise(
    n_celltypes = n_distinct(cell_type),
    pathway_description = dplyr::first(Description),
    .groups = "drop"
  ) %>%
  filter(n_celltypes >= 2) %>%
  arrange(desc(n_celltypes), ID)

shared_links <- enrich_sig %>%
  filter(ID %in% shared_pathways$ID)

links_all <- bind_rows(
  top2_per_ct %>% mutate(flag = "top2"),
  shared_links %>% mutate(flag = "shared")
) %>%
  distinct(cell_type, ID, .keep_all = TRUE) %>%
  transmute(
    from = cell_type,
    to = ID,
    value = score,
    flag = flag
  )

paper3_write_csv(shared_pathways, file.path(fig_paths$tables, "shared_pathways.csv"))
paper3_write_csv(links_all, file.path(fig_paths$tables, "chord_links.csv"))

pathways <- sort(unique(links_all$to))
path_ids <- setNames(seq_along(pathways), pathways)
paper3_write_csv(
  data.frame(
    pathway_id = seq_along(pathways),
    pathway = pathways,
    stringsAsFactors = FALSE
  ),
  file.path(fig_paths$tables, "pathway_id_legend.csv")
)

sector_order <- c(ct_all, pathways)
grid_col <- rep("grey80", length(sector_order))
names(grid_col) <- sector_order
grid_col[ct_all] <- cell_cols[ct_all]

lineage_ring_cols <- ifelse(
  sector_order %in% names(lineage_map),
  lineage_cols[lineage_map[sector_order]],
  NA_character_
)
names(lineage_ring_cols) <- sector_order

pathways_order <- c(
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_HYPOXIA",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "KEGG_ARGlNINE_AND_PROLINE_METABOLISM",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "KEGG_FOCAL_ADHESION",
  "KEGG_GLYCOSAMINOGLYCAN_BIOSYNTHESIS_KERATAN_SULFATE",
  "KEGG_LEISHMANIA_INFECTION",
  "KEGG_PHENYLALANINE_METABOLISM"
)

enrich_tmp <- enrich_tbl %>%
  filter(ID %in% pathways_order) %>%
  mutate(
    score = -log10(p.adjust),
    score = ifelse(!is.finite(score), max(score[is.finite(score)], na.rm = TRUE) + 1, score)
  ) %>%
  select(cell_type, ID, score)

heat_df <- enrich_tmp %>% complete(cell_type, ID, fill = list(score = 0))
heat_mat <- heat_df %>%
  mutate(ID = factor(ID, levels = pathways_order)) %>%
  arrange(ID, cell_type) %>%
  pivot_wider(names_from = cell_type, values_from = score) %>%
  column_to_rownames("ID") %>%
  as.matrix()

celltype_order <- intersect(
  c("Epithelial", "Myeloid", "T cell", "Mast cells", "B cells", "Fibroblast", "Endothelial", "Pericytes_SMC"),
  colnames(heat_mat)
)
heat_mat <- heat_mat[, celltype_order, drop = FALSE]
max_score <- quantile(heat_mat, 0.95, na.rm = TRUE)
col_fun <- colorRamp2(c(0, max_score / 2, max_score), c("white", "gold", "red"))

cell2lineage <- c(
  "Endothelial" = "Vascular",
  "Myeloid" = "Immune",
  "Fibroblast" = "Stromal",
  "T cell" = "Immune",
  "Pericytes_SMC" = "Vascular",
  "Mast cells" = "Immune",
  "Epithelial" = "Epithelial",
  "B cells" = "Immune"
)

ha_lineage <- HeatmapAnnotation(
  Lineage = factor(cell2lineage[celltype_order], levels = lineage_order),
  col = list(Lineage = lineage_cols),
  annotation_name_side = "left"
)

ht <- Heatmap(
  heat_mat,
  name = "-log10(FDR)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = ha_lineage,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  row_names_max_width = unit(6, "cm"),
  width = unit(6, "cm"),
  column_names_side = "bottom",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  border = TRUE,
  rect_gp = gpar(col = "grey70", lwd = 0.5)
)

lgd_celltype <- Legend(
  title = "Cell type",
  labels = celltype_order,
  legend_gp = gpar(fill = cell_cols[celltype_order])
)

heatmap_grob <- grid::grid.grabExpr({
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    annotation_legend_list = list(lgd_celltype)
  )
}, wrap = TRUE)

lineage_to_cells <- list(
  Epithelial = c("Epithelial"),
  Immune = c("Myeloid", "T cell", "Mast cells", "B cells"),
  Stromal = c("Fibroblast"),
  Vascular = c("Endothelial", "Pericytes_SMC")
)

computed_genes_by_lineage <- lapply(lineage_to_cells, function(cts) {
  cts <- intersect(cts, names(gene_lists))
  if (length(cts) == 0) {
    return(character(0))
  }
  unique(unlist(gene_lists[cts], use.names = FALSE))
})
computed_genes_by_lineage <- computed_genes_by_lineage[lineage_order]

legacy_genes_by_lineage_path <- paper3_legacy_path("genes_by_lineage.rds", paths = paths, required = FALSE)
genes_by_lineage_dest <- file.path(fig_paths$rds, "genes_by_lineage.rds")
if (file.exists(legacy_genes_by_lineage_path)) {
  paper3_stage_reference_file(
    source_path = legacy_genes_by_lineage_path,
    destination_path = genes_by_lineage_dest,
    label = "legacy genes_by_lineage.rds"
  )
  genes_by_lineage <- readRDS(genes_by_lineage_dest)
  genes_by_lineage_source <- legacy_genes_by_lineage_path
} else {
  genes_by_lineage <- computed_genes_by_lineage
  saveRDS(genes_by_lineage, file = genes_by_lineage_dest)
  genes_by_lineage_source <- paste("computed_from", gene_list_path)
}

if (file.exists(legacy_genes_by_lineage_path)) {
  paper3_write_csv(
    tibble(
      lineage = lineage_order,
      legacy_n = lengths(readRDS(legacy_genes_by_lineage_path)[lineage_order]),
      computed_n = lengths(computed_genes_by_lineage[lineage_order])
    ),
    file.path(fig_paths$tables, "lineage_gene_reconciliation.csv")
  )
}

paper3_write_csv(
  data.frame(
    lineage = names(genes_by_lineage),
    n_genes = lengths(genes_by_lineage),
    stringsAsFactors = FALSE
  ),
  file.path(fig_paths$tables, "gene_counts_by_lineage.csv")
)

pairwise_pairs <- combn(lineage_order, 2)
pairwise_df <- apply(pairwise_pairs, 2, function(pair) {
  tibble(
    lineage_1 = pair[[1]],
    lineage_2 = pair[[2]],
    n_overlap = length(intersect(genes_by_lineage[[pair[[1]]]], genes_by_lineage[[pair[[2]]]]))
  )
}) %>%
  bind_rows()
paper3_write_csv(pairwise_df, file.path(fig_paths$tables, "pairwise_overlap_counts.csv"))

pairwise_plot_df <- expand.grid(
  lineage_1 = lineage_order,
  lineage_2 = lineage_order,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    n_overlap = if (lineage_1 == lineage_2) {
      length(genes_by_lineage[[lineage_1]])
    } else {
      length(intersect(genes_by_lineage[[lineage_1]], genes_by_lineage[[lineage_2]]))
    }
  ) %>%
  ungroup()

pC <- ggplot(
  pairwise_plot_df,
  aes(x = factor(lineage_1, levels = lineage_order), y = factor(lineage_2, levels = rev(lineage_order)))
) +
  geom_tile(aes(fill = n_overlap), colour = "white") +
  geom_text(aes(label = n_overlap), size = 4) +
  scale_fill_gradient(low = "white", high = "#2c7fb8", name = "# shared genes") +
  labs(
    title = "(C) Lineage overlap summary",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0)
  )

