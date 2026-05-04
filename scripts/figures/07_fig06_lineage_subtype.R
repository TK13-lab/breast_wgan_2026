source(file.path(getwd(), "scripts", "00_project_paths.R"))

paths <- paper3_init()
fig_paths <- paper3_init_fig("fig06_lineage_subtype_synthesis", paths = paths)

paper3_require_packages(c(
  "dplyr", "tidyr", "purrr", "stringr", "forcats", "ggplot2", "scales",
  "ComplexHeatmap", "grid", "gridExtra", "patchwork", "msigdbr", "tibble"
))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(ComplexHeatmap)
  library(grid)
  library(gridExtra)
  library(patchwork)
  library(msigdbr)
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

sig_path <- paper3_fig05_subtype_gene_path(paths = paths, required = TRUE)
gene_list_path <- paper3_fig03_gene_list_path(paths = paths, required = TRUE)

sig_tbl <- readRDS(sig_path)
gene_lists_raw <- readRDS(gene_list_path)

subtype_order <- c("Basal", "Her2", "LumA", "LumB")
lineage_order <- c("Epithelial", "Immune", "Stromal", "Vascular")

lineage_cols <- c(
  "Epithelial" = "#33A02C",
  "Immune" = "#FF7F00",
  "Stromal" = "#1F78B4",
  "Vascular" = "#A6CEE3"
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
ct_by_lineage <- split(names(lineage_map), lineage_map)
subtype_cols <- c(
  "Basal" = "#E41A1C",
  "Her2" = "#377EB8",
  "LumA" = "#4DAF4A",
  "LumB" = "#984EA3"
)
pathways_order <- c(
  "HALLMARK_APOPTOSIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "KEGG_ALLOGRAFT_REJECTION",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "KEGG_ASTHMA",
  "KEGG_AUTOIMMUNE_THYROID_DISEASE",
  "KEGG_CELL_ADHESION_MOLECULES_CAMS",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "KEGG_FOCAL_ADHESION",
  "KEGG_GRAFT_VERSUS_HOST_DISEASE",
  "KEGG_LEISHMANIA_INFECTION",
  "KEGG_TYPE_I_DIABETES_MELLITUS",
  "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_BIOSYNTHESIS",
  "KEGG_VIRAL_MYOCARDITIS"
)

`%||%` <- function(a, b) if (!is.null(a)) a else b

parse_genes <- function(x) {
  if (is.null(x)) {
    return(character(0))
  }
  if (is.list(x)) {
    x <- unlist(x, use.names = FALSE)
  }
  if (is.character(x) && length(x) > 1) {
    return(unique(trimws(x[nzchar(x)])))
  }
  if (is.character(x) && length(x) == 1) {
    y <- unlist(str_split(x, pattern = "[,;|\\s]+"), use.names = FALSE)
    y <- trimws(y)
    return(unique(y[nzchar(y)]))
  }
  character(0)
}

upper_unique <- function(v) {
  unique(toupper(as.character(v)))
}

pretty_pathway_label <- function(x) {
  x %>%
    str_remove("^HALLMARK_") %>%
    str_remove("^KEGG_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

norm_ct <- function(nm) {
  x <- tolower(nm)
  x <- gsub("[._]+", " ", x)
  x <- trimws(x)
  if (grepl("endothel", x)) return("Endothelial")
  if (grepl("fibro", x)) return("Fibroblast")
  if (grepl("pericy", x) || grepl("smc", x) || grepl("smooth", x)) return("Pericytes_SMC")
  if (grepl("epithel", x) || grepl("^epi\\b", x)) return("Epithelial")
  if (grepl("myelo", x) || grepl("macro", x) || grepl("mono", x) || grepl("neutro", x)) return("Myeloid")
  if (grepl("mast", x)) return("Mast cells")
  if (grepl("\\bb\\b", x) || grepl("b cell", x) || grepl("bcell", x)) return("B cells")
  if (grepl("\\bt\\b", x) || grepl("t cell", x) || grepl("tcell", x)) return("T cell")
  NA_character_
}

canon_gene_lists <- function(gene_lists_raw, cell_cols) {
  old_nm <- names(gene_lists_raw)
  new_nm <- vapply(old_nm, norm_ct, character(1))
  keep <- !is.na(new_nm)
  tmp <- gene_lists_raw[keep]
  tmp_nm <- new_nm[keep]
  idx_split <- split(seq_along(tmp), tmp_nm)
  gene_lists <- imap(idx_split, function(idx, ct) {
    upper_unique(unlist(tmp[idx], use.names = FALSE))
  })
  gene_lists[intersect(names(cell_cols), names(gene_lists))]
}

build_sig_long <- function(sig_tbl) {
  sig_tbl %>%
    as_tibble() %>%
    transmute(
      subtype = as.character(.data$subtype),
      lineage = as.character(.data$lineage),
      genes_v = purrr::map(.data$genes, parse_genes)
    ) %>%
    tidyr::unnest_longer(genes_v, values_to = "gene") %>%
    filter(!is.na(gene), gene != "") %>%
    mutate(
      gene = toupper(gene),
      subtype = factor(subtype, levels = subtype_order),
      lineage = factor(lineage, levels = lineage_order)
    ) %>%
    filter(!is.na(subtype), !is.na(lineage)) %>%
    distinct(subtype, lineage, gene)
}

sig_long <- build_sig_long(sig_tbl)
stopifnot(nrow(sig_long) > 0)
gene_lists <- canon_gene_lists(gene_lists_raw, cell_cols)
stopifnot(length(gene_lists) >= 2)

paper3_write_csv(
  sig_long %>% count(subtype, lineage, name = "n_genes"),
  file.path(fig_paths$tables, "subtype_lineage_gene_counts.csv")
)
paper3_write_csv(
  data.frame(cell_type = names(gene_lists), n_genes = lengths(gene_lists), stringsAsFactors = FALSE),
  file.path(fig_paths$tables, "celltype_gene_list_sizes.csv")
)

rank_maps <- lapply(gene_lists, function(v) {
  r <- seq_along(v)
  names(r) <- v
  r
})

term2gene <- bind_rows(
  msigdbr(species = "Homo sapiens", collection = "H") %>% transmute(pathway = gs_name, gene = toupper(gene_symbol)),
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>% transmute(pathway = gs_name, gene = toupper(gene_symbol))
) %>%
  filter(pathway %in% pathways_order) %>%
  distinct()

keep_top_combinations <- function(cm, topN = 8) {
  cs <- ComplexHeatmap::comb_size(cm)
  idx <- which(cs > 0)
  if (length(idx) == 0) {
    return(cm)
  }
  idx <- idx[order(cs[idx], decreasing = TRUE)]
  idx <- idx[seq_len(min(topN, length(idx)))]
  cm[, idx]
}

make_upset_one_subtype_grob <- function(sub, topN = 8, pt_mm = 2.2, axis_fs = 11, title_fs = 13) {
  df_sub <- sig_long %>% filter(subtype == sub)
  sets <- setNames(vector("list", length(lineage_order)), lineage_order)
  for (lin in lineage_order) {
    sets[[lin]] <- df_sub %>% filter(lineage == lin) %>% pull(gene) %>% unique()
  }
  if (sum(lengths(sets) > 0) < 2) {
    return(grid::textGrob(
      paste0(as.character(sub), "\n(no enough lineage sets)"),
      x = 0,
      y = 1,
      just = c("left", "top"),
      gp = gpar(fontsize = title_fs, fontface = "bold")
    ))
  }

  cm <- make_comb_mat(sets)
  cm2 <- keep_top_combinations(cm, topN = topN)

  top_anno <- HeatmapAnnotation(
    "Intersection size" = anno_barplot(
      comb_size(cm2),
      border = FALSE,
      gp = gpar(fill = "grey60", col = NA),
      axis_param = list(side = "left", labels_rot = 0, gp = gpar(fontsize = axis_fs))
    ),
    annotation_name_side = "left",
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = axis_fs + 1, fontface = "bold")
  )

  right_anno <- rowAnnotation(
    "Set size" = anno_barplot(
      set_size(cm2),
      border = FALSE,
      gp = gpar(fill = unname(lineage_cols[lineage_order]), col = NA),
      axis_param = list(side = "bottom", labels_rot = 0, gp = gpar(fontsize = axis_fs))
    ),
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = axis_fs + 1, fontface = "bold"),
    width = unit(2.5, "cm")
  )

  ht <- UpSet(
    cm2,
    set_order = lineage_order,
    top_annotation = top_anno,
    right_annotation = right_anno,
    pt_size = unit(pt_mm, "mm")
  )

  grid.grabExpr({
    draw(ht)
    grid.text(
      label = as.character(sub),
      x = unit(0, "npc") + unit(1.5, "mm"),
      y = unit(1, "npc") - unit(1.5, "mm"),
      just = c("left", "top"),
      gp = gpar(fontsize = title_fs, fontface = "bold")
    )
  })
}

gA <- make_upset_one_subtype_grob("Basal", topN = 8)
gB <- make_upset_one_subtype_grob("Her2", topN = 8)
gC <- make_upset_one_subtype_grob("LumA", topN = 8)
gD <- make_upset_one_subtype_grob("LumB", topN = 8)

fig6A_core <- gridExtra::arrangeGrob(gA, gB, gC, gD, ncol = 2)
fig6A <- gridExtra::arrangeGrob(
  fig6A_core,
  top = grid::textGrob("(A)", x = 0, just = "left", gp = gpar(fontsize = 16, fontface = "bold"))
)

map_df <- sig_long %>%
  inner_join(term2gene, by = "gene", relationship = "many-to-many") %>%
  distinct(subtype, lineage, pathway, gene)

supp_table_dest <- file.path(fig_paths$tables, "supp_table_fig6b_by_lineage.csv")
legacy_supp_table_path <- paper3_fig06_supp_table_path(paths = paths, required = FALSE)
if (file.exists(legacy_supp_table_path)) {
  paper3_stage_reference_file(
    source_path = legacy_supp_table_path,
    destination_path = supp_table_dest,
    label = "legacy SuppTable_Fig6B_by_lineage.csv"
  )
  pw_sum_raw <- read.csv(supp_table_dest, check.names = FALSE) %>% as_tibble()
  gene_col <- if ("genes_all" %in% names(pw_sum_raw)) pw_sum_raw$genes_all else pw_sum_raw$genes
  pw_sum <- pw_sum_raw %>%
    transmute(
      subtype = factor(as.character(subtype), levels = subtype_order),
      lineage = factor(as.character(lineage), levels = lineage_order),
      pathway = as.character(pathway),
      n_genes = as.numeric(n_genes),
      genes = purrr::map(gene_col, parse_genes)
    ) %>%
    filter(!is.na(subtype), !is.na(lineage))
  supp_table_source <- legacy_supp_table_path
} else {
  pw_sum <- map_df %>%
    group_by(subtype, lineage, pathway) %>%
    summarise(n_genes = n_distinct(gene), genes = list(sort(unique(gene))), .groups = "drop")

  paper3_write_csv(
    pw_sum %>% mutate(genes_collapsed = vapply(genes, paste, collapse = "; ", character(1))),
    supp_table_dest
  )
  supp_table_source <- "computed_from_sig_long"
}

pw_keep <- pw_sum %>%
  group_by(pathway) %>%
  summarise(total = sum(n_genes), .groups = "drop") %>%
  filter(total > 0) %>%
  pull(pathway)
pw_sum <- pw_sum %>% filter(pathway %in% pw_keep)

top_pw_per_sub <- pw_sum %>%
  group_by(subtype, pathway) %>%
  summarise(total = sum(n_genes), .groups = "drop") %>%
  group_by(subtype) %>%
  slice_max(total, n = 8, with_ties = TRUE) %>%
  pull(pathway) %>%
  unique()

df_plot_B <- pw_sum %>%
  filter(pathway %in% top_pw_per_sub) %>%
  mutate(
    subtype = factor(as.character(subtype), levels = subtype_order),
    lineage = factor(as.character(lineage), levels = lineage_order),
    pathway = as.character(pathway)
  )

pw_levels <- df_plot_B %>%
  group_by(pathway) %>%
  summarise(total = sum(n_genes, na.rm = TRUE), .groups = "drop") %>%
  arrange(total) %>%
  pull(pathway)

df_plot_B <- df_plot_B %>%
  mutate(pathway = factor(pathway, levels = pw_levels)) %>%
  tidyr::complete(
    subtype = factor(subtype_order, levels = subtype_order),
    lineage = factor(lineage_order, levels = lineage_order),
    pathway = factor(pw_levels, levels = pw_levels),
    fill = list(n_genes = 0)
  )

x_max <- max(df_plot_B$n_genes, na.rm = TRUE)

make_one_sub_bar <- function(df, sub, show_y = FALSE, show_legend = FALSE) {
  dsub <- df %>% filter(subtype == sub)
  ggplot(dsub, aes(x = n_genes, y = fct_rev(pathway), fill = lineage)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = lineage_cols, drop = FALSE, name = "Lineage") +
    scale_x_continuous(limits = c(0, x_max * 1.05), expand = expansion(mult = c(0, 0.02))) +
    labs(title = sub, x = NULL, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position = if (show_legend) "right" else "none",
      axis.text.x = element_text(size = 11),
      axis.text.y = if (show_y) element_text(size = 10) else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      plot.margin = margin(4, 6, 4, 4)
    )
}

p_basal <- make_one_sub_bar(df_plot_B, "Basal", show_y = TRUE, show_legend = TRUE)
p_her2 <- make_one_sub_bar(df_plot_B, "Her2", show_y = FALSE, show_legend = FALSE)
p_luma <- make_one_sub_bar(df_plot_B, "LumA", show_y = FALSE, show_legend = FALSE)
p_lumb <- make_one_sub_bar(df_plot_B, "LumB", show_y = FALSE, show_legend = FALSE)

fig6B_core <- (p_basal | p_her2 | p_luma | p_lumb) +
  plot_layout(nrow = 1, widths = c(1.25, 1, 1, 1), guides = "collect") &
  theme(legend.position = "right")

fig6B <- gridExtra::arrangeGrob(
  grid::textGrob("(B)", x = 0, just = "left", gp = gpar(fontsize = 14, fontface = "bold")),
  patchwork::patchworkGrob(fig6B_core),
  ncol = 1,
  heights = c(0.10, 1)
)

ct_order <- names(cell_cols)
sig_sets <- sig_long %>%
  group_by(subtype, lineage) %>%
  summarise(sig_genes = list(unique(gene)), n_sig = n_distinct(gene), .groups = "drop")

ct_meta <- tibble(cell_type = ct_order, lineage_ct = unname(lineage_map[ct_order])) %>%
  filter(cell_type %in% names(gene_lists))

plot_df_C <- tidyr::crossing(subtype = subtype_order, ct_meta) %>%
  left_join(sig_sets, by = c("subtype" = "subtype", "lineage_ct" = "lineage")) %>%
  mutate(
    sig_genes = purrr::map2(sig_genes, n_sig, function(genes, n) {
      if (is.na(n)) {
        character(0)
      } else {
        genes
      }
    }),
    n_sig = ifelse(is.na(n_sig), 0L, n_sig),
    topk = map(cell_type, \(ct) gene_lists[[ct]]),
    n_hit = map2_int(topk, sig_genes, \(a, b) length(intersect(a, b))),
    pct_hit = ifelse(n_sig > 0, n_hit / n_sig, NA_real_)
  )

fig6C <- ggplot(
  plot_df_C,
  aes(x = factor(subtype, levels = subtype_order), y = factor(cell_type, levels = rev(ct_order)))
) +
  geom_tile(aes(fill = pct_hit), colour = "grey88", linewidth = 0.35, na.rm = FALSE) +
  geom_text(aes(label = ifelse(is.na(pct_hit), "", n_hit)), size = 3.3) +
  facet_grid(lineage_ct ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient(
    low = "white",
    high = "#2c7fb8",
    limits = c(0, 1),
    oob = squish,
    na.value = "grey85",
    labels = percent_format(accuracy = 1),
    name = "Recovered\n(% of sig genes)"
  ) +
  theme_classic(base_size = 12) +
  labs(title = "(C)") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    strip.background = element_blank(),
    strip.text.y.left = element_text(face = "bold", size = 12),
    legend.position = "right",
    plot.margin = margin(6, 8, 6, 6)
  )

pick_best_ct <- function(gene, lin, gene_lists, rank_maps) {
  allowed <- ct_by_lineage[[lin]] %||% character(0)
  cands <- intersect(allowed, names(gene_lists))
  if (length(cands) == 0) {
    cands <- names(gene_lists)
  }
  ranks <- map_dbl(cands, function(ct) {
    rm <- rank_maps[[ct]]
    v <- unname(rm[gene])
    ifelse(is.na(v), Inf, as.numeric(v))
  })
  if (all(!is.finite(ranks))) {
    return(NA_character_)
  }
  cands[which.min(ranks)]
}

build_flows_one_subtype <- function(sub, top_pathways = 10, max_genes = 40) {
  sg <- sig_long %>% filter(subtype == sub)
  gene_ct <- sg %>%
    distinct(lineage, gene) %>%
    mutate(cell_type = map2_chr(gene, lineage, ~ pick_best_ct(.x, .y, gene_lists, rank_maps))) %>%
    filter(!is.na(cell_type)) %>%
    mutate(
      rank_in_ct = map2_dbl(gene, cell_type, \(g, ct) {
        rm <- rank_maps[[ct]]
        v <- unname(rm[g])
        ifelse(is.na(v), Inf, as.numeric(v))
      })
    )

  gene_pw <- sg %>% distinct(gene) %>% inner_join(term2gene, by = "gene")
  if (nrow(gene_pw) == 0 || nrow(gene_ct) == 0) {
    return(tibble(subtype = sub, cell_type = character(0), gene = character(0), pathway = character(0), weight = numeric(0)))
  }

  pw_counts <- gene_pw %>% count(pathway, sort = TRUE)
  keep_pw <- pw_counts$pathway[seq_len(min(top_pathways, nrow(pw_counts)))]
  gene_pw <- gene_pw %>% filter(pathway %in% keep_pw)
  pw_rank <- setNames(seq_along(pathways_order), pathways_order)

  gene_pw1 <- gene_pw %>%
    left_join(pw_counts, by = "pathway") %>%
    mutate(ord = pw_rank[pathway]) %>%
    arrange(gene, desc(n), ord) %>%
    group_by(gene) %>%
    slice(1) %>%
    ungroup() %>%
    select(gene, pathway)

  flows <- gene_ct %>%
    inner_join(gene_pw1, by = "gene") %>%
    transmute(subtype = sub, cell_type = cell_type, gene = gene, pathway = pathway, weight = 1, rank_in_ct = rank_in_ct)

  if (n_distinct(flows$gene) > max_genes) {
    keep_genes <- flows %>%
      arrange(rank_in_ct) %>%
      distinct(gene, .keep_all = TRUE) %>%
      slice_head(n = max_genes) %>%
      pull(gene)
    flows <- flows %>% filter(gene %in% keep_genes)
  }
  flows %>% select(subtype, cell_type, gene, pathway, weight)
}

flows_sub <- map_dfr(subtype_order, build_flows_one_subtype)
paper3_write_csv(flows_sub, file.path(fig_paths$tables, "sankey_flows.csv"))

panel_d_summary <- flows_sub %>%
  group_by(subtype, cell_type, pathway) %>%
  summarise(
    n_genes = n_distinct(gene),
    genes = paste(sort(unique(gene)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    subtype = factor(subtype, levels = subtype_order),
    cell_type = factor(cell_type, levels = rev(ct_order)),
    lineage = unname(lineage_map[as.character(cell_type)])
  ) %>%
  filter(!is.na(subtype), !is.na(cell_type), !is.na(lineage))

paper3_write_csv(
  panel_d_summary %>%
    transmute(
      subtype = as.character(subtype),
      cell_type = as.character(cell_type),
      lineage = lineage,
      pathway = as.character(pathway),
      n_genes = n_genes,
      genes = genes
    ),
  file.path(fig_paths$tables, "panel_d_celltype_pathway_support.csv")
)

panel_d_pathways <- unique(c(
  intersect(pathways_order, unique(as.character(panel_d_summary$pathway))),
  setdiff(unique(as.character(panel_d_summary$pathway)), pathways_order)
))

panel_d_plot <- tidyr::expand_grid(
  subtype = factor(subtype_order, levels = subtype_order),
  cell_type = factor(rev(ct_order), levels = rev(ct_order)),
  pathway = factor(panel_d_pathways, levels = panel_d_pathways)
) %>%
  left_join(panel_d_summary, by = c("subtype", "cell_type", "pathway")) %>%
  mutate(
    lineage = coalesce(lineage, unname(lineage_map[as.character(cell_type)])),
    n_genes = coalesce(n_genes, 0L),
    genes = coalesce(genes, ""),
    pathway_label = pretty_pathway_label(as.character(pathway))
  ) %>%
  select(subtype, cell_type, pathway, pathway_label, lineage, n_genes, genes)

fig6D <- ggplot(panel_d_plot, aes(x = pathway, y = cell_type)) +
  geom_point(
    shape = 21,
    size = 2.6,
    stroke = 0.2,
    fill = "white",
    colour = "grey88"
  ) +
  geom_point(
    data = panel_d_plot %>% filter(n_genes > 0),
    aes(size = n_genes, fill = cell_type),
    shape = 21,
    stroke = 0.3,
    colour = "grey30",
    alpha = 0.95
  ) +
  geom_text(
    data = panel_d_plot %>% filter(n_genes > 0),
    aes(label = n_genes),
    size = 3,
    fontface = "bold"
  ) +
  facet_wrap(~ subtype, nrow = 1) +
  scale_fill_manual(values = cell_cols, name = "Cell type") +
  scale_size_continuous(
    range = c(5.2, 8.6),
    breaks = sort(unique(panel_d_summary$n_genes)),
    limits = c(min(panel_d_summary$n_genes), max(panel_d_summary$n_genes)),
    name = "# linked genes"
  ) +
  scale_x_discrete(labels = function(x) str_wrap(pretty_pathway_label(x), width = 18)) +
  coord_cartesian(clip = "off") +
  labs(title = "(D)", x = NULL, y = NULL) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 5.5)),
    size = guide_legend(order = 2)
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.x = element_text(size = 9.5, angle = 38, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "grey97", colour = "grey85"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right",
    panel.border = element_rect(colour = "grey88", fill = NA),
    panel.spacing.x = unit(6, "mm"),
    plot.margin = margin(6, 8, 6, 6)
  )

save_plot_bundle(wrap_elements(full = fig6A), "fig06_panel_A_upset", width = 11, height = 8.5)
save_plot_bundle(wrap_elements(full = fig6B), "fig06_panel_B_pathway_barplots", width = 16, height = 6)
save_plot_bundle(fig6C, "fig06_panel_C_backmapping_heatmap", width = 8, height = 7)
save_plot_bundle(fig6D, "fig06_panel_D_celltype_pathway_support", width = 18, height = 7.2)
# Keep the legacy stem so downstream references do not point to a stale Sankey export.
save_plot_bundle(fig6D, "fig06_panel_D_sankey", width = 18, height = 7.2)

paper3_write_lines(
  c(
    "Fig. 6 workspace direct-run summary",
    sprintf("Subtype-lineage signatures: %s", sig_path),
    sprintf("Fig. 3 top-k gene lists: %s", gene_list_path),
    sprintf("Supplementary pathway table source: %s", supp_table_source),
    "All panels were regenerated directly in the workspace, with the pathway summary table reconciled to the legacy canonical export when that reference was available.",
    "The public release keeps panel-level exports only and omits the combined figure layout."
  ),
  file.path(fig_paths$logs, "provenance_notes.txt")
)
paper3_capture_session_info(file.path(fig_paths$logs, "session_info.txt"))

message("DONE: Fig. 6 regenerated in workspace outputs/figures/fig06_lineage_subtype_synthesis/")
