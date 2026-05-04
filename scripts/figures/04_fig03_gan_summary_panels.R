source(file.path(getwd(), "scripts", "00_project_paths.R"))

paths <- paper3_paths()
paper3_mkdir(paths$outputs, paths$logs, paths$cache, paths$output_figures, paths$output_checkpoints)
fig_paths <- paper3_init_fig("fig03_gan_signature_discovery", paths = paths)

source_raw_dir <- {
  args <- commandArgs(trailingOnly = TRUE)
  arg_raw <- if (length(args) >= 1) args[[1]] else ""
  env_raw <- Sys.getenv("PAPER3_FIG03_RAW_DIR", unset = "")
  workspace_raw <- fig_paths$raw
  chosen <- if (nzchar(arg_raw)) {
    arg_raw
  } else if (nzchar(env_raw)) {
    env_raw
  } else if (paper3_fig03_raw_complete(workspace_raw)) {
    workspace_raw
  } else {
    ""
  }
  if (!nzchar(chosen)) {
    stop(
      "No Fig. 3 raw directory was found. Pass the raw directory as the first argument or set PAPER3_FIG03_RAW_DIR.",
      call. = FALSE
    )
  }
  normalizePath(chosen, winslash = "/", mustWork = TRUE)
}

paper3_require_packages(c("readr", "dplyr", "purrr", "tidyr", "tibble", "ggplot2", "patchwork"))
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
has_enrich <- requireNamespace("clusterProfiler", quietly = TRUE) &&
  requireNamespace("msigdbr", quietly = TRUE)

cell_specs <- tibble::tribble(
  ~output_stub, ~cell_type,
  "Endothelial", "Endothelial",
  "Fibroblast", "Fibroblast",
  "Pericytes", "Pericytes_SMC",
  "Epithelial", "Epithelial",
  "Myeloid", "Myeloid",
  "T_cells", "T cell",
  "Mast_cells", "Mast cells",
  "B_cells", "B cells"
)

cell_type_levels <- cell_specs$cell_type
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

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Missing required file: %s", path), call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

save_plot2 <- function(plot_obj, filename, width = 11, height = 7) {
  pdf_path <- file.path(fig_paths$plots, paste0(filename, ".pdf"))
  png_path <- file.path(fig_paths$plots, paste0(filename, ".png"))

  if ("cairo_pdf" %in% ls(getNamespace("grDevices"))) {
    ggsave(pdf_path, plot = plot_obj, width = width, height = height, device = cairo_pdf)
  } else {
    ggsave(pdf_path, plot = plot_obj, width = width, height = height, device = "pdf")
  }
  ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = 300)
}

rescale_safe <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(1, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

copy_raw_inputs <- function(source_dir, target_dir) {
  pattern <- paste(
    c(
      "^training_history_.*\\.csv$",
      "^ranked_genes_.*\\.csv$",
      "^k_curve_knn_.*\\.csv$",
      "^k_curve_lr_.*\\.csv$"
    ),
    collapse = "|"
  )
  files <- list.files(source_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop(sprintf("No Fig. 3 raw CSV files found under %s", source_dir), call. = FALSE)
  }
  invisible(file.copy(files, target_dir, overwrite = TRUE))
  normalizePath(target_dir, winslash = "/", mustWork = TRUE)
}

analysis_raw_dir <- if (normalizePath(source_raw_dir, winslash = "/", mustWork = TRUE) == normalizePath(fig_paths$raw, winslash = "/", mustWork = FALSE)) {
  fig_paths$raw
} else {
  copy_raw_inputs(source_raw_dir, fig_paths$raw)
}

read_kcurves <- function(raw_dir) {
  purrr::map2_dfr(cell_specs$output_stub, cell_specs$cell_type, function(output_stub, cell_type) {
    path <- file.path(raw_dir, paste0("k_curve_knn_", output_stub, ".csv"))
    safe_read_csv(path) %>%
      rename(top_k = any_of(c("k", "top_k"))) %>%
      arrange(top_k) %>%
      distinct(top_k, .keep_all = TRUE) %>%
      transmute(
        top_k = as.numeric(top_k),
        AUROC = as.numeric(AUROC),
        AUPRC = as.numeric(AUPRC),
        F1 = as.numeric(F1),
        ACC = as.numeric(ACC),
        cell_type = cell_type
      )
  }) %>%
    mutate(cell_type = factor(cell_type, levels = cell_type_levels))
}

read_histories <- function(raw_dir) {
  purrr::map2_dfr(cell_specs$output_stub, cell_specs$cell_type, function(output_stub, cell_type) {
    path <- file.path(raw_dir, paste0("training_history_", output_stub, ".csv"))
    safe_read_csv(path) %>% mutate(cell_type = cell_type)
  }) %>%
    group_by(cell_type, epoch) %>%
    summarise(
      train_loss = mean(train_loss, na.rm = TRUE),
      val_loss = mean(val_loss, na.rm = TRUE),
      val_acc = mean(val_acc, na.rm = TRUE),
      val_auroc = mean(val_auroc, na.rm = TRUE),
      val_auprc = mean(val_auprc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(cell_type = factor(cell_type, levels = cell_type_levels))
}

read_rankings <- function(raw_dir) {
  purrr::map2_dfr(cell_specs$output_stub, cell_specs$cell_type, function(output_stub, cell_type) {
    path <- file.path(raw_dir, paste0("ranked_genes_", output_stub, ".csv"))
    safe_read_csv(path) %>%
      transmute(
        gene = as.character(gene),
        score = as.numeric(score),
        rank = as.integer(rank),
        cell_type = cell_type
      )
  }) %>%
    mutate(cell_type = factor(cell_type, levels = cell_type_levels))
}

df_metrics <- read_kcurves(analysis_raw_dir)
hist_sum <- read_histories(analysis_raw_dir)
rank_all <- read_rankings(analysis_raw_dir)

best_k_tbl <- df_metrics %>%
  group_by(cell_type) %>%
  mutate(
    AUROC_n = rescale_safe(AUROC),
    AUPRC_n = rescale_safe(AUPRC),
    F1_n = rescale_safe(F1),
    ACC_n = rescale_safe(ACC),
    composite = rowMeans(across(c(AUROC_n, AUPRC_n, F1_n, ACC_n)), na.rm = TRUE)
  ) %>%
  arrange(desc(composite), top_k) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    cell_type,
    best_k = top_k,
    composite,
    AUROC,
    AUPRC,
    F1,
    ACC
  )

best_metric_tbl <- hist_sum %>%
  group_by(cell_type) %>%
  summarise(
    best_val_acc = max(val_acc, na.rm = TRUE),
    best_val_auroc = max(val_auroc, na.rm = TRUE),
    best_val_auprc = max(val_auprc, na.rm = TRUE),
    .groups = "drop"
  )

rank_best <- rank_all %>%
  inner_join(best_k_tbl %>% select(cell_type, best_k), by = "cell_type") %>%
  filter(rank <= best_k)

gene_lists_tbl <- rank_best %>%
  arrange(cell_type, rank) %>%
  group_by(cell_type) %>%
  summarise(
    gene_list = list(gene),
    genes_collapsed = paste(gene, collapse = "; "),
    .groups = "drop"
  )

gene_lists <- setNames(gene_lists_tbl$gene_list, gene_lists_tbl$cell_type)

paper3_write_csv(best_k_tbl, file.path(fig_paths$tables, "best_k_by_celltype.csv"))
paper3_write_csv(best_metric_tbl, file.path(fig_paths$tables, "best_validation_metrics_by_celltype.csv"))
paper3_write_csv(rank_best, file.path(fig_paths$tables, "topk_gene_lists_long.csv"))
paper3_write_csv(gene_lists_tbl %>% select(cell_type, genes_collapsed), file.path(fig_paths$tables, "topk_gene_lists_by_celltype.csv"))
saveRDS(gene_lists, file = file.path(fig_paths$rds, "topk_gene_lists_by_celltype_new.rds"))
paper3_write_csv(best_k_tbl, file.path(fig_paths$root, "best_k_by_celltype.csv"))
saveRDS(gene_lists, file = file.path(fig_paths$root, "topk_gene_lists_by_celltype_new.rds"))

plateau_df <- df_metrics %>%
  group_by(cell_type) %>%
  summarise(
    max_auroc = max(AUROC, na.rm = TRUE),
    threshold = 0.95 * max_auroc,
    xmin = min(top_k[AUROC >= threshold], na.rm = TRUE),
    xmax = max(top_k[AUROC >= threshold], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(best_k_tbl %>% select(cell_type, best_k), by = "cell_type")

if (has_ggh4x) {
  strip_obj <- ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(
      fill = unname(cell_cols[cell_type_levels]),
      colour = "black",
      size = 0.5
    ),
    text_x = ggh4x::elem_list_text(colour = "black", face = "bold")
  )

  pA <- ggplot(df_metrics, aes(top_k, AUROC, colour = cell_type, group = cell_type)) +
    geom_rect(
      data = plateau_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey92",
      colour = NA
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.4) +
    geom_vline(
      data = best_k_tbl,
      aes(xintercept = best_k),
      inherit.aes = FALSE,
      linetype = 2,
      linewidth = 0.7,
      colour = "black"
    ) +
    scale_colour_manual(values = cell_cols) +
    ggh4x::facet_wrap2(~ cell_type, ncol = 4, strip = strip_obj) +
    labs(
      x = "Top-k genes",
      y = "AUROC",
      title = "Fig. 3A: K-sensitivity of GAN-derived signatures",
      subtitle = "Grey band marks AUROC >= 95% of the within-cell-type maximum"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
} else {
  pA <- ggplot(df_metrics, aes(top_k, AUROC, colour = cell_type, group = cell_type)) +
    geom_rect(
      data = plateau_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey92",
      colour = NA
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.4) +
    geom_vline(
      data = best_k_tbl,
      aes(xintercept = best_k),
      inherit.aes = FALSE,
      linetype = 2,
      linewidth = 0.7,
      colour = "black"
    ) +
    scale_colour_manual(values = cell_cols) +
    facet_wrap(~ cell_type, ncol = 4) +
    labs(
      x = "Top-k genes",
      y = "AUROC",
      title = "Fig. 3A: K-sensitivity of GAN-derived signatures",
      subtitle = "Grey band marks AUROC >= 95% of the within-cell-type maximum"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

pB <- best_k_tbl %>%
  mutate(
    cell_type = factor(cell_type, levels = rev(cell_type_levels)),
    label = sprintf("k=%d; comp=%.2f", best_k, composite)
  ) %>%
  ggplot(aes(cell_type, best_k, fill = cell_type)) +
  geom_col(width = 0.75, colour = "grey20", linewidth = 0.4) +
  geom_text(aes(label = label), hjust = -0.05, size = 3.1) +
  scale_fill_manual(values = cell_cols, guide = "none") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Selected best_k",
    title = "Fig. 3B: Selected best_k per cell type",
    subtitle = "Composite score uses normalized AUROC, AUPRC, F1, and ACC"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

loss_long <- hist_sum %>%
  select(cell_type, epoch, train_loss, val_loss) %>%
  pivot_longer(
    cols = c(train_loss, val_loss),
    names_to = "split",
    values_to = "loss"
  ) %>%
  mutate(
    split = recode(split, train_loss = "Train loss", val_loss = "Validation loss")
  )

pC <- ggplot(loss_long, aes(epoch, loss, colour = cell_type)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ split, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = cell_cols, name = "Cell type") +
  labs(
    x = "Epoch",
    y = "Loss",
    title = "Fig. 3C: Training loss summary"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

metric_long <- hist_sum %>%
  select(cell_type, epoch, val_acc, val_auroc, val_auprc) %>%
  pivot_longer(
    cols = c(val_acc, val_auroc, val_auprc),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      val_acc = "Validation ACC",
      val_auroc = "Validation AUROC",
      val_auprc = "Validation AUPRC"
    )
  )

pD <- ggplot(metric_long, aes(epoch, value, colour = cell_type)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = cell_cols, name = "Cell type") +
  labs(
    x = "Epoch",
    y = "Metric value",
    title = "Fig. 3D: Validation metric summary"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

label_best <- rank_best %>%
  group_by(cell_type) %>%
  slice_max(score, n = 5, with_ties = FALSE) %>%
  ungroup()

pE <- ggplot(rank_best, aes(rank, score, colour = cell_type, group = cell_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = cell_cols) +
  facet_wrap(~ cell_type, ncol = 4, scales = "free_y") +
  labs(
    x = "Rank within selected top-k genes",
    y = "GAN score",
    title = "Fig. 3E: Ranked genes retained at selected best_k"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

if (has_ggrepel) {
  pE <- pE +
    ggrepel::geom_text_repel(
      data = label_best,
      aes(label = gene),
      size = 3,
      show.legend = FALSE,
      max.overlaps = Inf,
      min.segment.length = 0
    )
}

heatmap_tbl <- best_k_tbl %>%
  select(cell_type, AUROC, AUPRC, F1, ACC) %>%
  pivot_longer(
    cols = c(AUROC, AUPRC, F1, ACC),
    names_to = "metric",
    values_to = "value"
  )

pF <- ggplot(heatmap_tbl, aes(metric, factor(cell_type, levels = rev(cell_type_levels)), fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3.2) +
  scale_fill_gradient(low = "#f7fbff", high = "#084594") +
  labs(
    x = NULL,
    y = NULL,
    title = "Fig. 3F: Holdout metrics at selected best_k"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

benchmark_script <- file.path(getwd(), "scripts", "figures", "08_fig03_revision_benchmark_panels.R")
if (!file.exists(benchmark_script)) {
  stop(sprintf("Required benchmark panel script is missing: %s", benchmark_script), call. = FALSE)
}

benchmark_env <- new.env(parent = globalenv())
sys.source(benchmark_script, envir = benchmark_env)

pG <- benchmark_env$pG
pH <- benchmark_env$pH

if (!inherits(pG, "ggplot") || !inherits(pH, "ggplot")) {
  stop("Failed to assemble Fig. 3G/H panels from scripts/figures/08_fig03_revision_benchmark_panels.R", call. = FALSE)
}

save_plot2(pA, "Fig3_A_topk_sensitivity", width = 12, height = 7)
save_plot2(pB, "Fig3_B_bestk_summary", width = 10, height = 5)
save_plot2(pC, "Fig3_C_loss_summary", width = 11, height = 5)
save_plot2(pD, "Fig3_D_metric_summary", width = 11, height = 5)
save_plot2(pE, "Fig3_E_ranked_genes", width = 12, height = 8)
save_plot2(pF, "Fig3_F_bestval_heatmap", width = 7, height = 5)
save_plot2(pG, "Fig3_G_revision_patient_pair_AUROC", width = 13, height = 7.4)
save_plot2(pH, "Fig3_H_revision_patient_pair_AUPRC", width = 13, height = 7.4)

panel_inventory <- tibble::tribble(
  ~panel, ~description,
  "A", "K-sensitivity of GAN-derived signatures",
  "B", "Selected best_k per cell type",
  "C", "Training loss summary",
  "D", "Validation metric summary",
  "E", "Ranked genes retained at selected best_k",
  "F", "Holdout metrics at selected best_k",
  "G", "Patient-pair equal-budget benchmark by AUROC",
  "H", "Patient-pair equal-budget benchmark by AUPRC"
)
paper3_write_csv(panel_inventory, file.path(fig_paths$tables, "panel_inventory.csv"))
paper3_write_lines(
  c(
    "Fig. 3 workspace one-run summary",
    sprintf("Raw GAN source directory: %s", source_raw_dir),
    sprintf("Benchmark panel script: %s", benchmark_script),
    "Panels A-F come from the GAN signature discovery summary.",
    "Panels G-H are sourced from the revision patient-pair benchmark workflow.",
    "The public release keeps panel-level exports only and omits the combined A-H layout."
  ),
  file.path(fig_paths$logs, "provenance_notes.txt")
)

if (has_enrich) {
  msig_all <- bind_rows(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
  ) %>%
    select(gs_name, gene_symbol) %>%
    distinct()

  enrich_tbl <- purrr::imap_dfr(gene_lists, function(genes, cell_type) {
    enrich_obj <- tryCatch(
      clusterProfiler::enricher(
        genes,
        TERM2GENE = msig_all,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
      ),
      error = function(e) NULL
    )

    if (is.null(enrich_obj)) {
      return(tibble())
    }

    enrich_df <- as_tibble(as.data.frame(enrich_obj))
    if (nrow(enrich_df) == 0) {
      return(tibble())
    }

    enrich_df %>% mutate(cell_type = cell_type)
  })

  paper3_write_csv(enrich_tbl, file.path(fig_paths$tables, "Fig3_ORA_enrichment_table.csv"))

  if (nrow(enrich_tbl) > 0) {
    pathway_long <- enrich_tbl %>%
      filter(qvalue < 0.05) %>%
      mutate(score = pmin(-log10(qvalue), 10)) %>%
      select(Description, cell_type, score)

    pathway_long_sel <- pathway_long %>%
      group_by(Description) %>%
      filter(n_distinct(cell_type) >= 2) %>%
      ungroup()

    if (nrow(pathway_long_sel) == 0) {
      pathway_long_sel <- pathway_long %>%
        group_by(Description) %>%
        summarise(score = max(score, na.rm = TRUE), .groups = "drop") %>%
        slice_max(score, n = 20, with_ties = FALSE) %>%
        inner_join(pathway_long, by = "Description")
    }

    p_ora <- pathway_long_sel %>%
      mutate(cell_type = factor(cell_type, levels = cell_type_levels)) %>%
      ggplot(aes(cell_type, reorder(Description, score), fill = score)) +
      geom_tile(colour = "white") +
      scale_fill_gradient(low = "#fff5f0", high = "#99000d") +
      labs(
        x = NULL,
        y = NULL,
        title = "Fig. 3 supplementary: ORA heatmap of selected top-k gene lists",
        fill = "-log10(qvalue)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"),
        panel.grid = element_blank()
      )

    save_plot2(p_ora, "Fig3_ORA_heatmap", width = 11, height = 9)
  }
}

paper3_capture_session_info(file.path(fig_paths$logs, "session_info.txt"))
writeLines(
  c(
    sprintf("source_raw_dir=%s", source_raw_dir),
    sprintf("analysis_raw_dir=%s", analysis_raw_dir),
    sprintf("generated_at=%s", format(Sys.time(), tz = "UTC", usetz = TRUE))
  ),
  con = file.path(fig_paths$logs, "run_metadata.txt")
)

message("DONE: Fig. 3 panel exports regenerated in workspace outputs/figures/fig03_gan_signature_discovery/")
