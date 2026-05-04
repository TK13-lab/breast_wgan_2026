source(file.path(getwd(), "scripts", "00_project_paths.R"))

paths <- paper3_paths()
paper3_mkdir(paths$outputs, paths$logs, paths$cache, paths$output_figures, paths$output_checkpoints)
fig_paths <- paper3_init_fig("fig03_gan_signature_discovery", paths = paths)

paper3_require_packages(c("readr", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "forcats"))
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(forcats)
})

has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)

cell_specs <- tibble::tribble(
  ~cell_type_key, ~cell_type,
  "Endothelial", "Endothelial",
  "Fibroblast", "Fibroblast",
  "Pericytes_SMC", "Pericytes_SMC",
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

method_levels <- c(
  "GAN fixed top-k + LR",
  "scVI latent + LR",
  "PCA latent + LR",
  "All candidate genes + LR",
  "DEG fixed top-k + LR"
)

method_labels_plot <- c(
  "GAN fixed top-k + LR" = "GAN top-k",
  "scVI latent + LR" = "scVI latent",
  "PCA latent + LR" = "PCA latent",
  "All candidate genes + LR" = "All genes",
  "DEG fixed top-k + LR" = "DEG top-k"
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

benchmark_dir <- file.path(paths$outputs, "revision", "patient_pair_equal_budget", "tables")
latent_dir <- file.path(paths$outputs, "revision", "latent_patient_pair", "tables")

benchmark_by_cell <- safe_read_csv(file.path(benchmark_dir, "patient_pair_equal_budget_by_celltype.csv")) %>%
  filter(method_label %in% c(
    "GAN fixed top-k + LR",
    "PCA latent + LR",
    "All candidate genes + LR",
    "DEG fixed top-k + LR"
  )) %>%
  select(cell_type, method_label, mean_AUROC, mean_AUPRC, mean_F1, mean_ACC, n_folds)

scvi_by_cell <- safe_read_csv(file.path(latent_dir, "scvi_latent_patient_pair_by_celltype.csv")) %>%
  transmute(
    cell_type,
    method_label = "scVI latent + LR",
    mean_AUROC,
    mean_AUPRC,
    mean_F1,
    mean_ACC,
    n_folds
  )

plot_df <- bind_rows(benchmark_by_cell, scvi_by_cell) %>%
  mutate(
    cell_type = factor(cell_type, levels = cell_type_levels),
    method_label = factor(method_label, levels = method_levels)
  ) %>%
  arrange(cell_type, method_label)

paper3_write_csv(plot_df, file.path(fig_paths$tables, "Fig3_revision_patient_pair_benchmark_long.csv"))

auroc_winners <- plot_df %>%
  group_by(cell_type) %>%
  slice_max(mean_AUROC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  count(method_label, name = "auroc_win_count")

auprc_winners <- plot_df %>%
  group_by(cell_type) %>%
  slice_max(mean_AUPRC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  count(method_label, name = "auprc_win_count")

method_summary <- plot_df %>%
  group_by(method_label) %>%
  summarise(
    mean_AUROC = mean(mean_AUROC, na.rm = TRUE),
    mean_AUPRC = mean(mean_AUPRC, na.rm = TRUE),
    mean_F1 = mean(mean_F1, na.rm = TRUE),
    mean_ACC = mean(mean_ACC, na.rm = TRUE),
    n_cell_types = dplyr::n(),
    .groups = "drop"
  ) %>%
  left_join(auroc_winners, by = "method_label") %>%
  left_join(auprc_winners, by = "method_label") %>%
  mutate(
    auroc_win_count = tidyr::replace_na(auroc_win_count, 0L),
    auprc_win_count = tidyr::replace_na(auprc_win_count, 0L),
    method_label = factor(method_label, levels = method_levels)
  ) %>%
  arrange(method_label)

paper3_write_csv(method_summary, file.path(fig_paths$tables, "Fig3_revision_patient_pair_benchmark_summary.csv"))

main_table <- plot_df %>%
  group_by(cell_type) %>%
  summarise(
    gan_AUROC = mean_AUROC[method_label == "GAN fixed top-k + LR"],
    gan_AUPRC = mean_AUPRC[method_label == "GAN fixed top-k + LR"],
    scvi_AUROC = mean_AUROC[method_label == "scVI latent + LR"],
    scvi_AUPRC = mean_AUPRC[method_label == "scVI latent + LR"],
    best_non_gan_AUROC = max(mean_AUROC[method_label != "GAN fixed top-k + LR"], na.rm = TRUE),
    best_non_gan_AUPRC = max(mean_AUPRC[method_label != "GAN fixed top-k + LR"], na.rm = TRUE),
    best_non_gan_auroc_method = as.character(method_label[which.max(ifelse(method_label != "GAN fixed top-k + LR", mean_AUROC, -Inf))]),
    best_non_gan_auprc_method = as.character(method_label[which.max(ifelse(method_label != "GAN fixed top-k + LR", mean_AUPRC, -Inf))]),
    n_folds = first(n_folds),
    .groups = "drop"
  ) %>%
  mutate(
    delta_AUROC_vs_scvi = gan_AUROC - scvi_AUROC,
    delta_AUPRC_vs_scvi = gan_AUPRC - scvi_AUPRC,
    delta_AUROC_vs_best_non_gan = gan_AUROC - best_non_gan_AUROC,
    delta_AUPRC_vs_best_non_gan = gan_AUPRC - best_non_gan_AUPRC,
    gan_is_top_auroc = gan_AUROC >= best_non_gan_AUROC,
    gan_is_top_auprc = gan_AUPRC >= best_non_gan_AUPRC,
    cell_type = factor(cell_type, levels = cell_type_levels)
  ) %>%
  arrange(cell_type)

paper3_write_csv(main_table, file.path(fig_paths$tables, "Fig3_revision_patient_pair_benchmark_main.csv"))

round_metric_cols <- function(df) {
  df %>%
    mutate(
      across(
        .cols = where(is.numeric) & !matches("^n_"),
        .fns = ~ round(.x, digits = 3)
      )
    )
}

paper3_write_csv(
  round_metric_cols(main_table),
  file.path(fig_paths$tables, "Fig3_revision_patient_pair_benchmark_main_display.csv")
)
paper3_write_csv(
  round_metric_cols(method_summary),
  file.path(fig_paths$tables, "Fig3_revision_patient_pair_benchmark_summary_display.csv")
)

plot_table <- plot_df %>%
  mutate(
    method_plot = recode(as.character(method_label), !!!method_labels_plot),
    method_plot = factor(method_plot, levels = rev(unname(method_labels_plot))),
    is_gan = method_label == "GAN fixed top-k + LR",
    facet_fill = cell_cols[as.character(cell_type)]
  )

strip_obj <- NULL
if (has_ggh4x) {
  strip_obj <- ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(
      fill = unname(cell_cols[cell_type_levels]),
      colour = "black",
      size = 0.45
    ),
    text_x = ggh4x::elem_list_text(colour = "black", face = "bold")
  )
}

build_metric_panel <- function(metric_col, metric_label, title_text, filename) {
  panel_df <- plot_table %>%
    mutate(
      metric_value = .data[[metric_col]],
      point_fill = ifelse(is_gan, facet_fill, "#FFFFFF"),
      point_stroke = ifelse(is_gan, 1.0, 0.65),
      point_alpha = ifelse(is_gan, 1.0, 0.9)
    )

  axis_floor <- floor(min(panel_df$metric_value, na.rm = TRUE) * 10) / 10
  axis_floor <- max(0.2, axis_floor - 0.05)
  axis_ceiling <- min(0.95, ceiling(max(panel_df$metric_value, na.rm = TRUE) * 10) / 10 + 0.02)

  p <- ggplot(panel_df, aes(x = metric_value, y = method_plot)) +
    geom_segment(
      aes(x = axis_floor, xend = metric_value, yend = method_plot),
      linewidth = 0.55,
      colour = "grey80"
    ) +
    geom_vline(
      xintercept = 0.5,
      linewidth = 0.35,
      linetype = "dashed",
      colour = "grey60"
    ) +
    geom_point(
      aes(fill = point_fill, alpha = point_alpha, stroke = point_stroke),
      shape = 21,
      size = 2.8,
      colour = "black"
    ) +
    scale_fill_identity() +
    scale_alpha_identity() +
    scale_x_continuous(
      limits = c(axis_floor, axis_ceiling),
      breaks = seq(0.2, 1.0, by = 0.1),
      expand = c(0, 0)
    ) +
    labs(
      x = metric_label,
      y = NULL,
      title = title_text,
      subtitle = "Patient-pair holdout with equal signature budget across methods"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      panel.spacing = unit(0.85, "lines"),
      strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.45),
      strip.text = element_text(face = "bold", size = 10)
    )

  if (has_ggh4x) {
    p <- p + ggh4x::facet_wrap2(~ cell_type, ncol = 4, strip = strip_obj)
  } else {
    p <- p + facet_wrap(~ cell_type, ncol = 4)
  }

  save_plot2(p, filename, width = 13, height = 7.4)
  p
}

pG <- build_metric_panel(
  metric_col = "mean_AUROC",
  metric_label = "Mean AUROC across held-out patient pairs",
  title_text = "Fig. 3G. Patient-pair benchmark by AUROC",
  filename = "Fig3_G_revision_patient_pair_AUROC"
)

pH <- build_metric_panel(
  metric_col = "mean_AUPRC",
  metric_label = "Mean AUPRC across held-out patient pairs",
  title_text = "Fig. 3H. Patient-pair benchmark by AUPRC",
  filename = "Fig3_H_revision_patient_pair_AUPRC"
)

pSummary <- method_summary %>%
  mutate(
    method_plot = recode(as.character(method_label), !!!method_labels_plot),
    method_plot = factor(method_plot, levels = rev(unname(method_labels_plot)))
  ) %>%
  ggplot(aes(x = method_plot)) +
  geom_col(aes(y = mean_AUROC), fill = "grey85", colour = "grey20", linewidth = 0.35) +
  geom_point(aes(y = mean_AUPRC), colour = "#08519c", size = 2.5) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Overall mean metric",
    title = "Revision benchmark summary",
    subtitle = "Bars show mean AUROC and blue points show mean AUPRC"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 9)
  )


 
