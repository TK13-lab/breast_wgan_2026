source(file.path(getwd(), "scripts", "00_project_paths.R"))

paths <- paper3_init()
fig_paths <- paper3_init_fig("fig05_bulk_projection_survival", paths = paths)

paper3_require_packages(c(
  "dplyr", "tidyr", "tibble", "ggplot2", "scales", "survival",
  "survminer", "patchwork", "ggpubr", "purrr"
))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(survival)
  library(survminer)
  library(patchwork)
  library(ggpubr)
  library(purrr)
})

save_plot_bundle <- function(plot_obj, stem, width, height, dpi = 300) {
  pdf_path <- file.path(fig_paths$plots, paste0(stem, ".pdf"))
  png_path <- file.path(fig_paths$plots, paste0(stem, ".png"))
  tiff_path <- file.path(fig_paths$plots, paste0(stem, ".tiff"))
  ggsave(pdf_path, plot = plot_obj, width = width, height = height, useDingbats = FALSE)
  ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = dpi)
  ggsave(tiff_path, plot = plot_obj, width = width, height = height, dpi = 600, compression = "lzw")
}

time_col <- "OS.time"
event_col <- "OS"
lineages <- c("Epithelial", "Immune", "Stromal", "Vascular")
subtypes <- c("Basal", "Her2", "LumA", "LumB")

expr_path <- paper3_legacy_file("expr_merged_sub_TCGA_tumor.rds", paths = paths)
bulk_path <- paper3_legacy_file("bulk_all.with_PAM50.rds", paths = paths)
genes_by_lineage_path <- paper3_fig04_lineage_gene_path(paths = paths, required = TRUE)

expr_all <- as.matrix(readRDS(expr_path))
bulk_all_sub <- readRDS(bulk_path)
genes_by_lineage <- readRDS(genes_by_lineage_path)

stopifnot(all(lineages %in% names(genes_by_lineage)))

to_event01 <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }
  if (is.numeric(x)) {
    return(as.integer(x != 0))
  }
  x0 <- tolower(trimws(as.character(x)))
  event_tokens <- c("1", "true", "dead", "deceased", "death", "event", "yes", "recurrence", "progression")
  as.integer(x0 %in% event_tokens)
}

winsorize_sd <- function(z, k = 3) {
  pmax(pmin(z, k), -k)
}

meta_A <- bulk_all_sub %>%
  mutate(
    group = case_when(
      sample_type_simple == "Normal" ~ "Normal",
      sample_type_simple == "Tumor" ~ as.character(PAM50_4),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(group %in% c("Normal", subtypes)) %>%
  filter(sample_id %in% colnames(expr_all))

ids_A <- meta_A$sample_id
expr_A <- expr_all[, ids_A, drop = FALSE]
expr_A_z <- t(scale(t(expr_A)))
expr_A_z[!is.finite(expr_A_z)] <- NA_real_

score_A <- lapply(lineages, function(lin) {
  genes <- intersect(genes_by_lineage[[lin]], rownames(expr_A_z))
  if (length(genes) < 2) {
    return(NULL)
  }
  tibble(
    sample_id = ids_A,
    lineage = lin,
    score = colMeans(expr_A_z[genes, , drop = FALSE], na.rm = TRUE)
  )
}) %>%
  bind_rows() %>%
  left_join(meta_A %>% select(sample_id, group), by = "sample_id") %>%
  mutate(
    group = factor(group, levels = c("Normal", subtypes)),
    lineage = factor(lineage, levels = lineages)
  )

paper3_write_csv(score_A, file.path(fig_paths$tables, "lineage_scores_across_groups.csv"))

comps <- list(
  c("Normal", "Basal"),
  c("Normal", "Her2"),
  c("Normal", "LumA"),
  c("Normal", "LumB")
)

pA <- ggplot(score_A, aes(x = group, y = score, fill = group)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(color = group), width = 0.15, size = 0.7, alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~ lineage, nrow = 1, scales = "free_y") +
  ggpubr::stat_compare_means(
    comparisons = comps,
    method = "wilcox.test",
    label = "p.format",
    hide.ns = FALSE,
    step.increase = 0.08,
    size = 2.8
  ) +
  scale_fill_manual(values = c(
    "Normal" = "#1f78b4",
    "Basal" = "#e31a1c",
    "Her2" = "#33a02c",
    "LumA" = "#ff7f00",
    "LumB" = "#6a3d9a"
  )) +
  scale_color_manual(values = c(
    "Normal" = "#1f78b4",
    "Basal" = "#e31a1c",
    "Her2" = "#33a02c",
    "LumA" = "#ff7f00",
    "LumB" = "#6a3d9a"
  )) +
  labs(x = NULL, y = "Mean gene z-score") +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 10),
    legend.position = "right",
    plot.margin = margin(6, 6, 6, 6)
  )

tumor_ids <- bulk_all_sub %>%
  filter(sample_type_simple == "Tumor") %>%
  pull(sample_id) %>%
  intersect(colnames(expr_all))

expr_tumor <- expr_all[, tumor_ids, drop = FALSE]
bulk_surv <- bulk_all_sub %>%
  filter(
    sample_id %in% tumor_ids,
    PAM50_4 %in% subtypes,
    !is.na(.data[[time_col]]),
    !is.na(.data[[event_col]])
  ) %>%
  transmute(
    sample_id,
    patient_id,
    PAM50_4,
    time = .data[[time_col]],
    event = .data[[event_col]]
  )

fit_cox_lineage_subtype <- function(lineage, subtype, sig_list, expr_tumor, bulk_surv, min_n = 10, p_filter = 0.1) {
  meta_sub <- bulk_surv %>%
    filter(PAM50_4 == subtype) %>%
    mutate(event01 = to_event01(event))
  if (nrow(meta_sub) < min_n) {
    return(NULL)
  }

  genes0 <- intersect(sig_list[[lineage]], rownames(expr_tumor))
  if (length(genes0) < 2) {
    return(NULL)
  }

  expr_sub0 <- expr_tumor[genes0, meta_sub$sample_id, drop = FALSE]
  sd_g <- apply(expr_sub0, 1, sd, na.rm = TRUE)
  keep <- is.finite(sd_g) & sd_g > 0
  genes0 <- genes0[keep]
  expr_sub0 <- expr_sub0[keep, , drop = FALSE]
  if (length(genes0) < 2) {
    return(NULL)
  }

  uni_tbl <- lapply(genes0, function(gene_name) {
    df_g <- data.frame(
      time = meta_sub$time,
      event = meta_sub$event01,
      gene = as.numeric(expr_sub0[gene_name, meta_sub$sample_id])
    )
    fit_g <- tryCatch(coxph(Surv(time, event) ~ gene, data = df_g), error = function(e) NULL)
    if (is.null(fit_g)) {
      return(NULL)
    }
    s <- summary(fit_g)
    tibble(gene = gene_name, pval = s$coef[1, "Pr(>|z|)"])
  }) %>%
    bind_rows()

  if (nrow(uni_tbl) < 2) {
    return(NULL)
  }

  genes <- uni_tbl %>% filter(pval < p_filter) %>% pull(gene)
  if (length(genes) < 2) {
    return(NULL)
  }

  expr_sub <- expr_tumor[genes, meta_sub$sample_id, drop = FALSE]
  safe <- make.names(genes)
  X <- t(expr_sub)
  colnames(X) <- safe

  df_multi <- data.frame(time = meta_sub$time, event = meta_sub$event01, X, check.names = FALSE)
  form <- as.formula(paste0("Surv(time, event) ~ ", paste(safe, collapse = " + ")))

  fit <- tryCatch(coxph(form, data = df_multi), error = function(e) NULL)
  if (is.null(fit)) {
    return(NULL)
  }

  coefs <- coef(fit)
  safe_used <- names(coefs)
  orig_genes <- genes[match(safe_used, safe)]
  tibble(lineage = lineage, subtype = subtype, gene = orig_genes, coef = as.numeric(coefs))
}

coef_path <- file.path(fig_paths$rds, "lineage_subtype_coefs.coxph.rds")
legacy_coef_path <- paper3_fig05_coef_path(paths = paths, required = FALSE)
if (file.exists(legacy_coef_path)) {
  paper3_stage_reference_file(
    source_path = legacy_coef_path,
    destination_path = coef_path,
    label = "legacy lineage_subtype_coefs.coxph.rds"
  )
  coef_tbl <- readRDS(coef_path)
  coef_source <- legacy_coef_path
} else if (file.exists(coef_path)) {
  coef_tbl <- readRDS(coef_path)
  coef_source <- coef_path
} else {
  coef_tbl <- purrr::map_dfr(lineages, function(lin) {
    purrr::map_dfr(subtypes, function(sub) {
      fit_cox_lineage_subtype(
        lineage = lin,
        subtype = sub,
        sig_list = genes_by_lineage,
        expr_tumor = expr_tumor,
        bulk_surv = bulk_surv,
        min_n = 10,
        p_filter = 0.1
      )
    })
  })
  saveRDS(coef_tbl, coef_path)
  coef_source <- "computed_in_workspace"
}

paper3_write_csv(coef_tbl, file.path(fig_paths$tables, "lineage_subtype_coefs.coxph.csv"))

genes_by_lineage_sig_path <- file.path(fig_paths$rds, "genes_by_lineage_sig.rds")
legacy_genes_by_lineage_sig_path <- paper3_fig05_lineage_sig_path(paths = paths, required = FALSE)
if (file.exists(legacy_genes_by_lineage_sig_path)) {
  paper3_stage_reference_file(
    source_path = legacy_genes_by_lineage_sig_path,
    destination_path = genes_by_lineage_sig_path,
    label = "legacy genes_by_lineage_sig.rds"
  )
  genes_by_lineage_sig <- readRDS(genes_by_lineage_sig_path)
  genes_by_lineage_sig_source <- legacy_genes_by_lineage_sig_path
} else {
  genes_by_lineage_sig <- coef_tbl %>%
    group_by(lineage) %>%
    summarise(
      genes = list(sort(unique(gene))),
      n_genes = dplyr::n_distinct(gene),
      .groups = "drop"
    )
  saveRDS(genes_by_lineage_sig, genes_by_lineage_sig_path)
  genes_by_lineage_sig_source <- "derived_from_coef_tbl"
}
paper3_write_csv(
  genes_by_lineage_sig %>% mutate(genes_collapsed = vapply(genes, paste, collapse = "; ", character(1))),
  file.path(fig_paths$tables, "genes_by_lineage_sig.csv")
)

genes_by_lineage_subtype_sig_path <- file.path(fig_paths$rds, "genes_by_lineage_subtype_sig.rds")
legacy_genes_by_lineage_subtype_sig_path <- paper3_fig05_subtype_gene_path(paths = paths, required = FALSE)
if (file.exists(legacy_genes_by_lineage_subtype_sig_path)) {
  paper3_stage_reference_file(
    source_path = legacy_genes_by_lineage_subtype_sig_path,
    destination_path = genes_by_lineage_subtype_sig_path,
    label = "legacy genes_by_lineage_subtype_sig.rds"
  )
  genes_by_lineage_subtype_sig <- readRDS(genes_by_lineage_subtype_sig_path)
  genes_by_lineage_subtype_sig_source <- legacy_genes_by_lineage_subtype_sig_path
} else {
  genes_by_lineage_subtype_sig <- coef_tbl %>%
    group_by(lineage, subtype) %>%
    summarise(
      genes = list(sort(unique(gene))),
      n_genes = dplyr::n_distinct(gene),
      .groups = "drop"
    )
  saveRDS(genes_by_lineage_subtype_sig, genes_by_lineage_subtype_sig_path)
  genes_by_lineage_subtype_sig_source <- "derived_from_coef_tbl"
}
paper3_write_csv(
  genes_by_lineage_subtype_sig %>% mutate(genes_collapsed = vapply(genes, paste, collapse = "; ", character(1))),
  file.path(fig_paths$tables, "genes_by_lineage_subtype_sig.csv")
)

risk_long <- purrr::map_dfr(subtypes, function(sub) {
  meta_sub <- bulk_surv %>% filter(PAM50_4 == sub)
  if (nrow(meta_sub) == 0) {
    return(NULL)
  }
  purrr::map_dfr(lineages, function(lin) {
    coef_sub <- coef_tbl %>% filter(lineage == lin, subtype == sub)
    if (nrow(coef_sub) == 0) {
      return(NULL)
    }
    genes <- intersect(coef_sub$gene, rownames(expr_tumor))
    if (length(genes) < 1) {
      return(NULL)
    }
    beta <- coef_sub$coef[match(genes, coef_sub$gene)]
    expr_sub <- expr_tumor[genes, meta_sub$sample_id, drop = FALSE]
    scores <- as.numeric(t(expr_sub) %*% beta)
    tibble(
      sample_id = meta_sub$sample_id,
      PAM50_4 = sub,
      lineage = lin,
      score = scores
    )
  })
})
paper3_write_csv(risk_long, file.path(fig_paths$tables, "risk_scores_long.csv"))

risk_wide <- risk_long %>%
  select(sample_id, lineage, score) %>%
  pivot_wider(names_from = lineage, values_from = score)

bulk_surv2 <- bulk_surv %>%
  left_join(risk_wide, by = "sample_id")
saveRDS(bulk_surv2, file.path(fig_paths$rds, "bulk_surv2.rds"))
paper3_write_csv(bulk_surv2, file.path(fig_paths$tables, "bulk_surv2.csv"))

get_surv_plot_median <- function(df, lineage, subtype, tick_size = 11, axis_title_size = 11, title_size = 9) {
  df_sub <- df %>%
    filter(PAM50_4 == subtype) %>%
    mutate(score = .data[[lineage]], event01 = to_event01(event)) %>%
    filter(!is.na(score), !is.na(time), !is.na(event))

  if (nrow(df_sub) < 10) {
    return(NULL)
  }

  med <- median(df_sub$score, na.rm = TRUE)
  df_sub <- df_sub %>%
    mutate(group = ifelse(score >= med, "High risk", "Low risk"), group = factor(group, levels = c("Low risk", "High risk")))

  if (min(table(df_sub$group)) < 10) {
    return(NULL)
  }

  fit <- survival::survfit(survival::Surv(time, event01) ~ group, data = df_sub)
  survminer::ggsurvplot(
    fit,
    data = df_sub,
    risk.table = FALSE,
    pval = TRUE,
    conf.int = FALSE,
    legend.position = "bottom",
    legend.title = NULL,
    legend.labs = c("Low risk", "High risk"),
    palette = c("#1f78b4", "#e31a1c"),
    font.x = c(axis_title_size),
    font.y = c(axis_title_size),
    font.tickslab = c(tick_size),
    font.legend = c(tick_size),
    ggtheme = ggplot2::theme_classic(base_size = 12)
  )$plot +
    ggtitle(paste0(subtype, " - ", lineage, " lineage signature")) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.margin = ggplot2::margin(6, 6, 6, 6),
      legend.title = ggplot2::element_blank()
    )
}

plots_B <- purrr::map(subtypes, function(sub) {
  purrr::map(lineages, function(lin) get_surv_plot_median(bulk_surv2, lin, sub, tick_size = 10, axis_title_size = 10))
}) %>%
  unlist(recursive = FALSE)
plots_B <- plots_B[!vapply(plots_B, is.null, logical(1))]

pB <- wrap_plots(plots_B, nrow = length(subtypes), ncol = length(lineages)) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.size = grid::unit(6, "mm"),
    legend.spacing.x = grid::unit(4, "mm")
  )

extract_cox_stats <- function(fit) {
  sm <- summary(fit)
  cm <- sm$coefficients
  if (is.null(cm) || nrow(cm) < 1) {
    return(list(beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  list(
    beta = unname(cm[1, "coef"]),
    se = unname(cm[1, "se(coef)"]),
    p = unname(cm[1, "Pr(>|z|)"])
  )
}

fit_cox_robust_one <- function(d, winsor_k = 3, unstable_beta = 10) {
  dd <- d %>%
    transmute(time = time, event = to_event01(event), score = score) %>%
    filter(!is.na(time), !is.na(event), !is.na(score))

  n <- nrow(dd)
  events <- sum(dd$event == 1, na.rm = TRUE)
  if (n < 10 || events < 5) {
    return(list(beta = NA_real_, se = NA_real_, p = NA_real_, n = n, events = events, method = "skip"))
  }

  dd$score_z <- as.numeric(scale(dd$score))
  dd$score_z <- winsorize_sd(dd$score_z, k = winsor_k)

  warn_flag <- FALSE
  fit1 <- withCallingHandlers(
    try(coxph(Surv(time, event) ~ score_z, data = dd), silent = TRUE),
    warning = function(w) {
      warn_flag <<- TRUE
      invokeRestart("muffleWarning")
    }
  )

  if (!inherits(fit1, "try-error")) {
    st <- extract_cox_stats(fit1)
    ok <- is.finite(st$beta) && is.finite(st$se) && is.finite(st$p) && !warn_flag && abs(st$beta) <= unstable_beta
    if (ok) {
      return(list(beta = st$beta, se = st$se, p = st$p, n = n, events = events, method = "coxph"))
    }
  }

  if (requireNamespace("coxphf", quietly = TRUE)) {
    fit2 <- try(coxphf::coxphf(Surv(time, event) ~ score_z, data = dd), silent = TRUE)
    if (!inherits(fit2, "try-error")) {
      beta2 <- unname(fit2$coefficients[1])
      se2 <- sqrt(unname(fit2$var[1, 1]))
      p2 <- unname(fit2$prob[1])
      return(list(beta = beta2, se = se2, p = p2, n = n, events = events, method = "coxphf"))
    }
  }

  list(beta = NA_real_, se = NA_real_, p = NA_real_, n = n, events = events, method = "fail")
}

fit_cox_grid_robust <- function(df, lineages, subtypes, min_n = 40, min_events = 10, winsor_k = 3, unstable_beta = 10) {
  purrr::map_dfr(subtypes, function(sub) {
    df_sub <- df %>% filter(PAM50_4 == sub) %>% filter(!is.na(time), !is.na(event))
    purrr::map_dfr(lineages, function(lin) {
      d <- df_sub %>% transmute(time = time, event = event, score = .data[[lin]]) %>% filter(!is.na(score))
      n <- nrow(d)
      events <- sum(to_event01(d$event) == 1, na.rm = TRUE)
      if (n < min_n || events < min_events) {
        return(tibble(subtype = sub, lineage = lin, n = n, events = events, beta = NA_real_, se = NA_real_, p = NA_real_, method = "skip"))
      }
      ft <- fit_cox_robust_one(d, winsor_k = winsor_k, unstable_beta = unstable_beta)
      tibble(subtype = sub, lineage = lin, n = ft$n, events = ft$events, beta = ft$beta, se = ft$se, p = ft$p, method = ft$method)
    })
  }) %>%
    mutate(
      logHR = beta,
      p2 = ifelse(is.na(p), NA_real_, pmax(p, 1e-300)),
      FDR = p.adjust(p2, method = "BH"),
      sig = !is.na(FDR) & FDR < 0.05
    )
}

resC <- fit_cox_grid_robust(
  df = bulk_surv2,
  lineages = lineages,
  subtypes = subtypes,
  min_n = 40,
  min_events = 10,
  winsor_k = 3,
  unstable_beta = 10
)
paper3_write_csv(resC, file.path(fig_paths$tables, "cox_grid_robust.csv"))

cap_logHR <- 2
pC_df <- resC %>%
  mutate(
    subtype = factor(subtype, levels = rev(subtypes)),
    lineage = factor(lineage, levels = lineages),
    logHR_cap = pmax(pmin(logHR, cap_logHR), -cap_logHR),
    star = ifelse(sig, "*", "")
  )

pC <- ggplot(pC_df, aes(x = lineage, y = subtype)) +
  geom_tile(aes(fill = logHR_cap), colour = "grey85", linewidth = 0.5) +
  geom_text(aes(label = star), size = 5, colour = "black") +
  scale_fill_gradient2(
    low = "#1f78b4",
    mid = "white",
    high = "#e31a1c",
    limits = c(-cap_logHR, cap_logHR),
    oob = squish,
    na.value = "grey60",
    name = "log(HR)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.margin = margin(6, 6, 6, 6)
  ) +
  labs(caption = "* FDR < 0.05")

fit_beta_robust <- function(train_df, winsor_k = 3, unstable_beta = 10) {
  dd <- train_df %>%
    transmute(time = time, event = to_event01(event), score = score) %>%
    filter(!is.na(time), !is.na(event), !is.na(score))
  if (nrow(dd) < 10 || sum(dd$event == 1) < 5) {
    return(list(beta = NA_real_, method = "skip", mu = NA_real_, sdv = NA_real_))
  }
  mu <- mean(dd$score, na.rm = TRUE)
  sdv <- sd(dd$score, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) {
    sdv <- 1
  }
  dd$score_z <- winsorize_sd((dd$score - mu) / sdv, k = winsor_k)

  warn_flag <- FALSE
  fit1 <- withCallingHandlers(
    try(coxph(Surv(time, event) ~ score_z, data = dd), silent = TRUE),
    warning = function(w) {
      warn_flag <<- TRUE
      invokeRestart("muffleWarning")
    }
  )

  if (!inherits(fit1, "try-error")) {
    beta <- unname(summary(fit1)$coef[1, "coef"])
    if (is.finite(beta) && !warn_flag && abs(beta) <= unstable_beta) {
      return(list(beta = beta, method = "coxph", mu = mu, sdv = sdv))
    }
  }

  if (requireNamespace("coxphf", quietly = TRUE)) {
    fit2 <- try(coxphf::coxphf(Surv(time, event) ~ score_z, data = dd), silent = TRUE)
    if (!inherits(fit2, "try-error")) {
      beta2 <- unname(fit2$coefficients[1])
      return(list(beta = beta2, method = "coxphf", mu = mu, sdv = sdv))
    }
  }

  list(beta = NA_real_, method = "fail", mu = mu, sdv = sdv)
}

cindex_from_lp <- function(test_df, lp) {
  dd <- test_df %>%
    transmute(time = time, event = to_event01(event), lp = lp) %>%
    filter(!is.na(time), !is.na(event), !is.na(lp))
  if (nrow(dd) < 10 || sum(dd$event == 1) < 3) {
    return(NA_real_)
  }
  sc <- survival::survConcordance(Surv(time, event) ~ lp, data = dd)
  as.numeric(sc$concordance)
}

cv_cindex_one <- function(df_sub, k = 5, repeats = 20, seed = 1, winsor_k = 3, unstable_beta = 10) {
  set.seed(seed)
  n <- nrow(df_sub)
  if (n < 20) {
    return(list(c_mean = NA_real_, c_sd = NA_real_, method = NA_character_))
  }

  cvals <- c()
  methods <- c()
  for (r in seq_len(repeats)) {
    folds <- sample(rep(seq_len(k), length.out = n))
    for (f in seq_len(k)) {
      train <- df_sub[folds != f, , drop = FALSE]
      test <- df_sub[folds == f, , drop = FALSE]
      fb <- fit_beta_robust(train, winsor_k = winsor_k, unstable_beta = unstable_beta)
      methods <- c(methods, fb$method)
      if (!is.finite(fb$beta)) {
        cvals <- c(cvals, NA_real_)
        next
      }
      score_z_test <- winsorize_sd((test$score - fb$mu) / fb$sdv, k = winsor_k)
      lp_test <- fb$beta * score_z_test
      cvals <- c(cvals, cindex_from_lp(test, lp_test))
    }
  }
  cvals <- cvals[is.finite(cvals)]
  if (length(cvals) < 5) {
    return(list(c_mean = NA_real_, c_sd = NA_real_, method = NA_character_))
  }
  list(c_mean = mean(cvals), c_sd = sd(cvals), method = paste(sort(unique(methods)), collapse = "+"))
}

calc_cindex_grid <- function(df, lineages, subtypes, min_n = 40, min_events = 10, k = 5, repeats = 20, seed = 1, winsor_k = 3, unstable_beta = 10) {
  purrr::map_dfr(subtypes, function(sub) {
    df_s <- df %>% filter(PAM50_4 == sub) %>% filter(!is.na(time), !is.na(event))
    purrr::map_dfr(lineages, function(lin) {
      d <- df_s %>% transmute(time = time, event = event, score = .data[[lin]]) %>% filter(!is.na(score))
      n <- nrow(d)
      events <- sum(to_event01(d$event) == 1, na.rm = TRUE)
      if (n < min_n || events < min_events) {
        return(tibble(subtype = sub, lineage = lin, n = n, events = events, cindex = NA_real_, method = "skip"))
      }
      cv <- cv_cindex_one(d, k = k, repeats = repeats, seed = seed, winsor_k = winsor_k, unstable_beta = unstable_beta)
      tibble(subtype = sub, lineage = lin, n = n, events = events, cindex = cv$c_mean, method = cv$method)
    })
  })
}

resD <- calc_cindex_grid(
  df = bulk_surv2,
  lineages = lineages,
  subtypes = subtypes,
  min_n = 40,
  min_events = 10,
  k = 5,
  repeats = 20,
  seed = 1,
  winsor_k = 3,
  unstable_beta = 10
)
paper3_write_csv(resD, file.path(fig_paths$tables, "cindex_grid.csv"))

pD_df <- resD %>%
  mutate(
    subtype = factor(subtype, levels = rev(subtypes)),
    lineage = factor(lineage, levels = lineages),
    label = ifelse(is.na(cindex), "NA", sprintf("%.2f", cindex))
  )

pD <- ggplot(pD_df, aes(x = lineage, y = subtype)) +
  geom_tile(aes(fill = cindex), colour = "grey85", linewidth = 0.5, na.rm = FALSE) +
  geom_text(aes(label = label), size = 3.8) +
  scale_fill_gradient(
    low = "white",
    high = "#2c7fb8",
    limits = c(0.50, 1.00),
    oob = squish,
    na.value = "grey80",
    name = "C-index\n(repeated 5-fold CV)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.margin = margin(6, 6, 6, 6)
  )

