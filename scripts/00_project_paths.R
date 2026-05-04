`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

paper3_project_root <- function() {
  root <- Sys.getenv("PAPER3_PROJECT_ROOT", unset = getwd())
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

paper3_paths <- function(root = paper3_project_root()) {
  list(
    root = root,
    data_links = file.path(root, "data_links"),
    docs = file.path(root, "docs"),
    scripts = file.path(root, "scripts"),
    outputs = file.path(root, "outputs"),
    logs = file.path(root, "logs"),
    cache = file.path(root, "cache"),
    legacy_project = file.path(root, "data_links", "paper3_breast_wgan_win"),
    legacy_gan_output_v1 = file.path(root, "data_links", "paper3_breast_wgan_win", "results", "gan_output_v1", "output"),
    legacy_gan_output_v2 = file.path(root, "data_links", "paper3_breast_wgan_win", "results", "gan_output_v2", "output2"),
    output_figures = file.path(root, "outputs", "figures"),
    output_checkpoints = file.path(root, "outputs", "checkpoints"),
    python_venv = file.path(root, ".venv", "bin", "python")
  )
}

paper3_require_paths <- function(paths = paper3_paths()) {
  must_exist <- c("legacy_project")

  for (nm in must_exist) {
    if (!file.exists(paths[[nm]]) && !dir.exists(paths[[nm]])) {
      stop(sprintf("Required path is missing: %s -> %s", nm, paths[[nm]]), call. = FALSE)
    }
  }

  invisible(paths)
}

paper3_mkdir <- function(...) {
  dirs <- c(...)
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }
  invisible(dirs)
}

paper3_assert_file <- function(path, label = NULL) {
  if (!file.exists(path)) {
    label <- label %||% basename(path)
    stop(sprintf("Required file is missing: %s -> %s", label, path), call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

paper3_assert_dir <- function(path, label = NULL) {
  if (!dir.exists(path)) {
    label <- label %||% basename(path)
    stop(sprintf("Required directory is missing: %s -> %s", label, path), call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

paper3_legacy_candidates <- function(relative_path, paths = paper3_paths()) {
  legacy_root <- normalizePath(paths$legacy_project, winslash = "/", mustWork = TRUE)
  relative_path <- gsub("\\\\", "/", relative_path)
  direct <- gsub("\\\\", "/", file.path(legacy_root, relative_path))

  entries <- list.files(
    legacy_root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = TRUE,
    no.. = TRUE
  )
  entries <- gsub("\\\\", "/", entries)
  suffix <- paste0("/", relative_path)
  tail_name <- basename(relative_path)

  matches <- entries[
    endsWith(entries, suffix) |
      basename(entries) == tail_name
  ]

  unique(c(direct, matches[order(nchar(matches))]))
}

paper3_legacy_path <- function(..., paths = paper3_paths(), required = FALSE, label = NULL) {
  relative_path <- do.call(file.path, as.list(c(...)))
  candidates <- paper3_legacy_candidates(relative_path, paths = paths)

  for (candidate in candidates) {
    if (file.exists(candidate) || dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  if (required) {
    label <- label %||% basename(relative_path)
    stop(
      sprintf(
        "Required legacy path is missing: %s -> %s",
        label,
        paste(candidates, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  gsub("\\\\", "/", file.path(paths$legacy_project, relative_path))
}

paper3_legacy_windows_root <- function(paths = paper3_paths()) {
  dirname(normalizePath(paths$legacy_project, winslash = "/", mustWork = TRUE))
}

paper3_fig_paths <- function(fig_name, paths = paper3_paths()) {
  fig_root <- file.path(paths$output_figures, fig_name)
  list(
    root = fig_root,
    raw = file.path(fig_root, "raw"),
    plots = file.path(fig_root, "plots"),
    tables = file.path(fig_root, "tables"),
    rds = file.path(fig_root, "rds"),
    logs = file.path(fig_root, "logs"),
    checkpoints = file.path(paths$output_checkpoints, fig_name)
  )
}

paper3_init_fig <- function(fig_name, paths = paper3_init()) {
  fig_paths <- paper3_fig_paths(fig_name, paths = paths)
  paper3_mkdir(
    fig_paths$root,
    fig_paths$raw,
    fig_paths$plots,
    fig_paths$tables,
    fig_paths$rds,
    fig_paths$logs,
    fig_paths$checkpoints
  )
  fig_paths
}

paper3_prepare_for_csv <- function(x) {
  if (!is.data.frame(x)) {
    return(x)
  }

  out <- x
  list_cols <- names(out)[vapply(out, is.list, logical(1))]
  for (nm in list_cols) {
    out[[nm]] <- vapply(out[[nm]], function(value) {
      if (is.null(value) || length(value) == 0) {
        return("")
      }
      paste(as.character(unlist(value, use.names = FALSE)), collapse = "; ")
    }, character(1))
  }
  out
}

paper3_write_csv <- function(x, path, row.names = FALSE) {
  utils::write.csv(paper3_prepare_for_csv(x), file = path, row.names = row.names)
  invisible(path)
}

paper3_write_lines <- function(x, path) {
  writeLines(text = x, con = path, useBytes = TRUE)
  invisible(path)
}

paper3_stage_reference_file <- function(source_path, destination_path, label = NULL, required = TRUE) {
  if (!file.exists(source_path)) {
    if (required) {
      label <- label %||% basename(source_path)
      stop(sprintf("Required reference file is missing: %s -> %s", label, source_path), call. = FALSE)
    }
    return("")
  }

  paper3_mkdir(dirname(destination_path))
  ok <- file.copy(source_path, destination_path, overwrite = TRUE, copy.date = TRUE)
  if (!ok) {
    stop(sprintf("Failed to stage reference file: %s -> %s", source_path, destination_path), call. = FALSE)
  }
  normalizePath(destination_path, winslash = "/", mustWork = TRUE)
}

paper3_preferred_file <- function(candidates, label = "file", required = TRUE) {
  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  if (required) {
    stop(
      sprintf(
        "Required %s was not found. Checked: %s",
        label,
        paste(gsub("\\\\", "/", path.expand(candidates)), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ""
}

paper3_legacy_file <- function(..., paths = paper3_paths(), label = NULL) {
  target <- paper3_legacy_path(..., paths = paths, required = TRUE, label = label)
  paper3_assert_file(target, label = label %||% basename(target))
}

paper3_capture_session_info <- function(path) {
  info <- utils::capture.output(utils::sessionInfo())
  paper3_write_lines(info, path)
}

paper3_require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    stop(
      sprintf(
        "Missing required packages: %s. Install packages first or run scripts/check_package_status.R.",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

paper3_python_bin <- function(paths = paper3_paths(), required = FALSE) {
  candidates <- c(
    Sys.getenv("PAPER3_PYTHON_BIN", unset = ""),
    paths$python_venv,
    Sys.which("python3")
  )

  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(gsub("\\\\", "/", path.expand(candidate)))
    }
  }

  if (required) {
    stop(
      "No Python interpreter was found for breast_WGAN_model_2026. ",
      "Create .venv with scripts/setup/install_python_requirements.sh or set PAPER3_PYTHON_BIN.",
      call. = FALSE
    )
  }

  ""
}

paper3_fig03_output_stubs <- function() {
  c(
    "B_cells",
    "Endothelial",
    "Epithelial",
    "Fibroblast",
    "Mast_cells",
    "Myeloid",
    "Pericytes",
    "T_cells"
  )
}

paper3_fig03_expected_raw_files <- function(raw_dir) {
  output_stubs <- paper3_fig03_output_stubs()
  c(
    file.path(raw_dir, paste0("training_history_", output_stubs, ".csv")),
    file.path(raw_dir, paste0("ranked_genes_", output_stubs, ".csv")),
    file.path(raw_dir, paste0("k_curve_knn_", output_stubs, ".csv")),
    file.path(raw_dir, paste0("k_curve_lr_", output_stubs, ".csv"))
  )
}

paper3_fig03_raw_complete <- function(raw_dir) {
  dir.exists(raw_dir) && all(file.exists(paper3_fig03_expected_raw_files(raw_dir)))
}

paper3_fig03_preferred_raw_source <- function(paths = paper3_paths()) {
  fig_paths <- paper3_fig_paths("fig03_gan_signature_discovery", paths = paths)
  legacy_windows_root <- paper3_legacy_windows_root(paths = paths)
  candidates <- c(
    Sys.getenv("PAPER3_FIG03_RAW_DIR", unset = ""),
    file.path(paths$legacy_project, "results", "gan_output_v2", "output2"),
    file.path(paths$legacy_project, "results", "gan_output_v1", "output"),
    file.path(legacy_windows_root, "output-GAN2", "output2"),
    file.path(legacy_windows_root, "output-GAN", "output"),
    paths$legacy_gan_output_v2,
    paths$legacy_gan_output_v1,
    fig_paths$raw
  )

  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    if (paper3_fig03_raw_complete(candidate)) {
      return(gsub("\\\\", "/", path.expand(candidate)))
    }
  }

  fig_paths$raw
}

paper3_fig03_gene_list_path <- function(paths = paper3_paths(), required = TRUE) {
  fig3_paths <- paper3_fig_paths("fig03_gan_signature_discovery", paths = paths)
  paper3_preferred_file(
    c(
      file.path(fig3_paths$rds, "topk_gene_lists_by_celltype_new.rds"),
      file.path(fig3_paths$root, "topk_gene_lists_by_celltype_new.rds"),
      paper3_legacy_path("topk_gene_lists_by_celltype_new.rds", paths = paths, required = FALSE)
    ),
    label = "Fig. 3 top-k gene list RDS",
    required = required
  )
}

paper3_fig04_lineage_gene_path <- function(paths = paper3_paths(), required = TRUE) {
  fig4_paths <- paper3_fig_paths("fig04_pathway_lineage_analysis", paths = paths)
  paper3_preferred_file(
    c(
      paper3_legacy_path("genes_by_lineage.rds", paths = paths, required = FALSE),
      file.path(fig4_paths$rds, "genes_by_lineage.rds")
    ),
    label = "Fig. 4 lineage gene RDS",
    required = required
  )
}

paper3_fig05_subtype_gene_path <- function(paths = paper3_paths(), required = TRUE) {
  fig5_paths <- paper3_fig_paths("fig05_bulk_projection_survival", paths = paths)
  paper3_preferred_file(
    c(
      paper3_legacy_path("genes_by_lineage_subtype_sig.rds", paths = paths, required = FALSE),
      file.path(fig5_paths$rds, "genes_by_lineage_subtype_sig.rds")
    ),
    label = "Fig. 5 subtype-lineage signature RDS",
    required = required
  )
}

paper3_fig05_coef_path <- function(paths = paper3_paths(), required = TRUE) {
  fig5_paths <- paper3_fig_paths("fig05_bulk_projection_survival", paths = paths)
  paper3_preferred_file(
    c(
      paper3_legacy_path("lineage_subtype_coefs.coxph.rds", paths = paths, required = FALSE),
      file.path(fig5_paths$rds, "lineage_subtype_coefs.coxph.rds")
    ),
    label = "Fig. 5 subtype-lineage coefficient RDS",
    required = required
  )
}

paper3_fig05_lineage_sig_path <- function(paths = paper3_paths(), required = TRUE) {
  fig5_paths <- paper3_fig_paths("fig05_bulk_projection_survival", paths = paths)
  paper3_preferred_file(
    c(
      paper3_legacy_path("genes_by_lineage_sig.rds", paths = paths, required = FALSE),
      file.path(fig5_paths$rds, "genes_by_lineage_sig.rds")
    ),
    label = "Fig. 5 lineage signature RDS",
    required = required
  )
}

paper3_fig06_supp_table_path <- function(paths = paper3_paths(), required = TRUE) {
  fig6_paths <- paper3_fig_paths("fig06_lineage_subtype_synthesis", paths = paths)
  paper3_preferred_file(
    c(
      paper3_legacy_path("SuppTable_Fig6B_by_lineage.csv", paths = paths, required = FALSE),
      file.path(fig6_paths$tables, "supp_table_fig6b_by_lineage.csv")
    ),
    label = "Fig. 6 supplementary pathway table CSV",
    required = required
  )
}

paper3_init <- function() {
  paths <- paper3_paths()
  paper3_require_paths(paths)
  paper3_mkdir(paths$outputs, paths$logs, paths$cache, paths$output_figures, paths$output_checkpoints)
  paths
}

paper3_copy_asset <- function(source_path, destination_path) {
  src <- normalizePath(source_path, winslash = "/", mustWork = TRUE)
  paper3_mkdir(dirname(destination_path))
  ok <- file.copy(from = src, to = destination_path, overwrite = TRUE, copy.date = TRUE)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to copy asset: %s -> %s", src, destination_path), call. = FALSE)
  }
  normalizePath(destination_path, winslash = "/", mustWork = TRUE)
}

paper3_reference_exists <- function(path) {
  file.exists(path) || dir.exists(path)
}

paper3_stage_legacy_figure <- function(
    fig_name,
    asset_rows,
    legacy_scripts = character(),
    reference_objects = character(),
    notes = character(),
    paths = paper3_init()) {
  fig_paths <- paper3_init_fig(fig_name, paths = paths)

  required_cols <- c("source", "stage", "description")
  missing_cols <- setdiff(required_cols, names(asset_rows))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("asset_rows is missing required columns: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  asset_rows <- as.data.frame(asset_rows, stringsAsFactors = FALSE)
  if (!("dest_name" %in% names(asset_rows))) {
    asset_rows$dest_name <- basename(asset_rows$source)
  }

  asset_rows$source <- vapply(asset_rows$source, normalizePath, character(1), winslash = "/", mustWork = TRUE)
  asset_rows$destination_path <- vapply(seq_len(nrow(asset_rows)), function(i) {
    stage_dir <- file.path(fig_paths$plots, paste0(asset_rows$stage[[i]], "_exports"))
    paper3_copy_asset(
      source_path = asset_rows$source[[i]],
      destination_path = file.path(stage_dir, asset_rows$dest_name[[i]])
    )
  }, character(1))

  asset_manifest <- asset_rows[, c("stage", "description", "source", "destination_path")]
  names(asset_manifest) <- c("stage", "description", "source_path", "destination_path")
  asset_manifest$source_exists <- TRUE
  paper3_write_csv(asset_manifest, file.path(fig_paths$tables, "legacy_asset_manifest.csv"))

  reference_frames <- list()
  if (length(legacy_scripts) > 0) {
    reference_frames[[length(reference_frames) + 1]] <- data.frame(
      reference_type = "legacy_script",
      description = basename(legacy_scripts),
      path = legacy_scripts,
      exists = vapply(legacy_scripts, paper3_reference_exists, logical(1)),
      stringsAsFactors = FALSE
    )
  }
  if (length(reference_objects) > 0) {
    reference_frames[[length(reference_frames) + 1]] <- data.frame(
      reference_type = "reference_object",
      description = basename(reference_objects),
      path = reference_objects,
      exists = vapply(reference_objects, paper3_reference_exists, logical(1)),
      stringsAsFactors = FALSE
    )
  }
  reference_manifest <- if (length(reference_frames) > 0) {
    do.call(rbind, reference_frames)
  } else {
    data.frame(
      reference_type = character(0),
      description = character(0),
      path = character(0),
      exists = logical(0),
      stringsAsFactors = FALSE
    )
  }

  if (nrow(reference_manifest) > 0) {
    reference_manifest$path <- vapply(reference_manifest$path, function(path) {
      if (!paper3_reference_exists(path)) {
        return(gsub("\\\\", "/", path.expand(path)))
      }
      normalizePath(path, winslash = "/", mustWork = TRUE)
    }, character(1))
    paper3_write_csv(reference_manifest, file.path(fig_paths$tables, "reference_manifest.csv"))
  }

  note_lines <- c(
    sprintf("Figure workspace: %s", fig_name),
    sprintf("Generated at: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "Mode: reuse legacy Windows assets and objects into standardized workspace outputs."
  )

  if (length(notes) > 0) {
    note_lines <- c(note_lines, "", notes)
  }

  paper3_write_lines(note_lines, file.path(fig_paths$logs, "provenance_notes.txt"))
  paper3_capture_session_info(file.path(fig_paths$logs, "session_info.txt"))

  invisible(fig_paths)
}
