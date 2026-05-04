options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

install_if_missing <- function(pkgs, installer, extra_args = list()) {
  for (pkg in pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      message("[skip] ", pkg)
      next
    }

    message("[install] ", pkg)
    do.call(installer, c(list(pkg), extra_args))
  }
}

cran_pkgs <- c(
  "Seurat",
  "cowplot",
  "dplyr",
  "forcats",
  "ggh4x",
  "ggalluvial",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "gridExtra",
  "msigdbr",
  "patchwork",
  "png",
  "purrr",
  "RColorBrewer",
  "readr",
  "rstatix",
  "scales",
  "stringr",
  "survival",
  "survminer",
  "tibble",
  "tidyr",
  "yaml",
  "circlize"
)

bioc_pkgs <- c(
  "clusterProfiler",
  "ComplexHeatmap"
)

install_if_missing(cran_pkgs, install.packages)
install_if_missing(
  bioc_pkgs,
  BiocManager::install,
  extra_args = list(update = FALSE, ask = FALSE)
)

message("R package installation script finished.")
