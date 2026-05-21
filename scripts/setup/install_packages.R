options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

native_r_version <- "4.5.3"
bioc_version <- "3.22"

cran_specs <- c(
  BiocManager = "1.30.27",
  remotes = "2.5.0",
  devtools = "2.5.0",
  Seurat = "5.4.0",
  cowplot = "1.2.0",
  dplyr = "1.2.0",
  forcats = "1.0.1",
  ggh4x = "0.3.1",
  ggalluvial = "0.12.6",
  ggplot2 = "4.0.2",
  ggpubr = "0.6.3",
  ggrepel = "0.9.8",
  gridExtra = "2.3",
  msigdbr = "26.1.0",
  patchwork = "1.3.2",
  png = "0.1-9",
  purrr = "1.2.1",
  RColorBrewer = "1.1-3",
  readr = "2.2.0",
  rstatix = "0.7.3",
  scales = "1.4.0",
  stringr = "1.6.0",
  survival = "3.5-8",
  survminer = "0.5.2",
  tibble = "3.3.1",
  tidyr = "1.3.2",
  yaml = "2.3.12",
  circlize = "0.4.17"
)

bioc_specs <- c(
  clusterProfiler = "4.18.4",
  ComplexHeatmap = "2.26.1"
)

package_version_string <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(NA_character_)
  }
  as.character(utils::packageVersion(pkg))
}

canonical_version_string <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(NA_character_)
  }
  as.character(numeric_version(gsub("-", ".", x, fixed = TRUE)))
}

version_matches <- function(pkg, expected) {
  identical(
    canonical_version_string(package_version_string(pkg)),
    canonical_version_string(expected)
  )
}

bootstrap_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("[bootstrap] ", pkg)
    utils::install.packages(pkg)
  }
}

install_cran_version <- function(pkg, expected) {
  current <- package_version_string(pkg)
  if (version_matches(pkg, expected)) {
    message("[skip] ", pkg, " ", expected)
    return(invisible(TRUE))
  }

  note <- if (is.na(current)) "" else paste0(" (current: ", current, ")")
  message("[install] ", pkg, " ", expected, note)
  remotes::install_version(
    package = pkg,
    version = expected,
    repos = getOption("repos"),
    upgrade = "never",
    dependencies = TRUE
  )

  installed <- package_version_string(pkg)
  if (!identical(canonical_version_string(installed), canonical_version_string(expected))) {
    stop(sprintf("Expected %s %s but found %s", pkg, expected, installed), call. = FALSE)
  }
}

install_bioc_version <- function(pkg, expected, bioc_version) {
  current <- package_version_string(pkg)
  if (version_matches(pkg, expected)) {
    message("[skip] ", pkg, " ", expected)
    return(invisible(TRUE))
  }

  url <- sprintf(
    "https://bioconductor.org/packages/%s/bioc/src/contrib/%s_%s.tar.gz",
    bioc_version, pkg, expected
  )
  note <- if (is.na(current)) "" else paste0(" (current: ", current, ")")
  message("[install] ", pkg, " ", expected, note)
  remotes::install_url(
    url = url,
    repos = BiocManager::repositories(version = bioc_version),
    upgrade = "never",
    dependencies = TRUE
  )

  installed <- package_version_string(pkg)
  if (!identical(canonical_version_string(installed), canonical_version_string(expected))) {
    stop(sprintf("Expected %s %s but found %s", pkg, expected, installed), call. = FALSE)
  }
}

message("Expecting native R version ", native_r_version)
message("Expecting Bioconductor release ", bioc_version)
message("Detected R version ", as.character(getRversion()))

bootstrap_package("BiocManager")
bootstrap_package("remotes")

if (!version_matches("BiocManager", cran_specs[["BiocManager"]])) {
  install_cran_version("BiocManager", cran_specs[["BiocManager"]])
}

if (!version_matches("remotes", cran_specs[["remotes"]])) {
  install_cran_version("remotes", cran_specs[["remotes"]])
}

if (!identical(as.character(BiocManager::version()), bioc_version)) {
  BiocManager::install(version = bioc_version, ask = FALSE, update = FALSE)
}

for (pkg in names(cran_specs)) {
  install_cran_version(pkg, cran_specs[[pkg]])
}

for (pkg in names(bioc_specs)) {
  install_bioc_version(pkg, bioc_specs[[pkg]], bioc_version = bioc_version)
}

message("R package installation script finished.")
