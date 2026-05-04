source(file.path(getwd(), "scripts", "00_project_paths.R"))

paths <- paper3_init()
fig_paths <- paper3_stage_legacy_figure(
  fig_name = "fig02_archival_assets",
  asset_rows = data.frame(
    source = c(
      paper3_legacy_path("tiff_final", "fig2.tif", paths = paths, required = FALSE),
      paper3_legacy_path("tiff_final", "fig2_rationale.png", paths = paths, required = FALSE),
      paper3_legacy_path("tiff_final", "fig2_rationale.tif", paths = paths, required = FALSE)
    ),
    stage = c("final", "supporting", "supporting"),
    description = c(
      "Legacy final Fig. 2 TIFF export",
      "Legacy rationale panel PNG export",
      "Legacy rationale panel TIFF export"
    ),
    stringsAsFactors = FALSE
  ),
  legacy_scripts = character(),
  reference_objects = c(
    paper3_legacy_path("ppt", "fig2_rationale.pptx", paths = paths, required = FALSE)
  ),
  notes = c(
    "No standalone legacy Fig. 2 R script was found in the Windows project tree.",
    "This step standardizes the existing final figure assets and records their provenance so the scripts directory still contains a Fig. 2 entry."
  ),
  paths = paths
)

asset_note <- data.frame(
  item = c("figure_type", "legacy_script_status", "recommended_next_step"),
  value = c(
    "archival_asset_reuse",
    "not_found_in_windows_project_root",
    "port Fig. 2 generation later if upstream source code becomes available"
  ),
  stringsAsFactors = FALSE
)
paper3_write_csv(asset_note, file.path(fig_paths$tables, "fig02_asset_notes.csv"))
