script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grepl("--file=", commandArgs(FALSE))]), mustWork = TRUE))

run_script <- function(path) {
  message("Running: ", path)
  status <- system2("Rscript", path)
  if (!identical(status, 0L)) stop("Script failed: ", path)
}

run_script(file.path(script_dir, "01_oncology_analysis_pipeline.R"))
run_script(file.path(script_dir, "08_oncology_figures.R"))
run_script(file.path(script_dir, "08_restore_previous_drug_figures_keep_new_versions.R"))
run_script(file.path(script_dir, "09_oncology_qc_reports.R"))

message("Oncology therapeutic death landscape pipeline completed.")
