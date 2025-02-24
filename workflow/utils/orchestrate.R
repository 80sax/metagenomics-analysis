# ------------------------------------------------------
# File: orchestrate.R
# Authors: Abraham Sotelo
# Date: 2025-02-23
#
# Description: Orchestrate the execution of stages
# ------------------------------------------------------

# Load state management --------------------------------
source("workflow/utils/state.R")
update_state <- update_state
check_sample_stage <- check_sample_stage

apply_stage_to_project <- function(input_dir, output_dir = NULL, stage, func, ...) {
  if (is.null(output_dir)) {
    output_dir <- input_dir
  }

  for (full_dir in list.dirs(input_dir, full.names = TRUE, recursive = FALSE)) {
    cat("Processing directory:", full_dir, "\n")
    dir <- basename(full_dir)

    if (check_sample_stage(sample = dir, stage = stage)) {
      cat("The stage:", stage, "is already applied on sample:", dir, "\n")
      next
    }
    data <- func(full_dir, output_dir, ...)

    update_state(sample = dir, stage = stage, data = data)
  }
}