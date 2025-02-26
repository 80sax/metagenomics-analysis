# ------------------------------------------------------
# File: orchestrate.R
# Authors: Abraham Sotelo
# Date: 2025-02-25
#
# Description: Orchestrate the execution of stages
# ------------------------------------------------------

# Load state management --------------------------------
source("workflow/utils/state.R")
update_state <- update_state
check_sample_stage <- check_sample_stage

#' Apply a Processing Stage to Project Samples
#'
#' This function applies a specified processing stage to all sample directories within a project.
#' It tracks the processing state and avoids reprocessing samples that have already completed the stage.
#'
#' @param input_dir Character string specifying the input directory containing sample subdirectories
#' @param output_dir Character string specifying the output directory. If NULL, defaults to input_dir
#' @param stage Character string identifying the processing stage to be applied
#' @param func Function to be applied to each sample directory
#' @param ... Additional arguments passed to func
#'
#' @details
#' The function iterates through all immediate subdirectories in input_dir, checking if the specified
#' stage has already been applied to each sample. For samples that haven't been processed, it applies
#' the provided function and updates their processing state.
#'
#' @return None (called for side effects)
#'
#' @examples
#' \dontrun{
#' apply_stage_to_project("raw_data", "processed_data", "qc", quality_check_func)
#' }
#'
#' @export
apply_stage_to_project <- function(input_dir, output_dir = NULL, stage, func, state_file, ...) {
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