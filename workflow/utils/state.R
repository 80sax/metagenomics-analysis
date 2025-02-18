# ------------------------------------------------------
# File: state.R
# Authors: Abraham Sotelo
# Date: 2025-02-17
#
# Description: Data pipeline status management
# ------------------------------------------------------

# Load configuration ------------------------------------
source("workflow/utils/config.R")
config <- get_config("state")
path <- config$state_path
file_prefix <- config$file_prefix

# ------------------------------------------------------
# Internal functions
# ------------------------------------------------------

#' Creates a new state file at the specified path
#'
#' @param state_path Path where the state file will be created. Defaults to 'path'
#'
#' @return None
#'
#' @details This internal function is responsible for initializing a new state file
#' at the specified location. The state file is used to track and maintain
#' the application's state information.
#'
#' @keywords internal
internal_create_state_file <- function(state_path = path) {
  cat("Creating state file", "\n")
  if (!dir.exists(state_path)) {
    dir.create(state_path, recursive = TRUE)
  }
  state_files <- list.files(state_path, pattern = paste0("^", file_prefix, ".*\\.json$"), full.names = TRUE)
  if (length(state_files) > 0) {
    cat("Found existing state files:\n")
    for (file in state_files) {
      cat("  -", basename(file), "\n")
    }
  }
  if (length(state_files) > 0) {
    versions <- gsub(paste0(file_prefix, "_v(.*)\\.json"), "\\1", basename(state_files))
    max_version <- max(versions)
    version_parts <- strsplit(max_version, "\\.")[[1]]
    version_parts[2] <- as.character(as.numeric(version_parts[2]) + 1)
    version <- paste(version_parts, collapse = ".")
  } else {
    version <- "1.0"
  }

  file <- paste0(state_path, "/", file_prefix, "_v", version, ".json")

  state <- list(
    metadata = list(
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      last_modified = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      version = version,
      number_of_samples = 0
    ),
    samples = c()
  )
  jsonlite::write_json(state, file, pretty = TRUE, auto_unbox = TRUE)
  cat("New state file created:", file, "\n")
}


# ------------------------------------------------------
# External functions
# ------------------------------------------------------


#' Update Workflow State
#'
#' Updates the state tracking file for a sample at a specific processing stage
#'
#' @param state_file Path to the state tracking file. If NULL, uses latest version
#' @param sample Character string specifying the sample identifier
#' @param stage Character string indicating the processing stage
#' @param data List or data frame containing the state information to be saved
#'
#' @return None (called for side effects)
#'
#' @details This function manages the state tracking system by updating or creating
#' state files that track the progress of sample processing through various workflow stages.
#'
#' @export
#'
update_state <- function(state_file = NULL, sample, stage, data) {
  # Load state
  # If no state file is provided, find the latest one
  cat("Updating sample:", sample, "\n")
  if (is.null(state_file)) {
    state_file <- list.files(path, pattern = paste0("^", file_prefix, ".*\\.json$"), full.names = TRUE)
    if (length(state_file) == 0) {
      internal_create_state_file()
      state_file <- list.files(path, pattern = paste0("^", file_prefix, ".*\\.json$"), full.names = TRUE)
    }
    latest_state_file <- max(state_file)
  } else {
    latest_state_file <- state_file
  }
  cat("State file loaded:", latest_state_file, "\n")
  state <- jsonlite::read_json(latest_state_file)

  # Getting the sample object
  if (!sample %in% names(state$samples)) {
    state$samples[[sample]] <- c()
    state$metadata$number_of_samples <- state$metadata$number_of_samples + 1
  }
  sample_object <- state$samples[[sample]]

  # Update state
  sample_object$current_stage <- stage
  if (!(stage %in% sample_object$stages)) {
    sample_object$stages <- append(sample_object$stages, stage)
  }
  sample_object <- append(sample_object, data)

  cat("Stage updated:", sample_object$current_stage, "\n")
  cat("Actions conducted so far:", toString(sample_object$stages), "\n")
  cat("Data updated:", jsonlite::toJSON(data, pretty = TRUE), "\n")
  state$samples[[sample]] <- sample_object
  jsonlite::write_json(state, latest_state_file, pretty = TRUE, auto_unbox = TRUE)
}
