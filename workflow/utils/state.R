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
    samples = list()
  )
  jsonlite::write_json(state, file, pretty = TRUE, auto_unbox = TRUE)
  cat("New state file created:", file, "\n")
}


# ------------------------------------------------------
# External functions
# ------------------------------------------------------


write_state <- function() {

}
