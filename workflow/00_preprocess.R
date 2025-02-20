# ------------------------------------------------------
# File: load_DNA_samples.R
# Authors: Abraham Sotelo
# Date: 2025-02-20
#
# Description: Load and clean DNA samples
#
# Inputs:
#  - Directory containing DNA fastq files

# Outputs:
#
# ------------------------------------------------------
suppressPackageStartupMessages({
                                library(magrittr)
                                library(R.utils)})

# Load configuration ------------------------------------
source("workflow/utils/config.R")
config <- get_config("preprocessing")
raw_samples_path <- config$raw_samples_path
dna_sequences_path <- config$dna_sequences_path
raw_data_dictionary <- config$raw_data_dictionary

# Load state management ----------------------------------
source("workflow/utils/state.R")
update_state <- update_state
check_sample_stage <- check_sample_stage

# ------------------------------------------------------
# Internal functions
# ------------------------------------------------------

#' Clean DNA sequence filename
#'
#' @param filename A string containing the filename to clean
#' @return A string with "DNA-" prefix and ".gz" suffix removed, and fq extension replaced with fastq
#' @keywords internal
internal_clean_filename <- function(filename) {
  . <- NULL # Placeholder for filename, hack to avoid no visible binding for global variable ‘.’
  sub("DNA-", "", filename) %>% sub("[.]gz$", "", .) %>% sub("[.]fq$", ".fastq", .)
}


# ------------------------------------------------------
# External functions
# ------------------------------------------------------

#' Decompress raw DNA sample files
#'
#' @param raw_samples_path Path to directory containing compressed raw samples
#' @param dna_sequences_path Path to directory where decompressed files will be saved
#' @param raw_data_dictionary Dictionary mapping source directories to destination directories
#' @return None
#' @export
decompress_raw_samples <- function(raw_samples_path, dna_sequences_path, raw_data_dictionary) {
  start_time <- Sys.time()
  cat("Decompressing raw samples", "\n")
  if (!dir.exists(raw_samples_path)) {
    stop("Raw samples directory does not exist")
    return("Error")
  }
  directories <- list.dirs(raw_samples_path, full.names = FALSE, recursive = FALSE)
  for (dir in directories) {
    sample <- raw_data_dictionary[[dir]]
    new_dir <- file.path(dna_sequences_path, sample)
    new_dir_fwd <- file.path(new_dir, "fwd")
    new_dir_rev <- file.path(new_dir, "rev")

    # Check if these samples have been decompressed
    if (Sys.getenv("TESTING") != "TRUE") {
      if (check_sample_stage(sample = sample, stage = "decompress")) {
        cat("Sample", sample, "already decompressed", "\n")
        next
      }
    }

    # Actual decompression
    cat("Decompressing files in directory:", dir, "\n")
    files <- list.files(file.path(raw_samples_path, dir), pattern = "\\.gz$")
    decompressed_files_names <- internal_clean_filename(files)
    mapply(function(file, new_file) {
                                     dest_dir <- ifelse(grepl("1[.]fastq", new_file), new_dir_fwd, new_dir_rev)
                                     dest <- file.path(dest_dir, new_file)
                                     cat("Decompressing file:", file, "->", dest, "\n")
                                     gunzip(file, destname = dest, remove = FALSE, overwrite = TRUE)},
    file.path(raw_samples_path, dir, files), decompressed_files_names)

    # Update state
    if (Sys.getenv("TESTING") != "TRUE") {
      decompressed_files <- list.files(new_dir, recursive = TRUE)
      decompressed_files_paths <- list.files(new_dir, full.names = TRUE, recursive = TRUE)
      print(decompressed_files)
      data <- list(
        raw_directory = file.path(raw_samples_path, dir),
        raw_files = files,
        raw_files_number = length(files),
        raw_size = paste(round(sum(file.size(file.path(raw_samples_path, dir, files))) / 1024^3, 2), "GB"),
        decompress_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        decompress_files = decompressed_files,
        decompress_files_number = length(list.files(new_dir)),
        decompress_directory = new_dir,
        decompress_size = paste(round(sum(file.size(decompressed_files_paths)) / 1024^3, 2), "GB"),
        decompress_time = paste(round(as.numeric(difftime(Sys.time(), start_time, units = "sec")), 2), "sec")
      )
      update_state(sample = sample, stage = "decompress", data = data)
    }
  }
}


# ------------------------------------------------------
# Workflow
# ------------------------------------------------------
preprocessing <- function() {
  cat("Preprocessing DNA samples", "\n")
  decompress_raw_samples(raw_samples_path, dna_sequences_path, raw_data_dictionary)
}