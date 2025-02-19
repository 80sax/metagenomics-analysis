# ------------------------------------------------------
# File: load_DNA_samples.R
# Authors: Abraham Sotelo
# Date: 2025-02-18
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
  cat("Decompressing raw samples", "\n")
  if (!dir.exists(raw_samples_path)) {
    stop("Raw samples directory does not exist")
    return("Error")
  }
  directories <- list.dirs(raw_samples_path, full.names = FALSE, recursive = FALSE)
  for (dir in directories) {
    cat("Decompressing files in directory:", dir, "\n")
    files <- list.files(file.path(raw_samples_path, dir), pattern = "\\.gz$")
    mapply(function(file, dest) {
                                 cat("Decompressing file:", file, "->", dest, "\n")
                                 gunzip(file, destname = dest, remove = FALSE, overwrite = TRUE)},
    file.path(raw_samples_path, dir, files),
    file.path(dna_sequences_path, raw_data_dictionary[[dir]], internal_clean_filename(files)))
  }
}



#decompress_raw_samples(raw_samples_path, dna_sequences_path, raw_data_dictionary)