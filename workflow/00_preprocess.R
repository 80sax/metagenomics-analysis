# ------------------------------------------------------
# File: load_DNA_samples.R
# Authors: Abraham Sotelo
# Date: 2025-02-13
#
# Description: Load and clean DNA samples
#
# Inputs:
#  - Directory containing DNA fastq files

# Outputs:
#
# ------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(R.utils))

# Load configuration ------------------------------------
source("workflow/utils/config.R")
config <- get_config("preprocessing")
raw_samples_path <- config$raw_samples_path
dna_sequences_path <- config$dna_sequences_path
raw_data_dictionary <- config$raw_data_dictionary


# ------------------------------------------------------
# Internal functions
# ------------------------------------------------------
internal_clean_filename <- function(filename) {
  sub("DNA-", "", filename) %>% sub("[.]gz$", "", .)
}


# ------------------------------------------------------
# External functions
# ------------------------------------------------------

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


decompress_raw_samples(raw_samples_path, dna_sequences_path, raw_data_dictionary)