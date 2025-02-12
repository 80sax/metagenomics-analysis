# ------------------------------------------------------
# File: load_DNA_samples.R
# Authors: Abraham Sotelo
# Date: 2025-02-12
#
# Description: Load and clean DNA samples
#
# Inputs:
#  - Directory containing DNA fastq files

# Outputs:
#
# ------------------------------------------------------

suppressPackageStartupMessages(library(R.utils))

source("workflow/utils/config.R")

print(config)

#load_dna_samples <- function(raw_samples_path, dna_samples_path) {
#  # Load DNA samples
#  cat("Loading samples", "\n")
#  files <- list.files(raw_samples_path)
#  newname <- sub("DNA-", "", files) ## new name
#  file.rename(file.path(raw_samples_path,files), file.path(raw_samples_path, newname)) ##rename files
#
#  if (!dir.exists(dna_samples_path)) {
#    dir.create(dna_samples_path, recursive = TRUE)
#  }
#
#  extracted_files <- list()
#  for (file in list.files(path = raw_samples_path, pattern = ".fq.gz", full.names = TRUE)) {
#    cat("Extracting file: ", file, "\n")
#    extracted_file_path <- gunzip(file, remove = FALSE) ## unzip all files
#    new_file_path <- file.path(dna_samples_path, basename(extracted_file_path))
#    file.rename(extracted_file_path, new_file_path)
#    extracted_files <- c(extracted_files, new_file_path)
#  }
#  cat("\nExtracted files:\n")
#  cat(paste0(extracted_files, collapse = "\n"))
#  #for (file in extracted_files) {
#  #  cat(file, "\n")
#  }
#}
