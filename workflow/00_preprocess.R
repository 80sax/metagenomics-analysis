# ------------------------------------------------------
# File: load_DNA_samples.R
# Authors: Abraham Sotelo
# Date: 2025-02-27
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
  library(R.utils)
})

# Load configuration ------------------------------------
source("workflow/utils/config.R")
config <- get_config("preprocessing")
raw_samples_path <- config$raw_samples_path
dna_sequences_path <- config$dna_sequences_path
raw_data_dictionary <- config$raw_data_dictionary

# Load state management and orchestration --------------
source("workflow/utils/state.R")
update_state <- update_state
check_sample_stage <- check_sample_stage
source("workflow/utils/orchestrate.R")
apply_stage_to_project <- apply_stage_to_project

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
    mapply(
      function(file, new_file) {
        dest_dir <- ifelse(grepl("1[.]fastq", new_file), new_dir_fwd, new_dir_rev)
        dest <- file.path(dest_dir, new_file)
        cat("Decompressing file:", file, "->", dest, "\n")
        gunzip(file, destname = dest, remove = FALSE, overwrite = TRUE)
      }, file.path(raw_samples_path, dir, files), decompressed_files_names
    )

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


#' Identify Lines with Mismatched Sequence and Quality Lengths in FASTQ Files
#'
#' This function reads a FASTQ file and identifies lines where the sequence length
#' does not match its corresponding quality score length, which indicates potential
#' data corruption or format issues.
#'
#' @param file character, Path to the FASTQ file to be analyzed
#'
#' @return Returns either:
#'   - A numeric vector containing the indices of faulty lines where sequence and quality lengths don't match
#'   - NULL if no faulty lines are found
#'
#' @details
#' The function works by:
#' 1. Reading the file line by line
#' 2. Comparing the lengths of sequence lines (every 2nd line) with their corresponding quality score lines (every 4th line)
#' 3. Identifying any mismatches between these lengths
#'
#' @examples
#' \dontrun{
#' faulty_lines <- identify_faulty_lines("path/to/fastq/file.fastq")
#' }
#'
#' @keywords file-processing internal
identify_faulty_lines_file <- function(file) {
  cat("Analysing file:", file, "\n")
  lines <- readLines(file)
  cat("Number of lines:", length(lines), "\n")
  seqlens <- sapply(lines[seq(2, length(lines), 4)], nchar, USE.NAMES = FALSE)
  quallens <- sapply(lines[seq(4, length(lines), 4)], nchar, USE.NAMES = FALSE)
  faulty_lines <- which(seqlens != quallens)
  if (length(faulty_lines) > 0) {
    cat("Faulty lines:", faulty_lines, "\n")
    return(setNames(list(faulty_lines), file))
  } else {
    cat("No faulty lines found", "\n")
    return(NULL)
  }
}


#' Identify Faulty Lines in Sample Directory
#'
#' This function examines all forward and reverse sequence files in a sample directory
#' to identify faulty lines in each file.
#'
#' @param sample_dir Character string specifying the path to the sample directory
#'   containing 'fwd' and 'rev' subdirectories with sequence files
#'
#' @return A list with one element:
#'   \describe{
#'     \item{faulty_files}{A list of identified faulty lines from all files,
#'       with NULL results filtered out}
#'   }
#'
#' @details
#' The function looks for sequence files in both 'fwd' and 'rev' subdirectories
#' of the provided sample directory. It processes each file to identify faulty lines
#' and combines the results into a single list.
#'
#' @seealso identify_faulty_lines_file
#'
#' @examples
#' \dontrun{
#' faulty_lines <- identify_faulty_lines_sample("path/to/sample/directory")
#' }
identify_faulty_lines_sample <- function(sample_dir, ...) {
  cat("Identifying faulty lines in sample:", sample_dir, "\n")
  fwd_files <- list.files(file.path(sample_dir, "fwd"), full.names = TRUE)
  rev_files <- list.files(file.path(sample_dir, "rev"), full.names = TRUE)
  files <- c(fwd_files, rev_files)
  faulty_lines <- unlist(lapply(files, identify_faulty_lines_file), recursive = FALSE)
  faulty_lines <- faulty_lines[!sapply(faulty_lines, is.null)]
  return(list(faulty_files = faulty_lines))
}

#' Identifies faulty lines in DNA sequence data across samples
#'
#' This function processes DNA sequence data files to identify problematic or faulty lines
#' by applying a sample-level fault detection function across all samples in a project.
#'
#' @param dna_sequences_path Character string specifying the path to the DNA sequences data
#'
#' @return List of identified faulty lines per sample
#'
#' @details
#' The function utilizes the `apply_stage_to_project` framework to process multiple samples,
#' applying the `identify_faulty_lines_sample` function to each individual sample file.
#'
#' @seealso
#' \code{\link{identify_faulty_lines_sample}}
#' \code{\link{apply_stage_to_project}}
#'
#' @export
identify_faulty_lines <- function(dna_sequences_path) {
  apply_stage_to_project(dna_sequences_path, stage = "faulty_lines", func = identify_faulty_lines_sample)
}


#Quality profile



quality_profile_sample <- function(sample_dir) {
  start_time <- Sys.time()
  cat("Plotting quality profile for sample:", sample_dir, "\n")
  fwd_files <- list.files(file.path(sample_dir, "fwd"), full.names = TRUE)
  rev_files <- list.files(file.path(sample_dir, "rev"), full.names = TRUE)

  number_of_files <- length(fwd_files)
  number_of_samples <- min(number_of_files, 12)

  fwd_samples <- fwd_files[sample(1:number_of_files, number_of_samples, replace = FALSE)]
  rev_samples <- rev_files[sample(1:number_of_files, number_of_samples, replace = FALSE)]

  cat("Forward samples:", fwd_samples, sep = "\n")
  cat("Reverse samples:", rev_samples, sep = "\n")

  p <- dada2:::plotQualityProfile(fwd_samples)
  #q <- dada2:::plotQualityProfile(rev_samples)

  #ggplot2::ggsave("quality_profile.png", plot = p, width = 8, height = 6, dpi = 300)
  #ggplot2::ggsave("quality_profile.pdf", plot = q, width = 8, height = 6)

  saveRDS(p, "quality_profile_plot.rds")  # Save to file
  print(difftime(Sys.time(), start_time, units = "sec"))
}



# ------------------------------------------------------
# Workflow
# ------------------------------------------------------
preprocessing <- function() {
  cat("Preprocessing DNA samples", "\n")
  decompress_raw_samples(raw_samples_path, dna_sequences_path, raw_data_dictionary)
  identify_faulty_lines(dna_sequences_path)
  quality_profile(dna_sequences_path)
}