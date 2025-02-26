# ------------------------------------------------------
# File: test_00_preprocess.R
# Authors: Abraham Sotelo
# Date: 2025-02-25
#
# Description: Testing preprocessing
# ------------------------------------------------------
library(testthat)

Sys.setenv(TESTING = "TRUE")
setwd("..")
source("workflow/00_preprocess.R")

raw_files <- c("DNA-12345_1.fq.gz", "67890_2.fq.gz", "file.txt")

test_that("Clean raw samples file names", {
  expect_equal(internal_clean_filename(raw_files), c("12345_1.fastq", "67890_2.fastq", "file.txt"))
})

test_that("Decompress raw samples", {
  # Setup
  raw     <- "raw"
  dest    <- "dest"
  dirs    <- c("raw", "dest")
  samples <- c("sample1", "sample2")
  temp    <- tempdir()
  sapply(file.path(temp, dirs), dir.create)
  sapply(file.path(temp, raw, samples), dir.create)
  sapply(file.path(temp, raw, samples[1], raw_files), file.create)
  sapply(file.path(temp, raw, samples[2], raw_files), file.create)

  # Test
  utils::capture.output({
    decompress_raw_samples(
      raw_samples_path = file.path(temp, raw),
      dna_sequences_path = file.path(temp, dest),
      raw_data_dictionary = list(
        "sample1" = samples[1],
        "sample2" = samples[2]
      )
    )
  })

  # Assertions
  expect_true(file.exists(file.path(temp, dest, samples[1], "fwd", "12345_1.fastq")))
  expect_true(file.exists(file.path(temp, dest, samples[1], "rev", "67890_2.fastq")))
  expect_true(file.exists(file.path(temp, dest, samples[2], "fwd", "12345_1.fastq")))
  expect_true(file.exists(file.path(temp, dest, samples[2], "rev", "67890_2.fastq")))
  expect_false(file.exists(file.path(temp, dest, samples[1], "file.txt")))
  expect_false(file.exists(file.path(temp, dest, samples[2], "file.txt")))

  unlink(temp, recursive = TRUE)
})

correct_fastq <- c(
  "@SEQ_ID_1",
  "AGCTTAGCTAGCTACG",  # Sequence line
  "+",
  "!!''**))%%%%%%^^",   # Quality line
  "@SEQ_ID_2",
  "CGTAGCTAGCTA",      # Sequence line
  "+",
  "!!**%%$$$###"      # Quality line
)

faulty_fastq <- c(
  "@SEQ_ID_1",
  "AGCTTAGCTAGCTACG",  # Sequence line
  "+",
  "!!''**))%%%%%%^^",   # Quality line
  "@SEQ_ID_2",
  "CGTAGCTAGCTA",      # Sequence line
  "+",
  "!!**%%$$$###@",      # Quality line
  "@SEQ_ID_3",
  "AGCTTAGCTAGCTACG",  # Sequence line
  "+",
  "!!''**))%%%%%%^^",   # Quality line
  "@SEQ_ID_4",
  "CGTAGCTAGCTA",      # Sequence line
  "+",
  "!!**%%$$$###@"      # Quality line
)

test_that("Getting faulty lines in file", {
  # Setup
  temp    <- tempdir()
  correct_file <- file.path(temp, "correct.fastq")
  faulty_file <- file.path(temp, "faulty.fastq")
  writeLines(correct_fastq, correct_file)
  writeLines(faulty_fastq, faulty_file)

  #Testing
  utils::capture.output({
    expect_equal(identify_faulty_lines_file(correct_file), NULL)
    expect_equal(identify_faulty_lines_file(faulty_file), setNames(list(c(2, 4)), faulty_file))
  })
})


test_that("Getting faulty lines in sample", {
  # Setup
  temp <- tempdir()
  sample_dir <- file.path(temp, "sample")
  dirs <- c("fwd", "rev")

  sapply(dirs, function(d) dir.create(file.path(sample_dir, d), recursive = TRUE))
  files <- list(
    correct = file.path(sample_dir, "fwd", "correct.fastq"),
    faulty1 = file.path(sample_dir, "rev", "faulty1.fastq"),
    faulty2 = file.path(sample_dir, "rev", "faulty2.fastq")
  )
  mapply(writeLines, list(correct_fastq, faulty_fastq, faulty_fastq), files)

  # Testing
  utils::capture.output({
    expect_equal(identify_faulty_lines_sample(sample_dir),
      list(faulty_files = setNames(
        list(c(2, 4), c(2, 4)),          # List of faulty lines for both files
        c(files$faulty1, files$faulty2)  # Corresponding file names
      ))
    )
  })

  # Cleanup
  unlink(temp, recursive = TRUE)
})


Sys.unsetenv("TESTING")