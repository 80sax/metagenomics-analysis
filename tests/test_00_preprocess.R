# ------------------------------------------------------
# File: test_00_preprocess.R
# Authors: Abraham Sotelo
# Date: 2025-02-14
#
# Description: Testing preprocessing
# ------------------------------------------------------
library(testthat)

setwd("..")
source("workflow/00_preprocess.R")

test_that("Clean raw samples file names", {
  files <- c("DNA-12345.fq.gz", "67890.fq.gz", "abcde.fq")
  clean_files <- internal_clean_filename(files)
  expect_equal(internal_clean_filename(files), c("12345.fq", "67890.fq", "abcde.fq"))
})
