# ------------------------------------------------------
# File: test_script.R
# Authors: Abraham Sotelo
# Date: 2025-02-12
#
# Description: Testing framework
# ------------------------------------------------------

library(testthat)
source("../workflow/utils/status.R")

test_that("empty status file is created correctly", {
  # Setup
  temp_dir <- tempdir()
  temp_file <- file.path(temp_dir, "test_status.csv")

  # Test
  create_empty_status_file(temp_file)

  # Read the created file
  result <- read.csv(temp_file)

  # Assertions
  expect_true(file.exists(temp_file))
  expect_equal(ncol(result), 3)
  expect_equal(names(result), c("Raw.Sample", "Number.of.samples", "Stage"))
  expect_equal(nrow(result), 0)

  # Cleanup
  unlink(temp_file)
})