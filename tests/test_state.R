# ------------------------------------------------------
# File: test_state.R
# Authors: Abraham Sotelo
# Date: 2025-02-14
#
# Description: Testing state management
# ------------------------------------------------------
library(testthat)

setwd("..")
source("workflow/utils/state.R")

test_that("empty state file is created correctly", {
  # Setup
  temp_dir <- tempdir()
  file_1_0 <- paste0(temp_dir, "/", file_prefix, "_v1.0.json")
  file_1_1 <- paste0(temp_dir, "/", file_prefix, "_v1.1.json")

  # Test
  utils::capture.output({
    internal_create_state_file(state_path = temp_dir)
    internal_create_state_file(state_path = temp_dir)
  })
  json_data <- jsonlite::fromJSON(file_1_1)

  # Assertions
  expect_true(file.exists(file_1_0))  # Initial version should exist
  expect_true(file.exists(file_1_1))  # Subsequent version should exist
  expect_true("metadata" %in% names(json_data))  # Metadata list should be present
  expect_true(all(c("created", "last_modified", "version", "number_of_samples") %in% names(json_data$metadata)))

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})