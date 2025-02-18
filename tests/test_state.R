# ------------------------------------------------------
# File: test_state.R
# Authors: Abraham Sotelo
# Date: 2025-02-17
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
  expect_true("samples" %in% names(json_data))  # Metadata list should be present
  expect_true(all(c("created", "last_modified", "version", "number_of_samples") %in% names(json_data$metadata)))

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})


test_that("state file is updated correctly", {
  # Setup
  temp_dir <- tempdir()
  file_1_0 <- paste0(temp_dir, "/", file_prefix, "_v1.0.json")

  # Test
  utils::capture.output({
    internal_create_state_file(state_path = temp_dir)
    data <- list(texts = c("text1", "text2"), number = 0)
    update_state(state_file = file_1_0, sample = "sample1", stage = "stage1", data = data)
  })
  json_data <- jsonlite::fromJSON(file_1_0)

  # Assertions
  expect_true(file.exists(file_1_0))
  expect_true("sample1" %in% names(json_data$samples))
  expect_equal("stage1", json_data$samples$sample1$current_stage)
  expect_true("stage1" %in% json_data$samples$sample1$stages)
  expect_equal(json_data$samples$sample1$texts, data$texts)
  expect_equal(json_data$samples$sample1$number, data$number)

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})
