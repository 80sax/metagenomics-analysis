# ------------------------------------------------------
# File: load_DNA_samples.R
# Authors: Abraham Sotelo
# Date: 2025-02-13
#
# Description: Workflow configuration loader
# ------------------------------------------------------

get_config <- function(module = c("state", "preprocessing")) {
  module <- match.arg(module)
  config <- yaml::read_yaml("config.yaml")

  if (!module %in% names(config)) {
    stop(sprintf("Configuration module '%s' not found", module))
  }

  config[[module]]
}
