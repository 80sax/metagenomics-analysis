# ------------------------------------------------------
# File: setup.R
# Author: Abraham Sotelo
# Date: 2024-12-26
# Description: This script is used to setup the environment
# ------------------------------------------------------

restart_session <- function(argument) {
  cat("Restarting R session and executing setup.R with argument:", argument, "\n")
  system(paste("Rscript setup.R ", argument))
  quit(save = "no")
}

args <- commandArgs(trailingOnly = TRUE)

# Prepare virtual environment
if ("--prepare" %in% args) {
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
  }

  # Deactivate the current environment
  cat("Deactivating the current environment\n")
  renv::status()
  renv::deactivate(clean = TRUE)
  restart_session("--init")
}

# Initialize the renv project
if ("--init" %in% args) {
  # Changing the library path and terminal in VSCode settings
  cat("Project's library path:", renv::paths$library(), "\n")
  vscode_settings_path <- ".vscode/settings.json"
  cat("Changing VSCode settings:", vscode_settings_path, "\n")
  if (file.exists(vscode_settings_path)) {
    vscode_settings <- jsonlite::fromJSON(vscode_settings_path)
    vscode_settings$r.libPaths <- list(renv::paths$library())
    vscode_settings$r.rterm.linux <- Sys.getenv("RADIAN_TERMINAL_PATH")
    jsonlite::write_json(vscode_settings, vscode_settings_path, pretty = TRUE, auto_unbox = TRUE)
  } else {
    cat(".vscode/settings.json file does not exist\n")
  }

  renv::init(bioconductor = TRUE, restart = TRUE)
  restart_session("--post-init")
}

# Install packages
if ("--post-init" %in% args) {
  options(repos = c(getOption("repos"), vscDebugger = "https://manuelhentschel.r-universe.dev"))

  # Install dependencies
  description_file <- read.dcf("DESCRIPTION")
  suggested_packages <- strsplit(description_file[1, "Suggests"], ",\\s*")[[1]]
  # Suggested packages
  for (package in suggested_packages) {
    cat("Checking package:", package, "\n")
    renv::install(package)
  }
}