# ------------------------------------------------------
# File: status.R
# Authors: Abraham Sotelo
# Date: 2025-02-12
#
# Description: Data pipeline status management
# ------------------------------------------------------

create_empty_status_file <- function(output_file = "data/status/pipeline_status.csv") {
  titles <- c("Raw Sample", "Number of samples", "Stage")
  write.table(data.frame(matrix(ncol = length(titles), nrow = 0)), output_file,
              sep = ",", col.names = titles, row.names = FALSE)
}
