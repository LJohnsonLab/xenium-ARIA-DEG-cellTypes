# Export per-cell-type DEG results (byRegion RDS) to JSON
# RDS structure: named list by brain region, each element is a dataframe with gene rownames
# Mirrors the "save-results" chunk in 20260304-DEG_allCellTypes_brainRegions.qmd

library(tidyverse)
library(jsonlite)

rds_to_json <- function(rds_path) {
  json_path <- sub("\\.rds$", ".json", rds_path)
  readRDS(rds_path) |>
    map(\(df) rownames_to_column(df, var = "gene")) |>
    toJSON(pretty = TRUE) |>
    writeLines(json_path)
  message("Written: ", json_path)
}

rds_to_json("results/20260223-DEG_astrocytes_byRegion.rds")
rds_to_json("results/20260303-DEG_MG_byRegion.rds")
