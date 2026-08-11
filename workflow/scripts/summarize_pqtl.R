library(data.table)
library(dplyr)
library(stringr)

result_files <- snakemake@input[["results"]]

collect_snp_row <- function(f) {
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) < 2) return(NULL)
  snp_row <- dt[2, ]
  basename_no_ext <- tools::file_path_sans_ext(basename(f))
  parts <- str_match(basename_no_ext, "^(.+?)_(.+)_RINT$")
  snp_row$protein <- parts[, 2]
  snp_row$variant <- parts[, 3]
  snp_row
}

results <- rbindlist(lapply(result_files, collect_snp_row), fill = TRUE)
results <- results %>%
  select(protein, variant, everything()) %>%
  arrange(`Pr(>|t|)`)

fwrite(results, snakemake@output[["summary"]])
cat(sprintf("Summarized %d pQTL results\n", nrow(results)))
