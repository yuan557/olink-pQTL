library(data.table)
library(dplyr)

df_olink <- fread(snakemake@input[["olink"]])

ancestry_field <- snakemake@params[["ancestry_field"]]
ancestry_code  <- snakemake@params[["ancestry_code"]]
exclusion_field <- snakemake@params[["exclusion_field"]]

df_filtered <- df_olink %>%
  filter(.data[[ancestry_field]] == ancestry_code) %>%
  filter(is.na(.data[[exclusion_field]])) %>%
  select(olink_instance_0.eid) %>%
  rename(IID = olink_instance_0.eid) %>%
  mutate(`#FID` = IID) %>%
  select(`#FID`, IID)

fwrite(df_filtered, snakemake@output[["sample_ids"]], sep = "\t")
cat(sprintf("Selected %d samples\n", nrow(df_filtered)))
