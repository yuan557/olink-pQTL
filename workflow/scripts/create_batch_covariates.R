library(data.table)
library(dplyr)

df_dosage <- fread(snakemake@input[["dosage"]], select = c("FID", "IID"))
valid_eids <- df_dosage$IID

main_df <- fread(snakemake@input[["batch_main"]]) %>%
  select(eid, PlateID, batch, date_blood_collection) %>%
  filter(eid %in% valid_eids)

df_process_time <- fread(snakemake@input[["processing_dates"]])
df_panel <- fread(snakemake@input[["panel_map"]])

df_full <- left_join(df_process_time, df_panel, by = "Panel", relationship = "many-to-many") %>%
  left_join(main_df, by = "PlateID", relationship = "many-to-many") %>%
  mutate(
    time_to_analysis = as.numeric(
      difftime(as.Date(Processing_StartDate), as.Date(date_blood_collection), units = "days")
    ) / 365.25
  ) %>%
  select(eid, Protein, batch, time_to_analysis) %>%
  unique()

fwrite(df_full, snakemake@output[["batch_csv"]])
cat(sprintf("Batch covariates: %d rows\n", nrow(df_full)))
