library(data.table)
library(dplyr)

birth_year_baseline <- snakemake@params[["birth_year_baseline"]]

cat("[1/6] Loading cis_pQTLs manifest...\n")
cis_pQTLs <- fread(snakemake@input[["manifest"]]) %>%
  mutate(
    geno_col    = paste0(ID, "_", ref),
    protein_col = paste0("olink_instance_0.", tolower(protein))
  )
saveRDS(cis_pQTLs, snakemake@output[["cis_pqtls_rds"]])

cat("[2/6] Loading genotype dosage (selecting relevant columns)...\n")
geno_cols_needed <- c("FID", "IID", cis_pQTLs$geno_col) %>% unique()
df_geno <- fread(snakemake@input[["dosage"]], select = geno_cols_needed)
saveRDS(df_geno, snakemake@output[["geno_rds"]])

valid_eids <- df_geno$FID

cat("[3/6] Loading Olink data...\n")
df_olink <- fread(snakemake@input[["olink"]]) %>%
  filter(olink_instance_0.eid %in% valid_eids) %>%
  select("olink_instance_0.eid", contains("olink_instance_0")) %>%
  rename(eid = olink_instance_0.eid)

cat("[4/6] Loading covariates...\n")
df_covar <- fread(
  snakemake@input[["covariates"]],
  select = c("eid", "54-0.0", "22000-0.0", "3166-0.0", "34-0.0", "31-0.0")
) %>%
  filter(eid %in% valid_eids) %>%
  mutate(age = birth_year_baseline - `34-0.0`) %>%
  rename(
    sex            = `31-0.0`,
    UKB_center     = `54-0.0`,
    genotype_batch = `22000-0.0`,
    time_sample    = `3166-0.0`
  )

cat("[5/6] Loading PCA...\n")
df_pca <- fread(snakemake@input[["eigenvec"]]) %>%
  filter(IID %in% df_geno$IID) %>%
  rename(eid = IID) %>%
  select(-`#FID`)

cat("[6/6] Merging base data frame...\n")
df_base <- df_olink %>%
  inner_join(df_covar, by = "eid") %>%
  inner_join(df_pca,   by = "eid") %>%
  mutate(
    sex            = as.factor(sex),
    UKB_center     = as.factor(UKB_center),
    genotype_batch = as.factor(genotype_batch)
  )
saveRDS(df_base, snakemake@output[["base_rds"]])

cat("[batch] Loading batch data...\n")
df_batch <- fread(snakemake@input[["batch"]]) %>%
  mutate(Protein = ifelse(Protein == "WARS1", "WARS", Protein))
saveRDS(df_batch, snakemake@output[["batch_rds"]])

cat("Preprocessing complete.\n")
