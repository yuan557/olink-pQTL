library(data.table)
library(dplyr)

protein_name <- snakemake@params[["protein"]]
variant_id   <- snakemake@params[["variant"]]
use_rint     <- snakemake@params[["use_rint"]]
n_pcs        <- snakemake@params[["n_pcs"]]

cis_pQTLs  <- readRDS(snakemake@input[["cis_pqtls_rds"]])
df_base    <- readRDS(snakemake@input[["base_rds"]])
df_batch   <- readRDS(snakemake@input[["batch_rds"]])
df_geno_all <- readRDS(snakemake@input[["geno_rds"]])

row_idx <- which(cis_pQTLs$protein == protein_name & cis_pQTLs$ID == variant_id)
if (length(row_idx) == 0) stop(sprintf("No manifest entry for %s / %s", protein_name, variant_id))

geno_col    <- cis_pQTLs$geno_col[row_idx[1]]
protein_col <- cis_pQTLs$protein_col[row_idx[1]]

cat(sprintf("Protein: %s | SNP: %s\n", protein_name, geno_col))

if (!geno_col %in% names(df_geno_all)) {
  cat(sprintf("Skipping: geno_col '%s' not in genotype data.\n", geno_col))
  fwrite(data.table(term = character(), Estimate = numeric(), `Std. Error` = numeric(),
                    `t value` = numeric(), `Pr(>|t|)` = numeric()),
         snakemake@output[["result"]])
  quit(status = 0)
}

df_geno <- df_geno_all[, .(eid = FID, geno = get(geno_col))]
setnames(df_geno, "geno", geno_col)

df_protein_batch <- df_batch %>%
  filter(Protein == protein_name) %>%
  select(eid, batch, time_to_analysis)

setDT(df_base);   setkey(df_base,  eid)
setDT(df_geno);   setkey(df_geno,  eid)
df_all <- df_base[df_geno, nomatch = 0]

setDT(df_protein_batch); setkey(df_protein_batch, eid)
setkey(df_all, eid)
df_all <- df_protein_batch[df_all, on = "eid"]
df_all[, batch := as.factor(batch)]

if (!protein_col %in% names(df_all)) {
  cat(sprintf("Skipping: protein_col '%s' not found.\n", protein_col))
  fwrite(data.table(term = character(), Estimate = numeric(), `Std. Error` = numeric(),
                    `t value` = numeric(), `Pr(>|t|)` = numeric()),
         snakemake@output[["result"]])
  quit(status = 0)
}

if (use_rint) {
  rint <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
  df_all[[protein_col]] <- rint(df_all[[protein_col]])
}

pc_terms <- paste0("PC", seq_len(n_pcs), collapse = " + ")
formula <- as.formula(paste0(
  "`", protein_col, "` ~ `", geno_col, "`",
  " + age + I(age^2) + sex + age:sex + I(age^2):sex",
  " + UKB_center + genotype_batch + batch + time_to_analysis",
  " + ", pc_terms
))

df_all <- droplevels(df_all)
fit  <- lm(formula, data = df_all)
coef <- summary(fit)$coefficients

fwrite(data.table(term = rownames(coef), coef), snakemake@output[["result"]])
cat("Done.\n")
