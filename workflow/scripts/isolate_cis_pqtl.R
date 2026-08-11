library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

target_genes <- snakemake@params[["genes"]]
chr_ids_dir  <- snakemake@output[["chr_ids"]]
dir.create(chr_ids_dir, recursive = TRUE, showWarnings = FALSE)

df_sumstats <- fread(snakemake@input[["sumstats"]]) %>%
  rename(bp_x = `GENPOS (hg38)`)

df_clumped <- fread(snakemake@input[["clumped"]])

if (length(target_genes) > 0) {
  df_clumped <- df_clumped %>% filter(protein %in% target_genes)
}

df_clumped <- df_clumped %>% select(protein, chr, bp_x)

df_join <- left_join(df_clumped, df_sumstats, by = "bp_x") %>%
  select(`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, protein) %>%
  unique() %>%
  rename(hg37 = `Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`) %>%
  mutate(hg37 = str_remove(hg37, ":imp:v1")) %>%
  separate(hg37, into = c("chr", "bp", "ref", "alt"), sep = ":") %>%
  mutate(
    ID = paste0(chr, ":", bp, "_", ref, "_", alt),
    geno_col = paste0(ID, "_", ref),
    protein_col = paste0("olink_instance_0.", tolower(protein))
  )

for (chr_val in unique(df_join$chr)) {
  df_chr <- df_join %>% filter(chr == chr_val) %>% select(ID)
  fwrite(
    df_chr,
    file.path(chr_ids_dir, paste0("chr", chr_val, "_IDs.txt")),
    col.names = FALSE, quote = FALSE, sep = "\t"
  )
}

fwrite(df_join, snakemake@output[["manifest"]])
cat(sprintf("Identified %d cis-pQTL variant-protein pairs across %d chromosomes\n",
            nrow(df_join), length(unique(df_join$chr))))
