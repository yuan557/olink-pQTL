# olink-pQTL

Snakemake pipeline that tests whether a single-nucleotide variant (SNV) is
associated with changes in plasma protein levels measured by the Olink Explore
platform. For each variant-protein pair the pipeline fits a covariate-adjusted
linear model and outputs the effect size (beta), standard error, and p-value.

Designed for UK Biobank Olink Explore data on SLURM-based HPC clusters
(tested on the Digital Research Alliance of Canada).

## Purpose

Given a list of target genes and their known cis-pQTL variants (from published
summary statistics such as Sun et al. 2023 / UKB-PPP), the pipeline asks:

> **Is this SNV associated with a change in plasma protein level, after
> adjusting for age, sex, ancestry PCs, genotyping batch, Olink plate batch,
> and processing time?**

For every variant-protein pair it reports beta, SE, and p-value from:

```
RINT(NPX) ~ dosage + age + age^2 + sex + age:sex + age^2:sex
           + assessment_centre + genotype_array + plate_batch
           + processing_time + PC1 + ... + PC20
```

## Pipeline steps

1. **Sample selection** -- filter participants by genetic ancestry and exclude
   disease cases
2. **Genotype QC** -- extract genotyped variants per chromosome, merge,
   apply MAF / missingness / HWE filters, LD-prune, remove related
   individuals (KING), compute principal components
3. **cis-pQTL identification** -- look up cis-variants for each target gene
   from external summary statistics (Snakemake checkpoint; the number of
   downstream jobs is determined at runtime)
4. **Dosage extraction** -- extract cis-variants from imputed data and export
   additive dosages
5. **Covariate construction** -- build per-protein batch and processing-time
   covariates
6. **Association testing** -- fit one linear model per variant-protein pair
   (rank-inverse normal transformed NPX as outcome, additive dosage as
   predictor, covariates as above)
7. **Summary** -- collect all per-pair results into a single CSV

Adding a new gene is a single YAML entry -- no code changes needed.

## Project structure

```
olink-pQTL/
  config/
    config.yaml.example     # template -- copy to config.yaml and fill paths
  envs/
    containers/
      r_pqtl.def            # Apptainer definition for R 4.4 + packages
      plink2.def            # Apptainer definition for PLINK2
  profiles/
    slurm/
      config.yaml           # SLURM executor settings (account, mem, time)
  workflow/
    Snakefile               # single entry point
    rules/
      qc.smk                # sample selection, genotype QC, LD pruning, PCA
      pqtl.smk              # cis-pQTL extraction, dosage export, association
    scripts/
      select_samples.R
      isolate_cis_pqtl.R
      create_batch_covariates.R
      preprocess_cache.R
      run_pqtl_association.R
      summarize_pqtl.R
```

## Requirements

- Snakemake >= 8.0
- PLINK2
- R >= 4.4 with packages: `data.table`, `dplyr`, `stringr`, `tidyr`

On Alliance Canada clusters:

```bash
module load plink/2.00a5.8
module load r/4.4.0
pip install snakemake
```

Apptainer container definitions are provided in `envs/containers/` as an
alternative to `module load`.

## Quick start

```bash
# 1. Clone
git clone https://github.com/yuan557/olink-pQTL.git
cd olink-pQTL

# 2. Configure
cp config/config.yaml.example config/config.yaml
# Edit config/config.yaml -- fill in your data paths and target genes

# 3. Dry run (resolves the full DAG including checkpoint)
snakemake -n

# 4. Run on SLURM
snakemake --profile profiles/slurm/
```

`config/config.yaml` is gitignored -- only the example template is committed.

## Configuration

All paths, QC parameters, and analysis targets live in `config/config.yaml`.

| Section | What it controls |
|-----------|-----------------|
| `paths` | Input file locations (genotypes, Olink NPX, covariates, summary stats) |
| `cohort` | Genetic ancestry filter and disease exclusion field |
| `qc` | MAF, missingness, HWE, LD pruning, KING cutoff, number of PCs |
| `targets` | List of genes to analyse |
| `pqtl` | RINT toggle, covariate list, birth year baseline |
| `resources` | PLINK2 thread and memory defaults |

## Adding a target gene

Append to `config/config.yaml`:

```yaml
targets:
  - gene: "CD40"
  - gene: "IL7R"
```

Re-run `snakemake --profile profiles/slurm/`. The checkpoint in `pqtl.smk`
discovers the new gene's cis-variants from the summary statistics and
schedules association jobs automatically.

## Output

The final output is `results/pqtl/pqtl_summary.csv` with one row per
variant-protein pair containing:

- Protein name and variant ID
- Effect size (beta) and standard error
- p-value
- Sample size

