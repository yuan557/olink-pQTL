# olink-pQTL

Snakemake workflow that runs covariate-adjusted linear models on genotype
dosages and Olink NPX values to produce cis-pQTL summary statistics
(beta, SE, p-value).

Designed for UK Biobank Olink Explore data on SLURM-based HPC clusters
(tested on the Digital Research Alliance of Canada).

## What it does

Given a set of target genes, the pipeline:

1. Selects participants by genetic ancestry and excludes disease cases
2. Extracts and QCs genotyped variants (per chromosome), merges, LD-prunes,
   removes related individuals (KING), and computes principal components
3. Identifies cis-pQTL variants for each target gene from external summary
   statistics (e.g. Sun et al. 2023, UKB-PPP)
4. Extracts those variants from imputed data and exports additive dosages
5. Builds per-protein batch and processing-time covariates
6. Fits a linear model per variant-protein pair:

```
RINT(NPX) ~ dosage + age + age^2 + sex + age:sex + age^2:sex
           + assessment_centre + genotype_array + plate_batch
           + processing_time + PC1 + ... + PC20
```

7. Collects all results into a single summary table

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
      select_samples.R      # filter participants by ancestry / exclusion
      isolate_cis_pqtl.R    # identify cis-variants from summary stats
      create_batch_covariates.R
      preprocess_cache.R    # merge and cache data for association jobs
      run_pqtl_association.R  # RINT + linear model per variant-protein pair
      summarize_pqtl.R      # collect results into one table
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

# 3. Dry run
snakemake -n

# 4. Run on SLURM
snakemake --profile profiles/slurm/
```

`config/config.yaml` is gitignored -- only the example template is committed.

## Configuration

All paths, QC parameters, and analysis targets are in `config/config.yaml`.

Key sections:

| Section | What it controls |
|---------|-----------------|
| `paths` | Input file locations (genotypes, Olink data, covariates, summary stats) |
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
  - gene: "IL7R"   # new
```

Re-run `snakemake --profile profiles/slurm/`. The checkpoint in `pqtl.smk`
discovers the new gene's cis-variants from the summary statistics and
schedules association jobs automatically.

## Test data

This repository does not include real data. To test with public data:

1. Download 1000 Genomes phase 3 pgen files from
   https://www.cog-genomics.org/plink/2.0/resources
2. Simulate NPX values (random normal) for a subset of samples
3. Create a mock `cis_pQTL.csv` manifest pointing to a few variants
4. Point `config.yaml` paths at the simulated files
