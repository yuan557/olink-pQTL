# pQTL Pipeline

Snakemake workflow for cis-pQTL association analysis using UK Biobank Olink
proteomics and imputed genotype data. Runs RINT-transformed linear models
regressing NPX protein levels on genotype dosage, adjusting for standard
covariates (age, sex, PCs, batch, assessment centre, genotype array).

## Project structure

```
pqtl-pipeline/
  config/
    config.yaml.example   # template -- copy to config.yaml and fill paths
  envs/
    containers/
      r_pqtl.def          # Apptainer definition for R 4.4 + packages
      plink2.def           # Apptainer definition for PLINK2
  profiles/
    slurm/
      config.yaml          # SLURM executor settings (account, mem, time)
  workflow/
    Snakefile              # entry point
    rules/
      qc.smk              # sample selection, genotype QC, LD pruning, PCA
      pqtl.smk             # cis-pQTL extraction, dosage export, association
    scripts/
      select_samples.R
      isolate_cis_pqtl.R
      create_batch_covariates.R
      preprocess_cache.R
      run_pqtl_association.R
      summarize_pqtl.R
```

## Prerequisites

On an HPC cluster with SLURM, load modules directly:

```bash
module load plink/2.00a5.8
module load r/4.4.0
```

Or build Apptainer containers from the definitions in `envs/containers/`:

```bash
apptainer build envs/containers/plink2.sif envs/containers/plink2.def
apptainer build envs/containers/r_pqtl.sif envs/containers/r_pqtl.def
```

R packages required: `data.table`, `dplyr`, `stringr`, `tidyr`, `ggplot2`.

Install Snakemake (>= 8.0):

```bash
pip install snakemake
```

## Setup

1. Copy the config template and fill in your data paths:

```bash
cp config/config.yaml.example config/config.yaml
```

2. Edit `config/config.yaml` -- all input paths, QC thresholds, and target
   genes are defined there. To add a new gene, append an entry under `targets`.

3. `config/config.yaml` is gitignored. Only the example template is committed.

## Running

Dry run (check the DAG without executing):

```bash
snakemake -n
```

Run locally:

```bash
snakemake --cores 4
```

Run on SLURM:

```bash
snakemake --profile profiles/slurm/
```

To use Apptainer containers, make sure `use-singularity: true` is set in the
SLURM profile (it is by default), and that `.sif` files are built in
`envs/containers/`. Adjust `singularity-args` bind paths for your cluster.

To use `module load` instead, remove the `container:` directives from the rule
files or set `use-singularity: false` in the profile.

## Pipeline overview

1. **Sample selection** -- filter Olink participants by genetic ancestry and
   exclude disease cases.
2. **Genotype QC** -- extract per-chromosome genotyped variants for selected
   samples, merge, LD-prune, remove related individuals (KING), compute PCs.
3. **cis-pQTL identification** -- read external pQTL summary statistics, extract
   cis-variants for target genes from imputed data, merge, export dosage.
4. **Batch covariates** -- build per-protein plate batch and processing time
   covariates.
5. **Association** -- for each variant-protein pair, run a linear model:
   `RINT(NPX) ~ dosage + age + age^2 + sex + age:sex + age^2:sex + centre +
   array + batch + processing_time + PC1..PC20`. Output: beta, SE, p-value.
6. **Summary** -- collect all per-pair results into one table.

## Test data

This scaffold does not include real data. To test with public data:

1. Download 1000 Genomes phase 3 pgen files from
   https://www.cog-genomics.org/plink/2.0/resources
2. Simulate NPX values (random normal) for a subset of samples.
3. Create a mock `cis_pQTL.csv` manifest pointing to a few variants.
4. Point `config.yaml` paths at the simulated files.

## Adding a new target gene

Append to `config/config.yaml`:

```yaml
targets:
  - gene: "CD40"
  - gene: "IL7R"   # new
```

Then re-run `snakemake`. The checkpoint in `pqtl.smk` discovers the new gene's
cis-variants from the summary statistics and schedules association jobs
automatically.
