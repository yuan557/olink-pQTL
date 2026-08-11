checkpoint isolate_cis_pqtl:
    """Identify cis-pQTL variants from summary statistics for target genes.

    This is a checkpoint because the number of variant-protein pairs
    is discovered at runtime from the summary statistics files.
    """
    input:
        sumstats=config["paths"]["pqtl_sumstats"],
        clumped=config["paths"]["pqtl_clumped"],
    output:
        manifest=f"{OUTDIR}/pqtl/cis_pqtl_manifest.csv",
        chr_ids=directory(f"{OUTDIR}/pqtl/chr_ids"),
    params:
        genes=[t["gene"] for t in TARGETS],
    script:
        "../scripts/isolate_cis_pqtl.R"


def get_cis_chroms(wildcards):
    """Return chromosome ID files that actually have variants."""
    import glob
    ck = checkpoints.isolate_cis_pqtl.get(**wildcards)
    id_dir = ck.output.chr_ids
    return sorted(glob.glob(f"{id_dir}/chr*_IDs.txt"))


rule extract_cis_pqtl_chr:
    """Extract cis-pQTL variants from imputed data for one chromosome."""
    input:
        chr_ids=f"{OUTDIR}/pqtl/chr_ids/chr{{chr}}_IDs.txt",
        keep_ids=f"{OUTDIR}/qc/pca/unrelated.king.cutoff.in.id",
    output:
        pgen=f"{OUTDIR}/pqtl/extracted/chr{{chr}}.pgen",
        pvar=f"{OUTDIR}/pqtl/extracted/chr{{chr}}.pvar",
        psam=f"{OUTDIR}/pqtl/extracted/chr{{chr}}.psam",
    params:
        imputed_prefix=lambda wc: config["paths"]["imputed_pgen_prefix"].format(chr=wc.chr),
        out_prefix=lambda wc: f"{OUTDIR}/pqtl/extracted/chr{wc.chr}",
        threads=config["resources"]["plink2_threads"],
        memory=config["resources"]["plink2_memory_mb"],
    shell:
        """
        plink2 --threads {params.threads} --memory {params.memory} \
            --pfile {params.imputed_prefix} \
            --keep {input.keep_ids} \
            --extract {input.chr_ids} \
            --make-pgen \
            --out {params.out_prefix}
        """


rule merge_cis_pqtl:
    """Merge per-chromosome extracted pgens into one file."""
    input:
        chr_ids=get_cis_chroms,
    output:
        pgen=f"{OUTDIR}/pqtl/merged/all_cis.pgen",
        pvar=f"{OUTDIR}/pqtl/merged/all_cis.pvar",
        psam=f"{OUTDIR}/pqtl/merged/all_cis.psam",
        merge_list=temp(f"{OUTDIR}/pqtl/merged/merge_list.txt"),
    params:
        out_prefix=f"{OUTDIR}/pqtl/merged/all_cis",
        memory=config["resources"]["plink2_memory_mb"],
    run:
        import os, re
        with open(output.merge_list, "w") as f:
            for id_file in input.chr_ids:
                chr_num = re.search(r"chr(\d+)_IDs", id_file).group(1)
                prefix = os.path.join(OUTDIR, "pqtl", "extracted", f"chr{chr_num}")
                f.write(prefix + "\n")
        shell(
            "plink2 --memory {params.memory} "
            "--pmerge-list {output.merge_list} "
            "--make-pgen "
            "--out {params.out_prefix}"
        )


rule export_dosage:
    """Export merged pgen to additive dosage format (.raw)."""
    input:
        pgen=f"{OUTDIR}/pqtl/merged/all_cis.pgen",
        pvar=f"{OUTDIR}/pqtl/merged/all_cis.pvar",
        psam=f"{OUTDIR}/pqtl/merged/all_cis.psam",
    output:
        raw=f"{OUTDIR}/pqtl/merged/all_cis_dosage.raw",
    params:
        pfile_prefix=f"{OUTDIR}/pqtl/merged/all_cis",
        out_prefix=f"{OUTDIR}/pqtl/merged/all_cis_dosage",
        memory=config["resources"]["plink2_memory_mb"],
    shell:
        """
        plink2 --memory {params.memory} \
            --pfile {params.pfile_prefix} \
            --export A \
            --out {params.out_prefix}
        """


rule create_batch_covariates:
    """Build protein-level batch and processing-time covariates."""
    input:
        batch_main=config["paths"]["batch_main_df"],
        processing_dates=config["paths"]["batch_processing_dates"],
        panel_map=config["paths"]["panel_protein_map"],
        dosage=f"{OUTDIR}/pqtl/merged/all_cis_dosage.raw",
    output:
        batch_csv=f"{OUTDIR}/pqtl/batch_covariates.csv",
    script:
        "../scripts/create_batch_covariates.R"


rule preprocess_cache:
    """Load and cache all data frames needed by the association step."""
    input:
        manifest=f"{OUTDIR}/pqtl/cis_pqtl_manifest.csv",
        dosage=f"{OUTDIR}/pqtl/merged/all_cis_dosage.raw",
        olink=config["paths"]["olink_data"],
        covariates=config["paths"]["covariates"],
        eigenvec=f"{OUTDIR}/qc/pca/pca.eigenvec",
        batch=f"{OUTDIR}/pqtl/batch_covariates.csv",
    output:
        cis_pqtls_rds=f"{OUTDIR}/pqtl/cache/cis_pqtls.rds",
        geno_rds=f"{OUTDIR}/pqtl/cache/geno_slim.rds",
        base_rds=f"{OUTDIR}/pqtl/cache/df_base.rds",
        batch_rds=f"{OUTDIR}/pqtl/cache/df_batch.rds",
    params:
        birth_year_baseline=config["pqtl"]["birth_year_baseline"],
    script:
        "../scripts/preprocess_cache.R"


def get_pqtl_pairs(wildcards):
    """Return the manifest CSV so Snakemake can scatter over pairs."""
    ck = checkpoints.isolate_cis_pqtl.get()
    return ck.output.manifest


def list_pqtl_result_files(wildcards):
    """Expand over all variant-protein pairs discovered by the checkpoint."""
    import pandas as pd
    ck = checkpoints.isolate_cis_pqtl.get()
    manifest = pd.read_csv(ck.output.manifest)
    return [
        f"{OUTDIR}/pqtl/results/{row.protein}_{row.ID}_RINT.csv"
        for row in manifest.itertuples()
    ]


rule run_pqtl_association:
    """Run RINT + linear model for one variant-protein pair."""
    input:
        cis_pqtls_rds=f"{OUTDIR}/pqtl/cache/cis_pqtls.rds",
        geno_rds=f"{OUTDIR}/pqtl/cache/geno_slim.rds",
        base_rds=f"{OUTDIR}/pqtl/cache/df_base.rds",
        batch_rds=f"{OUTDIR}/pqtl/cache/df_batch.rds",
    output:
        result=f"{OUTDIR}/pqtl/results/{{protein}}_{{variant}}_RINT.csv",
    params:
        protein=lambda wc: wc.protein,
        variant=lambda wc: wc.variant,
        use_rint=config["pqtl"]["use_rint"],
        n_pcs=config["qc"]["n_pcs"],
    script:
        "../scripts/run_pqtl_association.R"


rule summarize_pqtl:
    """Collect all per-pair pQTL results into a single summary table."""
    input:
        results=list_pqtl_result_files,
    output:
        summary=f"{OUTDIR}/pqtl/pqtl_summary.csv",
    script:
        "../scripts/summarize_pqtl.R"
