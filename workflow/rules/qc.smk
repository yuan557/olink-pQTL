rule select_samples:
    """Select Olink participants by genetic ancestry, excluding disease cases."""
    input:
        olink=config["paths"]["olink_data"],
    output:
        sample_ids=f"{OUTDIR}/qc/sample_ids.txt",
    params:
        ancestry_field=config["cohort"]["genetic_ancestry_field"],
        ancestry_code=config["cohort"]["genetic_ancestry_code"],
        exclusion_field=config["cohort"]["exclusion_field"],
    script:
        "../scripts/select_samples.R"


rule select_genotype_chr:
    """Extract genotyped variants for selected samples, per chromosome."""
    input:
        sample_ids=f"{OUTDIR}/qc/sample_ids.txt",
    output:
        pgen=f"{OUTDIR}/qc/genotype/chr{{chr}}.pgen",
        pvar=f"{OUTDIR}/qc/genotype/chr{{chr}}.pvar",
        psam=f"{OUTDIR}/qc/genotype/chr{{chr}}.psam",
    params:
        bfile=lambda wc: config["paths"]["genotype_calls_prefix"].format(chr=wc.chr),
        out_prefix=lambda wc: f"{OUTDIR}/qc/genotype/chr{wc.chr}",
        threads=config["resources"]["plink2_threads"],
        memory=config["resources"]["plink2_memory_mb"],
    shell:
        """
        plink2 --threads {params.threads} --memory {params.memory} \
            --bfile {params.bfile} \
            --keep {input.sample_ids} \
            --make-pgen \
            --out {params.out_prefix}
        """


rule merge_chromosomes:
    """Merge all per-chromosome pgen files into one."""
    input:
        pgens=expand(f"{OUTDIR}/qc/genotype/chr{{chr}}.pgen", chr=CHROMS),
    output:
        pgen=f"{OUTDIR}/qc/genotype/allchr.pgen",
        pvar=f"{OUTDIR}/qc/genotype/allchr.pvar",
        psam=f"{OUTDIR}/qc/genotype/allchr.psam",
        merge_list=temp(f"{OUTDIR}/qc/genotype/merge_list.txt"),
    params:
        out_prefix=f"{OUTDIR}/qc/genotype/allchr",
        memory=config["resources"]["plink2_memory_mb"],
    run:
        import os
        with open(output.merge_list, "w") as f:
            for chr in CHROMS:
                f.write(os.path.join(OUTDIR, "qc", "genotype", f"chr{chr}") + "\n")
        shell(
            "plink2 --memory {params.memory} "
            "--pmerge-list {output.merge_list} "
            "--make-pgen "
            "--out {params.out_prefix}"
        )


rule ld_prune:
    """LD-prune genotyped variants for PCA."""
    input:
        pgen=f"{OUTDIR}/qc/genotype/allchr.pgen",
        pvar=f"{OUTDIR}/qc/genotype/allchr.pvar",
        psam=f"{OUTDIR}/qc/genotype/allchr.psam",
        high_ld=config["paths"]["high_ld_regions"],
    output:
        prune_in=f"{OUTDIR}/qc/pca/ld_prune.prune.in",
        pgen=f"{OUTDIR}/qc/pca/pruned.pgen",
        pvar=f"{OUTDIR}/qc/pca/pruned.pvar",
        psam=f"{OUTDIR}/qc/pca/pruned.psam",
    params:
        allchr_prefix=f"{OUTDIR}/qc/genotype/allchr",
        prune_prefix=f"{OUTDIR}/qc/pca/ld_prune",
        pruned_prefix=f"{OUTDIR}/qc/pca/pruned",
        maf=config["qc"]["maf"],
        geno=config["qc"]["geno_missing"],
        hwe=config["qc"]["hwe"],
        window=config["qc"]["ld_window"],
        step=config["qc"]["ld_step"],
        r2=config["qc"]["ld_r2"],
        threads=config["resources"]["plink2_threads"],
        memory=config["resources"]["plink2_memory_mb"],
    shell:
        """
        plink2 --threads {params.threads} --memory {params.memory} \
            --pfile {params.allchr_prefix} \
            --maf {params.maf} --geno {params.geno} --hwe {params.hwe} \
            --exclude range {input.high_ld} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --out {params.prune_prefix}

        plink2 --threads {params.threads} --memory {params.memory} \
            --pfile {params.allchr_prefix} \
            --extract {output.prune_in} \
            --make-pgen \
            --out {params.pruned_prefix}
        """


rule remove_related:
    """Remove related individuals using KING kinship estimator."""
    input:
        pgen=f"{OUTDIR}/qc/pca/pruned.pgen",
        pvar=f"{OUTDIR}/qc/pca/pruned.pvar",
        psam=f"{OUTDIR}/qc/pca/pruned.psam",
    output:
        keep_ids=f"{OUTDIR}/qc/pca/unrelated.king.cutoff.in.id",
        remove_ids=f"{OUTDIR}/qc/pca/unrelated.king.cutoff.out.id",
    params:
        pfile_prefix=f"{OUTDIR}/qc/pca/pruned",
        out_prefix=f"{OUTDIR}/qc/pca/unrelated",
        cutoff=config["qc"]["king_cutoff"],
        threads=config["resources"]["plink2_threads"],
        memory=config["resources"]["plink2_memory_mb"],
    shell:
        """
        plink2 --threads {params.threads} --memory {params.memory} \
            --pfile {params.pfile_prefix} \
            --king-cutoff {params.cutoff} \
            --out {params.out_prefix}
        """


rule compute_pca:
    """Compute principal components on unrelated, LD-pruned samples."""
    input:
        pgen=f"{OUTDIR}/qc/pca/pruned.pgen",
        pvar=f"{OUTDIR}/qc/pca/pruned.pvar",
        psam=f"{OUTDIR}/qc/pca/pruned.psam",
        remove_ids=f"{OUTDIR}/qc/pca/unrelated.king.cutoff.out.id",
    output:
        eigenvec=f"{OUTDIR}/qc/pca/pca.eigenvec",
        eigenval=f"{OUTDIR}/qc/pca/pca.eigenval",
    params:
        pfile_prefix=f"{OUTDIR}/qc/pca/pruned",
        out_prefix=f"{OUTDIR}/qc/pca/pca",
        n_pcs=config["qc"]["n_pcs"],
        threads=config["resources"]["plink2_threads"],
        memory=config["resources"]["plink2_memory_mb"],
    shell:
        """
        plink2 --threads {params.threads} --memory {params.memory} \
            --pfile {params.pfile_prefix} \
            --remove {input.remove_ids} \
            --pca {params.n_pcs} approx biallelic-var-wts \
            --out {params.out_prefix}
        """
