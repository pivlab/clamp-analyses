R2_CFG = config["recount2"]
R2_DIR = R2_CFG["paths"]["production"]

# Builds the recount2 inputs used by the coverage comparator.


rule download_recount2:
    output:
        plier=R2_CFG["plier_rds"],
        rpkm=R2_CFG["rpkm_rds"],
    log:
        f"{R2_DIR}/download.log"
    params:
        url=R2_CFG["url"],
        out_dir=lambda wildcards, output: os.path.dirname(os.path.dirname(output.plier)),
    resources:
        mem_mb=8000,
        runtime=1440,
    shell:
        "mkdir -p {params.out_dir} && "
        "curl -L --fail -o {params.out_dir}/recount2.zip {params.url} && "
        "unzip -o {params.out_dir}/recount2.zip -d {params.out_dir} "
        "> {log} 2>&1"


rule preprocess_recount2:
    input:
        plier=rules.download_recount2.output.plier,
        rpkm=rules.download_recount2.output.rpkm,
        script="scripts/recount2/preprocess.R",
        common="scripts/recount2/common.R",
    output:
        fbm_filt=f"{R2_DIR}/recount2_fbm_filt.rds",
        genes=f"{R2_DIR}/recount2_genes.rds",
        samples=f"{R2_DIR}/recount2_samples.rds",
    log:
        f"{R2_DIR}/preprocess.log"
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.fbm_filt),
        mean_cutoff=R2_CFG["preprocess"]["mean_cutoff"],
        var_cutoff=R2_CFG["preprocess"]["var_cutoff"],
        block_size=R2_CFG["preprocess"]["block_size"],
        seed=R2_CFG["preprocess"]["seed"],
    resources:
        mem_mb=R2_CFG["resources"]["preprocess_mem_mb"],
        runtime=18000,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --plier-rds {input.plier} --rpkm-rds {input.rpkm} "
        "--out-dir {params.out_dir} --mean-cutoff {params.mean_cutoff} "
        "--var-cutoff {params.var_cutoff} --block-size {params.block_size} "
        "--seed {params.seed} > {log} 2>&1"


rule svd_recount2:
    input:
        fbm_filt=rules.preprocess_recount2.output.fbm_filt,
        genes=rules.preprocess_recount2.output.genes,
        samples=rules.preprocess_recount2.output.samples,
        script="scripts/recount2/svd.R",
        common="scripts/recount2/common.R",
    output:
        svd=f"{R2_DIR}/recount2_svdRes.rds",
        k=f"{R2_DIR}/CLAMP_K_recount2.rds",
    log:
        f"{R2_DIR}/svd.log"
    params:
        n_cores=R2_CFG["svd"]["n_cores"],
        seed=R2_CFG["svd"]["seed"],
    threads: R2_CFG["svd"]["n_cores"]
    resources:
        mem_mb=R2_CFG["resources"]["svd_mem_mb"],
        runtime=6000,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --fbm-rds {input.fbm_filt} --genes-rds {input.genes} "
        "--samples-rds {input.samples} --svd {output.svd} --k {output.k} "
        "--n-cores {params.n_cores} --seed {params.seed} > {log} 2>&1"


rule clampbase_recount2:
    input:
        fbm_filt=rules.preprocess_recount2.output.fbm_filt,
        genes=rules.preprocess_recount2.output.genes,
        samples=rules.preprocess_recount2.output.samples,
        svd=rules.svd_recount2.output.svd,
        k=rules.svd_recount2.output.k,
        script="scripts/recount2/clamp.R",
        common="scripts/recount2/common.R",
    output:
        rds=f"{R2_DIR}/CLAMPbase.rds",
        b=f"{R2_DIR}/CLAMPbase/B.csv",
        z=f"{R2_DIR}/CLAMPbase/Z.csv",
    log:
        f"{R2_DIR}/clampbase.log"
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.rds),
        seed=R2_CFG["clamp"]["seed"],
    resources:
        mem_mb=R2_CFG["resources"]["clampbase_mem_mb"],
        runtime=12000,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --fbm-rds {input.fbm_filt} --genes-rds {input.genes} "
        "--samples-rds {input.samples} --svd {input.svd} --k {input.k} "
        "--out-dir {params.out_dir} --seed {params.seed} > {log} 2>&1"


rule recount2_models:
    input:
        rules.clampbase_recount2.output.rds,
        rules.clampbase_recount2.output.z,


# Adopts existing recount2 outputs without recomputing them.
rule recount2_precomputed:
    input:
        rules.preprocess_recount2.output,
        rules.svd_recount2.output,
        rules.clampbase_recount2.output,
