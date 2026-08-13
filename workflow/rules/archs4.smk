A4_CFG = config["archs4"]
A4_PROD = A4_CFG["paths"]["production"]
A4_BIO = A4_CFG["paths"]["biology"]
A4_MODEL_NB = os.path.join(REPO_ROOT, A4_CFG["paths"]["model_notebooks"])
A4_BIO_NB = os.path.join(REPO_ROOT, A4_CFG["paths"]["biology_notebooks"])

A4_PREP = f"{A4_PROD}/00_preprocess"
A4_FINAL = A4_CFG["final_model"]["dir"]
A4_SEEDS = [int(s) for s in A4_CFG["final_model"]["seeds"]]
A4_SEED_PATTERN = "|".join(str(s) for s in A4_SEEDS)
A4_LEAF = A4_CFG["final_model"]["leaf"]
A4_MODEL_NAME = A4_CFG["final_model"]["model_name"]

# The prior is several GMT collections merged into one sparse matrix; keep the
# collection order fixed so the two --prior-* arguments stay aligned.
A4_PRIOR = A4_CFG["prior"]["gmts"]
A4_PRIOR_NAMES = list(A4_PRIOR)
A4_PRIOR_PATHS = [A4_PRIOR[name] for name in A4_PRIOR_NAMES]


A4_REF_SEED = int(A4_CFG["final_model"]["reference_seed"])


def a4_seed_dir(seed):
    """Directory holding one seed of the shipped ARCHS4 model."""
    return f"{A4_FINAL}/{A4_LEAF.format(seed=seed)}"


# At 100% of the compendium there is no subsampling, so the input matrix, its SVD
# and the resulting CLAMPbase fit are identical for every seed -- verified on
# disk: all three seed directories hold byte-identical svd.rds and CLAMPbase/Z.csv.
# Only the CLAMPfull optimisation actually varies. Computing those two stages once
# and sharing them saves ~60 GB and two multi-day jobs on a fresh run; the seed
# directories keep their redundant copies for now.
A4_REF_DIR = a4_seed_dir(A4_REF_SEED)


# ============================================================
# ARCHS4 base pipeline: the filtered expression matrix and the full-compendium
# CLAMP model every downstream ARCHS4 analysis projects into.
#
# These rules are expensive enough that they are normally not run at all: the
# results already exist and are adopted with
#   snakemake --touch archs4_precomputed
# See the ARCHS4 section of the README.  They are written out in full so the
# pipeline remains reproducible from the raw HDF5 for anyone starting fresh.
# ============================================================

rule preprocess_archs4:
    input:
        h5=A4_CFG["raw_h5"],
        gene_lengths=A4_CFG["gene_lengths"],
        script="scripts/archs4/preprocess.R",
        common="scripts/archs4/common.R",
    output:
        metadata=f"{A4_PREP}/metadata.rds",
        metadata_filtered=f"{A4_PREP}/metadata_filtered.rds",
        all_samples=f"{A4_PREP}/all_samples.rds",
        fbm=f"{A4_PREP}/fbm.bk",
        fbm_filtered=f"{A4_PREP}/fbm_filtered.bk",
    log:
        f"{A4_PREP}/preprocess.log"
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.metadata),
        mean_cutoff=A4_CFG["preprocess"]["mean_cutoff"],
        var_cutoff=A4_CFG["preprocess"]["var_cutoff"],
        sc_max=A4_CFG["preprocess"]["single_cell_probability_max"],
        block_size=A4_CFG["raw_data"]["block_size"],
        n_cores=A4_CFG["raw_data"]["n_cores"],
        seed=A4_CFG["preprocess"]["seed"],
    threads: A4_CFG["raw_data"]["n_cores"]
    resources:
        mem_mb=200000,
        runtime=18000,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --h5 {input.h5} --gene-lengths {input.gene_lengths} "
        "--out-dir {params.out_dir} --mean-cutoff {params.mean_cutoff} "
        "--var-cutoff {params.var_cutoff} "
        "--single-cell-probability-max {params.sc_max} "
        "--block-size {params.block_size} --n-cores {params.n_cores} --seed {params.seed} "
        "> {log} 2>&1"


rule svd_archs4:
    input:
        metadata=rules.preprocess_archs4.output.metadata_filtered,
        fbm=rules.preprocess_archs4.output.fbm_filtered,
        script="scripts/archs4/svd.R",
        common="scripts/archs4/common.R",
    output:
        # svd_full.rds (the pre-NaN-filter object) is deliberately not declared:
        # it is a 22 GB debugging artifact that was not retained for the shipped
        # model, and declaring it would make --touch fail on a file we do not have.
        svd=f"{A4_REF_DIR}/svd.rds",
        k=f"{A4_REF_DIR}/CLAMP_K.rds",
    log:
        f"{A4_REF_DIR}/svd.log"
    params:
        # bigstatsr addresses an FBM by its backing file without the extension.
        fbm_prefix=lambda wildcards, input: input.fbm[: -len(".bk")],
        n_cores=A4_CFG["svd"]["n_cores"],
        k_multiplier=A4_CFG["clamp"]["k_multiplier"],
        seed=A4_CFG["svd"]["seed"],
    threads: A4_CFG["svd"]["n_cores"]
    resources:
        mem_mb=200000,
        runtime=6000,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --metadata {input.metadata} --fbm {params.fbm_prefix} "
        "--svd {output.svd} --k {output.k} "
        "--n-cores {params.n_cores} --k-multiplier {params.k_multiplier} --seed {params.seed} "
        "> {log} 2>&1"


rule clampbase_archs4:
    input:
        metadata=rules.preprocess_archs4.output.metadata_filtered,
        samples=rules.preprocess_archs4.output.all_samples,
        fbm=rules.preprocess_archs4.output.fbm_filtered,
        svd=rules.svd_archs4.output.svd,
        k=rules.svd_archs4.output.k,
        script="scripts/archs4/clamp.R",
        common="scripts/archs4/common.R",
    output:
        rds=f"{A4_REF_DIR}/CLAMPbase.rds",
        b=f"{A4_REF_DIR}/CLAMPbase/B.csv",
        z=f"{A4_REF_DIR}/CLAMPbase/Z.csv",
    log:
        f"{A4_REF_DIR}/clampbase.log"
    params:
        fbm_prefix=lambda wildcards, input: input.fbm[: -len(".bk")],
        out_dir=A4_REF_DIR,
        seed=A4_CFG["clamp"]["seed"],
    resources:
        mem_mb=500000,
        runtime=12000,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --stage base --metadata {input.metadata} "
        "--samples {input.samples} --fbm {params.fbm_prefix} --svd {input.svd} "
        "--k {input.k} --out-dir {params.out_dir} --seed {params.seed} "
        "> {log} 2>&1"


rule clampfull_archs4:
    input:
        metadata=rules.preprocess_archs4.output.metadata_filtered,
        samples=rules.preprocess_archs4.output.all_samples,
        fbm=rules.preprocess_archs4.output.fbm_filtered,
        svd=rules.svd_archs4.output.svd,
        k=rules.svd_archs4.output.k,
        base=rules.clampbase_archs4.output.rds,
        priors=A4_PRIOR_PATHS,
        script="scripts/archs4/clamp.R",
        common="scripts/archs4/common.R",
    output:
        rds=f"{A4_FINAL}/{A4_LEAF}/{A4_MODEL_NAME}.rds",
        b=f"{A4_FINAL}/{A4_LEAF}/{A4_MODEL_NAME}/B.csv",
        z=f"{A4_FINAL}/{A4_LEAF}/{A4_MODEL_NAME}/Z.csv",
        summary=f"{A4_FINAL}/{A4_LEAF}/{A4_MODEL_NAME}/summary.csv",
    log:
        f"{A4_FINAL}/{A4_LEAF}/clampfull.log"
    params:
        fbm_prefix=lambda wildcards, input: input.fbm[: -len(".bk")],
        out_dir=lambda wildcards: a4_seed_dir(wildcards.seed),
        prior_names=",".join(A4_PRIOR_NAMES),
        prior_gmts=",".join(A4_PRIOR_PATHS),
        multiplier=A4_CFG["clamp"]["multiplier"],
        max_iter=A4_CFG["clamp"]["max_iter"],
        # CLAMPfull is the only stage that varies between the three model seeds,
        # so this is the one place the run index has to reach the RNG.
        seed=lambda wildcards: A4_CFG["clamp"]["seed"] + int(wildcards.seed),
    wildcard_constraints:
        seed=A4_SEED_PATTERN,
    resources:
        mem_mb=500000,
        runtime=18000,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --stage full --metadata {input.metadata} "
        "--samples {input.samples} --fbm {params.fbm_prefix} --svd {input.svd} "
        "--k {input.k} --base-model {input.base} --out-dir {params.out_dir} "
        "--prior-names {params.prior_names} --prior-gmts {params.prior_gmts} "
        "--multiplier {params.multiplier} --max-iter {params.max_iter} --seed {params.seed} "
        "> {log} 2>&1"


rule archs4_models:
    input:
        expand(f"{A4_FINAL}/{A4_LEAF}/{A4_MODEL_NAME}/Z.csv", seed=A4_SEEDS),
        expand(f"{A4_FINAL}/{A4_LEAF}/{A4_MODEL_NAME}.rds", seed=A4_SEEDS),
        f"{A4_REF_DIR}/CLAMPbase.rds",


rule model_building_qc_archs4:
    input:
        metadata=rules.preprocess_archs4.output.metadata_filtered,
        summaries=expand(
            f"{A4_FINAL}/{A4_LEAF}/{A4_MODEL_NAME}/summary.csv", seed=A4_SEEDS
        ),
        notebook=f"{A4_MODEL_NB}/00_model_building_qc.ipynb",
    output:
        matrix_qc=f"{A4_PROD}/qc/matrix_qc.csv",
        seed_qc=f"{A4_PROD}/qc/seed_qc.csv",
        seed_overlap=f"{A4_PROD}/qc/seed_overlap.csv",
        complete=touch(f"{A4_PROD}/qc/notebook.complete"),
    log:
        notebook=f"{A4_MODEL_NB}/00_model_building_qc.executed.ipynb",
    params:
        seeds=A4_SEEDS,
        fdr=A4_CFG["ora"]["fdr"][0],
        out_dir=f"{A4_PROD}/qc",
    conda: "clamp-analyses"
    notebook:
        f"{A4_MODEL_NB}/00_model_building_qc.ipynb"


# ============================================================
# Adoption handle for results that already exist.
#
# `snakemake --touch archs4_precomputed` marks the migrated outputs as current
# so the expensive rules above never fire.  --touch fails loudly on a missing
# file, which makes this double as a check that scripts/archs4/migrate_outputs.sh
# put everything where the rules expect it.  A fresh user simply skips this
# target and computes from scratch.
# ============================================================

rule archs4_precomputed:
    input:
        rules.preprocess_archs4.output,
        rules.archs4_models.input,
