import os

# ============================================================
# LV-trait association (GLS) via pivlab/phenoplier-cli, run against the CLAMP
# models the ARCHS4 coverage, saturation, final-model and canonical rules of
# this workflow produce.
#
# This module is the UPSTREAM of archs4_traits.smk: every rule here writes one
# per-model GLS summary into the archs4.yaml `traits.*` directory that
# aggregate_archs4_{coverage,saturation}_traits scan, named the way
# scripts/archs4/traits/aggregate_*_traits.R expect (cov_rs<f>_seed<s>,
# sat_rs<f>_k<k>_seed<s>, final_<dataset>; the directory, not the filename,
# says whether a file is CLAMPfull_bp or CLAMPbase).
#
# Every model path comes from the helpers/config of the rule file that builds
# that model -- archs4_coverage.smk (a4_cov_*), archs4_saturation.smk
# (a4_sat_*), archs4_canonical.smk (A4_CAN_*), archs4.yaml: final_models --
# so a model is fitted, validated and published before GLS runs on it, and
# nothing here names a model location of its own. Where the fitted model is a
# directory() output (CLAMPfull_bp), the rule depends on that model's
# validated.json exactly like the ORA rules do.
#
# phenoplier-cli ships its own Snakemake pipeline on snakemake>=9 /
# Python>=3.12, incompatible with envs/snakemake.yaml (snakemake=8 / 3.11), so
# it runs in its own conda env (phenoplier.yaml: conda_env) through
# scripts/phenoplier/run_gls.sh, one workspace project per model, rather than
# as rules of this DAG. run_gls.sh is the hardened path every model in
# scripts/phenoplier/RUN_SUMMARY.md went through: register + init, pin the
# step-6/7 pools, run (resume-safe), gate on a complete summary, publish.
#
# phenoplier-cli is pinned (phenoplier.yaml: version, 0.5.2 -- the release
# the as-run models used): every job checks the installed version and fails
# on a mismatch, so results never silently mix releases. The GLS
# standard-error collapse that emitted artifactual ~0 p-values for a handful
# of pathological LVs in the sub-sampling models is fixed in that release
# (pivlab/phenoplier-cli#85); the completeness gate still refuses any summary
# with a zero p-value.
# ============================================================

PHENOPLIER_CFG = config["phenoplier"]
PHENOPLIER_TARGET = config.get("phenoplier_target", PHENOPLIER_CFG["target"])
PHENOPLIER_CLUSTER = PHENOPLIER_CFG["clusters"][PHENOPLIER_TARGET]
PHENOPLIER_WORKERS = PHENOPLIER_CFG["workers"]
PHENOPLIER_GLS_RES = PHENOPLIER_CFG["resources"]["gls"]
PHENOPLIER_STORE_RES = PHENOPLIER_CFG["resources"]["store"]
PHENOPLIER_REPORT = PHENOPLIER_CFG["report_root"]
PHENOPLIER_TRAIT_FILTER = PHENOPLIER_CFG["trait_filter"]
PHENOPLIER_EXPECTED_PHENOTYPES = int(
    PHENOPLIER_CFG["expected_phenotypes"][PHENOPLIER_TRAIT_FILTER]
)

# Where the summaries land: the directories archs4_traits.smk consumes.
A4_TRAIT_DIRS = A4_CFG["traits"]
PHENOPLIER_COV_FULL_DIR = A4_TRAIT_DIRS["coverage"]["clampfull_bp_dir"]
PHENOPLIER_COV_BASE_DIR = A4_TRAIT_DIRS["coverage"]["clampbase_dir"]
PHENOPLIER_SAT_FULL_DIR = A4_TRAIT_DIRS["saturation"]["clampfull_bp_dir"]
PHENOPLIER_SAT_BASE_DIR = A4_TRAIT_DIRS["saturation"]["clampbase_dir"]
PHENOPLIER_FIN_FULL_DIR = A4_TRAIT_DIRS["finals"]["clampfull_bp_dir"]
PHENOPLIER_FIN_BASE_DIR = A4_TRAIT_DIRS["finals"]["clampbase_dir"]

# Saturation: only the ranks the trait-recovery report uses.
PHENOPLIER_SAT_KS = [int(k) for k in PHENOPLIER_CFG["saturation_k_values"]]
PHENOPLIER_SAT_K_PATTERN = "|".join(map(str, PHENOPLIER_SAT_KS))
PHENOPLIER_SAT_CELLS = [
    (fraction, k, seed)
    for (fraction, k, seed) in A4_SAT_CELLS
    if k in PHENOPLIER_SAT_KS
]

# Final models: the published full-data fits. CLAMPfull_bp is what
# publish_bp_model links into archs4.yaml: final_models.clampfull.bp from the
# coverage rs100 / reference-seed fit of each dataset. CLAMPbase comes from
# archs4.yaml: coverage.clampbase_source, i.e. the rule-produced CLAMPbase of
# each compendium (clamp_gtex, clampbase_recount2, or the archs4 coverage cell
# named by coverage.clampbase_seed).
#
# For archs4 both finals ARE coverage cells the rules above already run GLS
# on (rs100 / seed 1: publish_bp_model links that exact CLAMPfull_bp fit, and
# at 100% of the compendium CLAMPbase does not depend on the subsample seed).
# Their final summaries are therefore copies of the coverage ones, not a
# second ~20 h GLS run. GTEx / recount2 finals come from the comparator fits,
# which the coverage trait report never runs GLS on, so those run here.
A4_FINAL_MODELS = A4_CFG["final_models"]
A4_FINAL_DATASETS = list(A4_COV_BP_DATASETS)
A4_FINAL_DATASET_PATTERN = "|".join(A4_FINAL_DATASETS)
A4_FINAL_LINKED = ["archs4"]
A4_FINAL_LINKED_PATTERN = "|".join(A4_FINAL_LINKED)
A4_FINAL_RUN_DATASETS = [d for d in A4_FINAL_DATASETS if d not in A4_FINAL_LINKED]
A4_FINAL_RUN_PATTERN = "|".join(A4_FINAL_RUN_DATASETS)


def a4_final_bp_rds(dataset):
    spec = A4_FINAL_MODELS["clampfull"]["bp"][dataset]
    return f"{spec['root']}/{spec['rds']}"


def a4_final_base_seed(dataset):
    return int(A4_COV_CFG["clampbase_seed"][dataset])


def a4_final_base_rds(dataset):
    if dataset in A4_FINAL_LINKED:
        return f"{a4_cov_cell(100, a4_final_base_seed(dataset))}/CLAMPbase.rds"
    return A4_COV_CFG["clampbase_source"][dataset]["model"]


def phenoplier_final_name(dataset, model):
    """Workspace project that holds a final model's GLS run."""
    if dataset in A4_FINAL_LINKED:
        seed = A4_COV_REFERENCE_SEED if model == A4_COV_CFG["model_name"] else a4_final_base_seed(dataset)
        return f"cov_rs100_seed{seed}_{model}"
    return f"final_{dataset}_{model}"


def phenoplier_store_dir(summary_dir):
    """Stores sit beside their summaries: <root>/summaries -> <root>/stores."""
    return os.path.join(os.path.dirname(summary_dir), "stores")


# Every summary each family produces. These lists are the inputs of
# archs4_traits.smk's aggregate rules and of the convenience targets below.
PHENOPLIER_COV_FULL = [
    f"{PHENOPLIER_COV_FULL_DIR}/cov_rs{fraction}_seed{seed}.tsv.gz"
    for fraction in A4_COV_LEVELS
    for seed in a4_cov_seeds_for(fraction)
]
PHENOPLIER_COV_BASE = [
    f"{PHENOPLIER_COV_BASE_DIR}/cov_rs{fraction}_seed{seed}.tsv.gz"
    for fraction in A4_COV_LEVELS
    for seed in a4_cov_seeds_for(fraction)
]
PHENOPLIER_SAT_FULL = [
    f"{PHENOPLIER_SAT_FULL_DIR}/sat_rs{fraction}_k{k}_seed{seed}.tsv.gz"
    for (fraction, k, seed) in PHENOPLIER_SAT_CELLS
]
PHENOPLIER_SAT_BASE = [
    f"{PHENOPLIER_SAT_BASE_DIR}/sat_rs{fraction}_k{k}_seed{seed}.tsv.gz"
    for (fraction, k, seed) in PHENOPLIER_SAT_CELLS
]
PHENOPLIER_FIN_FULL = [
    f"{PHENOPLIER_FIN_FULL_DIR}/final_{dataset}.tsv.gz" for dataset in A4_FINAL_DATASETS
]
PHENOPLIER_FIN_BASE = [
    f"{PHENOPLIER_FIN_BASE_DIR}/final_{dataset}.tsv.gz" for dataset in A4_FINAL_DATASETS
]
PHENOPLIER_CANONICAL = [
    f"{A4_CAN_FINAL_ROOT}/{dataset}/traits/canon_{dataset}.tsv.gz"
    for dataset in A4_CAN_COMPENDIA
]
PHENOPLIER_FIN_STORES = [
    f"{phenoplier_store_dir(PHENOPLIER_FIN_FULL_DIR)}/final_{dataset}_{A4_COV_CFG['model_name']}.h5"
    for dataset in A4_FINAL_DATASETS
] + [
    f"{phenoplier_store_dir(PHENOPLIER_FIN_BASE_DIR)}/final_{dataset}_CLAMPbase.h5"
    for dataset in A4_FINAL_DATASETS
]

# One command for every GLS rule: the rules differ only in which model they
# point run_gls.sh at (params.rds / params.name) and where the summary goes.
# Everything else is workflow-wide config, baked in here once.
PHENOPLIER_GLS_SHELL = "".join([
    "bash {input.script} --rds {params.rds} --name {params.name} ",
    "--summary-out {output.summary} --n-jobs {threads} ",
    f"--conda-env {PHENOPLIER_CFG['conda_env']} ",
    f"--require-version {PHENOPLIER_CFG['version']} ",
    f"--namespace {PHENOPLIER_CFG['namespace']} ",
    f"--cohort {PHENOPLIER_CFG['cohort']} ",
    f"--lv-percentile {PHENOPLIER_CFG['lv_percentile']} ",
    f"--trait-filter {PHENOPLIER_TRAIT_FILTER} ",
    f"--expected-phenotypes {PHENOPLIER_EXPECTED_PHENOTYPES} ",
    f"--executor {PHENOPLIER_CLUSTER['executor']} ",
    f"--cluster '{PHENOPLIER_CLUSTER.get('cluster', '')}' ",
    f"--workers-step6 {PHENOPLIER_WORKERS['step6']} ",
    f"--blas-step6 {PHENOPLIER_WORKERS['blas_step6']} ",
    f"--workers-step7 {PHENOPLIER_WORKERS['step7']} ",
    f"--blas-step7 {PHENOPLIER_WORKERS['blas_step7']} ",
    "> {log} 2>&1",
])

PHENOPLIER_STORE_SHELL = "".join([
    "bash {input.script} --rds {params.rds} --name {params.name} ",
    "--store-out {output.store} ",
    f"--conda-env {PHENOPLIER_CFG['conda_env']} ",
    f"--require-version {PHENOPLIER_CFG['version']} ",
    f"--cohort {PHENOPLIER_CFG['cohort']} ",
    "> {log} 2>&1",
])


# ============================================================
# Step 1: one GLS run per model. Project names match the as-run pico
# workspace (scripts/phenoplier/*/), so a workspace that already holds a
# finished project is reused rather than recomputed.
# ============================================================

rule gls_bp_coverage_full:
    """GLS on one ARCHS4 coverage CLAMPfull_bp fit (native rank)."""
    input:
        validated=lambda wc: a4_cov_validated("archs4", wc.fraction, wc.seed),
        script="scripts/phenoplier/run_gls.sh",
    output:
        summary=f"{PHENOPLIER_COV_FULL_DIR}/cov_rs{{fraction}}_seed{{seed}}.tsv.gz",
    log:
        f"{PHENOPLIER_COV_FULL_DIR}/cov_rs{{fraction}}_seed{{seed}}.log"
    params:
        rds=lambda wc: (
            f"{a4_cov_model_dir('archs4', wc.fraction, wc.seed)}/{A4_COV_CFG['model_name']}.rds"
        ),
        name=lambda wc: f"cov_rs{wc.fraction}_seed{wc.seed}_{A4_COV_CFG['model_name']}",
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        fraction=A4_COV_LEVEL_PATTERN,
        seed=A4_COV_SEED_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


rule gls_bp_coverage_base:
    """GLS on the CLAMPbase of one ARCHS4 coverage cell."""
    input:
        rds=lambda wc: f"{a4_cov_cell(wc.fraction, wc.seed)}/CLAMPbase.rds",
        script="scripts/phenoplier/run_gls.sh",
    output:
        summary=f"{PHENOPLIER_COV_BASE_DIR}/cov_rs{{fraction}}_seed{{seed}}.tsv.gz",
    log:
        f"{PHENOPLIER_COV_BASE_DIR}/cov_rs{{fraction}}_seed{{seed}}.log"
    params:
        rds=lambda wc, input: input.rds,
        name=lambda wc: f"cov_rs{wc.fraction}_seed{wc.seed}_CLAMPbase",
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        fraction=A4_COV_LEVEL_PATTERN,
        seed=A4_COV_SEED_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


rule gls_bp_saturation_full:
    """GLS on one ARCHS4 saturation CLAMPfull_bp fit (forced K)."""
    input:
        validated=lambda wc: a4_sat_validated(wc.fraction, wc.k, wc.seed),
        script="scripts/phenoplier/run_gls.sh",
    output:
        summary=f"{PHENOPLIER_SAT_FULL_DIR}/sat_rs{{fraction}}_k{{k}}_seed{{seed}}.tsv.gz",
    log:
        f"{PHENOPLIER_SAT_FULL_DIR}/sat_rs{{fraction}}_k{{k}}_seed{{seed}}.log"
    params:
        rds=lambda wc: (
            f"{a4_sat_model_dir(wc.fraction, wc.k, wc.seed)}/{A4_SAT_CFG['model_name']}.rds"
        ),
        name=lambda wc: (
            f"sat_rs{wc.fraction}_k{wc.k}_seed{wc.seed}_{A4_SAT_CFG['model_name']}"
        ),
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        fraction=A4_SAT_FRACTION_PATTERN,
        k=PHENOPLIER_SAT_K_PATTERN,
        seed=A4_SAT_SEED_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


rule gls_bp_saturation_base:
    """GLS on one ARCHS4 saturation CLAMPbase (forced K)."""
    input:
        # ancient(): like the saturation ORA rules, an existing CLAMPbase from
        # the earlier sweep must not be refit just because it is old.
        rds=lambda wc: ancient(a4_sat_base_rds(wc.fraction, wc.k, wc.seed)),
        script="scripts/phenoplier/run_gls.sh",
    output:
        summary=f"{PHENOPLIER_SAT_BASE_DIR}/sat_rs{{fraction}}_k{{k}}_seed{{seed}}.tsv.gz",
    log:
        f"{PHENOPLIER_SAT_BASE_DIR}/sat_rs{{fraction}}_k{{k}}_seed{{seed}}.log"
    params:
        rds=lambda wc, input: input.rds,
        name=lambda wc: f"sat_rs{wc.fraction}_k{wc.k}_seed{wc.seed}_CLAMPbase",
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        fraction=A4_SAT_FRACTION_PATTERN,
        k=PHENOPLIER_SAT_K_PATTERN,
        seed=A4_SAT_SEED_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


rule gls_bp_final_full:
    """GLS on the published full-data CLAMPfull_bp of a comparator compendium."""
    input:
        manifest=rules.publish_bp_model.output.manifest,
        script="scripts/phenoplier/run_gls.sh",
    output:
        summary=f"{PHENOPLIER_FIN_FULL_DIR}/final_{{dataset}}.tsv.gz",
    log:
        f"{PHENOPLIER_FIN_FULL_DIR}/final_{{dataset}}.log"
    params:
        rds=lambda wc: a4_final_bp_rds(wc.dataset),
        name=lambda wc: phenoplier_final_name(wc.dataset, A4_COV_CFG["model_name"]),
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        dataset=A4_FINAL_RUN_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


rule gls_final_base:
    """GLS on the full-data CLAMPbase of a comparator compendium."""
    input:
        rds=lambda wc: a4_final_base_rds(wc.dataset),
        script="scripts/phenoplier/run_gls.sh",
    output:
        summary=f"{PHENOPLIER_FIN_BASE_DIR}/final_{{dataset}}.tsv.gz",
    log:
        f"{PHENOPLIER_FIN_BASE_DIR}/final_{{dataset}}.log"
    params:
        rds=lambda wc, input: input.rds,
        name=lambda wc: phenoplier_final_name(wc.dataset, "CLAMPbase"),
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        dataset=A4_FINAL_RUN_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


rule link_bp_final_full:
    """The archs4 final CLAMPfull_bp is the coverage rs100 / reference-seed
    fit publish_bp_model links; reuse that cell's GLS summary."""
    input:
        summary=f"{PHENOPLIER_COV_FULL_DIR}/cov_rs100_seed{A4_COV_REFERENCE_SEED}.tsv.gz",
        manifest=rules.publish_bp_model.output.manifest,
    output:
        summary=f"{PHENOPLIER_FIN_FULL_DIR}/final_{{dataset}}.tsv.gz",
    wildcard_constraints:
        dataset=A4_FINAL_LINKED_PATTERN,
    shell:
        "cp -f {input.summary} {output.summary}"


rule link_final_base:
    """The archs4 final CLAMPbase is the coverage rs100 CLAMPbase (seed from
    archs4.yaml: coverage.clampbase_seed); reuse that cell's GLS summary."""
    input:
        summary=lambda wc: (
            f"{PHENOPLIER_COV_BASE_DIR}/cov_rs100_seed{a4_final_base_seed(wc.dataset)}.tsv.gz"
        ),
    output:
        summary=f"{PHENOPLIER_FIN_BASE_DIR}/final_{{dataset}}.tsv.gz",
    wildcard_constraints:
        dataset=A4_FINAL_LINKED_PATTERN,
    shell:
        "cp -f {input.summary} {output.summary}"


rule gls_canonical:
    """GLS on the published canonical-prior CLAMPfull of one compendium.

    Written beside the model like its ORA results (archs4_canonical.smk).
    """
    input:
        rds=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/{A4_CAN_MODEL_NAME}.rds",
        script="scripts/phenoplier/run_gls.sh",
    output:
        summary=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/traits/canon_{{dataset}}.tsv.gz",
    log:
        f"{A4_CAN_FINAL_ROOT}/{{dataset}}/traits/canon_{{dataset}}.log"
    params:
        rds=lambda wc, input: input.rds,
        name=lambda wc: f"canon_{wc.dataset}_{A4_CAN_MODEL_NAME}",
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        dataset=A4_CAN_COMPENDIUM_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


# ============================================================
# Step 2: one-file HDF5 composite study store per final model: the CLAMP .rds
# (Z, hyperparameters, provenance) plus that model's GLS results, via
# `phenoplier store build`. Depends on the model's published summary so the
# run has finished and passed the completeness gate.
# ============================================================

rule store_bp_final_full:
    input:
        summary=f"{PHENOPLIER_FIN_FULL_DIR}/final_{{dataset}}.tsv.gz",
        script="scripts/phenoplier/store_build.sh",
    output:
        store=f"{phenoplier_store_dir(PHENOPLIER_FIN_FULL_DIR)}/final_{{dataset}}_{A4_COV_CFG['model_name']}.h5",
    log:
        f"{phenoplier_store_dir(PHENOPLIER_FIN_FULL_DIR)}/final_{{dataset}}_{A4_COV_CFG['model_name']}.log"
    params:
        rds=lambda wc: a4_final_bp_rds(wc.dataset),
        name=lambda wc: phenoplier_final_name(wc.dataset, A4_COV_CFG["model_name"]),
    resources:
        mem_mb=int(PHENOPLIER_STORE_RES["mem_mb"]),
        runtime=int(PHENOPLIER_STORE_RES["runtime"]),
    wildcard_constraints:
        dataset=A4_FINAL_DATASET_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_STORE_SHELL


rule store_final_base:
    input:
        summary=f"{PHENOPLIER_FIN_BASE_DIR}/final_{{dataset}}.tsv.gz",
        script="scripts/phenoplier/store_build.sh",
    output:
        store=f"{phenoplier_store_dir(PHENOPLIER_FIN_BASE_DIR)}/final_{{dataset}}_CLAMPbase.h5",
    log:
        f"{phenoplier_store_dir(PHENOPLIER_FIN_BASE_DIR)}/final_{{dataset}}_CLAMPbase.log"
    params:
        rds=lambda wc: a4_final_base_rds(wc.dataset),
        name=lambda wc: phenoplier_final_name(wc.dataset, "CLAMPbase"),
    resources:
        mem_mb=int(PHENOPLIER_STORE_RES["mem_mb"]),
        runtime=int(PHENOPLIER_STORE_RES["runtime"]),
    wildcard_constraints:
        dataset=A4_FINAL_DATASET_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_STORE_SHELL


# ============================================================
# Step 3: cross-model long table (every LV x trait row of every summary).
# The per-family trait-RECOVERY tables and reports are archs4_traits.smk's.
# A summary is ~4 M rows (1,728 LVs x 2,366 traits), so the table is streamed
# to a gzipped CSV and is tens of GB for the full set; it is its own target
# (phenoplier_long_table), not part of phenoplier_traits.
# ============================================================

rule aggregate_phenoplier_traits:
    input:
        clampfull_bp=PHENOPLIER_COV_FULL + PHENOPLIER_SAT_FULL + PHENOPLIER_FIN_FULL,
        clampbase=PHENOPLIER_COV_BASE + PHENOPLIER_SAT_BASE + PHENOPLIER_FIN_BASE,
        canonical=PHENOPLIER_CANONICAL,
        script="scripts/phenoplier/aggregate_traits.py",
        wrapper="scripts/phenoplier/aggregate_traits.sh",
    output:
        long=f"{PHENOPLIER_REPORT}/phenoplier_traits_long.csv.gz",
    log:
        f"{PHENOPLIER_REPORT}/aggregate.log"
    params:
        conda_env=PHENOPLIER_CFG["conda_env"],
    resources:
        mem_mb=8000,
        runtime=720,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        "bash {input.wrapper} {params.conda_env} --out {output.long} "
        "--clampfull-bp {input.clampfull_bp} --clampbase {input.clampbase} "
        "--canonical {input.canonical} > {log} 2>&1"


# ============================================================
# Step 4: convenience targets
# ============================================================

rule phenoplier_coverage:
    input:
        PHENOPLIER_COV_FULL + PHENOPLIER_COV_BASE,


rule phenoplier_saturation:
    input:
        PHENOPLIER_SAT_FULL + PHENOPLIER_SAT_BASE,


rule phenoplier_finals:
    input:
        PHENOPLIER_FIN_FULL + PHENOPLIER_FIN_BASE,


rule phenoplier_canonical:
    input:
        PHENOPLIER_CANONICAL,


rule phenoplier_final_stores:
    input:
        PHENOPLIER_FIN_STORES,


rule phenoplier_traits:
    """Every per-model GLS summary this workflow produces."""
    input:
        PHENOPLIER_COV_FULL + PHENOPLIER_COV_BASE
        + PHENOPLIER_SAT_FULL + PHENOPLIER_SAT_BASE
        + PHENOPLIER_FIN_FULL + PHENOPLIER_FIN_BASE
        + PHENOPLIER_CANONICAL,


rule phenoplier_long_table:
    input:
        rules.aggregate_phenoplier_traits.output.long,
