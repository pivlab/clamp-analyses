import re

# ============================================================
# LV-trait association (GLS) via pivlab/phenoplier-cli, run against every
# CLAMP model this repo already produces. phenoplier-cli ships its own
# Snakemake pipeline on snakemake>=9/Python>=3.12, incompatible with this
# repo's snakemake=8.*/Python 3.11 (envs/snakemake.yaml) -- it runs in its own
# conda env (workflow/config/phenoplier.yaml: conda_env) via
# scripts/phenoplier/run_gls.sh, one `phenoplier shortcut gls` invocation per
# model, rather than as part of this Snakefile's own DAG.
#
# Model paths below are read from config/helpers each dataset's own rule file
# already defines (archs4.smk, archs4_coverage.smk, archs4_saturation.smk,
# gtex.smk) -- nothing here duplicates a path that file owns.
#
# Depends on pivlab/phenoplier-cli#85 (open at the time this was written):
# that PR fixes a GLS standard-error collapse that emits artifactual ~0
# p-values for a handful of pathological LVs, specifically in the
# `hall_coverage_rs*`-style sub-sampling models -- i.e. exactly the
# archs4_coverage/archs4_saturation models this rule file runs against.
# Until #85 merges, the coverage/saturation "traits recovered" curves this
# produces can show the same spurious 10%/50% spikes that PR documents. See
# scripts/phenoplier/setup_env.sh, which pins the install ref accordingly.
# ============================================================

PHENOPLIER_CFG = config["phenoplier"]
PHENOPLIER_OUT = PHENOPLIER_CFG["output_root"]
PHENOPLIER_TARGET = config.get("phenoplier_target", PHENOPLIER_CFG["target"])
PHENOPLIER_CLUSTER_CFG = PHENOPLIER_CFG["clusters"][PHENOPLIER_TARGET]

PHENOPLIER_MODELS = {}

# ARCHS4 final model: CLAMPfull per seed, CLAMPbase once (reference seed only
# -- see archs4.smk's A4_REF_DIR comment: CLAMPbase is identical across seeds
# at 100% of the compendium, so only one copy is kept).
for _seed in A4_SEEDS:
    PHENOPLIER_MODELS[f"archs4_final/seed{_seed}/CLAMPfull"] = (
        f"{a4_seed_dir(_seed)}/{A4_MODEL_NAME}.rds"
    )
PHENOPLIER_MODELS["archs4_final/CLAMPbase"] = f"{A4_REF_DIR}/CLAMPbase.rds"

# ARCHS4 coverage: CLAMPfull_bp + CLAMPbase per (fraction, seed).
for _fraction in A4_COV_LEVELS:
    for _seed in a4_cov_seeds_for(_fraction):
        PHENOPLIER_MODELS[f"archs4_coverage/archs4/rs{_fraction}/seed{_seed}/CLAMPfull"] = (
            f"{a4_cov_model_dir('archs4', _fraction, _seed)}/{A4_COV_CFG['model_name']}.rds"
        )
        PHENOPLIER_MODELS[f"archs4_coverage/archs4/rs{_fraction}/seed{_seed}/CLAMPbase"] = (
            f"{a4_cov_cell(_fraction, _seed)}/CLAMPbase.rds"
        )

# ARCHS4 coverage comparators (GTEx, Recount2 -- fit at rs100 only).
for _dataset in A4_COV_COMPARATORS:
    for _seed in a4_cov_seeds_for(100):
        PHENOPLIER_MODELS[f"archs4_coverage/{_dataset}/seed{_seed}/CLAMPfull"] = (
            f"{a4_cov_model_dir(_dataset, 100, _seed)}/{A4_COV_CFG['model_name']}.rds"
        )
    # CLAMPbase is a single pre-existing artifact per comparator dataset, not
    # seed-varied -- see archs4.yaml: coverage.comparators.<dataset>.base.
    # For recount2 this is this repo's only CLAMPbase model: there is no
    # recount2.smk building it from raw data (flagged in the PR description).
    PHENOPLIER_MODELS[f"archs4_coverage/{_dataset}/CLAMPbase"] = (
        A4_COV_CFG["comparators"][_dataset]["base"]
    )

# ARCHS4 saturation: CLAMPfull_bp + CLAMPbase per (fraction, K, seed).
for _fraction, _k, _seed in A4_SAT_CELLS:
    PHENOPLIER_MODELS[f"archs4_saturation/rs{_fraction}/k{_k}/seed{_seed}/CLAMPfull"] = (
        f"{a4_sat_model_dir(_fraction, _k, _seed)}/{A4_SAT_CFG['model_name']}.rds"
    )
    PHENOPLIER_MODELS[f"archs4_saturation/rs{_fraction}/k{_k}/seed{_seed}/CLAMPbase"] = (
        a4_sat_base_rds(_fraction, _k, _seed)
    )

# GTEx production model (clamp_gtex in gtex.smk).
PHENOPLIER_MODELS["gtex/CLAMPfull"] = f"{GTEX_PROD}/CLAMPfull.rds"
PHENOPLIER_MODELS["gtex/CLAMPbase"] = f"{GTEX_PROD}/CLAMPbase.rds"

PHENOPLIER_MODEL_KEYS = list(PHENOPLIER_MODELS)
PHENOPLIER_KEY_PATTERN = "|".join(re.escape(key) for key in PHENOPLIER_MODEL_KEYS)


def phenoplier_model_rds(wildcards):
    return PHENOPLIER_MODELS[wildcards.model_key]


def phenoplier_keys(prefix):
    return [key for key in PHENOPLIER_MODEL_KEYS if key.startswith(prefix)]


rule phenoplier_gls:
    input:
        rds=phenoplier_model_rds,
        script="scripts/phenoplier/run_gls.sh",
    output:
        summary=f"{PHENOPLIER_OUT}/{{model_key}}/gls-summary-phenomexcan.tsv.gz",
    log:
        f"{PHENOPLIER_OUT}/{{model_key}}/phenoplier.log"
    params:
        conda_env=PHENOPLIER_CFG["conda_env"],
        namespace=PHENOPLIER_CFG["namespace"],
        lv_percentile=PHENOPLIER_CFG["lv_percentile"],
        executor=PHENOPLIER_CLUSTER_CFG["executor"],
        cluster=PHENOPLIER_CLUSTER_CFG.get("cluster", ""),
        trait_filter=PHENOPLIER_CFG.get("trait_filter", "biomedical"),
    resources:
        # Matches phenoplier-cli's own documented per-job defaults
        # (scripts/slurm/submit.py: 15 CPUs / 45 GB / 7-day wall clock).
        mem_mb=46000,
        runtime=10080,
    threads: 15
    wildcard_constraints:
        model_key=PHENOPLIER_KEY_PATTERN,
    shell:
        "bash {input.script} {input.rds} {wildcards.model_key} "
        "{params.conda_env} {params.namespace} {params.lv_percentile} "
        "{params.executor} '{params.cluster}' {threads} {params.trait_filter} "
        "{output.summary} > {log} 2>&1"


# Traits are filtered by phenoplier-cli itself, at run time (see
# phenoplier.yaml: trait_filter), so the excluded ones never reach these
# summaries and there is nothing left to filter here -- this rule is now a
# plain concatenation. Each model's own exclusion log (copied next to its
# summary by run_gls.sh as trait_filter_excluded.tsv) records what was
# dropped and why.
rule aggregate_phenoplier_traits:
    input:
        summaries=expand(
            f"{PHENOPLIER_OUT}/{{model_key}}/gls-summary-phenomexcan.tsv.gz",
            model_key=PHENOPLIER_MODEL_KEYS,
        ),
        aggregate_script="scripts/phenoplier/aggregate_traits.py",
        wrapper="scripts/phenoplier/aggregate_traits.sh",
    output:
        long=f"{PHENOPLIER_OUT}/phenoplier_traits_long.csv",
    log:
        f"{PHENOPLIER_OUT}/aggregate.log"
    params:
        conda_env=PHENOPLIER_CFG["conda_env"],
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        "bash {input.wrapper} {PHENOPLIER_OUT} {output.long} "
        "{params.conda_env} > {log} 2>&1"


rule phenoplier_traits:
    input:
        rules.aggregate_phenoplier_traits.output.long,


rule phenoplier_archs4_final:
    input:
        expand(
            f"{PHENOPLIER_OUT}/{{model_key}}/gls-summary-phenomexcan.tsv.gz",
            model_key=phenoplier_keys("archs4_final/"),
        ),


rule phenoplier_archs4_coverage:
    input:
        expand(
            f"{PHENOPLIER_OUT}/{{model_key}}/gls-summary-phenomexcan.tsv.gz",
            model_key=phenoplier_keys("archs4_coverage/"),
        ),


rule phenoplier_archs4_saturation:
    input:
        expand(
            f"{PHENOPLIER_OUT}/{{model_key}}/gls-summary-phenomexcan.tsv.gz",
            model_key=phenoplier_keys("archs4_saturation/"),
        ),


rule phenoplier_gtex:
    input:
        expand(
            f"{PHENOPLIER_OUT}/{{model_key}}/gls-summary-phenomexcan.tsv.gz",
            model_key=phenoplier_keys("gtex/"),
        ),
