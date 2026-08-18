RUNTIME_GTEX_CFG = config["runtime_benchmark_gtex"]
RUNTIME_GTEX_ROOT = RUNTIME_GTEX_CFG["output_root"]
RUNTIME_GTEX_DATASET = RUNTIME_GTEX_CFG["dataset"]
RUNTIME_GTEX_METHODS = list(RUNTIME_GTEX_CFG["methods"])
RUNTIME_GTEX_SEEDS = list(RUNTIME_GTEX_CFG["seeds"])
RUNTIME_GTEX_THREADS = int(RUNTIME_GTEX_CFG["threads"])
RUNTIME_GTEX_METHOD_PATTERN = "|".join(RUNTIME_GTEX_METHODS)
RUNTIME_GTEX_SEED_PATTERN = "|".join(str(seed) for seed in RUNTIME_GTEX_SEEDS)
RUNTIME_GTEX_NB = os.path.join(REPO_ROOT, RUNTIME_GTEX_CFG["notebook"])

RUNTIME_GTEX_TIMING_FILES = expand(
    f"{RUNTIME_GTEX_ROOT}/timings/{{method}}/seed_{{seed}}.csv",
    method=RUNTIME_GTEX_METHODS,
    seed=RUNTIME_GTEX_SEEDS,
)

# Advisory only: resources: timing_slot=1 already guarantees a single fit at a
# time, so nothing competes for RAM on this 125 GB machine.  These are the
# expected peaks, recorded so the intent is reviewable and so the rules stay
# portable if the benchmark is ever moved to a cluster.  Flashier materialises
# t(as.matrix(df)) and its production model rds alone is 5.5 GB; PLIER's is 3.0 GB.
RUNTIME_GTEX_MEM_MB = {
    "Flashier": 110000,
    "PLIER": 96000,
    "CLAMPfull": 80000,
    "CLAMPbase": 80000,
    "GSSig": 72000,
    "MOFA-FLEX": 64000,
    "NMF": 48000,
    "ICA": 48000,
    "PCA": 48000,
}


# ============================================================
# Controlled computational-timing analysis over the GTEx production inputs
#
# 9 methods x 3 seeds = 27 fits on one dataset.  Every fit reuses the preprocessed
# GTEx artifacts produced by rule clamp_gtex, so -- exactly like the pseudobulk
# benchmark, which reuses rules.preprocess_pseudobulk.output -- the reported
# runtimes measure factorization only and exclude the ~48 min GCT parse / FBM
# build / preprocess / z-score / rsvd / K-estimation stage.
#
# Fits are strictly serialized: resources timing_slot=1 plus `--resources
# timing_slot=1` on the command line.  This is what makes the wall-clock numbers
# uncontended, and it also makes OOM structurally impossible.  Snakemake treats a
# resource that is NOT named on --resources as unlimited, so omitting the flag
# silently disables serialization; the report notebook re-derives non-overlap from
# started_epoch/finished_epoch as a backstop and hard-fails if it was skipped.
#
# CoGAPS is deliberately absent from RUNTIME_GTEX_METHODS.  rule cogaps_gtex was
# run on GTEx and killed at its 7-day budget without converging
# (output/01_model_building/01_gtex/CoGAPS/NOT_CONVERGED.txt), which is also why it
# is excluded from rule full_models_gtex.  It has no finite runtime to report and
# is documented as a censored observation rather than dropped.  This exclusion and
# its reason live in workflow/config/runtime_benchmark_gtex.yaml and are printed
# by the report notebook.
# ============================================================

rule runtime_gtex_matrix_shape:
    # df_gtex_fbm_filt.csv is 6.8 GB.  Probe its shape once for the whole benchmark
    # rather than once per fit, which would read 185 GB and flush the page cache
    # immediately after each timed interval.
    input:
        df_csv=rules.clamp_gtex.output.df_csv,
    output:
        shape=f"{RUNTIME_GTEX_ROOT}/matrix_shape.csv",
    resources:
        mem_mb=2000,
        runtime=60,
    conda: "clamp-analyses"
    shell:
        "python scripts/gtex/run_timing_benchmark.py --probe-shape "
        "--matrix-csv {input.df_csv} --shape-csv {output.shape}"


rule runtime_fit_gtex:
    input:
        fbm_filt=rules.clamp_gtex.output.fbm_filt,
        # The FBM rds is a 102 KB handle; the data lives in this backing file, whose
        # absolute path is baked into the rds.  Declared so the dependency is visible
        # even though no rule declares it as an output.
        fbm_bk=f"{GTEX_PROD}/FBMgtex_preproc_filtered.bk",
        svd_res=rules.clamp_gtex.output.svd_res,
        genes=rules.clamp_gtex.output.genes,
        samples=rules.clamp_gtex.output.samples,
        df_rds=rules.clamp_gtex.output.df_rds,
        df_csv=rules.clamp_gtex.output.df_csv,
        k_rds=rules.clamp_gtex.output.k_rds,
        k_csv=rules.clamp_gtex.output.k_csv,
        gmt=rules.pathway_prior.output,
        shape=rules.runtime_gtex_matrix_shape.output.shape,
    output:
        model_dir=directory(f"{RUNTIME_GTEX_ROOT}/models/{{method}}/seed_{{seed}}"),
        timing=f"{RUNTIME_GTEX_ROOT}/timings/{{method}}/seed_{{seed}}.csv",
    log:
        f"{RUNTIME_GTEX_ROOT}/logs/{{method}}/seed_{{seed}}.log"
    params:
        dataset=RUNTIME_GTEX_DATASET,
        # GTEx production hyperparameters, so the timings reflect the published GTEx
        # models rather than the pseudobulk settings.  These differ: flashier 20 vs
        # 10 backfit iterations, MOFA-FLEX 1000 vs 200 epochs.
        clamp_max_iter=config["gtex"]["clamp"]["max_iter"],
        flashier_backfit=config["gtex"]["flashier"]["backfit_maxiter"],
        gss_d_cluster=config["gtex"]["gss"]["d_cluster"],
        mofa_max_epochs=config["gtex"]["mofa_flex_prior"]["max_epochs"],
        min_fraction=config["gtex"]["mofa_flex_prior"]["min_fraction"],
        min_count=config["gtex"]["mofa_flex_prior"]["min_count"],
        max_count=config["gtex"]["mofa_flex_prior"]["max_count"],
        similarity_threshold=config["gtex"]["mofa_flex_prior"]["similarity_threshold"],
    threads: RUNTIME_GTEX_THREADS
    resources:
        timing_slot=1,
        mem_mb=lambda wildcards: RUNTIME_GTEX_MEM_MB[wildcards.method],
        runtime=7 * 24 * 60,
    wildcard_constraints:
        method=RUNTIME_GTEX_METHOD_PATTERN,
        seed=RUNTIME_GTEX_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python scripts/gtex/run_timing_benchmark.py "
        "--dataset {params.dataset} --method {wildcards.method} "
        "--seed {wildcards.seed} --threads {threads} "
        "--fbm-filt {input.fbm_filt} --svd-res {input.svd_res} "
        "--genes {input.genes} --samples {input.samples} "
        "--df-rds {input.df_rds} --df-csv {input.df_csv} "
        "--k-rds {input.k_rds} --k-csv {input.k_csv} --gmt {input.gmt} "
        "--shape-csv {input.shape} "
        "--output-dir {output.model_dir} --timing-csv {output.timing} --log-file {log} "
        "--clamp-max-iter {params.clamp_max_iter} "
        "--flashier-backfit {params.flashier_backfit} "
        "--gss-d-cluster {params.gss_d_cluster} "
        "--mofa-max-epochs {params.mofa_max_epochs} "
        "--min-fraction {params.min_fraction} --min-count {params.min_count} "
        "--max-count {params.max_count} "
        "--similarity-threshold {params.similarity_threshold}"


rule computational_timing_report_gtex:
    input:
        timings=RUNTIME_GTEX_TIMING_FILES,
        shape=rules.runtime_gtex_matrix_shape.output.shape,
        # The pinned GO:BP prior every model is trained against; the report re-checks
        # its sha256 against the pin in pseudobulk.yaml.
        gmt=rules.pathway_prior.output,
        # The evidence for the CoGAPS exclusion, tracked by the DAG so the claim in
        # the report cannot drift from what is on disk.
        cogaps_evidence=RUNTIME_GTEX_CFG["cogaps_evidence"],
        gtex_config="workflow/config/gtex.yaml",
        runtime_config="workflow/config/runtime_benchmark_gtex.yaml",
        figure_config="config.yaml",
        notebook=RUNTIME_GTEX_NB,
    output:
        long=f"{RUNTIME_GTEX_ROOT}/runtime_long.csv",
        # Named per_fit rather than seed_totals: the pseudobulk benchmark sums one
        # seed across six datasets, GTEx is a single dataset so each record is one
        # fit and there is nothing to sum.
        per_fit=f"{RUNTIME_GTEX_ROOT}/runtime_per_fit.csv",
        summary=f"{RUNTIME_GTEX_ROOT}/runtime_summary.csv",
        complete=touch(f"{RUNTIME_GTEX_ROOT}/notebook.complete"),
    log:
        notebook=f"{RUNTIME_GTEX_ROOT}/00_computational_timing.executed.ipynb"
    params:
        dataset=RUNTIME_GTEX_DATASET,
        methods=RUNTIME_GTEX_METHODS,
        seeds=RUNTIME_GTEX_SEEDS,
        threads=RUNTIME_GTEX_THREADS,
        out_dir=RUNTIME_GTEX_ROOT,
    conda: "clamp-analyses"
    notebook:
        RUNTIME_GTEX_NB


rule computational_timing_analysis_gtex:
    input:
        rules.computational_timing_report_gtex.output.complete
