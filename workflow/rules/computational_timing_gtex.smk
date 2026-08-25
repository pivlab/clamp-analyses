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
# GTEx computational-timing analysis 
#
# 9 methods x 3 seeds = 27 fits on one dataset
# ============================================================

rule runtime_gtex_matrix_shape:
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
        gmt=rules.pathway_prior.output,
        cogaps_evidence=RUNTIME_GTEX_CFG["cogaps_evidence"],
        gtex_config="workflow/config/gtex.yaml",
        runtime_config="workflow/config/runtime_benchmark_gtex.yaml",
        figure_config="config.yaml",
        notebook=RUNTIME_GTEX_NB,
    output:
        long=f"{RUNTIME_GTEX_ROOT}/runtime_long.csv",
        per_fit=f"{RUNTIME_GTEX_ROOT}/runtime_per_fit.csv",
        summary=f"{RUNTIME_GTEX_ROOT}/runtime_summary.csv",
        complete=touch(f"{RUNTIME_GTEX_ROOT}/notebook.complete"),
    log:
        notebook=os.path.join(os.path.dirname(RUNTIME_GTEX_NB),
                              "00_computational_timing.executed.ipynb")
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
