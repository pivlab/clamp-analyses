RUNTIME_CFG = config["runtime_benchmark"]
RUNTIME_ROOT = RUNTIME_CFG["output_root"]
RUNTIME_METHODS = list(RUNTIME_CFG["methods"])
RUNTIME_SEEDS = list(RUNTIME_CFG["seeds"])
RUNTIME_THREADS = int(RUNTIME_CFG["threads"])
RUNTIME_METHOD_PATTERN = "|".join(RUNTIME_METHODS)
RUNTIME_SEED_PATTERN = "|".join(str(seed) for seed in RUNTIME_SEEDS)
RUNTIME_NB = os.path.join(REPO_ROOT, RUNTIME_CFG["notebook"])

RUNTIME_TIMING_FILES = expand(
    f"{RUNTIME_ROOT}/timings/{{dataset}}/{{method}}/seed_{{seed}}.csv",
    dataset=DATASETS,
    method=RUNTIME_METHODS,
    seed=RUNTIME_SEEDS,
)


# ============================================================
# Controlled computational-timing analysis over pseudobulk model inputs
#
# One job owns one dataset/method/seed model directory.  Snakemake can
# therefore resume an interrupted 180-fit benchmark without touching the
# production models under output/01_model_building.
# ============================================================

rule runtime_fit_pseudobulk:
    input:
        norm=rules.preprocess_pseudobulk.output.norm,
        cpm_filt=rules.preprocess_pseudobulk.output.cpm_filt,
        k=rules.preprocess_pseudobulk.output.k,
        gmt=rules.pathway_prior.output,
    output:
        model_dir=directory(
            f"{RUNTIME_ROOT}/models/{{dataset}}/{{method}}/seed_{{seed}}"
        ),
        timing=f"{RUNTIME_ROOT}/timings/{{dataset}}/{{method}}/seed_{{seed}}.csv",
    log:
        f"{RUNTIME_ROOT}/logs/{{dataset}}/{{method}}/seed_{{seed}}.log"
    params:
        clamp_max_iter=config["clamp"]["max_iter"],
        flashier_backfit=config["flashier"]["backfit_maxiter"],
        mofa_max_epochs=config["mofa_flex"]["max_epochs"],
        cogaps_iterations=config["cogaps"]["n_iterations"],
    threads: RUNTIME_THREADS
    resources:
        mem_mb=48000,
        runtime=7 * 24 * 60,
    wildcard_constraints:
        dataset=DATASET_PATTERN,
        method=RUNTIME_METHOD_PATTERN,
        seed=RUNTIME_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python scripts/pseudobulk/run_timing_benchmark.py "
        "--dataset {wildcards.dataset} --method {wildcards.method} "
        "--seed {wildcards.seed} --threads {threads} "
        "--norm {input.norm} --cpm-filt {input.cpm_filt} --k {input.k} --gmt {input.gmt} "
        "--output-dir {output.model_dir} --timing-csv {output.timing} --log-file {log} "
        "--clamp-max-iter {params.clamp_max_iter} "
        "--flashier-backfit {params.flashier_backfit} "
        "--mofa-max-epochs {params.mofa_max_epochs} "
        "--cogaps-iterations {params.cogaps_iterations}"


rule computational_timing_report_pseudobulk:
    input:
        timings=RUNTIME_TIMING_FILES,
        pseudobulk_config="workflow/config/pseudobulk.yaml",
        runtime_config="workflow/config/runtime_benchmark.yaml",
        figure_config="config.yaml",
        notebook=RUNTIME_NB,
    output:
        long=f"{RUNTIME_ROOT}/runtime_long.csv",
        seed_totals=f"{RUNTIME_ROOT}/runtime_seed_totals.csv",
        summary=f"{RUNTIME_ROOT}/runtime_summary.csv",
        png=f"{RUNTIME_ROOT}/computational_timing.png",
        pdf=f"{RUNTIME_ROOT}/computational_timing.pdf",
        svg=f"{RUNTIME_ROOT}/computational_timing.svg",
        complete=touch(f"{RUNTIME_ROOT}/notebook.complete"),
    log:
        notebook=f"{RUNTIME_ROOT}/00_computational_timing.executed.ipynb"
    params:
        datasets=DATASETS,
        methods=RUNTIME_METHODS,
        seeds=RUNTIME_SEEDS,
        out_dir=RUNTIME_ROOT,
    conda: "clamp-analyses"
    notebook:
        RUNTIME_NB


rule computational_timing_analysis:
    input:
        rules.computational_timing_report_pseudobulk.output.complete
