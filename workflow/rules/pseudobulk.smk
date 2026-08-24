DATASETS = list(config["datasets"])
DATASET_PATTERN = "|".join(DATASETS)
BUILD_DATASETS = [
    dataset for dataset in DATASETS
    if config["datasets"][dataset].get("build_pseudobulk", False)
]
BUILD_DATASET_PATTERN = "|".join(BUILD_DATASETS)
PROD = config["paths"]["production"]
BIO = config["paths"]["biology"]
MODEL_NB = os.path.join(REPO_ROOT, config["paths"]["model_notebooks"])
BIO_NB = os.path.join(REPO_ROOT, config["paths"]["biology_notebooks"])
GMT = config["references"]["go_bp_file"]
CELL_MARKER_FILE = config["references"]["cell_marker_file"]
ALLEN_BRAIN_GMT_FILE = config["references"]["allen_brain_gmt_file"]
AZIMUTH_FILE = config["references"]["azimuth_file"]
CV_FOLDS = list(range(1, config["grouped_cv"]["n_folds"] + 1))
CV_FOLD_PATTERN = "|".join(str(fold) for fold in CV_FOLDS)

METHOD_FILES = config["methods"]
METHOD_B = list(METHOD_FILES.values())
METHOD_LABELS = {name: f"models/{path}" for name, path in METHOD_FILES.items()}

def dataset_raw_inputs(wildcards):
    cfg = config["datasets"][wildcards.dataset]
    keys = ["raw"]
    if cfg["kind"] == "matrix_market":
        keys.extend(["features", "barcodes", "metadata"])
    return [cfg[key] for key in keys]


def pseudobulk_path(dataset, filename):
    cfg = config["datasets"][dataset]
    if cfg.get("build_pseudobulk", False):
        return f"{PROD}/{dataset}/pseudobulk/{filename}"
    return f"{cfg['pseudobulk_dir']}/{filename}"


def pseudobulk_input(filename):
    return lambda wildcards: pseudobulk_path(wildcards.dataset, filename)


def all_pseudobulk_inputs(filename):
    return [pseudobulk_path(dataset, filename) for dataset in DATASETS]


# ============================================================
# Reference gene set for pathway priors
# ============================================================
# Downloads and checksums the fixed GO:BP gene-set file used as the shared
# pathway prior by every method that needs one (CLAMPfull, PLIER, MOFA-FLEX).

rule pathway_prior:
    output:
        GMT
    params:
        url=config["references"]["go_bp_url"],
        sha256=config["references"]["go_bp_sha256"],
    conda: "clamp-analyses"
    shell:
        "mkdir -p $(dirname {output}) && "
        "curl --fail --location --retry 3 '{params.url}' --output {output}.download && "
        "echo '{params.sha256}  {output}.download' | sha256sum --check --status && "
        "mv {output}.download {output}"


# ============================================================
# Step 1: build pseudobulk data
# ============================================================
# For each dataset, produces
# the pseudobulk expression matrix: which samples qualify, the CPM-averaged
# expression matrix built from raw single cells, and the matching
# ground-truth cell-type composition table used to score every model later.

rule cohort_manifest:
    input:
        raw=dataset_raw_inputs,
        config="workflow/config/pseudobulk.yaml",
    output:
        manifest=f"{PROD}/{{dataset}}/pseudobulk/cohort_manifest.csv",
    wildcard_constraints:
        dataset=BUILD_DATASET_PATTERN
    resources:
        mem_mb=8000,
        runtime=60,
    conda: "clamp-analyses"
    shell:
        "python scripts/pseudobulk/build_cohort_manifest.py "
        "--config {input.config} --dataset {wildcards.dataset} --output {output.manifest}"


rule build_pseudobulk:
    input:
        raw=dataset_raw_inputs,
        manifest=rules.cohort_manifest.output.manifest,
        config="workflow/config/pseudobulk.yaml",
    output:
        bulk=f"{PROD}/{{dataset}}/pseudobulk/bulk_expr.csv",
        info=f"{PROD}/{{dataset}}/pseudobulk/patient_info.csv",
        truth0=f"{PROD}/{{dataset}}/pseudobulk/truthFrac_v0.csv",
        truth1=f"{PROD}/{{dataset}}/pseudobulk/truthFrac_v1.csv",
        summary=f"{PROD}/{{dataset}}/pseudobulk/build_summary.csv",
    wildcard_constraints:
        dataset=BUILD_DATASET_PATTERN
    resources:
        mem_mb=32000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "python scripts/pseudobulk/build_pseudobulk.py "
        "--config {input.config} --dataset {wildcards.dataset} "
        "--cohort-manifest {input.manifest} "
        "--bulk {output.bulk} --patient-info {output.info} "
        "--truth-v0 {output.truth0} --truth-v1 {output.truth1} --summary {output.summary}"

rule pseudobulk_data:
    input:
        bulk=all_pseudobulk_inputs("bulk_expr.csv"),
        info=all_pseudobulk_inputs("patient_info.csv"),
        truth=all_pseudobulk_inputs("truthFrac_v0.csv"),
        builds=expand(f"{PROD}/{{dataset}}/pseudobulk/build_summary.csv", dataset=BUILD_DATASETS)


# ============================================================
# Step 2: preprocess pseudobulk data
# ============================================================
# filters out low-signal genes, z-scores and
# picks the model rank (k) via SVD elbow selection

rule preprocess_pseudobulk:
    input:
        bulk=pseudobulk_input("bulk_expr.csv")
    output:
        norm=f"{PROD}/{{dataset}}/preprocessing/norm.csv",
        cpm_filt=f"{PROD}/{{dataset}}/preprocessing/cpm_filt.csv",
        row_stats=f"{PROD}/{{dataset}}/preprocessing/row_stats.csv",
        k=f"{PROD}/{{dataset}}/preprocessing/k.csv",
        diagnostics=f"{PROD}/{{dataset}}/preprocessing/rank_diagnostics.csv",
    params:
        mean_cutoff=config["preprocess"]["mean_cutoff"],
        var_cutoff=config["preprocess"]["var_cutoff"],
        seed=config["preprocess"]["seed"],
    wildcard_constraints:
        dataset=DATASET_PATTERN
    resources:
        mem_mb=32000,
        runtime=240,
    conda: "clamp-analyses"
    shell:
        "MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "Rscript scripts/pseudobulk/preprocess.R --input {input.bulk} --input-scale cpm "
        "--norm {output.norm} --cpm-filt {output.cpm_filt} --row-stats {output.row_stats} --k {output.k} "
        "--diagnostics {output.diagnostics} "
        "--mean-cutoff {params.mean_cutoff} --var-cutoff {params.var_cutoff} --seed {params.seed}"


# ============================================================
# Step 3: decomposition pseudobulk models
# ============================================================
# CLAMPbase/CLAMPfull, PLIER, PCA/NMF/ICA,
# Flashier, MOFA-FLEX, GSSig, CoGAPS models with pseudobulk data

rule clamp_pseudobulk:
    input:
        norm=rules.preprocess_pseudobulk.output.norm,
        k=rules.preprocess_pseudobulk.output.k,
        gmt=rules.pathway_prior.output,
    output:
        base_B=f"{PROD}/{{dataset}}/models/CLAMPbase/B.csv",
        base_Z=f"{PROD}/{{dataset}}/models/CLAMPbase/Z.csv",
        base_rds=f"{PROD}/{{dataset}}/models/CLAMPbase/CLAMPbase.rds",
        full_B=f"{PROD}/{{dataset}}/models/CLAMPfull/B.csv",
        full_Z=f"{PROD}/{{dataset}}/models/CLAMPfull/Z.csv",
        full_rds=f"{PROD}/{{dataset}}/models/CLAMPfull/CLAMPfull.rds",
        l2=f"{PROD}/{{dataset}}/models/CLAMPfull/L2.csv",
        summary=f"{PROD}/{{dataset}}/models/CLAMPfull/summary.csv",
    params:
        base_dir=lambda wc: f"{PROD}/{wc.dataset}/models/CLAMPbase",
        full_dir=lambda wc: f"{PROD}/{wc.dataset}/models/CLAMPfull",
        max_iter=config["clamp"]["max_iter"],
        seed=config["clamp"]["seed"],
    wildcard_constraints:
        dataset=DATASET_PATTERN
    resources:
        mem_mb=48000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "Rscript scripts/pseudobulk/clamp.R --norm {input.norm} --k {input.k} --gmt {input.gmt} "
        "--base-dir {params.base_dir} --full-dir {params.full_dir} "
        "--max-iter {params.max_iter} --seed {params.seed}"


rule plier_pseudobulk:
    input:
        norm=rules.preprocess_pseudobulk.output.norm,
        k=rules.preprocess_pseudobulk.output.k,
        gmt=rules.pathway_prior.output,
    output:
        B=f"{PROD}/{{dataset}}/models/PLIER/B.csv",
        Z=f"{PROD}/{{dataset}}/models/PLIER/Z.csv",
        rds=f"{PROD}/{{dataset}}/models/PLIER/PLIER.rds",
        summary=f"{PROD}/{{dataset}}/models/PLIER/summary.csv",
    params:
        out_dir=lambda wc: f"{PROD}/{wc.dataset}/models/PLIER",
        seed=config["preprocess"]["seed"],
    wildcard_constraints: dataset=DATASET_PATTERN
    resources: mem_mb=48000, runtime=1440
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/pseudobulk/plier.R --norm {input.norm} --k {input.k} --gmt {input.gmt} "
        "--out-dir {params.out_dir} --seed {params.seed}"


rule pca_nmf_ica_pseudobulk:
    input:
        norm=rules.preprocess_pseudobulk.output.norm,
        k=rules.preprocess_pseudobulk.output.k,
    output:
        pca_B=f"{PROD}/{{dataset}}/models/PCA/B.csv",
        pca_Z=f"{PROD}/{{dataset}}/models/PCA/Z.csv",
        pca_model=f"{PROD}/{{dataset}}/models/PCA/pca_model.pkl",
        nmf_B=f"{PROD}/{{dataset}}/models/NMF/B.csv",
        nmf_Z=f"{PROD}/{{dataset}}/models/NMF/Z.csv",
        nmf_model=f"{PROD}/{{dataset}}/models/NMF/nmf_model.pkl",
        ica_B=f"{PROD}/{{dataset}}/models/ICA/B.csv",
        ica_Z=f"{PROD}/{{dataset}}/models/ICA/Z.csv",
        ica_model=f"{PROD}/{{dataset}}/models/ICA/ica_model.pkl",
    params:
        models_dir=lambda wc: f"{PROD}/{wc.dataset}/models",
        seed=config["preprocess"]["seed"],
    wildcard_constraints: dataset=DATASET_PATTERN
    resources: mem_mb=48000, runtime=720
    conda: "clamp-analyses"
    shell:
        "python scripts/pseudobulk/pca_nmf_ica.py --norm {input.norm} --k {input.k} "
        "--models-dir {params.models_dir} --seed {params.seed}"


rule flashier_pseudobulk:
    input:
        norm=rules.preprocess_pseudobulk.output.norm,
        k=rules.preprocess_pseudobulk.output.k,
    output:
        B=f"{PROD}/{{dataset}}/models/flashier/B.csv",
        Z=f"{PROD}/{{dataset}}/models/flashier/Z.csv",
        rds=f"{PROD}/{{dataset}}/models/flashier/flashier_model.rds",
    params:
        out_dir=lambda wc: f"{PROD}/{wc.dataset}/models/flashier",
        backfit=config["flashier"]["backfit_maxiter"],
        seed=config["preprocess"]["seed"],
    wildcard_constraints: dataset=DATASET_PATTERN
    resources: mem_mb=48000, runtime=1440
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/pseudobulk/flashier.R --norm {input.norm} --k {input.k} "
        "--out-dir {params.out_dir} --backfit-maxiter {params.backfit} --seed {params.seed}"


rule mofa_flex_pseudobulk:
    input:
        norm=rules.preprocess_pseudobulk.output.norm,
        k=rules.preprocess_pseudobulk.output.k,
        gmt=rules.pathway_prior.output,
    output:
        B=f"{PROD}/{{dataset}}/models/MOFA_FLEX_PRIOR/B.csv",
        Z=f"{PROD}/{{dataset}}/models/MOFA_FLEX_PRIOR/Z.csv",
        model=f"{PROD}/{{dataset}}/models/MOFA_FLEX_PRIOR/model.pkl",
    params:
        out_dir=lambda wc: f"{PROD}/{wc.dataset}/models/MOFA_FLEX_PRIOR",
        epochs=config["mofa_flex"]["max_epochs"],
        seed=config["mofa_flex"]["seed"],
    wildcard_constraints: dataset=DATASET_PATTERN
    resources: mem_mb=48000, runtime=1440
    conda: "clamp-analyses"
    shell:
        "python scripts/pseudobulk/mofa_flex.py --norm {input.norm} --k {input.k} --gmt {input.gmt} "
        "--out-dir {params.out_dir} --max-epochs {params.epochs} --seed {params.seed}"


rule gss_pseudobulk:
    input:
        norm=rules.preprocess_pseudobulk.output.norm,
        k=rules.preprocess_pseudobulk.output.k,
    output:
        B=f"{PROD}/{{dataset}}/models/GSS/B.csv",
        Z=f"{PROD}/{{dataset}}/models/GSS/Z.csv",
        rds=f"{PROD}/{{dataset}}/models/GSS/RAVmodel.rds",
    params:
        out_dir=lambda wc: f"{PROD}/{wc.dataset}/models/GSS",
        seed=config["preprocess"]["seed"],
    wildcard_constraints: dataset=DATASET_PATTERN
    resources: mem_mb=48000, runtime=720
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/pseudobulk/gss.R --dataset {wildcards.dataset} --norm {input.norm} "
        "--k {input.k} --out-dir {params.out_dir} --seed {params.seed}"


rule cogaps_pseudobulk:
    input:
        cpm_filt=rules.preprocess_pseudobulk.output.cpm_filt,
        k=rules.preprocess_pseudobulk.output.k,
    output:
        B=f"{PROD}/{{dataset}}/models/CoGAPS/B.csv",
        Z=f"{PROD}/{{dataset}}/models/CoGAPS/Z.csv",
        rds=f"{PROD}/{{dataset}}/models/CoGAPS/cogaps_model.rds",
    params:
        out_dir=lambda wc: f"{PROD}/{wc.dataset}/models/CoGAPS",
        n_iterations=config["cogaps"]["n_iterations"],
        n_threads=config["cogaps"]["n_threads"],
        seed=config["cogaps"]["seed"],
    wildcard_constraints: dataset=DATASET_PATTERN
    resources: mem_mb=48000, runtime=1440
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/pseudobulk/cogaps.R --cpm-filt {input.cpm_filt} --k {input.k} "
        "--out-dir {params.out_dir} --n-iterations {params.n_iterations} "
        "--n-threads {params.n_threads} --seed {params.seed}"

rule full_models_pseudobulk:
    input:
        expand(f"{PROD}/{{dataset}}/models/{{method}}", dataset=DATASETS, method=METHOD_B)


rule model_building_qc_pseudobulk:
    input:
        pseudobulk=all_pseudobulk_inputs("bulk_expr.csv"),
        builds=expand(f"{PROD}/{{dataset}}/pseudobulk/build_summary.csv", dataset=BUILD_DATASETS),
        ranks=expand(f"{PROD}/{{dataset}}/preprocessing/rank_diagnostics.csv", dataset=DATASETS),
        models=expand(f"{PROD}/{{dataset}}/models/{{method}}", dataset=DATASETS, method=METHOD_B),
        projections=expand(f"{PROD}/{{dataset}}/single_cell_projection/projection_summary.csv", dataset=DATASETS),
        notebook=f"{MODEL_NB}/00_model_building_qc.ipynb",
    output:
        summary=f"{PROD}/qc/model_building_qc.csv",
        matrices=f"{PROD}/qc/model_matrix_qc.csv",
        complete=touch(f"{PROD}/qc/notebook.complete"),
    log:
        notebook=f"{MODEL_NB}/00_model_building_qc.executed.ipynb",
    params:
        datasets=DATASETS,
        methods=METHOD_FILES,
        mod_root=PROD,
        out_dir=f"{PROD}/qc",
    conda: "clamp-analyses"
    notebook:
        f"{MODEL_NB}/00_model_building_qc.ipynb"


# ============================================================
# Step 4: grouped cross-validation for CLAMPfull
# ============================================================
# Splits each dataset into 5 donor-preserving folds, refits CLAMPfull per
# fold on the training samples, projects the held-out samples onto the
# resulting latent variables, and pools all folds into out-of-fold
# predictions and metrics

rule grouped_cv_folds_pseudobulk:
    input:
        bulk=pseudobulk_input("bulk_expr.csv"),
        truth=pseudobulk_input("truthFrac_v0.csv"),
        info=pseudobulk_input("patient_info.csv"),
    output:
        membership=f"{PROD}/{{dataset}}/grouped_cv/fold_membership.csv",
        summary=f"{PROD}/{{dataset}}/grouped_cv/fold_summary.csv",
    params:
        n_folds=config["grouped_cv"]["n_folds"],
        group_col=config["grouped_cv"]["group_col"],
        seed=config["grouped_cv"]["split_seed"],
    wildcard_constraints: dataset=DATASET_PATTERN
    resources: mem_mb=8000, runtime=60
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/pseudobulk/make_grouped_cv_folds.R --dataset {wildcards.dataset} "
        "--bulk {input.bulk} --truth {input.truth} --patient-info {input.info} "
        "--output {output.membership} --summary {output.summary} "
        "--n-folds {params.n_folds} --group-col {params.group_col} --seed {params.seed}"


rule clampfull_grouped_cv_pseudobulk:
    input:
        bulk=pseudobulk_input("bulk_expr.csv"),
        truth=pseudobulk_input("truthFrac_v0.csv"),
        membership=rules.grouped_cv_folds_pseudobulk.output.membership,
        gmt=rules.pathway_prior.output,
    output:
        train_B=f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/train_B.csv",
        train_Z=f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/train_Z.csv",
        test_B=f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/test_B.csv",
        row_stats=f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/row_stats.csv",
        train_truth=f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/train_truth.csv",
        test_truth=f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/test_truth.csv",
        summary=f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/summary.csv",
        model=f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/CLAMPfull.rds",
    params:
        out_dir=lambda wc: f"{PROD}/{wc.dataset}/grouped_cv/fold{wc.fold}/CLAMPfull",
        seed=config["clamp"]["seed"],
        max_iter=config["clamp"]["max_iter"],
        mean_cutoff=config["preprocess"]["mean_cutoff"],
        var_cutoff=config["preprocess"]["var_cutoff"],
    wildcard_constraints:
        dataset=DATASET_PATTERN,
        fold=CV_FOLD_PATTERN,
    resources: mem_mb=48000, runtime=1440
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/pseudobulk/clampfull_grouped_cv.R --dataset {wildcards.dataset} "
        "--fold {wildcards.fold} --bulk {input.bulk} --input-scale cpm --truth {input.truth} "
        "--membership {input.membership} --gmt {input.gmt} "
        "--out-dir {params.out_dir} --seed {params.seed} "
        "--max-iter {params.max_iter} --mean-cutoff {params.mean_cutoff} --var-cutoff {params.var_cutoff}"

rule grouped_cv_models_pseudobulk:
    input:
        expand(
            f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/CLAMPfull.rds",
            dataset=DATASETS,
            fold=CV_FOLDS,
        )

rule grouped_cv_analysis_pseudobulk:
    input:
        memberships=expand(f"{PROD}/{{dataset}}/grouped_cv/fold_membership.csv", dataset=DATASETS),
        models=expand(
            f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/CLAMPfull.rds",
            dataset=DATASETS,
            fold=CV_FOLDS,
        ),
        train_B=expand(
            f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/train_B.csv",
            dataset=DATASETS,
            fold=CV_FOLDS,
        ),
        train_Z=expand(
            f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/train_Z.csv",
            dataset=DATASETS,
            fold=CV_FOLDS,
        ),
        test_B=expand(
            f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/test_B.csv",
            dataset=DATASETS,
            fold=CV_FOLDS,
        ),
        train_truth=expand(
            f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/train_truth.csv",
            dataset=DATASETS,
            fold=CV_FOLDS,
        ),
        test_truth=expand(
            f"{PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/test_truth.csv",
            dataset=DATASETS,
            fold=CV_FOLDS,
        ),
        analysis="scripts/pseudobulk/analyze_clampfull_grouped_cv.R",
        common="scripts/pseudobulk/common.R",
        config="workflow/config/pseudobulk.yaml",
        analysis_config="workflow/config/cell_type_analysis.yaml",
    output:
        membership=f"{PROD}/grouped_cv_analysis/fold_membership.csv",
        calibrations=f"{PROD}/grouped_cv_analysis/fold_calibrations.csv",
        predictions=f"{PROD}/grouped_cv_analysis/oof_predictions.csv",
        metrics=f"{PROD}/grouped_cv_analysis/oof_metrics.csv",
        thresholded_metrics=f"{PROD}/grouped_cv_analysis/thresholded_metrics.csv",
        thresholded_summary=f"{PROD}/grouped_cv_analysis/thresholded_summary.csv",
        thresholded_dataset_summary=f"{PROD}/grouped_cv_analysis/thresholded_dataset_summary.csv",
        stability=f"{PROD}/grouped_cv_analysis/lv_selection_stability.csv",
        summary=f"{PROD}/grouped_cv_analysis/method_summary.csv",
    resources: mem_mb=8000, runtime=60
    conda: "clamp-analyses"
    shell:
        "Rscript {input.analysis}"


# ============================================================
# Step 5: project single cells onto the CLAMP latent variables
# ============================================================
# Projects every raw single cell from each
# dataset onto the full-data CLAMPfull models
# so each individual cell's LV activity can later be checked against its annotated
# cell type

rule single_cell_projection_pseudobulk:
    input:
        raw=dataset_raw_inputs,
        z=rules.clamp_pseudobulk.output.full_Z,
        l2=rules.clamp_pseudobulk.output.l2,
        row_stats=rules.preprocess_pseudobulk.output.row_stats,
        config="workflow/config/pseudobulk.yaml",
    output:
        scores=f"{PROD}/{{dataset}}/single_cell_projection/single_cell_lv_scores.h5",
        summary=f"{PROD}/{{dataset}}/single_cell_projection/projection_summary.csv",
    params:
        chunk_cells=config["projection"]["chunk_cells"],
        chunk_nnz=config["projection"]["chunk_nnz"],
    wildcard_constraints: dataset=DATASET_PATTERN
    resources: mem_mb=64000, runtime=4320
    conda: "clamp-analyses"
    shell:
        "python scripts/pseudobulk/single_cell_projection.py --config {input.config} "
        "--dataset {wildcards.dataset} --z {input.z} --l2 {input.l2} "
        "--row-stats {input.row_stats} --output {output.scores} --summary {output.summary} "
        "--chunk-cells {params.chunk_cells} --chunk-nnz {params.chunk_nnz}"


rule single_cell_projections:
    input:
        expand(f"{PROD}/{{dataset}}/single_cell_projection/single_cell_lv_scores.h5", dataset=DATASETS)


# ============================================================
# Step 6: biology analyses
# ============================================================
# Turns model outputs into biological results: benchmarking every method
# against the ground truth, the held-out cross-validation report,
# disentangling latent variables against marker genes, single-cell recovery,
# and the hard-to-distinguish cell-type-pair panels.

rule benchmark_pseudobulk:
    input:
        models=expand(f"{PROD}/{{dataset}}/models/{{method}}", dataset=DATASETS, method=METHOD_B),
        truths=all_pseudobulk_inputs("truthFrac_v0.csv"),
        notebook=f"{BIO_NB}/00_benchmark.ipynb",
    output:
        long=f"{BIO}/00_benchmark/benchmark_long.csv",
        summary=f"{BIO}/00_benchmark/benchmark_summary.csv",
        assignments=f"{BIO}/00_benchmark/clampfull_lv_assignments.csv",
        corr=f"{BIO}/00_benchmark/clampfull_lv_ct_corr_full.csv",
        bootstrap=f"{BIO}/00_benchmark/bootstrap_winrate.csv",
        effects=f"{BIO}/00_benchmark/paired_wilcoxon_effect_sizes.csv",
        qvalues=f"{BIO}/00_benchmark/paired_top3_bracket_qvalues.csv",
        complete=touch(f"{BIO}/00_benchmark/notebook.complete"),
    log:
        notebook=f"{BIO_NB}/00_benchmark.executed.ipynb",
    params:
        methods=METHOD_LABELS,
        datasets=DATASETS,
        mod_root=PROD,
        out_dir=f"{BIO}/00_benchmark",
    conda: "clamp-analyses"
    notebook:
        f"{BIO_NB}/00_benchmark.ipynb"


rule holdout_report_pseudobulk:
    input:
        predictions=rules.grouped_cv_analysis_pseudobulk.output.predictions,
        thresholded_metrics=rules.grouped_cv_analysis_pseudobulk.output.thresholded_metrics,
        thresholded_summary=rules.grouped_cv_analysis_pseudobulk.output.thresholded_summary,
        thresholded_dataset_summary=rules.grouped_cv_analysis_pseudobulk.output.thresholded_dataset_summary,
        notebook=f"{BIO_NB}/01_holdout80.ipynb",
    output:
        complete=touch(f"{BIO}/01_holdout80/notebook.complete"),
    log:
        notebook=f"{BIO_NB}/01_holdout80.executed.ipynb",
    params:
        datasets=DATASETS,
        out_dir=f"{BIO}/01_holdout80",
    conda: "clamp-analyses"
    notebook:
        f"{BIO_NB}/01_holdout80.ipynb"


rule disentangle_pseudobulk:
    input:
        models=expand(f"{PROD}/{{dataset}}/models/CLAMPfull/CLAMPfull.rds", dataset=DATASETS),
        assignments=rules.benchmark_pseudobulk.output.assignments,
        corr=rules.benchmark_pseudobulk.output.corr,
        truths=all_pseudobulk_inputs("truthFrac_v0.csv"),
        cell_marker_file=CELL_MARKER_FILE,
        allen_brain_gmt_file=ALLEN_BRAIN_GMT_FILE,
        azimuth_file=AZIMUTH_FILE,
        notebook=f"{BIO_NB}/02_disentangle.ipynb",
    output:
        top=f"{BIO}/02_disentangle/top_lvs_per_celltype.csv",
        margins=f"{BIO}/02_disentangle/lv_specificity_margins.csv",
        quality=f"{BIO}/02_disentangle/lv_assignment_quality.csv",
        quality_by_dataset=f"{BIO}/02_disentangle/lv_assignment_quality_by_dataset.csv",
        selected=f"{BIO}/02_disentangle/module_selected_lvs.csv",
        pathways=f"{BIO}/02_disentangle/module_top_pathways.csv",
        enrichment=f"{BIO}/02_disentangle/marker_enrichment_long.csv",
        effects=f"{BIO}/02_disentangle/lv_marker_effects.csv",
        recovery=f"{BIO}/02_disentangle/marker_pathway_recovery.csv",
        module_panel_ready=f"{BIO}/02_disentangle/module_panel_ready.csv",
        complete=touch(f"{BIO}/02_disentangle/notebook.complete"),
    log:
        notebook=f"{BIO_NB}/02_disentangle.executed.ipynb",
    params:
        datasets=DATASETS,
        mod_root=PROD,
        out_dir=f"{BIO}/02_disentangle",
    conda: "clamp-analyses"
    notebook:
        f"{BIO_NB}/02_disentangle.ipynb"


rule single_cell_recovery_pseudobulk:
    input:
        scores=expand(f"{PROD}/{{dataset}}/single_cell_projection/single_cell_lv_scores.h5", dataset=DATASETS),
        assignments=rules.benchmark_pseudobulk.output.assignments,
        notebook=f"{BIO_NB}/03_b_matrix_singlecell.ipynb",
    output:
        heatmap=f"{BIO}/03_b_matrix_singlecell/heatmap_long.csv",
        recovery=f"{BIO}/03_b_matrix_singlecell/recovery_summary.csv",
        purity=f"{BIO}/03_b_matrix_singlecell/purity_corrected_summary.csv",
        assignments=f"{BIO}/03_b_matrix_singlecell/lv_assignments.csv",
        overall=f"{BIO}/03_b_matrix_singlecell/recovery_overall.csv",
        complete=touch(f"{BIO}/03_b_matrix_singlecell/notebook.complete"),
    log:
        notebook=f"{BIO_NB}/03_b_matrix_singlecell.executed.ipynb",
    params:
        datasets=DATASETS,
        mod_root=PROD,
        out_dir=f"{BIO}/03_b_matrix_singlecell",
    conda: "clamp-analyses"
    notebook:
        f"{BIO_NB}/03_b_matrix_singlecell.ipynb"


rule hard_cell_types_pseudobulk:
    input:
        top=rules.disentangle_pseudobulk.output.top,
        margins=rules.disentangle_pseudobulk.output.margins,
        selected=rules.disentangle_pseudobulk.output.selected,
        pathways=rules.disentangle_pseudobulk.output.pathways,
        enrichment=rules.disentangle_pseudobulk.output.enrichment,
        effects=rules.disentangle_pseudobulk.output.effects,
        pathway_recovery=rules.disentangle_pseudobulk.output.recovery,
        corr=rules.benchmark_pseudobulk.output.corr,
        recovery=rules.single_cell_recovery_pseudobulk.output.recovery,
        notebook=f"{BIO_NB}/04_hard_cell_types.ipynb",
    output:
        table=f"{BIO}/04_hard_cell_types/hard_cell_types.csv",
        correlations=f"{BIO}/04_hard_cell_types/related_lv_correlations.csv",
        enrichment=f"{BIO}/04_hard_cell_types/hard_marker_enrichment.csv",
        group_dot_ready=f"{BIO}/04_hard_cell_types/hard_group_dot_ready.csv",
        complete=touch(f"{BIO}/04_hard_cell_types/notebook.complete"),
    log:
        notebook=f"{BIO_NB}/04_hard_cell_types.executed.ipynb",
    params:
        out_dir=f"{BIO}/04_hard_cell_types",
    conda: "clamp-analyses"
    notebook:
        f"{BIO_NB}/04_hard_cell_types.ipynb"


rule hard_pair_loadings_pseudobulk:
    input:
        models=expand(f"{PROD}/{{dataset}}/models/CLAMPfull/Z.csv", dataset=DATASETS),
        top=rules.disentangle_pseudobulk.output.top,
        corr=rules.benchmark_pseudobulk.output.corr,
        cell_marker_file=CELL_MARKER_FILE,
        allen_brain_gmt_file=ALLEN_BRAIN_GMT_FILE,
        azimuth_file=AZIMUTH_FILE,
        notebook=f"{BIO_NB}/05_hard_pair_loadings.ipynb",
    output:
        panel_ready=f"{BIO}/05_hard_pair_loadings/hard_pair_loadings_panel_ready.csv",
        stats=f"{BIO}/05_hard_pair_loadings/hard_pair_loading_stats.csv",
        sweep=f"{BIO}/05_hard_pair_loadings/hard_pair_marker_resource_sweep.csv",
        complete=touch(f"{BIO}/05_hard_pair_loadings/notebook.complete"),
    log:
        notebook=f"{BIO_NB}/05_hard_pair_loadings.executed.ipynb",
    params:
        out_dir=f"{BIO}/05_hard_pair_loadings",
    conda: "clamp-analyses"
    notebook:
        f"{BIO_NB}/05_hard_pair_loadings.ipynb"


rule biology_pseudobulk:
    input:
        rules.benchmark_pseudobulk.output.complete,
        rules.holdout_report_pseudobulk.output.complete,
        rules.disentangle_pseudobulk.output.complete,
        rules.single_cell_recovery_pseudobulk.output.complete,
        rules.hard_cell_types_pseudobulk.output.complete,
        rules.hard_pair_loadings_pseudobulk.output.complete,
