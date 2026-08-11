DB_CFG = config["donor_bulk"]
DB_DATASETS = list(DB_CFG["datasets"])
DB_DATASET_PATTERN = "|".join(DB_DATASETS)
DB_PROD = DB_CFG["paths"]["production"]
DB_BIO = DB_CFG["paths"]["biology"]
DB_MODEL_NB = os.path.join(REPO_ROOT, DB_CFG["paths"]["model_notebooks"])
DB_BIO_NB = os.path.join(REPO_ROOT, DB_CFG["paths"]["biology_notebooks"])


def donor_bulk_raw_inputs(wildcards):
    return dataset_raw_inputs(wildcards)


def donor_bulk_manifest(wildcards):
    cfg = config["datasets"][wildcards.dataset]
    if cfg.get("cohort_manifest"):
        return cfg["cohort_manifest"]
    return pseudobulk_path(wildcards.dataset, "cohort_manifest.csv")


# ============================================================
# Donor-bulk libraries and independent expression UMAPs
# ============================================================

rule build_donor_bulk:
    input:
        raw=donor_bulk_raw_inputs,
        manifest=donor_bulk_manifest,
        config="workflow/config/pseudobulk.yaml",
        script="scripts/donor_bulk/build_donor_bulk.py",
    output:
        counts=f"{DB_PROD}/{{dataset}}/bulk/bulk_counts.csv",
        mean_cell_cpm=f"{DB_PROD}/{{dataset}}/bulk/bulk_mean_cell_cpm.csv",
        truth=f"{DB_PROD}/{{dataset}}/bulk/truthFrac_v0.csv",
        info=f"{DB_PROD}/{{dataset}}/bulk/patient_info.csv",
        summary=f"{DB_PROD}/{{dataset}}/bulk/aggregation_summary.csv",
        sample_counts=f"{DB_PROD}/{{dataset}}/single_cell_umap/sampled_counts.npz",
        sample_cells=f"{DB_PROD}/{{dataset}}/single_cell_umap/sampled_cells.csv",
        sample_genes=f"{DB_PROD}/{{dataset}}/single_cell_umap/sampled_genes.csv",
    params:
        chunk_cells=DB_CFG["aggregation"]["chunk_cells"],
        chunk_nnz=DB_CFG["aggregation"]["chunk_nnz"],
        sample_cells=DB_CFG["aggregation"]["sample_cells"],
        minimum_per_type=DB_CFG["aggregation"]["minimum_per_type"],
        seed=DB_CFG["aggregation"]["seed"],
    wildcard_constraints:
        dataset=DB_DATASET_PATTERN
    resources:
        mem_mb=64000,
        runtime=10080,
        donor_bulk_io=1,
    priority: 100
    conda: "clamp-analyses"
    shell:
        "python {input.script} --config {input.config} --dataset {wildcards.dataset} "
        "--manifest {input.manifest} --counts {output.counts} "
        "--mean-cell-cpm {output.mean_cell_cpm} --truth {output.truth} "
        "--patient-info {output.info} --summary {output.summary} "
        "--sample-counts {output.sample_counts} --sample-cells {output.sample_cells} "
        "--sample-genes {output.sample_genes} --chunk-cells {params.chunk_cells} "
        "--chunk-nnz {params.chunk_nnz} --sample-cells-target {params.sample_cells} "
        "--minimum-per-type {params.minimum_per_type} --seed {params.seed}"


rule donor_bulk_data:
    input:
        counts=expand(f"{DB_PROD}/{{dataset}}/bulk/bulk_counts.csv", dataset=DB_DATASETS),
        controls=expand(f"{DB_PROD}/{{dataset}}/bulk/bulk_mean_cell_cpm.csv", dataset=DB_DATASETS),
        truth=expand(f"{DB_PROD}/{{dataset}}/bulk/truthFrac_v0.csv", dataset=DB_DATASETS),
        samples=expand(f"{DB_PROD}/{{dataset}}/single_cell_umap/sampled_counts.npz", dataset=DB_DATASETS),


rule donor_bulk_expression_umap:
    input:
        counts=rules.build_donor_bulk.output.sample_counts,
        cells=rules.build_donor_bulk.output.sample_cells,
        genes=rules.build_donor_bulk.output.sample_genes,
        script="scripts/donor_bulk/build_expression_umap.py",
    output:
        points=f"{DB_PROD}/{{dataset}}/single_cell_umap/umap_points.csv",
        summary=f"{DB_PROD}/{{dataset}}/single_cell_umap/umap_summary.csv",
    params:
        target_sum=DB_CFG["umap"]["target_sum"],
        n_hvg=DB_CFG["umap"]["n_hvg"],
        n_pcs=DB_CFG["umap"]["n_pcs"],
        n_neighbors=DB_CFG["umap"]["n_neighbors"],
        min_dist=DB_CFG["umap"]["min_dist"],
        seed=DB_CFG["umap"]["seed"],
    wildcard_constraints:
        dataset=DB_DATASET_PATTERN
    resources:
        mem_mb=64000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "python {input.script} --dataset {wildcards.dataset} --counts {input.counts} "
        "--cells {input.cells} --genes {input.genes} --output {output.points} "
        "--summary {output.summary} --target-sum {params.target_sum} "
        "--n-hvg {params.n_hvg} --n-pcs {params.n_pcs} "
        "--n-neighbors {params.n_neighbors} --min-dist {params.min_dist} "
        "--seed {params.seed}"


rule donor_bulk_umaps:
    input:
        expand(f"{DB_PROD}/{{dataset}}/single_cell_umap/umap_points.csv", dataset=DB_DATASETS)


# ============================================================
# Full-data donor-bulk models
# ============================================================

rule preprocess_donor_bulk:
    input:
        bulk=rules.build_donor_bulk.output.counts,
    output:
        norm=f"{DB_PROD}/{{dataset}}/preprocessing/norm.csv",
        cpm_filt=f"{DB_PROD}/{{dataset}}/preprocessing/cpm_filt.csv",
        row_stats=f"{DB_PROD}/{{dataset}}/preprocessing/row_stats.csv",
        k=f"{DB_PROD}/{{dataset}}/preprocessing/k.csv",
        diagnostics=f"{DB_PROD}/{{dataset}}/preprocessing/rank_diagnostics.csv",
    wildcard_constraints:
        dataset=DB_DATASET_PATTERN
    resources:
        mem_mb=32000,
        runtime=240,
    conda: "clamp-analyses"
    shell:
        "MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "Rscript scripts/pseudobulk/preprocess.R --input {input.bulk} "
        "--norm {output.norm} --cpm-filt {output.cpm_filt} "
        "--row-stats {output.row_stats} --k {output.k} "
        "--diagnostics {output.diagnostics} --mean-cutoff {config[preprocess][mean_cutoff]} "
        "--var-cutoff {config[preprocess][var_cutoff]} --seed {config[preprocess][seed]}"


rule clampfull_donor_bulk:
    input:
        norm=rules.preprocess_donor_bulk.output.norm,
        k=rules.preprocess_donor_bulk.output.k,
        gmt=rules.pathway_prior.output,
    output:
        B=f"{DB_PROD}/{{dataset}}/models/CLAMPfull/B.csv",
        Z=f"{DB_PROD}/{{dataset}}/models/CLAMPfull/Z.csv",
        model=f"{DB_PROD}/{{dataset}}/models/CLAMPfull/CLAMPfull.rds",
        l2=f"{DB_PROD}/{{dataset}}/models/CLAMPfull/L2.csv",
        summary=f"{DB_PROD}/{{dataset}}/models/CLAMPfull/summary.csv",
    params:
        base_dir=lambda wc: f"{DB_PROD}/{wc.dataset}/models/CLAMPbase_internal",
        full_dir=lambda wc: f"{DB_PROD}/{wc.dataset}/models/CLAMPfull",
    wildcard_constraints:
        dataset=DB_DATASET_PATTERN
    resources:
        mem_mb=48000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "Rscript scripts/pseudobulk/clamp.R --norm {input.norm} --k {input.k} "
        "--gmt {input.gmt} --base-dir {params.base_dir} --full-dir {params.full_dir} "
        "--model full --max-iter {config[clamp][max_iter]} --seed {config[clamp][seed]}"


rule donor_bulk_full_models:
    input:
        expand(f"{DB_PROD}/{{dataset}}/models/CLAMPfull/CLAMPfull.rds", dataset=DB_DATASETS)


# ============================================================
# Leakage-free donor-bulk and matched-control grouped CV
# ============================================================

rule donor_bulk_grouped_cv_folds:
    input:
        bulk=rules.build_donor_bulk.output.counts,
        truth=rules.build_donor_bulk.output.truth,
        info=rules.build_donor_bulk.output.info,
    output:
        membership=f"{DB_PROD}/{{dataset}}/grouped_cv/fold_membership.csv",
        summary=f"{DB_PROD}/{{dataset}}/grouped_cv/fold_summary.csv",
    wildcard_constraints:
        dataset=DB_DATASET_PATTERN
    resources:
        mem_mb=8000,
        runtime=60,
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/pseudobulk/make_grouped_cv_folds.R --dataset {wildcards.dataset} "
        "--bulk {input.bulk} --truth {input.truth} --patient-info {input.info} "
        "--output {output.membership} --summary {output.summary} "
        "--n-folds {config[grouped_cv][n_folds]} --group-col {config[grouped_cv][group_col]} "
        "--seed {config[grouped_cv][split_seed]}"


rule clampfull_grouped_cv_donor_bulk:
    input:
        bulk=rules.build_donor_bulk.output.counts,
        truth=rules.build_donor_bulk.output.truth,
        membership=rules.donor_bulk_grouped_cv_folds.output.membership,
        gmt=rules.pathway_prior.output,
    output:
        train_B=f"{DB_PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/train_B.csv",
        train_Z=f"{DB_PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/train_Z.csv",
        test_B=f"{DB_PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/test_B.csv",
        row_stats=f"{DB_PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/row_stats.csv",
        train_truth=f"{DB_PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/train_truth.csv",
        test_truth=f"{DB_PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/test_truth.csv",
        summary=f"{DB_PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/summary.csv",
        model=f"{DB_PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/CLAMPfull.rds",
    params:
        out_dir=lambda wc: f"{DB_PROD}/{wc.dataset}/grouped_cv/fold{wc.fold}/CLAMPfull",
    wildcard_constraints:
        dataset=DB_DATASET_PATTERN,
        fold=CV_FOLD_PATTERN,
    resources:
        mem_mb=48000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "Rscript scripts/pseudobulk/clampfull_grouped_cv.R --dataset {wildcards.dataset} "
        "--fold {wildcards.fold} --bulk {input.bulk} --truth {input.truth} "
        "--membership {input.membership} --gmt {input.gmt} --out-dir {params.out_dir} "
        "--seed {config[clamp][seed]} --max-iter {config[clamp][max_iter]} "
        "--mean-cutoff {config[preprocess][mean_cutoff]} --var-cutoff {config[preprocess][var_cutoff]}"


rule clampfull_grouped_cv_donor_bulk_control:
    input:
        bulk=rules.build_donor_bulk.output.mean_cell_cpm,
        truth=rules.build_donor_bulk.output.truth,
        membership=rules.donor_bulk_grouped_cv_folds.output.membership,
        gmt=rules.pathway_prior.output,
    output:
        train_B=f"{DB_PROD}/{{dataset}}/grouped_cv_mean_cell_cpm/fold{{fold}}/CLAMPfull/train_B.csv",
        train_Z=f"{DB_PROD}/{{dataset}}/grouped_cv_mean_cell_cpm/fold{{fold}}/CLAMPfull/train_Z.csv",
        test_B=f"{DB_PROD}/{{dataset}}/grouped_cv_mean_cell_cpm/fold{{fold}}/CLAMPfull/test_B.csv",
        row_stats=f"{DB_PROD}/{{dataset}}/grouped_cv_mean_cell_cpm/fold{{fold}}/CLAMPfull/row_stats.csv",
        train_truth=f"{DB_PROD}/{{dataset}}/grouped_cv_mean_cell_cpm/fold{{fold}}/CLAMPfull/train_truth.csv",
        test_truth=f"{DB_PROD}/{{dataset}}/grouped_cv_mean_cell_cpm/fold{{fold}}/CLAMPfull/test_truth.csv",
        summary=f"{DB_PROD}/{{dataset}}/grouped_cv_mean_cell_cpm/fold{{fold}}/CLAMPfull/summary.csv",
        model=f"{DB_PROD}/{{dataset}}/grouped_cv_mean_cell_cpm/fold{{fold}}/CLAMPfull/CLAMPfull.rds",
    params:
        out_dir=lambda wc: f"{DB_PROD}/{wc.dataset}/grouped_cv_mean_cell_cpm/fold{wc.fold}/CLAMPfull",
    wildcard_constraints:
        dataset=DB_DATASET_PATTERN,
        fold=CV_FOLD_PATTERN,
    resources:
        mem_mb=48000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
        "Rscript scripts/pseudobulk/clampfull_grouped_cv.R --dataset {wildcards.dataset} "
        "--fold {wildcards.fold} --bulk {input.bulk} --truth {input.truth} "
        "--membership {input.membership} --gmt {input.gmt} --out-dir {params.out_dir} "
        "--seed {config[clamp][seed]} --max-iter {config[clamp][max_iter]} "
        "--mean-cutoff {config[preprocess][mean_cutoff]} --var-cutoff {config[preprocess][var_cutoff]}"


rule donor_bulk_grouped_cv_models:
    input:
        donor=expand(
            f"{DB_PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/CLAMPfull.rds",
            dataset=DB_DATASETS, fold=CV_FOLDS,
        ),
        control=expand(
            f"{DB_PROD}/{{dataset}}/grouped_cv_mean_cell_cpm/fold{{fold}}/CLAMPfull/CLAMPfull.rds",
            dataset=DB_DATASETS, fold=CV_FOLDS,
        )


rule donor_bulk_models:
    input:
        full=rules.donor_bulk_full_models.input,
        grouped_cv=expand(
            f"{DB_PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/CLAMPfull.rds",
            dataset=DB_DATASETS, fold=CV_FOLDS,
        )


rule donor_bulk_controls:
    input:
        expand(
            f"{DB_PROD}/{{dataset}}/grouped_cv_mean_cell_cpm/fold{{fold}}/CLAMPfull/CLAMPfull.rds",
            dataset=DB_DATASETS, fold=CV_FOLDS,
        )


# ============================================================
# Original-cell projection, analysis, validation, and notebooks
# ============================================================

rule single_cell_projection_donor_bulk:
    input:
        raw=donor_bulk_raw_inputs,
        z=rules.clampfull_donor_bulk.output.Z,
        l2=rules.clampfull_donor_bulk.output.l2,
        row_stats=rules.preprocess_donor_bulk.output.row_stats,
        config="workflow/config/pseudobulk.yaml",
    output:
        scores=f"{DB_PROD}/{{dataset}}/single_cell_projection/single_cell_lv_scores.h5",
        summary=f"{DB_PROD}/{{dataset}}/single_cell_projection/projection_summary.csv",
    params:
        chunk_cells=DB_CFG["projection"]["chunk_cells"],
        chunk_nnz=DB_CFG["projection"]["chunk_nnz"],
        duplicate_policy=DB_CFG["projection"]["duplicate_gene_policy"],
    wildcard_constraints:
        dataset=DB_DATASET_PATTERN
    resources:
        mem_mb=64000,
        runtime=4320,
        donor_bulk_io=1,
    conda: "clamp-analyses"
    shell:
        "python scripts/pseudobulk/single_cell_projection.py --config {input.config} "
        "--dataset {wildcards.dataset} --z {input.z} --l2 {input.l2} "
        "--row-stats {input.row_stats} --output {output.scores} --summary {output.summary} "
        "--chunk-cells {params.chunk_cells} --chunk-nnz {params.chunk_nnz} "
        "--duplicate-gene-policy {params.duplicate_policy} "
        "--normalization-label donor_bulk_training_statistics"


rule donor_bulk_single_cell_projections:
    input:
        expand(f"{DB_PROD}/{{dataset}}/single_cell_projection/single_cell_lv_scores.h5", dataset=DB_DATASETS)


rule donor_bulk_analysis_tables:
    input:
        full_models=expand(f"{DB_PROD}/{{dataset}}/models/CLAMPfull/CLAMPfull.rds", dataset=DB_DATASETS),
        donor_cv=expand(
            f"{DB_PROD}/{{dataset}}/grouped_cv/fold{{fold}}/CLAMPfull/CLAMPfull.rds",
            dataset=DB_DATASETS, fold=CV_FOLDS,
        ),
        control_cv=expand(
            f"{DB_PROD}/{{dataset}}/grouped_cv_mean_cell_cpm/fold{{fold}}/CLAMPfull/CLAMPfull.rds",
            dataset=DB_DATASETS, fold=CV_FOLDS,
        ),
        projections=expand(f"{DB_PROD}/{{dataset}}/single_cell_projection/single_cell_lv_scores.h5", dataset=DB_DATASETS),
        script="scripts/donor_bulk/analyze_donor_bulk.py",
    output:
        assignments=f"{DB_BIO}/analysis/full_lv_assignments.csv",
        correlations=f"{DB_BIO}/analysis/full_lv_celltype_correlations.csv",
        full_summary=f"{DB_BIO}/analysis/full_recovery_by_dataset.csv",
        predictions=f"{DB_BIO}/analysis/oof_predictions.csv",
        calibrations=f"{DB_BIO}/analysis/fold_calibrations.csv",
        metrics=f"{DB_BIO}/analysis/oof_metrics.csv",
        cv_dataset=f"{DB_BIO}/analysis/cv_recovery_by_dataset.csv",
        cv_overall=f"{DB_BIO}/analysis/cv_recovery_overall.csv",
        control_predictions=f"{DB_BIO}/analysis/control_oof_predictions.csv",
        control_calibrations=f"{DB_BIO}/analysis/control_fold_calibrations.csv",
        control_metrics=f"{DB_BIO}/analysis/control_oof_metrics.csv",
        control_dataset=f"{DB_BIO}/analysis/control_cv_recovery_by_dataset.csv",
        control_overall=f"{DB_BIO}/analysis/control_cv_recovery_overall.csv",
        comparison=f"{DB_BIO}/analysis/comparison_cv_samecell.csv",
        recovery=f"{DB_BIO}/analysis/single_cell_recovery.csv",
        recovery_dataset=f"{DB_BIO}/analysis/single_cell_recovery_by_dataset.csv",
        recovery_overall=f"{DB_BIO}/analysis/single_cell_recovery_overall.csv",
        specificity=f"{DB_BIO}/analysis/single_cell_specificity_matrix.csv",
        specificity_summary=f"{DB_BIO}/analysis/single_cell_specificity_summary.csv",
    resources:
        mem_mb=64000,
        runtime=1440,
    params:
        datasets=" ".join(DB_DATASETS),
    conda: "clamp-analyses"
    shell:
        "python {input.script} --production-root {DB_PROD} --output-dir {DB_BIO}/analysis "
        "--datasets {params.datasets}"


rule validate_donor_bulk:
    input:
        analysis=rules.donor_bulk_analysis_tables.output,
        umaps=expand(f"{DB_PROD}/{{dataset}}/single_cell_umap/umap_points.csv", dataset=DB_DATASETS),
        script="scripts/donor_bulk/validate_donor_bulk.py",
    output:
        summary=f"{DB_PROD}/qc/validation_summary.csv",
        json=f"{DB_PROD}/qc/validation.json",
    resources:
        mem_mb=16000,
        runtime=240,
    params:
        datasets=" ".join(DB_DATASETS),
        sample_target=DB_CFG["aggregation"]["sample_cells"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} --production-root {DB_PROD} --analysis-dir {DB_BIO}/analysis "
        "--output-summary {output.summary} --output-json {output.json} "
        "--sample-target {params.sample_target} --datasets {params.datasets}"


rule donor_bulk_model_qc:
    input:
        validation=rules.validate_donor_bulk.output.summary,
        aggregation=expand(f"{DB_PROD}/{{dataset}}/bulk/aggregation_summary.csv", dataset=DB_DATASETS),
        ranks=expand(f"{DB_PROD}/{{dataset}}/preprocessing/rank_diagnostics.csv", dataset=DB_DATASETS),
        projections=expand(f"{DB_PROD}/{{dataset}}/single_cell_projection/projection_summary.csv", dataset=DB_DATASETS),
        umaps=expand(f"{DB_PROD}/{{dataset}}/single_cell_umap/umap_summary.csv", dataset=DB_DATASETS),
        notebook=f"{DB_MODEL_NB}/01_donor_bulk_qc.ipynb",
    output:
        summary=f"{DB_PROD}/qc/donor_bulk_qc.csv",
        complete=touch(f"{DB_PROD}/qc/notebook.complete"),
    log:
        notebook=f"{DB_MODEL_NB}/01_donor_bulk_qc.executed.ipynb",
    params:
        datasets=DB_DATASETS,
        out_dir=f"{DB_PROD}/qc",
    conda: "clamp-analyses"
    notebook:
        f"{DB_MODEL_NB}/01_donor_bulk_qc.ipynb"


rule donor_bulk_recovery_report:
    input:
        assignments=rules.donor_bulk_analysis_tables.output.assignments,
        predictions=rules.donor_bulk_analysis_tables.output.predictions,
        metrics=rules.donor_bulk_analysis_tables.output.metrics,
        cv_dataset=rules.donor_bulk_analysis_tables.output.cv_dataset,
        comparison=rules.donor_bulk_analysis_tables.output.comparison,
        recovery=rules.donor_bulk_analysis_tables.output.recovery,
        specificity=rules.donor_bulk_analysis_tables.output.specificity,
        specificity_summary=rules.donor_bulk_analysis_tables.output.specificity_summary,
        validation=rules.validate_donor_bulk.output.summary,
        umaps=expand(f"{DB_PROD}/{{dataset}}/single_cell_umap/umap_points.csv", dataset=DB_DATASETS),
        projections=expand(f"{DB_PROD}/{{dataset}}/single_cell_projection/single_cell_lv_scores.h5", dataset=DB_DATASETS),
        notebook=f"{DB_BIO_NB}/06_donor_bulk_recovery.ipynb",
    output:
        fig2_predictions=f"{DB_BIO}/fig2_oof_predictions.csv",
        fig2_statistics=f"{DB_BIO}/fig2_oof_statistics.csv",
        perez_annotation=f"{DB_BIO}/fig2_perez_annotation.csv",
        perez_activity=f"{DB_BIO}/fig2_perez_activity.csv",
        perez_summary=f"{DB_BIO}/fig2_perez_summary.csv",
        complete=touch(f"{DB_BIO}/notebook.complete"),
    log:
        notebook=f"{DB_BIO_NB}/06_donor_bulk_recovery.executed.ipynb",
    params:
        datasets=DB_DATASETS,
        production_root=DB_PROD,
        out_dir=DB_BIO,
    conda: "clamp-analyses"
    notebook:
        f"{DB_BIO_NB}/06_donor_bulk_recovery.ipynb"


rule donor_bulk_report:
    input:
        rules.donor_bulk_model_qc.output.complete,
        rules.donor_bulk_recovery_report.output.complete,


rule donor_bulk_qc:
    input:
        rules.donor_bulk_model_qc.output.complete,


rule donor_bulk_biology:
    input:
        rules.donor_bulk_recovery_report.output.complete,
