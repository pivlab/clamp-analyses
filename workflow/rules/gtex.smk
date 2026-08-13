GTEX_PROD = config["gtex"]["paths"]["production"]
GTEX_BIO = config["gtex"]["paths"]["biology"]
GTEX_MODEL_NB = os.path.join(REPO_ROOT, config["gtex"]["paths"]["model_notebooks"])
GTEX_BIO_NB = os.path.join(REPO_ROOT, config["gtex"]["paths"]["biology_notebooks"])
GTEX_GMT = config["references"]["go_bp_file"]

import shlex

SUBTISSUE_CV_FOLDS = list(range(1, config["gtex"]["subtissue_inference"]["n_folds"] + 1))
SUBTISSUE_CV_FOLD_PATTERN = "|".join(str(f) for f in SUBTISSUE_CV_FOLDS)
SUBTISSUE_TISSUES_ARG = " ".join(shlex.quote(t) for t in config["gtex"]["subtissue_inference"]["tissues"])


# ============================================================
# Step 1: download + build the base FBM, SVD, CLAMPbase/full
# ============================================================

rule download_gtex_raw:
    output:
        config["gtex"]["raw_gz"]
    params:
        url=config["gtex"]["raw_url"],
    conda: "clamp-analyses"
    shell:
        "mkdir -p $(dirname {output}) && "
        "curl --fail --location --retry 3 '{params.url}' --output {output}.download && "
        "mv {output}.download {output}"


rule clamp_gtex:
    input:
        raw_gct=rules.download_gtex_raw.output,
        gmt=rules.pathway_prior.output,
    output:
        samples=f"{GTEX_PROD}/gtex_samples.rds",
        genes=f"{GTEX_PROD}/gtex_genes.rds",
        fbm_filt=f"{GTEX_PROD}/gtex_fbm_filt.rds",
        df_rds=f"{GTEX_PROD}/df_gtex_fbm_filt.rds",
        df_csv=f"{GTEX_PROD}/df_gtex_fbm_filt.csv",
        svd_res=f"{GTEX_PROD}/gtex_svdRes.rds",
        k_rds=f"{GTEX_PROD}/CLAMP_K_gtex.rds",
        k_csv=f"{GTEX_PROD}/CLAMP_K_gtex.csv",
        base_rds=f"{GTEX_PROD}/CLAMPbase.rds",
        base_B=f"{GTEX_PROD}/CLAMPbase/B.csv",
        base_Z=f"{GTEX_PROD}/CLAMPbase/Z.csv",
        full_rds=f"{GTEX_PROD}/CLAMPfull.rds",
        full_B=f"{GTEX_PROD}/CLAMPfull/B.csv",
        full_Z=f"{GTEX_PROD}/CLAMPfull/Z.csv",
        full_summary=f"{GTEX_PROD}/CLAMPfull/summary.csv",
    params:
        mean_cutoff=config["gtex"]["preprocess"]["mean_cutoff"],
        var_cutoff=config["gtex"]["preprocess"]["var_cutoff"],
        chunk_size=config["gtex"]["raw_data"]["chunk_size"],
        n_cores=config["gtex"]["raw_data"]["n_cores"],
        max_iter=config["gtex"]["clamp"]["max_iter"],
        seed=config["gtex"]["clamp"]["seed"],
    resources:
        mem_mb=96000,
        runtime=2880,
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/gtex/clamp.R --raw-gct {input.raw_gct} --gmt {input.gmt} "
        "--out-dir {GTEX_PROD} --mean-cutoff {params.mean_cutoff} --var-cutoff {params.var_cutoff} "
        "--chunk-size {params.chunk_size} --n-cores {params.n_cores} "
        "--max-iter {params.max_iter} --seed {params.seed}"


rule plier_gtex:
    input:
        fbm_filt=rules.clamp_gtex.output.fbm_filt,
        svd_res=rules.clamp_gtex.output.svd_res,
        k=rules.clamp_gtex.output.k_rds,
        genes=rules.clamp_gtex.output.genes,
        samples=rules.clamp_gtex.output.samples,
        gmt=rules.pathway_prior.output,
    output:
        rds=f"{GTEX_PROD}/PLIER.rds",
        B=f"{GTEX_PROD}/PLIER/B.csv",
        Z=f"{GTEX_PROD}/PLIER/Z.csv",
        summary=f"{GTEX_PROD}/PLIER/summary.csv",
    params:
        seed=config["gtex"]["preprocess"]["seed"],
    resources: mem_mb=64000, runtime=1440
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/gtex/plier.R --fbm-filt {input.fbm_filt} --svd-res {input.svd_res} "
        "--k {input.k} --genes {input.genes} --samples {input.samples} --gmt {input.gmt} "
        "--out-dir {GTEX_PROD} --seed {params.seed}"


rule pca_nmf_ica_gtex:
    input:
        df_gtex_fbm_filt=rules.clamp_gtex.output.df_csv,
        k=rules.clamp_gtex.output.k_csv,
    output:
        pca_B=f"{GTEX_PROD}/PCA/gtex_pca_B.pkl",
        pca_scores=f"{GTEX_PROD}/PCA/gtex_pca_scores.pkl",
        pca_loadings=f"{GTEX_PROD}/PCA/gtex_pca_loadings.pkl",
        ica_B=f"{GTEX_PROD}/ICA/gtex_ica_B.pkl",
        ica_scores=f"{GTEX_PROD}/ICA/gtex_ica_scores.pkl",
        ica_loadings=f"{GTEX_PROD}/ICA/gtex_ica_loadings.pkl",
        nmf_B=f"{GTEX_PROD}/NMF/gtex_nmf_B.pkl",
        nmf_scores=f"{GTEX_PROD}/NMF/gtex_nmf_scores.pkl",
        nmf_loadings=f"{GTEX_PROD}/NMF/gtex_nmf_loadings.pkl",
    params:
        seed=config["gtex"]["preprocess"]["seed"],
    resources: mem_mb=64000, runtime=720
    conda: "clamp-analyses"
    shell:
        "python scripts/gtex/pca_nmf_ica.py --df-gtex-fbm-filt {input.df_gtex_fbm_filt} --k {input.k} "
        "--out-dir {GTEX_PROD} --seed {params.seed}"


rule flashier_gtex:
    input:
        df_gtex_fbm_filt=rules.clamp_gtex.output.df_rds,
        k=rules.clamp_gtex.output.k_rds,
    output:
        B=f"{GTEX_PROD}/flashier/gtex_B.csv",
        rds=f"{GTEX_PROD}/flashier/flashier_model.rds",
    params:
        backfit_maxiter=config["gtex"]["flashier"]["backfit_maxiter"],
        seed=config["gtex"]["preprocess"]["seed"],
    resources: mem_mb=64000, runtime=1440
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/gtex/flashier.R --df-gtex-fbm-filt {input.df_gtex_fbm_filt} --k {input.k} "
        "--out-dir {GTEX_PROD}/flashier --backfit-maxiter {params.backfit_maxiter} --seed {params.seed}"


rule mofa_flex_prior_gtex:
    input:
        df_gtex_fbm_filt=rules.clamp_gtex.output.df_csv,
        k=rules.clamp_gtex.output.k_csv,
        gmt=rules.pathway_prior.output,
    output:
        B=f"{GTEX_PROD}/MOFA_FLEX_PRIOR/B_matrix.csv",
        Z=f"{GTEX_PROD}/MOFA_FLEX_PRIOR/Z_matrix.csv",
        model=f"{GTEX_PROD}/MOFA_FLEX_PRIOR/model.pkl",
    params:
        min_fraction=config["gtex"]["mofa_flex_prior"]["min_fraction"],
        min_count=config["gtex"]["mofa_flex_prior"]["min_count"],
        max_count=config["gtex"]["mofa_flex_prior"]["max_count"],
        similarity_threshold=config["gtex"]["mofa_flex_prior"]["similarity_threshold"],
        seed=config["gtex"]["mofa_flex_prior"]["seed"],
        max_epochs=config["gtex"]["mofa_flex_prior"]["max_epochs"],
    resources: mem_mb=64000, runtime=1440
    conda: "clamp-analyses"
    shell:
        "python scripts/gtex/mofa_flex_prior.py --df-gtex-fbm-filt {input.df_gtex_fbm_filt} --k {input.k} "
        "--out-dir {GTEX_PROD}/MOFA_FLEX_PRIOR --gmt {input.gmt} "
        "--min-fraction {params.min_fraction} "
        "--min-count {params.min_count} --max-count {params.max_count} "
        "--similarity-threshold {params.similarity_threshold} --seed {params.seed} "
        "--max-epochs {params.max_epochs}"


rule cogaps_gtex:
    # Distributed genome-wide CoGAPS is slow enough on this dataset that convergence
    # isn't guaranteed within any fixed budget. Rather than let a stalled run fail the
    # whole job (or hang indefinitely), it gets a hard wall-clock budget
    # (cogaps.timeout_days, default 7d): on convergence the usual gtex_B.csv/
    # cogaps_model.rds are written; on timeout NOT_CONVERGED.txt is written instead so
    # the rule still completes and the non-convergence is visible on disk rather than
    # only in a log. A genuine R/CoGAPS error (any other non-zero exit) still fails the
    # rule normally.
    input:
        df_gtex_fbm_filt=rules.clamp_gtex.output.df_rds,
        k=rules.clamp_gtex.output.k_rds,
    output:
        out_dir=directory(f"{GTEX_PROD}/CoGAPS"),
    params:
        n_iterations=config["gtex"]["cogaps"]["n_iterations"],
        n_threads=config["gtex"]["cogaps"]["n_threads"],
        seed=config["gtex"]["cogaps"]["seed"],
        timeout_days=config["gtex"]["cogaps"]["timeout_days"],
    resources: mem_mb=64000, runtime=7 * 24 * 60 + 30
    conda: "clamp-analyses"
    shell:
        """
        mkdir -p {output.out_dir}
        rm -f {output.out_dir}/NOT_CONVERGED.txt
        set +e
        timeout -k 60s {params.timeout_days}d Rscript scripts/gtex/cogaps.R \
            --df-gtex-fbm-filt {input.df_gtex_fbm_filt} --k {input.k} \
            --out-dir {output.out_dir} --n-iterations {params.n_iterations} \
            --n-threads {params.n_threads} --seed {params.seed}
        ec=$?
        set -e
        if [ $ec -eq 124 ] || [ $ec -eq 137 ]; then
            echo "CoGAPS did not converge within the {params.timeout_days}-day time budget (killed $(date -Iseconds))" > {output.out_dir}/NOT_CONVERGED.txt
        elif [ $ec -ne 0 ]; then
            exit $ec
        fi
        """


rule gss_gtex:
    input:
        df_gtex_fbm_filt=rules.clamp_gtex.output.df_rds,
        k=rules.clamp_gtex.output.k_rds,
    output:
        B=f"{GTEX_PROD}/GSS/gtex_B.csv",
    params:
        d_cluster=config["gtex"]["gss"]["d_cluster"],
        seed=config["gtex"]["gss"]["seed"],
    resources: mem_mb=64000, runtime=720
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/gtex/gss.R --dataset GTEx --df-gtex-fbm-filt {input.df_gtex_fbm_filt} "
        "--k {input.k} --d-cluster {params.d_cluster} --out-dir {GTEX_PROD}/GSS --seed {params.seed}"


rule full_models_gtex:
    input:
        rules.clamp_gtex.output.full_B,
        rules.plier_gtex.output.B,
        rules.pca_nmf_ica_gtex.output.pca_B,
        rules.flashier_gtex.output.B,
        rules.mofa_flex_prior_gtex.output.B,
        rules.gss_gtex.output.B,
        # cogaps_gtex intentionally excluded: distributed genome-wide CoGAPS with
        # nIterations=5000 does not converge in practical time on this dataset.
        # Rule kept and fully runnable on demand (`snakemake cogaps_gtex`), just
        # not part of the default pipeline.


rule model_building_qc_gtex:
    input:
        models=rules.full_models_gtex.input,
        k=rules.clamp_gtex.output.k_csv,
        df_csv=rules.clamp_gtex.output.df_csv,
        notebook=f"{GTEX_MODEL_NB}/00_model_building_qc.ipynb",
    output:
        summary=f"{GTEX_PROD}/qc/model_building_qc.csv",
        matrices=f"{GTEX_PROD}/qc/model_matrix_qc.csv",
        figure=f"{GTEX_PROD}/qc/model_building_qc.png",
        complete=touch(f"{GTEX_PROD}/qc/notebook.complete"),
    log:
        notebook=f"{GTEX_MODEL_NB}/00_model_building_qc.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{GTEX_MODEL_NB}/00_model_building_qc.ipynb"


# ============================================================
# Step 2: GPU k-means clustering ensemble
# ============================================================

rule kmeans_clustering_gtex:
    input:
        gtex_meta=config["gtex"]["metadata"],
        df_gtex_fbm_filt=rules.clamp_gtex.output.df_csv,
        clamp_base_rds=rules.clamp_gtex.output.base_rds,
        clamp_full_rds=rules.clamp_gtex.output.full_rds,
        plier_rds=rules.plier_gtex.output.rds,
        gss_b=rules.gss_gtex.output.B,
        flashier_b=rules.flashier_gtex.output.B,
        mofa_flex_prior_b=rules.mofa_flex_prior_gtex.output.B,
        pca_b=rules.pca_nmf_ica_gtex.output.pca_B,
        ica_b=rules.pca_nmf_ica_gtex.output.ica_B,
        nmf_b=rules.pca_nmf_ica_gtex.output.nmf_B,
    output:
        cache=directory(f"{GTEX_BIO}/00_kmeans_clustering/kmeans_results"),
        model=f"{GTEX_BIO}/00_kmeans_clustering/kmeans_models/gtex_CLAMPfull_kmeans_model.pkl",
    params:
        min_samples=config["gtex"]["kmeans_clustering"]["ari_min_samples"],
        base_seed=config["gtex"]["kmeans_clustering"]["base_seed"],
        n_reps_per_k=config["gtex"]["kmeans_clustering"]["n_reps_per_k"],
        gene_fractions=" ".join(str(x) for x in config["gtex"]["kmeans_clustering"]["gene_fractions"]),
    resources:
        mem_mb=128000,
        runtime=2880,
    conda: "gpu-kmeans"
    shell:
        "python scripts/gtex/kmeans_clustering.py --gtex-meta {input.gtex_meta} "
        "--df-gtex-fbm-filt {input.df_gtex_fbm_filt} --clamp-base-rds {input.clamp_base_rds} "
        "--clamp-full-rds {input.clamp_full_rds} --plier-rds {input.plier_rds} "
        "--gss-b {input.gss_b} --flashier-b {input.flashier_b} "
        "--mofa-flex-prior-b {input.mofa_flex_prior_b} --pca-b {input.pca_b} "
        "--ica-b {input.ica_b} --nmf-b {input.nmf_b} --min-samples {params.min_samples} "
        "--base-seed {params.base_seed} --n-reps-per-k {params.n_reps_per_k} "
        "--gene-fractions {params.gene_fractions} --cache-dir {output.cache} "
        "--model-out-dir {GTEX_BIO}/00_kmeans_clustering/kmeans_models"


rule kmeans_clustering_report_gtex:
    input:
        cache=rules.kmeans_clustering_gtex.output.cache,
        model=rules.kmeans_clustering_gtex.output.model,
        notebook=f"{GTEX_BIO_NB}/00_kmeans_clustering.ipynb",
    output:
        ari_data=f"{GTEX_BIO}/00_kmeans_clustering/ari_data.csv",
        ari_comparisons=f"{GTEX_BIO}/00_kmeans_clustering/ari_comparisons.csv",
        gene_fraction_ari_data=f"{GTEX_BIO}/00_kmeans_clustering/gene_fraction_ari_data.csv",
        gene_fraction_ari_comparisons=f"{GTEX_BIO}/00_kmeans_clustering/gene_fraction_ari_comparisons.csv",
        complete=touch(f"{GTEX_BIO}/00_kmeans_clustering/notebook.complete"),
    log:
        notebook=f"{GTEX_BIO_NB}/00_kmeans_clustering.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{GTEX_BIO_NB}/00_kmeans_clustering.ipynb"


# ============================================================
# Step 3: LV importance (RF + SHAP), two label sources
# ============================================================

rule lv_importance_rf_true_labels_gtex:
    input:
        clamp_full_rds=rules.clamp_gtex.output.full_rds,
        gtex_meta=config["gtex"]["metadata"],
    output:
        out_dir=directory(f"{GTEX_BIO}/01_LV_importance_rf_true_labels/gtex_feature_importance_true_labels_binary_shap"),
    params:
        min_samples=config["gtex"]["kmeans_clustering"]["min_samples"],
        global_seed=config["gtex"]["lv_importance"]["global_seed"],
        min_rf_accuracy=config["gtex"]["lv_importance"]["min_rf_accuracy"],
    resources:
        mem_mb=64000,
        runtime=2880,
    conda: "clamp-analyses"
    shell:
        "python scripts/gtex/lv_importance_true_labels.py --clamp-full-rds {input.clamp_full_rds} "
        "--gtex-meta {input.gtex_meta} --min-samples {params.min_samples} "
        "--global-seed {params.global_seed} --min-rf-accuracy {params.min_rf_accuracy} "
        "--out-dir {output.out_dir}"


rule lv_importance_rf_true_labels_report_gtex:
    input:
        shap_dir=rules.lv_importance_rf_true_labels_gtex.output.out_dir,
        notebook=f"{GTEX_BIO_NB}/01_LV_importance.ipynb",
    output:
        complete=touch(f"{GTEX_BIO}/01_LV_importance_rf_true_labels/notebook.complete"),
    log:
        notebook=f"{GTEX_BIO_NB}/01_LV_importance.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{GTEX_BIO_NB}/01_LV_importance.ipynb"


# ============================================================
# Step 4: biology / ORA / plotting notebooks (small analyses,
# kept as notebooks per the "nb are for plotting" convention)
# ============================================================

rule lv_importance_rf_true_labels_biology_gtex:
    input:
        shap_dir=rules.lv_importance_rf_true_labels_gtex.output.out_dir,
        notebook=f"{GTEX_BIO_NB}/02_b_matrix.ipynb",
    output:
        tissue_concordance=f"{GTEX_BIO}/02_LV_importance_rf_true_labels_biology/tissue_concordance.csv",
        tissue_group_heatmap=f"{GTEX_BIO}/02_LV_importance_rf_true_labels_biology/tissue_group_heatmap.csv",
        tissue_subtissue_heatmap=f"{GTEX_BIO}/02_LV_importance_rf_true_labels_biology/tissue_subtissue_heatmap.csv",
        complete=touch(f"{GTEX_BIO}/02_LV_importance_rf_true_labels_biology/notebook.complete"),
    log:
        notebook=f"{GTEX_BIO_NB}/02_b_matrix.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{GTEX_BIO_NB}/02_b_matrix.ipynb"


rule global_alignment_rf_true_labels_gtex:
    input:
        shap_dir=rules.lv_importance_rf_true_labels_gtex.output.out_dir,
        notebook=f"{GTEX_BIO_NB}/03_global_alignment.ipynb",
    output:
        grouped_summary=f"{GTEX_BIO}/07_global_alignment_rf_true_labels/gtex_global_alignment_grouped_summary.csv",
        subtissue_summary=f"{GTEX_BIO}/07_global_alignment_rf_true_labels/cumulative20/gtex_global_alignment_summary.csv",
        complete=touch(f"{GTEX_BIO}/07_global_alignment_rf_true_labels/notebook.complete"),
    log:
        notebook=f"{GTEX_BIO_NB}/03_global_alignment.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{GTEX_BIO_NB}/03_global_alignment.ipynb"


rule liver_disentangle_xcell_rf_true_labels_gtex:
    input:
        shap_dir=rules.lv_importance_rf_true_labels_gtex.output.out_dir,
        notebook=f"{GTEX_BIO_NB}/05_liver_disentangle_xcell.ipynb",
    output:
        xcell_scatter=f"{GTEX_BIO}/04_liver_disentangle_xcell_rf_true_labels/liver_xcell_scatter.csv",
        lv_pathways=f"{GTEX_BIO}/04_liver_disentangle_xcell_rf_true_labels/liver_lv_pathways.csv",
        complete=touch(f"{GTEX_BIO}/04_liver_disentangle_xcell_rf_true_labels/notebook.complete"),
    log:
        notebook=f"{GTEX_BIO_NB}/05_liver_disentangle_xcell.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{GTEX_BIO_NB}/05_liver_disentangle_xcell.ipynb"


rule multitissue_xcell_recovery_gtex:
    input:
        shap_dir=rules.lv_importance_rf_true_labels_gtex.output.out_dir,
        clamp_model=rules.clamp_gtex.output.full_rds,
        z_matrix=rules.clamp_gtex.output.full_Z,
        xcell=config["gtex"]["xcell_scores"],
        cell_marker_file=config["references"]["cell_marker_file"],
        notebook=f"{GTEX_BIO_NB}/06_multitissue_xcell_recovery.ipynb",
    output:
        panel_ready=f"{GTEX_BIO}/06_multitissue_xcell_recovery/multitissue_recovery_panel_ready.csv",
        shap_by_model=f"{GTEX_BIO}/06_multitissue_xcell_recovery/multitissue_shap_by_model.csv",
        ora_top1=f"{GTEX_BIO}/06_multitissue_xcell_recovery/cellmarker_ora_top1_per_lv.csv",
        summary=f"{GTEX_BIO}/06_multitissue_xcell_recovery/multitissue_recovery_summary.csv",
        complete=touch(f"{GTEX_BIO}/06_multitissue_xcell_recovery/notebook.complete"),
    log:
        notebook=f"{GTEX_BIO_NB}/06_multitissue_xcell_recovery.executed.ipynb",
    params:
        out_dir=f"{GTEX_BIO}/06_multitissue_xcell_recovery",
    conda: "clamp-analyses"
    notebook:
        f"{GTEX_BIO_NB}/06_multitissue_xcell_recovery.ipynb"


# Measures the shared GO:BP prior rather than anything GTEx-specific; it lives in
# this rule file because the GTEx supplement is what consumes the output.
rule geneset_orthogonality_gtex:
    input:
        go_bp_file=GTEX_GMT,
        cell_marker_file=config["references"]["cell_marker_file"],
        allen_brain_gmt_file=config["references"]["allen_brain_gmt_file"],
        azimuth_file=config["references"]["azimuth_file"],
        pseudobulk_marker_recovery=rules.disentangle_pseudobulk.output.recovery,
        gtex_tissues_pathmat=config["references"]["gtex_tissues_pathmat"],
        notebook=f"{GTEX_BIO_NB}/07_geneset_orthogonality.ipynb",
    output:
        per_term=f"{GTEX_BIO}/07_geneset_orthogonality/orthogonality_per_term.csv",
        summary=f"{GTEX_BIO}/07_geneset_orthogonality/orthogonality_summary.csv",
        complete=touch(f"{GTEX_BIO}/07_geneset_orthogonality/notebook.complete"),
    log:
        notebook=f"{GTEX_BIO_NB}/07_geneset_orthogonality.executed.ipynb",
    params:
        out_dir=f"{GTEX_BIO}/07_geneset_orthogonality",
    conda: "clamp-analyses"
    notebook:
        f"{GTEX_BIO_NB}/07_geneset_orthogonality.ipynb"


# ============================================================
# Step 5: subtissue recovery from out-of-fold RF/SHAP (independent
# of the true-labels chain above: this RF is always trained blind
# to subtissue, on broad SMTS only, regardless of Steps 3-4's label scheme)
# ============================================================

rule subtissue_cv_folds_gtex:
    input:
        clamp_full_rds=rules.clamp_gtex.output.full_rds,
        gtex_meta=config["gtex"]["metadata"],
    output:
        membership=f"{GTEX_BIO}/04_subtissues/fold_membership.tsv",
        summary=f"{GTEX_BIO}/04_subtissues/fold_summary.tsv",
    params:
        min_samples=config["gtex"]["kmeans_clustering"]["min_samples"],
        n_folds=config["gtex"]["subtissue_inference"]["n_folds"],
        global_seed=config["gtex"]["subtissue_inference"]["global_seed"],
    resources:
        mem_mb=16000,
        runtime=60,
    conda: "clamp-analyses"
    shell:
        "python scripts/gtex/subtissue_cv_folds.py --clamp-full-rds {input.clamp_full_rds} "
        "--gtex-meta {input.gtex_meta} --min-samples {params.min_samples} "
        "--n-folds {params.n_folds} --global-seed {params.global_seed} "
        "--out-membership {output.membership} --out-summary {output.summary}"


rule subtissue_fold_rf_shap_gtex:
    input:
        clamp_full_rds=rules.clamp_gtex.output.full_rds,
        gtex_meta=config["gtex"]["metadata"],
        membership=rules.subtissue_cv_folds_gtex.output.membership,
    output:
        oof_shap=f"{GTEX_BIO}/04_subtissues/fold{{fold}}/oof_shap.tsv",
        accuracy=f"{GTEX_BIO}/04_subtissues/fold{{fold}}/accuracy_summary.tsv",
    params:
        tissues=SUBTISSUE_TISSUES_ARG,
        global_seed=config["gtex"]["subtissue_inference"]["global_seed"],
        min_rf_accuracy=config["gtex"]["subtissue_inference"]["min_rf_accuracy"],
        out_dir=f"{GTEX_BIO}/04_subtissues/fold{{fold}}",
    wildcard_constraints:
        fold=SUBTISSUE_CV_FOLD_PATTERN,
    threads: config["gtex"]["subtissue_inference"]["n_jobs_per_fold"]
    resources:
        mem_mb=64000,
        runtime=2160,
    conda: "clamp-analyses"
    shell:
        "python scripts/gtex/subtissue_fold_rf_shap.py --clamp-full-rds {input.clamp_full_rds} "
        "--gtex-meta {input.gtex_meta} --membership {input.membership} --fold {wildcards.fold} "
        "--tissues {params.tissues} --global-seed {params.global_seed} "
        "--min-rf-accuracy {params.min_rf_accuracy} --n-jobs {threads} --out-dir {params.out_dir}"


rule subtissue_lr_eval_gtex:
    input:
        membership=rules.subtissue_cv_folds_gtex.output.membership,
        oof_shap=expand(f"{GTEX_BIO}/04_subtissues/fold{{fold}}/oof_shap.tsv", fold=SUBTISSUE_CV_FOLDS),
        accuracy=expand(f"{GTEX_BIO}/04_subtissues/fold{{fold}}/accuracy_summary.tsv", fold=SUBTISSUE_CV_FOLDS),
        clamp_full_rds=rules.clamp_gtex.output.full_rds,
        gtex_meta=config["gtex"]["metadata"],
        anatomical_spec=config["gtex"]["subtissue_inference"]["anatomical_spec"],
    output:
        anatomical=f"{GTEX_BIO}/04_subtissues/anatomical_subtissue_results.tsv",
        all_smtsd=f"{GTEX_BIO}/04_subtissues/all_smtsd_results.tsv",
        supplementary=f"{GTEX_BIO}/04_subtissues/subtissue_supplementary_table.tsv",
        canonical_predictions=f"{GTEX_BIO}/04_subtissues/subtissue_oof_predictions.tsv",
        canonical_confusion=f"{GTEX_BIO}/04_subtissues/subtissue_confusion_matrices.tsv",
        permutation_nulls=f"{GTEX_BIO}/04_subtissues/subtissue_permutation_nulls.tsv",
        bootstrap=f"{GTEX_BIO}/04_subtissues/subtissue_bootstrap_distributions.tsv",
        counts=f"{GTEX_BIO}/04_subtissues/subtissue_counts.tsv",
        exclusions=f"{GTEX_BIO}/04_subtissues/subtissue_exclusions.tsv",
        fold_audit=f"{GTEX_BIO}/04_subtissues/fold_smtsd_audit.tsv",
        summary=f"{GTEX_BIO}/04_subtissues/subtissue_lr_summary.tsv",
        predictions=f"{GTEX_BIO}/04_subtissues/subtissue_lr_oof_predictions.tsv",
        confusion=f"{GTEX_BIO}/04_subtissues/subtissue_lr_confusion_matrices.tsv",
        dropped=f"{GTEX_BIO}/04_subtissues/dropped_subtissues.tsv",
        rf_accuracy=f"{GTEX_BIO}/04_subtissues/rf_oof_accuracy_summary.tsv",
    params:
        fold_dir=f"{GTEX_BIO}/04_subtissues",
        n_folds=config["gtex"]["subtissue_inference"]["n_folds"],
        tissues=SUBTISSUE_TISSUES_ARG,
        min_subtissue_samples=config["gtex"]["subtissue_inference"]["min_subtissue_samples"],
        min_rf_accuracy=config["gtex"]["subtissue_inference"]["min_rf_accuracy"],
        global_seed=config["gtex"]["subtissue_inference"]["global_seed"],
        permutation_repeats=config["gtex"]["subtissue_inference"]["permutation_repeats"],
        bootstrap_repeats=config["gtex"]["subtissue_inference"]["bootstrap_repeats"],
        out_dir=f"{GTEX_BIO}/04_subtissues",
    threads: config["gtex"]["subtissue_inference"]["evaluation_jobs"]
    resources:
        mem_mb=64000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "python scripts/gtex/subtissue_lr_eval.py --membership {input.membership} "
        "--fold-dir {params.fold_dir} --n-folds {params.n_folds} "
        "--clamp-full-rds {input.clamp_full_rds} --gtex-meta {input.gtex_meta} "
        "--anatomical-spec {input.anatomical_spec} "
        "--tissues {params.tissues} --min-subtissue-samples {params.min_subtissue_samples} "
        "--min-rf-accuracy {params.min_rf_accuracy} --global-seed {params.global_seed} "
        "--permutation-repeats {params.permutation_repeats} "
        "--bootstrap-repeats {params.bootstrap_repeats} --n-jobs {threads} --out-dir {params.out_dir}"


rule subtissues_report_gtex:
    input:
        anatomical=rules.subtissue_lr_eval_gtex.output.anatomical,
        all_smtsd=rules.subtissue_lr_eval_gtex.output.all_smtsd,
        supplementary=rules.subtissue_lr_eval_gtex.output.supplementary,
        confusion=rules.subtissue_lr_eval_gtex.output.canonical_confusion,
        permutation_nulls=rules.subtissue_lr_eval_gtex.output.permutation_nulls,
        bootstrap=rules.subtissue_lr_eval_gtex.output.bootstrap,
        counts=rules.subtissue_lr_eval_gtex.output.counts,
        exclusions=rules.subtissue_lr_eval_gtex.output.exclusions,
        fold_audit=rules.subtissue_lr_eval_gtex.output.fold_audit,
        notebook=f"{GTEX_BIO_NB}/04_subtissues.ipynb",
    output:
        complete=touch(f"{GTEX_BIO}/04_subtissues/notebook.complete"),
    log:
        notebook=f"{GTEX_BIO_NB}/04_subtissues.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{GTEX_BIO_NB}/04_subtissues.ipynb"


rule biology_gtex:
    input:
        rules.kmeans_clustering_report_gtex.output.complete,
        rules.lv_importance_rf_true_labels_report_gtex.output.complete,
        rules.lv_importance_rf_true_labels_biology_gtex.output.complete,
        rules.global_alignment_rf_true_labels_gtex.output.complete,
        rules.liver_disentangle_xcell_rf_true_labels_gtex.output.complete,
        rules.multitissue_xcell_recovery_gtex.output.complete,
        rules.geneset_orthogonality_gtex.output.complete,
        rules.subtissues_report_gtex.output.complete,
