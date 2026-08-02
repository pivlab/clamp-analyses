PANELS = config["paths"]["panels"]
PANELS_NB = os.path.join(REPO_ROOT, config["paths"]["panels_notebooks"])

# ============================================================
# Publication figure panels: load CSV outputs from the
# biology_pseudobulk analysis notebooks directly and build native
# ggplot2/patchwork/cowplot panels here (no pre-rendered PNGs).
# ============================================================

rule fig2_panel:
    input:
        benchmark_long=rules.benchmark_pseudobulk.output.long,
        holdout_predictions=rules.grouped_cv_analysis_pseudobulk.output.predictions,
        holdout_thresholded_metrics=rules.grouped_cv_analysis_pseudobulk.output.thresholded_metrics,
        holdout_thresholded_summary=rules.grouped_cv_analysis_pseudobulk.output.thresholded_summary,
        marker=rules.disentangle_pseudobulk.output.module_panel_ready,
        hard=rules.hard_cell_types_pseudobulk.output.group_dot_ready,
        subtissue_lr_summary=rules.subtissue_lr_eval_gtex.output.summary,
        subtissues_complete=rules.subtissues_report_gtex.output.complete,
        liver_xcell_scatter=rules.liver_disentangle_xcell_rf_true_labels_gtex.output.xcell_scatter,
        notebook=f"{PANELS_NB}/fig2.ipynb",
    output:
        png=f"{PANELS}/fig2/fig2.png",
        pdf=f"{PANELS}/fig2/fig2.pdf",
        svg=f"{PANELS}/fig2/fig2.svg",
        complete=touch(f"{PANELS}/fig2/notebook.complete"),
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/fig2.ipynb"


rule supp1_panel:
    input:
        benchmark_long=rules.benchmark_pseudobulk.output.long,
        bootstrap=rules.benchmark_pseudobulk.output.bootstrap,
        heatmap_long=rules.single_cell_recovery_pseudobulk.output.heatmap,
        corr_full=rules.benchmark_pseudobulk.output.corr,
        assignments=rules.benchmark_pseudobulk.output.assignments,
        quality_by_dataset=rules.disentangle_pseudobulk.output.quality_by_dataset,
        related_corr=rules.hard_cell_types_pseudobulk.output.correlations,
        ari_data=rules.kmeans_clustering_report_gtex.output.ari_data,
        ari_comparisons=rules.kmeans_clustering_report_gtex.output.ari_comparisons,
        tissue_concordance=rules.lv_importance_rf_true_labels_biology_gtex.output.tissue_concordance,
        subtissue_confusion=rules.subtissue_lr_eval_gtex.output.confusion,
        notebook=f"{PANELS_NB}/supp1.ipynb",
    output:
        png=f"{PANELS}/supp1/supp1.png",
        pdf=f"{PANELS}/supp1/supp1.pdf",
        svg=f"{PANELS}/supp1/supp1.svg",
        complete=touch(f"{PANELS}/supp1/notebook.complete"),
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/supp1.ipynb"


rule panels_pseudobulk:
    input:
        rules.fig2_panel.output.complete,
        rules.supp1_panel.output.complete,
