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
        runtime_seed_totals=rules.computational_timing_report_pseudobulk.output.seed_totals,
        holdout_predictions=rules.grouped_cv_analysis_pseudobulk.output.predictions,
        holdout_thresholded_metrics=rules.grouped_cv_analysis_pseudobulk.output.thresholded_metrics,
        holdout_thresholded_summary=rules.grouped_cv_analysis_pseudobulk.output.thresholded_summary,
        marker=rules.disentangle_pseudobulk.output.module_panel_ready,
        hard=rules.hard_cell_types_pseudobulk.output.group_dot_ready,
        related_corr=rules.hard_cell_types_pseudobulk.output.correlations,
        hard_pair_loading_ready=rules.hard_pair_loadings_pseudobulk.output.panel_ready,
        hard_pair_loading_stats=rules.hard_pair_loadings_pseudobulk.output.stats,
        subtissue_confusion=rules.subtissue_lr_eval_gtex.output.canonical_confusion,
        subtissue_anatomical=rules.subtissue_lr_eval_gtex.output.anatomical,
        subtissues_complete=rules.subtissues_report_gtex.output.complete,
        liver_xcell_scatter=rules.liver_disentangle_xcell_rf_true_labels_gtex.output.xcell_scatter,
        liver_lv_pathways=rules.liver_disentangle_xcell_rf_true_labels_gtex.output.lv_pathways,
        multitissue_panel_ready=rules.multitissue_xcell_recovery_gtex.output.panel_ready,
        donor_bulk_predictions=rules.donor_bulk_recovery_report.output.fig2_predictions,
        donor_bulk_statistics=rules.donor_bulk_recovery_report.output.fig2_statistics,
        donor_bulk_perez_annotation=rules.donor_bulk_recovery_report.output.perez_annotation,
        donor_bulk_perez_activity=rules.donor_bulk_recovery_report.output.perez_activity,
        donor_bulk_perez_summary=rules.donor_bulk_recovery_report.output.perez_summary,
        notebook=f"{PANELS_NB}/fig2.ipynb",
    output:
        png=f"{PANELS}/fig2/fig2.png",
        pdf=f"{PANELS}/fig2/fig2.pdf",
        svg=f"{PANELS}/fig2/fig2.svg",
        complete=touch(f"{PANELS}/fig2/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/fig2.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/fig2.ipynb"


rule donor_bulk_figure2:
    input:
        rules.fig2_panel.output.complete,


rule supp1_panel:
    input:
        benchmark_long=rules.benchmark_pseudobulk.output.long,
        bootstrap=rules.benchmark_pseudobulk.output.bootstrap,
        heatmap_long=rules.single_cell_recovery_pseudobulk.output.heatmap,
        corr_full=rules.benchmark_pseudobulk.output.corr,
        assignments=rules.benchmark_pseudobulk.output.assignments,
        related_corr=rules.hard_cell_types_pseudobulk.output.correlations,
        ari_data=rules.kmeans_clustering_report_gtex.output.ari_data,
        ari_comparisons=rules.kmeans_clustering_report_gtex.output.ari_comparisons,
        gene_fraction_ari_data=rules.kmeans_clustering_report_gtex.output.gene_fraction_ari_data,
        gene_fraction_ari_comparisons=rules.kmeans_clustering_report_gtex.output.gene_fraction_ari_comparisons,
        tissue_subtissue_heatmap=rules.lv_importance_rf_true_labels_biology_gtex.output.tissue_subtissue_heatmap,
        z_matrix_subtissue=rules.global_alignment_rf_true_labels_gtex.output.subtissue_summary,
        marker=rules.disentangle_pseudobulk.output.module_panel_ready,
        hard=rules.hard_cell_types_pseudobulk.output.group_dot_ready,
        orthogonality_per_term=rules.geneset_orthogonality_gtex.output.per_term,
        orthogonality_summary=rules.geneset_orthogonality_gtex.output.summary,
        notebook=f"{PANELS_NB}/supp1.ipynb",
    output:
        png=f"{PANELS}/supp1/supp1.png",
        pdf=f"{PANELS}/supp1/supp1.pdf",
        svg=f"{PANELS}/supp1/supp1.svg",
        complete=touch(f"{PANELS}/supp1/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/supp1.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/supp1.ipynb"


rule supp2_panel:
    input:
        benchmark_long=rules.benchmark_pseudobulk.output.long,
        bootstrap=rules.benchmark_pseudobulk.output.bootstrap,
        heatmap_long=rules.single_cell_recovery_pseudobulk.output.heatmap,
        corr_full=rules.benchmark_pseudobulk.output.corr,
        assignments=rules.benchmark_pseudobulk.output.assignments,
        related_corr=rules.hard_cell_types_pseudobulk.output.correlations,
        ari_data=rules.kmeans_clustering_report_gtex.output.ari_data,
        ari_comparisons=rules.kmeans_clustering_report_gtex.output.ari_comparisons,
        gene_fraction_ari_data=rules.kmeans_clustering_report_gtex.output.gene_fraction_ari_data,
        gene_fraction_ari_comparisons=rules.kmeans_clustering_report_gtex.output.gene_fraction_ari_comparisons,
        tissue_subtissue_heatmap=rules.lv_importance_rf_true_labels_biology_gtex.output.tissue_subtissue_heatmap,
        z_matrix_subtissue=rules.global_alignment_rf_true_labels_gtex.output.subtissue_summary,
        marker=rules.disentangle_pseudobulk.output.module_panel_ready,
        hard=rules.hard_cell_types_pseudobulk.output.group_dot_ready,
        orthogonality_per_term=rules.geneset_orthogonality_gtex.output.per_term,
        orthogonality_summary=rules.geneset_orthogonality_gtex.output.summary,
        notebook=f"{PANELS_NB}/supp2.ipynb",
    output:
        png=f"{PANELS}/supp2/supp2.png",
        pdf=f"{PANELS}/supp2/supp2.pdf",
        svg=f"{PANELS}/supp2/supp2.svg",
        complete=touch(f"{PANELS}/supp2/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/supp2.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/supp2.ipynb"


rule panels_pseudobulk:
    input:
        rules.fig2_panel.output.complete,
        rules.supp1_panel.output.complete,
        rules.supp2_panel.output.complete,
