PANELS = config["paths"]["panels"]
PANELS_NB = os.path.join(REPO_ROOT, config["paths"]["panels_notebooks"])


def frozen_gtex(relative_path):
    """Reuse an existing GTEx result without adding its producer to this DAG."""
    return ancient(f"{GTEX_BIO}/{relative_path}")

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
        related_corr=rules.hard_cell_types_pseudobulk.output.correlations,
        hard_pair_loading_ready=rules.hard_pair_loadings_pseudobulk.output.panel_ready,
        hard_pair_loading_stats=rules.hard_pair_loadings_pseudobulk.output.stats,
        subtissue_confusion=frozen_gtex("04_subtissues/subtissue_confusion_matrices.tsv"),
        subtissue_anatomical=frozen_gtex("04_subtissues/anatomical_subtissue_results.tsv"),
        liver_xcell_scatter=frozen_gtex("04_liver_disentangle_xcell_rf_true_labels/liver_xcell_scatter.csv"),
        liver_lv_pathways=frozen_gtex("04_liver_disentangle_xcell_rf_true_labels/liver_lv_pathways.csv"),
        multitissue_panel_ready=frozen_gtex("06_multitissue_xcell_recovery/multitissue_recovery_panel_ready.csv"),
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
        runtime_seed_totals=rules.computational_timing_report_pseudobulk.output.seed_totals,
        # Provenance anchors: Supplementary Fig. 1e is the pseudobulk projection
        # benchmark (including Perez B cell = LV10), never the donor-bulk result.
        heatmap_long=rules.single_cell_recovery_pseudobulk.output.heatmap,
        pseudobulk_recovery=rules.single_cell_recovery_pseudobulk.output.recovery,
        pseudobulk_overall=rules.single_cell_recovery_pseudobulk.output.overall,
        corr_full=rules.benchmark_pseudobulk.output.corr,
        assignments=rules.benchmark_pseudobulk.output.assignments,
        related_corr=rules.hard_cell_types_pseudobulk.output.correlations,
        ari_data=frozen_gtex("00_kmeans_clustering/ari_data.csv"),
        ari_comparisons=frozen_gtex("00_kmeans_clustering/ari_comparisons.csv"),
        gene_fraction_ari_data=frozen_gtex("00_kmeans_clustering/gene_fraction_ari_data.csv"),
        gene_fraction_ari_comparisons=frozen_gtex("00_kmeans_clustering/gene_fraction_ari_comparisons.csv"),
        tissue_subtissue_heatmap=frozen_gtex("02_LV_importance_rf_true_labels_biology/tissue_subtissue_heatmap.csv"),
        z_matrix_subtissue=frozen_gtex("07_global_alignment_rf_true_labels/cumulative20/gtex_global_alignment_summary.csv"),
        marker=rules.disentangle_pseudobulk.output.module_panel_ready,
        hard=rules.hard_cell_types_pseudobulk.output.group_dot_ready,
        orthogonality_per_term=frozen_gtex("07_geneset_orthogonality/orthogonality_per_term.csv"),
        orthogonality_summary=frozen_gtex("07_geneset_orthogonality/orthogonality_summary.csv"),
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
        # Provenance anchors: Supplementary Fig. 2a is the donor-bulk projection
        # benchmark (including Perez B cell = LV32), never the pseudobulk result.
        purity=rules.donor_bulk_recovery_report.output.supp2_purity,
        umap_cells=rules.donor_bulk_recovery_report.output.supp2_umap_cells,
        umap_lvs=rules.donor_bulk_recovery_report.output.supp2_umap_lvs,
        donor_bulk_recovery=rules.donor_bulk_analysis_tables.output.recovery,
        donor_bulk_overall=rules.donor_bulk_analysis_tables.output.recovery_overall,
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


rule supp3_panel:
    input:
        benchmark_long=rules.benchmark_pseudobulk.output.long,
        bootstrap=rules.benchmark_pseudobulk.output.bootstrap,
        heatmap_long=rules.single_cell_recovery_pseudobulk.output.heatmap,
        corr_full=rules.benchmark_pseudobulk.output.corr,
        assignments=rules.benchmark_pseudobulk.output.assignments,
        related_corr=rules.hard_cell_types_pseudobulk.output.correlations,
        ari_data=frozen_gtex("00_kmeans_clustering/ari_data.csv"),
        ari_comparisons=frozen_gtex("00_kmeans_clustering/ari_comparisons.csv"),
        gene_fraction_ari_data=frozen_gtex("00_kmeans_clustering/gene_fraction_ari_data.csv"),
        gene_fraction_ari_comparisons=frozen_gtex("00_kmeans_clustering/gene_fraction_ari_comparisons.csv"),
        tissue_subtissue_heatmap=frozen_gtex("02_LV_importance_rf_true_labels_biology/tissue_subtissue_heatmap.csv"),
        z_matrix_subtissue=frozen_gtex("07_global_alignment_rf_true_labels/cumulative20/gtex_global_alignment_summary.csv"),
        marker=rules.disentangle_pseudobulk.output.module_panel_ready,
        hard=rules.hard_cell_types_pseudobulk.output.group_dot_ready,
        orthogonality_per_term=frozen_gtex("07_geneset_orthogonality/orthogonality_per_term.csv"),
        orthogonality_summary=frozen_gtex("07_geneset_orthogonality/orthogonality_summary.csv"),
        runtime_gtex_per_fit=rules.computational_timing_report_gtex.output.per_fit,
        notebook=f"{PANELS_NB}/supp3.ipynb",
    output:
        png=f"{PANELS}/supp3/supp3.png",
        pdf=f"{PANELS}/supp3/supp3.pdf",
        svg=f"{PANELS}/supp3/supp3.svg",
        complete=touch(f"{PANELS}/supp3/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/supp3.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/supp3.ipynb"


rule panels_pseudobulk:
    input:
        rules.fig2_panel.output.complete,
        rules.supp1_panel.output.complete,
        rules.supp2_panel.output.complete,
        rules.supp3_panel.output.complete,
