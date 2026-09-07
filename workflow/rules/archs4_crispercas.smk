import os

# CRISPR-Cas9 lipid-screen analysis of the canonical ARCHS4 model.
A4_CRISPERCAS_CFG = A4_CFG["crispercas"]
A4_CRISPERCAS_ROOT = A4_CRISPERCAS_CFG["compute_root"]
A4_CRISPERCAS_NB = os.path.join(REPO_ROOT, A4_CRISPERCAS_CFG["notebooks"])


# ============================================================
# Step 1: fgsea CRISPR-Cas9 screen scoring
# ============================================================
# Scores every LV of the canonical ARCHS4 model against the CRISPR-Cas9
# lipid screen gene sets and selects the LVs significant at FDR < 0.05.

rule fgsea_crispercas:
    input:
        z=f"{A4_CAN_FINAL_ROOT}/archs4/Z.csv",
        lipid_deg=A4_CRISPERCAS_CFG["lipid_deg"],
        script="scripts/archs4/crispercas/gene_enrichment_crisprcas.R",
    output:
        fgsea=f"{A4_CRISPERCAS_ROOT}/00_gene_enrichment_CRISPRCas9/fgsea_all_lvs.csv",
    log:
        f"{A4_CRISPERCAS_ROOT}/00_gene_enrichment_CRISPRCas9/fgsea.log"
    resources:
        mem_mb=16000,
        runtime=60,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --lipid-deg {input.lipid_deg} "
        "--out {output.fgsea} --n-reps {A4_CRISPERCAS_CFG[n_reps]} "
        "> {log} 2>&1"


rule gene_enrichment_crispercas_report:
    input:
        fgsea=rules.fgsea_crispercas.output.fgsea,
        notebook=f"{A4_CRISPERCAS_NB}/00_gene_enrichment_CRISPRCas9.ipynb",
    output:
        lv_sig=f"{A4_CRISPERCAS_ROOT}/00_gene_enrichment_CRISPRCas9/lv_sig_sub_ordered.csv",
        complete=touch(f"{A4_CRISPERCAS_ROOT}/00_gene_enrichment_CRISPRCas9/notebook.complete"),
    log:
        notebook=f"{A4_CRISPERCAS_NB}/00_gene_enrichment_CRISPRCas9.executed.ipynb",
    params:
        fdr=A4_CRISPERCAS_CFG["fdr"],
    conda: "clamp-analyses"
    notebook:
        f"{A4_CRISPERCAS_NB}/00_gene_enrichment_CRISPRCas9.ipynb"


# ============================================================
# Step 2: ORA against BP, canonical, CellMarker, ARCHS4 cell-lines, Azimuth
# ============================================================

rule fetch_azimuth_gmt:
    output:
        A4_CRISPERCAS_CFG["azimuth_gmt"],
    log:
        A4_CRISPERCAS_CFG["azimuth_gmt"] + ".log"
    shell:
        "curl --fail --location --retry 3 '{A4_CRISPERCAS_CFG[azimuth_url]}' "
        "--output {output} 2> {log}"


rule ora_azimuth_crispercas:
    input:
        z=f"{A4_CAN_FINAL_ROOT}/archs4/Z.csv",
        manifest=f"{A4_CAN_FINAL_ROOT}/archs4/manifest.json",
        database=rules.fetch_azimuth_gmt.output,
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(f"{A4_CAN_FINAL_ROOT}/archs4/ora/azimuth"),
    log:
        f"{A4_CAN_FINAL_ROOT}/archs4/ora/azimuth.log"
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --out-dir {output.ora_dir} "
        "--dataset archs4 --fraction 100 --seed 1 "
        "--model {A4_CAN_MODEL_NAME} --model-manifest {input.manifest} "
        "--database azimuth --database-label 'Azimuth' "
        "--database-type gmt --database-path {input.database} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {A4_CRISPERCAS_CFG[azimuth_min_size]} "
        "--max-size {A4_CRISPERCAS_CFG[azimuth_max_size]} --pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} "
        "--qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} --fdr {A4_CFG[ora][fdr][0]} "
        "> {log} 2>&1"


rule fetch_archs4_cell_lines_gmt:
    output:
        A4_CRISPERCAS_CFG["archs4_cell_lines_gmt"],
    log:
        A4_CRISPERCAS_CFG["archs4_cell_lines_gmt"] + ".log"
    shell:
        "curl --fail --location --retry 3 '{A4_CRISPERCAS_CFG[archs4_cell_lines_url]}' "
        "--output {output} 2> {log}"


rule ora_archs4_cell_lines_crispercas:
    input:
        z=f"{A4_CAN_FINAL_ROOT}/archs4/Z.csv",
        manifest=f"{A4_CAN_FINAL_ROOT}/archs4/manifest.json",
        database=rules.fetch_archs4_cell_lines_gmt.output,
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(f"{A4_CAN_FINAL_ROOT}/archs4/ora/archs4_cell_lines"),
    log:
        f"{A4_CAN_FINAL_ROOT}/archs4/ora/archs4_cell_lines.log"
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --out-dir {output.ora_dir} "
        "--dataset archs4 --fraction 100 --seed 1 "
        "--model {A4_CAN_MODEL_NAME} --model-manifest {input.manifest} "
        "--database archs4_cell_lines --database-label 'ARCHS4 Cell-lines' "
        "--database-type gmt --database-path {input.database} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {A4_CFG[ora][min_size]} "
        "--max-size {A4_CRISPERCAS_CFG[archs4_cell_lines_max_size]} --pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} "
        "--qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} --fdr {A4_CFG[ora][fdr][0]} "
        "> {log} 2>&1"


# ============================================================
# Step 3: top-1% LINCS L1000 perturbagens per significant LV
# ============================================================
# Ranking only (no p-value), same top_pct convention as ORA queries.

rule lincs_top_crispercas:
    input:
        lincs=f"{A4_CAN_FINAL_ROOT}/archs4/lincs-projection.pkl",
        lv_sig=rules.gene_enrichment_crispercas_report.output.lv_sig,
        script="scripts/archs4/crispercas/extract_lincs_top.py",
    output:
        f"{A4_CRISPERCAS_ROOT}/01_biology_LVs/lincs_top.csv",
    log:
        f"{A4_CRISPERCAS_ROOT}/01_biology_LVs/lincs_top.log"
    conda: "clamp-analyses"
    shell:
        "python {input.script} --lincs {input.lincs} --lv-sig {input.lv_sig} "
        "--top-pct {A4_CRISPERCAS_CFG[lincs_top_pct]} --out {output} > {log} 2>&1"


# ============================================================
# Step 4: report -- one table each for ORA, traits, LINCS
# ============================================================
# One row per CRISPR-significant LV per table.

rule biology_lvs_crispercas_report:
    input:
        lv_sig=rules.gene_enrichment_crispercas_report.output.lv_sig,
        ora_bp=expand(rules.ora_canonical.output.ora_dir, dataset="archs4", database="bp"),
        ora_canonical=expand(rules.ora_canonical.output.ora_dir, dataset="archs4", database="canonical"),
        ora_cellmarker=expand(rules.ora_canonical.output.ora_dir, dataset="archs4", database="cellmarker"),
        ora_azimuth=rules.ora_azimuth_crispercas.output.ora_dir,
        ora_archs4_cell_lines=rules.ora_archs4_cell_lines_crispercas.output.ora_dir,
        bp_gmt=A4_CFG["ora"]["databases"]["bp"]["path"],
        canonical_gls=A4_CRISPERCAS_CFG["canonical_gls"],
        lincs_top=rules.lincs_top_crispercas.output,
        drug_names=A4_CRISPERCAS_CFG["drug_names"],
        notebook=f"{A4_CRISPERCAS_NB}/01_biology_LVs.ipynb",
    output:
        ora_table=f"{A4_CRISPERCAS_ROOT}/01_biology_LVs/ora_table.csv",
        traits_table=f"{A4_CRISPERCAS_ROOT}/01_biology_LVs/traits_table.csv",
        lincs_table=f"{A4_CRISPERCAS_ROOT}/01_biology_LVs/lincs_table.csv",
        complete=touch(f"{A4_CRISPERCAS_ROOT}/01_biology_LVs/notebook.complete"),
    log:
        notebook=f"{A4_CRISPERCAS_NB}/01_biology_LVs.executed.ipynb",
    params:
        fdr=A4_CRISPERCAS_CFG["fdr"],
    conda: "clamp-analyses"
    notebook:
        f"{A4_CRISPERCAS_NB}/01_biology_LVs.ipynb"


# ============================================================
# Step 5: trait-centred stratagenic network and TWAS figure
# ============================================================

rule crispercas_symbol_ensembl_mapping:
    input:
        script="scripts/archs4/crispercas/extract_symbol_ensembl.R",
    output:
        mapping=f"{A4_CRISPERCAS_CFG['network_root']}/symbol_ensembl.csv",
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} {output.mapping}"


rule stratagenic_network_crispercas_report:
    input:
        z=f"{A4_CAN_FINAL_ROOT}/archs4/Z.csv",
        canonical_ora=f"{A4_CAN_FINAL_ROOT}/archs4/ora/canonical/enrichment.csv.gz",
        gls=A4_CRISPERCAS_CFG["canonical_gls"],
        smultixcan=A4_CRISPERCAS_CFG["smultixcan_zscores"],
        mapping=rules.crispercas_symbol_ensembl_mapping.output.mapping,
        notebook=f"{A4_CRISPERCAS_NB}/02_network.ipynb",
    output:
        community_edges=f"{A4_CRISPERCAS_CFG['network_root']}/community_edges.csv",
        community_genes=f"{A4_CRISPERCAS_CFG['network_root']}/community_genes.csv",
        trait_edges=f"{A4_CRISPERCAS_CFG['network_root']}/trait_supported_gene_edges.csv",
        shared_gene_edges=f"{A4_CRISPERCAS_CFG['network_root']}/inter_lv_shared_gene_edges.csv",
        lv_summary=f"{A4_CRISPERCAS_CFG['network_root']}/network_lv_summary.csv",
        lv_traits=f"{A4_CRISPERCAS_CFG['network_root']}/network_lv_significant_traits.csv",
        lv_pathways=f"{A4_CRISPERCAS_CFG['network_root']}/network_lv_significant_pathways.csv",
        twas_matrix=f"{A4_CRISPERCAS_CFG['network_root']}/poster_style_twas_matrix.csv",
        complete=touch(f"{A4_CRISPERCAS_ROOT}/02_network/notebook.complete"),
    log:
        notebook=f"{A4_CRISPERCAS_NB}/02_network.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{A4_CRISPERCAS_NB}/02_network.ipynb"


rule biology_archs4_crispercas:
    input:
        rules.gene_enrichment_crispercas_report.output.complete,
        rules.biology_lvs_crispercas_report.output.complete,
        rules.stratagenic_network_crispercas_report.output.complete,
