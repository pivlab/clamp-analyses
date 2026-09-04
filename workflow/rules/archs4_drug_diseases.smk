A4_DD_CFG = A4_CFG["drug_diseases"]
A4_DD_DATA = A4_DD_CFG["data_dir"]
A4_DD_SPREDIXCAN = A4_DD_CFG["spredixcan_dir"]
A4_DD_COMPUTE = A4_DD_CFG["compute_root"]
A4_DD_REPORT = A4_DD_CFG["report_root"]
A4_DD_NLVS = ",".join(map(str, A4_DD_CFG["n_top_lvs"]))
A4_DD_NGENES = ",".join(map(str, A4_DD_CFG["n_top_genes"]))

def dd_model_dir(compendium):
    return f"{A4_CAN_FINAL_ROOT}/{compendium}"

def dd_prediction_dir(compendium):
    return f"{A4_DD_COMPUTE}/predictions/module_based_{compendium}"

# ============================================================
# Step 1: project raw data into each compendium's canonical CLAMP model
# ============================================================

rule project_spredixcan_drug_diseases:
    input:
        model=f"{A4_CAN_FINAL_ROOT}/{{compendium}}/{A4_CAN_MODEL_NAME}.rds",
        provenance=f"{A4_CAN_FINAL_ROOT}/{{compendium}}/adoption_validation.json",
        spredixcan=A4_DD_SPREDIXCAN,
        script="scripts/archs4/drug_diseases/project_spredixcan.py"
    output:
        directory(f"{A4_DD_COMPUTE}/projections/{{compendium}}/spredixcan")
    wildcard_constraints:
        compendium=A4_CAN_COMPENDIUM_PATTERN
    resources:
        mem_mb=32000, runtime=720
    params:
        model_dir=lambda wc: dd_model_dir(wc.compendium)
    conda: "clamp-analyses"
    shell:
        "python {input.script} --model-dir {params.model_dir} --model-name {A4_CAN_MODEL_NAME} --spredixcan-dir {input.spredixcan} "
        "--compendium {wildcards.compendium} --output-dir {output}"

rule project_lincs_drug_diseases:
    input:
        model=f"{A4_CAN_FINAL_ROOT}/{{compendium}}/{A4_CAN_MODEL_NAME}.rds",
        provenance=f"{A4_CAN_FINAL_ROOT}/{{compendium}}/adoption_validation.json",
        lincs=f"{A4_DD_DATA}/lincs-data.pkl",
        script="scripts/archs4/drug_diseases/project_lincs.py"
    output:
        f"{A4_CAN_FINAL_ROOT}/{{compendium}}/lincs-projection.pkl"
    wildcard_constraints:
        compendium=A4_CAN_COMPENDIUM_PATTERN
    resources:
        mem_mb=16000, runtime=360
    params:
        model_dir=lambda wc: dd_model_dir(wc.compendium)
    conda: "clamp-analyses"
    shell:
        "python {input.script} --model-dir {params.model_dir} --model-name {A4_CAN_MODEL_NAME} --lincs-file {input.lincs} "
        "--compendium {wildcards.compendium} --output {output}"

# ============================================================
# Step 2: predict drug-disease associations (module-based per compendium,
# plus the gene-space baseline)
# ============================================================

rule predict_module_drug_diseases:
    input:
        lincs=f"{A4_CAN_FINAL_ROOT}/{{compendium}}/lincs-projection.pkl",
        spredixcan=f"{A4_DD_COMPUTE}/projections/{{compendium}}/spredixcan",
        gold=f"{A4_DD_DATA}/gold_standard.pkl", ukb=f"{A4_DD_DATA}/phenomexcan_traits_fullcode_to_efo.tsv",
        efo=f"{A4_DD_DATA}/term_id_xrefs.tsv.gz", do=f"{A4_DD_DATA}/xrefs-prop-slim.tsv",
        script="scripts/archs4/drug_diseases/predict_module_based.py"
    output:
        directory(f"{A4_DD_COMPUTE}/predictions/module_based_{{compendium}}")
    wildcard_constraints:
        compendium=A4_CAN_COMPENDIUM_PATTERN
    resources:
        mem_mb=16000, runtime=720
    conda: "clamp-analyses"
    shell:
        "python {input.script} --compendium {wildcards.compendium} --lincs-projection {input.lincs} --spredixcan-proj-dir {input.spredixcan}/proj "
        "--gold-standard {input.gold} --ukb-efo {input.ukb} --efo-xrefs {input.efo} --do-xrefs {input.do} --n-top-lvs {A4_DD_NLVS} --output-dir {output}"

rule predict_gene_drug_diseases:
    input:
        spredixcan_raw=f"{A4_DD_COMPUTE}/projections/archs4/spredixcan",
        lincs=f"{A4_DD_DATA}/lincs-data.pkl",
        gold=f"{A4_DD_DATA}/gold_standard.pkl", ukb=f"{A4_DD_DATA}/phenomexcan_traits_fullcode_to_efo.tsv",
        efo=f"{A4_DD_DATA}/term_id_xrefs.tsv.gz", do=f"{A4_DD_DATA}/xrefs-prop-slim.tsv",
        script="scripts/archs4/drug_diseases/predict_single_gene.py"
    output:
        directory(f"{A4_DD_COMPUTE}/predictions/gene_based")
    resources:
        mem_mb=16000, runtime=360
    conda: "clamp-analyses"
    shell:
        "python {input.script} --spredixcan-raw-dir {input.spredixcan_raw}/raw --lincs-file {input.lincs} "
        "--gold-standard {input.gold} --ukb-efo {input.ukb} --efo-xrefs {input.efo} --do-xrefs {input.do} --n-top-genes {A4_DD_NGENES} --output-dir {output}"

# ============================================================
# Step 3: aggregate predictions and compare all 4 methods
# ============================================================

rule aggregate_drug_diseases:
    input:
        arch=dd_prediction_dir("archs4"), recount=dd_prediction_dir("recount2"), gtex=dd_prediction_dir("gtex"),
        gene=f"{A4_DD_COMPUTE}/predictions/gene_based",
        gold=f"{A4_DD_DATA}/gold_standard.pkl", ukb=f"{A4_DD_DATA}/phenomexcan_traits_fullcode_to_efo.tsv",
        efo=f"{A4_DD_DATA}/term_id_xrefs.tsv.gz", do=f"{A4_DD_DATA}/xrefs-prop-slim.tsv",
        script="scripts/archs4/drug_diseases/aggregate_predictions.py"
    output:
        directory(f"{A4_DD_REPORT}/aggregate/max")
    resources:
        mem_mb=32000, runtime=720
    conda: "clamp-analyses"
    shell:
        "python {input.script} --module-archs4-dir {input.arch} --module-recount2-dir {input.recount} --module-gtex-dir {input.gtex} --gene-dir {input.gene} "
        "--gold-standard {input.gold} --ukb-efo {input.ukb} --efo-xrefs {input.efo} --do-xrefs {input.do} --n-top-lvs {A4_DD_NLVS} --n-top-genes {A4_DD_NGENES} "
        "--output-dir {output}"

rule prepare_three_model_comparison_data:
    input:
        arch=dd_prediction_dir('archs4'), recount=dd_prediction_dir('recount2'), gtex=dd_prediction_dir('gtex'),
        gene=f"{A4_DD_COMPUTE}/predictions/gene_based",
        gold=f"{A4_DD_DATA}/gold_standard.pkl",
        prep="scripts/archs4/drug_diseases/prepare_three_model_comparison.py"
    output:
        stats=f"{A4_DD_REPORT}/three_model_comparison/data/figure_statistics.txt"
    resources:
        mem_mb=16000, runtime=360
    conda: "clamp-analyses"
    shell:
        "python {input.prep} --archs4 {input.arch}/lincs/predictions/dotprod_neg --recount2 {input.recount}/lincs/predictions/dotprod_neg --gtex {input.gtex}/lincs/predictions/dotprod_neg --gene {input.gene}/lincs/predictions/dotprod_neg --gold-standard {input.gold} "
        "--n-top-lvs {A4_DD_CFG[figure_n_top_lvs]} --n-top-genes {A4_DD_CFG[figure_n_top_genes]} --n-boot {A4_DD_CFG[bootstrap][n_boot]} --seed {A4_DD_CFG[bootstrap][seed]} --output-dir {A4_DD_REPORT}/three_model_comparison/data"

# ============================================================
# Step 4: reports
# ============================================================

rule three_model_comparison_report:
    input:
        stats=rules.prepare_three_model_comparison_data.output.stats,
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/02_three_model_comparison.ipynb"
    output:
        complete=touch(f"{A4_DD_REPORT}/three_model_comparison/notebook.complete")
    log:
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/02_three_model_comparison.executed.ipynb"
    resources:
        mem_mb=4000, runtime=120
    conda: "clamp-analyses"
    notebook:
        f"{A4_BIO_NB}/02_drug_diseases_associations/02_three_model_comparison.ipynb"

rule drug_diseases_prediction_performance_report:
    input:
        aggregate=f"{A4_DD_REPORT}/aggregate/max",
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/00_prediction_performance.ipynb"
    output:
        complete=touch(f"{A4_DD_REPORT}/00_prediction_performance/notebook.complete")
    log:
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/00_prediction_performance.executed.ipynb"
    conda: "clamp-analyses"
    notebook:
        f"{A4_BIO_NB}/02_drug_diseases_associations/00_prediction_performance.ipynb"

rule drug_diseases_prediction_plots_report:
    input:
        aggregate=f"{A4_DD_REPORT}/aggregate/max",
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/01_prediction_performance_plots.ipynb"
    output:
        complete=touch(f"{A4_DD_REPORT}/01_prediction_performance_plots/notebook.complete")
    log:
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/01_prediction_performance_plots.executed.ipynb"
    conda: "clamp-analyses"
    notebook:
        f"{A4_BIO_NB}/02_drug_diseases_associations/01_prediction_performance_plots.ipynb"

rule archs4_drug_diseases:
    input:
        expand(f"{A4_CAN_FINAL_ROOT}/{{compendium}}/{A4_CAN_MODEL_NAME}.rds", compendium=A4_CAN_COMPENDIA),
        expand(f"{A4_CAN_FINAL_ROOT}/{{compendium}}/lincs-projection.pkl", compendium=A4_CAN_COMPENDIA),
        rules.aggregate_drug_diseases.output,
        rules.prepare_three_model_comparison_data.output,
        rules.three_model_comparison_report.output,
        rules.drug_diseases_prediction_performance_report.output,
        rules.drug_diseases_prediction_plots_report.output,
