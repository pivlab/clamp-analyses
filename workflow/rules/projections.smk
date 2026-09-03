import json
import os

PROJ_CFG = config["projections"]
PROJ_PROD = PROJ_CFG["paths"]["production"]
PROJ_BIO = PROJ_CFG["paths"]["biology"]
PROJ_NB = os.path.join(REPO_ROOT, PROJ_CFG["paths"]["notebooks"])
PROJ_RAW = PROJ_CFG["paths"]["raw_data"]

PROJ_A4_CFG = PROJ_CFG["archs4_model"]
PROJ_A4_DIR = PROJ_A4_CFG["root"]
PROJ_A4_RDS = f"{PROJ_A4_DIR}/{PROJ_A4_CFG['rds']}"
PROJ_A4_Z = f"{PROJ_A4_DIR}/{PROJ_A4_CFG['z']}"
PROJ_A4_SUMMARY = f"{PROJ_A4_DIR}/{PROJ_A4_CFG['summary']}"

PROJ_DATASETS = list(PROJ_CFG["datasets"])
PROJ_DB = list(PROJ_CFG["ora"]["databases"])
PROJ_MODELS = ["ARCHS4", "local"]
PROJ_THRESH = PROJ_CFG["thresholds"]
PROJ_ORA = PROJ_CFG["ora"]

PROJ_DATASET_PATTERN = "|".join(PROJ_DATASETS)
PROJ_DB_PATTERN = "|".join(PROJ_DB)
PROJ_MODEL_PATTERN = "|".join(PROJ_MODELS)


wildcard_constraints:
    dataset=PROJ_DATASET_PATTERN,
    database=PROJ_DB_PATTERN,
    model=PROJ_MODEL_PATTERN,
    fname=r"[^/.][^/]*",


def proj_ds(name):
    return PROJ_CFG["datasets"][name]


def proj_dir(dataset):
    return f"{PROJ_PROD}/{dataset}"


def proj_model_dir(dataset, model):
    return f"{proj_dir(dataset)}/models/{model}"


def proj_lv_dir(dataset, model):
    return f"{proj_dir(dataset)}/lv_stats/{model}"


def proj_ora_dir(dataset, model, database):
    return f"{proj_dir(dataset)}/ora/{model}/{database}"


def proj_raw_files(dataset):
    return [f"{PROJ_RAW}/{dataset}/raw/{d['name']}" for d in proj_ds(dataset)["downloads"]]


def proj_datasets_of_kind(kind):
    return [d for d in PROJ_DATASETS if proj_ds(d)["kind"] == kind]


def proj_model_input(dataset, model):
    return (f"{proj_model_dir(dataset, 'ARCHS4')}/B.csv" if model == "ARCHS4"
            else proj_model_dir(dataset, "local"))


def proj_z(dataset, model):
    return PROJ_A4_Z if model == "ARCHS4" else f"{proj_model_dir(dataset, 'local')}/CLAMPfull/Z.csv"


def proj_summary(dataset, model):
    return PROJ_A4_SUMMARY if model == "ARCHS4" else f"{proj_model_dir(dataset, 'local')}/CLAMPfull/summary.csv"


def proj_b(dataset, model):
    return (f"{proj_model_dir(dataset, 'ARCHS4')}/B.csv" if model == "ARCHS4"
            else f"{proj_model_dir(dataset, 'local')}/CLAMPfull/B.csv")


def proj_db_min_size(database):
    return int(PROJ_ORA["databases"][database].get("min_size", PROJ_ORA["min_size"]))


PROJ_ALL_MECH = [f"{proj_dir(d)}/mechanism_recovery" for d in PROJ_DATASETS]


# ============================================================
# Step 1: download and prepare projection datasets
# ============================================================

localrules: projection_raw_data, projection_unpack


def _download_entry(dataset, fname):
    for d in proj_ds(dataset)["downloads"]:
        if d["name"] == fname:
            return d
    raise KeyError(f"{fname} is not declared for {dataset}")


rule projection_raw_data:
    output:
        f"{PROJ_RAW}/{{dataset}}/raw/{{fname}}",
    params:
        url=lambda w: _download_entry(w.dataset, w.fname)["url"],
        sha256=lambda w: _download_entry(w.dataset, w.fname)["sha256"],
    conda: "clamp-analyses"
    shell:
        "mkdir -p $(dirname {output}) && "
        "curl --fail --location --retry 3 '{params.url}' --output {output}.download && "
        "echo '{params.sha256}  {output}.download' | sha256sum --check --status && "
        "mv {output}.download {output}"


rule projection_unpack:
    """GSE133218 ships its 30 per-sample GTFs inside one tar."""
    input:
        tar=lambda w: f"{PROJ_RAW}/{w.dataset}/raw/{proj_ds(w.dataset)['unpack']['tar']}",
    output:
        marker=touch(f"{PROJ_RAW}/{{dataset}}/raw/.unpacked"),
    conda: "clamp-analyses"
    shell:
        "tar -xf {input.tar} -C $(dirname {input.tar})"


rule prepare_geo_matrix:
    input:
        raw=lambda w: proj_raw_files(w.dataset),
        script="scripts/archs4/projections/prepare_geo_matrix.R",
        common="scripts/archs4/common.R",
    output:
        prepared=directory(f"{PROJ_PROD}/{{dataset}}/prepared"),
    params:
        counts=lambda w: proj_raw_files(w.dataset)[0],
        annotation=lambda w: json.dumps(proj_ds(w.dataset)["sample_annotation"]),
        id_type=lambda w: proj_ds(w.dataset)["id_type"],
    log:
        f"{PROJ_PROD}/{{dataset}}/prepared.log",
    conda: "clamp-analyses"
    resources: mem_mb=16000, runtime=60
    wildcard_constraints:
        dataset="|".join(proj_datasets_of_kind("geo_matrix")),
    shell:
        "Rscript {input.script} --counts-file {params.counts} "
        "--annotation '{params.annotation}' --id-type {params.id_type} "
        "--out-dir {output.prepared} > {log} 2>&1"


rule prepare_geo_gtf_tar:
    input:
        raw=lambda w: proj_raw_files(w.dataset),
        unpacked=f"{PROJ_RAW}/{{dataset}}/raw/.unpacked",
        script="scripts/archs4/projections/prepare_geo_gtf_tar.R",
        common="scripts/archs4/common.R",
    output:
        prepared=directory(f"{PROJ_PROD}/{{dataset}}/prepared"),
    params:
        raw_dir=lambda w: f"{PROJ_RAW}/{w.dataset}/raw",
        annotation=lambda w: json.dumps(proj_ds(w.dataset)["sample_annotation"]),
        id_type=lambda w: proj_ds(w.dataset)["id_type"],
    log:
        f"{PROJ_PROD}/{{dataset}}/prepared.log",
    conda: "clamp-analyses"
    resources: mem_mb=32000, runtime=240
    wildcard_constraints:
        dataset="|".join(proj_datasets_of_kind("geo_gtf_tar")),
    shell:
        "Rscript {input.script} --raw-dir {params.raw_dir} "
        "--annotation '{params.annotation}' --id-type {params.id_type} "
        "--out-dir {output.prepared} > {log} 2>&1"


rule prepare_sc_pseudobulk:
    input:
        raw=lambda w: proj_raw_files(w.dataset),
        script="scripts/archs4/projections/prepare_sc_pseudobulk.R",
        common="scripts/archs4/common.R",
    output:
        prepared=directory(f"{PROJ_PROD}/{{dataset}}/prepared"),
    params:
        raw_dir=lambda w: f"{PROJ_RAW}/{w.dataset}/raw",
        kind=lambda w: proj_ds(w.dataset)["kind"],
        spec=lambda w: json.dumps(proj_ds(w.dataset)),
    log:
        f"{PROJ_PROD}/{{dataset}}/prepared.log",
    conda: "clamp-analyses"
    resources: mem_mb=96000, runtime=480
    wildcard_constraints:
        dataset="|".join(proj_datasets_of_kind("ebi_sc_atlas")),
    shell:
        "Rscript {input.script} --raw-dir {params.raw_dir} --kind {params.kind} "
        "--spec '{params.spec}' --out-dir {output.prepared} > {log} 2>&1"


# ============================================================
# Step 2: fit local models and project into ARCHS4
# ============================================================

rule fit_local_clamp:
    input:
        prepared=f"{PROJ_PROD}/{{dataset}}/prepared",
        priors=list(PROJ_CFG["local_model"]["prior"]["gmts"].values()) +
               [PROJ_CFG["local_model"]["prior"]["cellmarker"]["path"]],
        script="scripts/archs4/projections/fit_local_clamp.R",
        common="scripts/archs4/common.R",
    output:
        model=directory(f"{PROJ_PROD}/{{dataset}}/models/local"),
    params:
        norm=f"{PROJ_PROD}/{{dataset}}/prepared/norm.csv",
        prior_names=",".join(PROJ_CFG["local_model"]["prior"]["gmts"]),
        prior_gmts=",".join(PROJ_CFG["local_model"]["prior"]["gmts"].values()),
        cellmarker=PROJ_CFG["local_model"]["prior"]["cellmarker"]["path"],
        cellmarker_sheet=PROJ_CFG["local_model"]["prior"]["cellmarker"]["sheet"],
        cellmarker_term=PROJ_CFG["local_model"]["prior"]["cellmarker"]["term_column"],
        cellmarker_gene=PROJ_CFG["local_model"]["prior"]["cellmarker"]["gene_column"],
        k_rule=PROJ_CFG["local_model"]["k_rule"],
        max_iter=PROJ_CFG["local_model"]["max_iter"],
        seed=PROJ_CFG["local_model"]["seed"],
    log:
        f"{PROJ_PROD}/{{dataset}}/models/local.log",
    conda: "clamp-analyses"
    resources: mem_mb=32000, runtime=240
    shell:
        "Rscript {input.script} --norm {params.norm} "
        "--prior-names {params.prior_names} --prior-gmts {params.prior_gmts} "
        "--cellmarker-file {params.cellmarker} --cellmarker-sheet {params.cellmarker_sheet} "
        "--cellmarker-term-column {params.cellmarker_term} --cellmarker-gene-column {params.cellmarker_gene} "
        "--k-rule {params.k_rule} --max-iter {params.max_iter} --seed {params.seed} "
        "--out-dir {output.model} > {log} 2>&1"


rule project_into_archs4:
    input:
        prepared=f"{PROJ_PROD}/{{dataset}}/prepared",
        # Deliberately not `ancient`: a published canonical-model update must
        # invalidate every projection and its downstream interpretation.
        model=PROJ_A4_RDS,
        provenance=PROJ_A4_CFG["provenance"],
        script="scripts/archs4/projections/project_clamp.R",
        common="scripts/archs4/common.R",
    output:
        b=f"{PROJ_PROD}/{{dataset}}/models/ARCHS4/B.csv",
        run=directory(f"{PROJ_PROD}/{{dataset}}/models/ARCHS4_run"),
    params:
        targets=f"{PROJ_PROD}/{{dataset}}/models/ARCHS4_targets.tsv",
        norm=f"{PROJ_PROD}/{{dataset}}/prepared/norm.csv",
    log:
        f"{PROJ_PROD}/{{dataset}}/models/ARCHS4.log",
    conda: "clamp-analyses"
    resources: mem_mb=24000, runtime=120
    shell:
        "mkdir -p $(dirname {params.targets}) && "
        "printf 'id\\tinput_path\\tinput_kind\\tout_path\\n{wildcards.dataset}\\t{params.norm}\\tcsv\\t{output.b}\\n' "
        "> {params.targets} && "
        "Rscript {input.script} --model {input.model} --targets {params.targets} "
        "--out-dir {output.run} > {log} 2>&1"


# ============================================================
# Step 3: test LVs, run ORA, and score expected mechanisms
# ============================================================

rule test_lvs:
    input:
        prepared=f"{PROJ_PROD}/{{dataset}}/prepared",
        model=lambda w: proj_model_input(w.dataset, w.model),
        script="scripts/archs4/projections/test_lvs.R",
        common="scripts/archs4/common.R",
    output:
        lv=directory(f"{PROJ_PROD}/{{dataset}}/lv_stats/{{model}}"),
    params:
        b=lambda w: proj_b(w.dataset, w.model),
        summary=lambda w: proj_summary(w.dataset, w.model),
        samples=f"{PROJ_PROD}/{{dataset}}/prepared/samples.csv",
        contrasts=lambda w: json.dumps(proj_ds(w.dataset)["contrasts"]),
        lv_fdr=PROJ_THRESH["lv_fdr"],
        lv_logfc=PROJ_THRESH["lv_logfc"],
        pathway_fdr=PROJ_THRESH["pathway_fdr"],
        pathway_auc=PROJ_THRESH["pathway_auc"],
    log:
        f"{PROJ_PROD}/{{dataset}}/lv_stats/{{model}}.log",
    conda: "clamp-analyses"
    resources: mem_mb=8000, runtime=30
    shell:
        "Rscript {input.script} --b {params.b} --samples {params.samples} "
        "--summary {params.summary} --contrasts '{params.contrasts}' "
        "--lv-fdr {params.lv_fdr} --lv-logfc {params.lv_logfc} "
        "--pathway-fdr {params.pathway_fdr} --pathway-auc {params.pathway_auc} "
        "--out-dir {output.lv} > {log} 2>&1"


rule ora_projection:
    input:
        lv=f"{PROJ_PROD}/{{dataset}}/lv_stats/{{model}}",
        model=lambda w: proj_model_input(w.dataset, w.model),
        database=lambda w: PROJ_ORA["databases"][w.database]["path"],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora=directory(f"{PROJ_PROD}/{{dataset}}/ora/{{model}}/{{database}}"),
    params:
        z=lambda w: proj_z(w.dataset, w.model),
        lv_subset=f"{PROJ_PROD}/{{dataset}}/lv_stats/{{model}}/significant_lvs.txt",
        db_type=lambda w: PROJ_ORA["databases"][w.database]["type"],
        min_size=lambda w: proj_db_min_size(w.database),
        max_size=PROJ_ORA["max_size"],
        top_pct=PROJ_ORA["top_pct"],
        fdr=PROJ_ORA["fdr"],
        n_samples=lambda w: 0,
        extra=lambda w: (
            " --sheet {sheet} --term-column {term} --gene-column {gene}".format(
                sheet=PROJ_ORA["databases"][w.database].get("sheet", "human"),
                term=PROJ_ORA["databases"][w.database].get("term_column", "cell_name"),
                gene=PROJ_ORA["databases"][w.database].get("gene_column", "Symbol"))
            if PROJ_ORA["databases"][w.database]["type"] == "cellmarker_xlsx" else ""),
    log:
        f"{PROJ_PROD}/{{dataset}}/ora/{{model}}/{{database}}.log",
    conda: "clamp-analyses"
    resources:
        mem_mb=PROJ_ORA["resources"]["mem_mb"],
        runtime=PROJ_ORA["resources"]["runtime"],
    shell:
        "Rscript {input.script} --z {params.z} --out-dir {output.ora} "
        "--dataset {wildcards.dataset} --model {wildcards.model} "
        "--n-samples {params.n_samples} --database {wildcards.database} "
        "--database-type {params.db_type} --database-path {input.database} "
        "--lv-subset {params.lv_subset} --top-pct {params.top_pct} "
        "--min-size {params.min_size} --max-size {params.max_size} --fdr {params.fdr} "
        "--pvalue-cutoff 1 --qvalue-cutoff 1{params.extra} > {log} 2>&1"


rule score_mechanisms:
    input:
        ora=lambda w: [proj_ora_dir(w.dataset, m, db) for m in PROJ_MODELS for db in PROJ_DB],
        lv=lambda w: [proj_lv_dir(w.dataset, m) for m in PROJ_MODELS],
        gene_set_gmt=PROJ_ORA["databases"]["canonical"]["path"],
        hallmark_gmt=PROJ_ORA["databases"]["hallmark"]["path"],
        cellmarker=PROJ_ORA["databases"]["cellmarker"]["path"],
        script="scripts/archs4/projections/score_mechanisms.R",
        common="scripts/archs4/common.R",
    output:
        mech=directory(f"{PROJ_PROD}/{{dataset}}/mechanism_recovery"),
    params:
        models_tsv=f"{PROJ_PROD}/{{dataset}}/mechanism_models.tsv",
        ora_tsv=f"{PROJ_PROD}/{{dataset}}/mechanism_ora.tsv",
        models=lambda w: "\\n".join(
            f"{m}\\t{proj_z(w.dataset, m)}\\t{proj_lv_dir(w.dataset, m)}/lv_stats.csv"
            f"\\t{proj_lv_dir(w.dataset, m)}/significant_lvs.txt" for m in PROJ_MODELS),
        ora_rows=lambda w: "\\n".join(
            f"{m}\\t{db}\\t{proj_ora_dir(w.dataset, m, db)}"
            for m in PROJ_MODELS for db in PROJ_DB),
        mechanisms=lambda w: json.dumps(PROJ_CFG["mechanisms"][w.dataset]),
        top_pct=PROJ_ORA["top_pct"],
        min_count=PROJ_ORA["min_count"],
        fdr=PROJ_ORA["fdr"],
        lv_fdr=PROJ_THRESH["lv_fdr"],
        lv_logfc=PROJ_THRESH["lv_logfc"],
        cellmarker_sheet=PROJ_ORA["databases"]["cellmarker"]["sheet"],
        cellmarker_term=PROJ_ORA["databases"]["cellmarker"]["term_column"],
        cellmarker_gene=PROJ_ORA["databases"]["cellmarker"]["gene_column"],
    log:
        f"{PROJ_PROD}/{{dataset}}/mechanism_recovery.log",
    conda: "clamp-analyses"
    resources: mem_mb=16000, runtime=60
    shell:
        "printf 'model\\tz\\tlv_stats\\tsignificant_lvs\\n{params.models}\\n' > {params.models_tsv} && "
        "printf 'model\\tdatabase\\tdir\\n{params.ora_rows}\\n' > {params.ora_tsv} && "
        "Rscript {input.script} --models {params.models_tsv} --ora {params.ora_tsv} "
        "--mechanisms '{params.mechanisms}' --top-pct {params.top_pct} "
        "--min-count {params.min_count} --fdr {params.fdr} "
        "--lv-fdr {params.lv_fdr} --lv-logfc {params.lv_logfc} "
        "--gene-set-gmt {input.gene_set_gmt} --hallmark-gmt {input.hallmark_gmt} "
        "--cellmarker-file {input.cellmarker} "
        "--cellmarker-sheet {params.cellmarker_sheet} --cellmarker-term-column {params.cellmarker_term} "
        "--cellmarker-gene-column {params.cellmarker_gene} "
        "--out-dir {output.mech} > {log} 2>&1"


# ============================================================
# Step 4: report the ARCHS4 projection benchmarks
# ============================================================

PROJ_BENCH = PROJ_CFG["benchmarks"]
PROJ_BENCH_PB_RESULTS = list(PROJ_BENCH["pseudobulk"].values())
PROJ_BENCH_GTEX_RESULTS = list(PROJ_BENCH["gtex"].values())
PROJ_BENCH_PB_DIR = os.path.dirname(PROJ_BENCH["pseudobulk"]["projection"])
PROJ_BENCH_GTEX_DIR = os.path.dirname(PROJ_BENCH["gtex"]["projection"])


rule projection_report_pseudobulk_recovery:
    input:
        results=PROJ_BENCH_PB_RESULTS,
        full_null=ancient(f"{PROJ_BIO}/pseudobulk_recovery/pseudobulk_recovery_long.csv"),
        notebook=f"{PROJ_NB}/03_pseudobulk_recovery.ipynb",
    output:
        complete=touch(f"{PROJ_BENCH_PB_DIR}/notebook.complete"),
    log:
        notebook=f"{PROJ_NB}/03_pseudobulk_recovery.executed.ipynb",
    conda: "clamp-analyses"
    resources: mem_mb=16000, runtime=30
    notebook:
        f"{PROJ_NB}/03_pseudobulk_recovery.ipynb"


rule projection_report_gtex_ari:
    input:
        results=PROJ_BENCH_GTEX_RESULTS,
        local_ari=ancient(f"{GTEX_BIO}/00_kmeans_clustering/ari_data.csv"),
        notebook=f"{PROJ_NB}/04_gtex_ari.ipynb",
    output:
        complete=touch(f"{PROJ_BENCH_GTEX_DIR}/notebook.complete"),
    log:
        notebook=f"{PROJ_NB}/04_gtex_ari.executed.ipynb",
    conda: "clamp-analyses"
    resources: mem_mb=16000, runtime=30
    notebook:
        f"{PROJ_NB}/04_gtex_ari.ipynb"


# ============================================================
# Step 5: aggregate mechanism recovery and render reports
# ============================================================

rule aggregate_projections:
    input:
        mech=PROJ_ALL_MECH,
        lv=[proj_lv_dir(d, m) for d in PROJ_DATASETS for m in PROJ_MODELS],
        script="scripts/archs4/projections/aggregate.R",
        common="scripts/archs4/common.R",
    output:
        agg=directory(f"{PROJ_BIO}/aggregate"),
    params:
        prod_root=PROJ_PROD,
        datasets=",".join(PROJ_DATASETS),
        groups=",".join(proj_ds(d)["group"] for d in PROJ_DATASETS),
    log:
        f"{PROJ_BIO}/aggregate.log",
    conda: "clamp-analyses"
    resources: mem_mb=16000, runtime=60
    shell:
        "Rscript {input.script} --prod-root {params.prod_root} "
        "--datasets {params.datasets} --groups {params.groups} "
        "--out-dir {output.agg} > {log} 2>&1"


rule projection_report_cytokines:
    input:
        agg=rules.aggregate_projections.output.agg,
        notebook=f"{PROJ_NB}/00_cytokines.ipynb",
    output:
        complete=touch(f"{PROJ_BIO}/01_cytokines/notebook.complete"),
    log:
        notebook=f"{PROJ_NB}/00_cytokines.executed.ipynb",
    params:
        group_label="cytokines",
        agg_dir=f"{PROJ_BIO}/aggregate",
        prod_root=PROJ_PROD,
    conda: "clamp-analyses"
    resources: mem_mb=16000, runtime=60
    notebook:
        f"{PROJ_NB}/00_cytokines.ipynb"


rule projection_report_monocyte:
    input:
        agg=rules.aggregate_projections.output.agg,
        notebook=f"{PROJ_NB}/01_monocyte.ipynb",
    output:
        complete=touch(f"{PROJ_BIO}/02_monocyte/notebook.complete"),
    log:
        notebook=f"{PROJ_NB}/01_monocyte.executed.ipynb",
    params:
        group_label="monocyte",
        agg_dir=f"{PROJ_BIO}/aggregate",
        prod_root=PROJ_PROD,
    conda: "clamp-analyses"
    resources: mem_mb=16000, runtime=60
    notebook:
        f"{PROJ_NB}/01_monocyte.ipynb"


rule projection_report_placenta:
    input:
        agg=rules.aggregate_projections.output.agg,
        notebook=f"{PROJ_NB}/02_placenta.ipynb",
    output:
        complete=touch(f"{PROJ_BIO}/03_placenta/notebook.complete"),
    log:
        notebook=f"{PROJ_NB}/02_placenta.executed.ipynb",
    params:
        group_label="placenta",
        agg_dir=f"{PROJ_BIO}/aggregate",
        prod_root=PROJ_PROD,
    conda: "clamp-analyses"
    resources: mem_mb=16000, runtime=60
    notebook:
        f"{PROJ_NB}/02_placenta.ipynb"

PROJ_REPORTS = [
    f"{PROJ_BIO}/01_cytokines/notebook.complete",
    f"{PROJ_BIO}/02_monocyte/notebook.complete",
    f"{PROJ_BIO}/03_placenta/notebook.complete",
]


rule archs4_projections:
    input:
        PROJ_ALL_MECH,
        rules.aggregate_projections.output.agg,
        PROJ_REPORTS,
        rules.projection_report_pseudobulk_recovery.output.complete,
        rules.projection_report_gtex_ari.output.complete,
