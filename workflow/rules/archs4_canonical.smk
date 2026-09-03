# Canonical-prior model evaluation.  The three model artifacts are published
# inputs; this branch deliberately does not fit, adopt, or republish models.

A4_CAN_CFG = A4_CFG["canonical"]
A4_CAN_MODEL_NAME = A4_CAN_CFG["model_name"]
A4_CAN_FINAL_ROOT = A4_CAN_CFG["final_root"]
A4_CAN_COMPENDIA = list(A4_CAN_CFG["compendia"])
A4_CAN_DATABASES = list(A4_CAN_CFG["ora"]["databases"])
A4_CAN_COMPENDIUM_PATTERN = "|".join(A4_CAN_COMPENDIA)
A4_CAN_DATABASE_PATTERN = "|".join(A4_CAN_DATABASES)


def a4_can_database(wildcards):
    return A4_CFG["ora"]["databases"][wildcards.database]


rule ora_canonical:
    input:
        model=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/{A4_CAN_MODEL_NAME}.rds",
        z=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/Z.csv",
        manifest=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/manifest.json",
        provenance=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/adoption_validation.json",
        database=lambda wc: a4_can_database(wc)["path"],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora/{{database}}"),
    log:
        f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora/{{database}}.log"
    params:
        database_type=lambda wc: a4_can_database(wc)["type"],
        database_label=lambda wc: a4_can_database(wc)["label"],
        exclude=lambda wc: (
            f"--exclude-term-regex '{a4_can_database(wc)['exclude_term_regex']}'"
            if "exclude_term_regex" in a4_can_database(wc) else ""
        ),
        sheet=lambda wc: f"--sheet {a4_can_database(wc)['sheet']}" if "sheet" in a4_can_database(wc) else "",
        columns=lambda wc: (
            f"--term-column {a4_can_database(wc)['term_column']} "
            f"--gene-column {a4_can_database(wc)['gene_column']}"
            if "term_column" in a4_can_database(wc) else ""
        ),
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
    wildcard_constraints:
        dataset=A4_CAN_COMPENDIUM_PATTERN,
        database=A4_CAN_DATABASE_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --out-dir {output.ora_dir} "
        "--dataset {wildcards.dataset} --fraction 100 --seed 1 "
        "--model {A4_CAN_MODEL_NAME} --model-manifest {input.manifest} "
        "--database {wildcards.database} --database-label '{params.database_label}' "
        "--database-type {params.database_type} --database-path {input.database} "
        "{params.exclude} {params.sheet} {params.columns} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {A4_CFG[ora][min_size]} "
        "--max-size {A4_CFG[ora][max_size]} --pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} "
        "--qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} --fdr {A4_CFG[ora][fdr][0]} "
        "> {log} 2>&1"


rule aggregate_canonical_ora:
    input:
        ora=lambda wc: [
            f"{A4_CAN_FINAL_ROOT}/{wc.dataset}/ora/{database}"
            for database in A4_CAN_DATABASES
        ],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        coverage_long=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora.csv",
        cross_dataset=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora_cross.csv",
        panel_ready=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora_panel.csv",
    log:
        f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora_aggregate.log"
    params:
        databases=",".join(A4_CAN_DATABASES)
    resources:
        mem_mb=16000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --ora-root {A4_CAN_FINAL_ROOT}/{wildcards.dataset} "
        "--coverage-out {output.coverage_long} --cross-out {output.cross_dataset} "
        "--panel-out {output.panel_ready} --databases {params.databases} "
        "> {log} 2>&1"


rule archs4_canonical_ora:
    input:
        expand(rules.aggregate_canonical_ora.output, dataset=A4_CAN_COMPENDIA)
