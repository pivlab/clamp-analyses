A4_SAT_CFG = A4_CFG["saturation"]
A4_SAT_BIO = f"{A4_BIO}/01_saturation"
A4_SAT_NB = f"{A4_BIO_NB}/01_saturation"
A4_SAT_MODEL_ROOT = A4_SAT_CFG["model_root"]
A4_SAT_ORA_ROOT = A4_SAT_CFG["ora_root"]
A4_SAT_INPUT_ROOT = A4_SAT_CFG["input_root"]
A4_SAT_BASE_ROOT = A4_SAT_CFG["base_root"]
A4_SAT_FRACTIONS = [int(f) for f in A4_SAT_CFG["fractions"]]
A4_SAT_KS = [int(k) for k in A4_SAT_CFG["k_values"]]
A4_SAT_SEEDS = [int(s) for s in A4_SAT_CFG["seeds"]]
A4_SAT_RNG_SEEDS = [int(s) for s in A4_SAT_CFG["rng_seeds"]]
A4_SAT_DATABASES = list(A4_CFG["ora"]["databases"])
A4_SAT_FRACTION_PATTERN = "|".join(map(str, A4_SAT_FRACTIONS))
A4_SAT_K_PATTERN = "|".join(map(str, A4_SAT_KS))
A4_SAT_SEED_PATTERN = "|".join(map(str, A4_SAT_SEEDS))
A4_SAT_DATABASE_PATTERN = "|".join(A4_SAT_DATABASES)
A4_SAT_RETRY_FACTORS = [float(v) for v in A4_SAT_CFG["routing"]["retry_mem_factors"]]
A4_SAT_RETRY_CAP = int(A4_SAT_CFG["routing"]["retry_mem_cap_mb"])
A4_SAT_THREADS = int(A4_SAT_CFG["routing"]["threads_per_fit"])


# ============================================================
# Pathway recovery against model rank K, at several study-subsampling levels.
#
# Coverage asks how much data CLAMP needs; saturation asks how much model
# capacity, and whether the two interact.  Both are study-subsampled off the
# same draws, which is what makes reusing coverage's SVD/FBM per
# (fraction, seed) correct: K changes the model, not the data.  Nothing is
# re-subsampled and no SVD is recomputed.
#
# CLAMPbase is unsupervised and therefore prior-independent, so the fits the
# earlier Hallmark sweep left behind are valid here unchanged.  Only CLAMPfull
# is refit, with the same GO:BP prior coverage trains on, and scored against
# the same three databases.
# ============================================================


def a4_sat_input_cell(fraction, seed):
    """Coverage cell holding the subsample/SVD/FBM for this (fraction, seed).

    One cell backs all six K values.
    """
    level = A4_SAT_CFG["input_levels"][str(fraction)]
    leaf = A4_SAT_CFG["input_leaf"].format(fraction=fraction, seed=seed)
    return f"{A4_SAT_INPUT_ROOT}/{level}/{leaf}"


def a4_sat_base_cell(fraction, k, seed):
    leaf = A4_SAT_CFG["base_leaf"].format(fraction=fraction, k=k, seed=seed)
    return f"{A4_SAT_BASE_ROOT}/{leaf}"


def a4_sat_base_rds(fraction, k, seed):
    return f"{a4_sat_base_cell(fraction, k, seed)}/CLAMPbase.rds"


def a4_sat_base_z(fraction, k, seed):
    return f"{a4_sat_base_cell(fraction, k, seed)}/CLAMPbase/Z.csv"


def a4_sat_model_dir(fraction, k, seed):
    return (
        f"{A4_SAT_MODEL_ROOT}/rs{fraction}/k{k}/seed{seed}/"
        f"{A4_SAT_CFG['model_name']}"
    )


def a4_sat_validated(fraction, k, seed):
    return f"{A4_SAT_MODEL_ROOT}/rs{fraction}/k{k}/seed{seed}/validated.json"


def a4_sat_ora_dir(fraction, k, seed, model, database):
    return f"{A4_SAT_ORA_ROOT}/rs{fraction}/k{k}/seed{seed}/{model}/{database}"


def a4_sat_rng_seed(seed):
    return A4_SAT_RNG_SEEDS[int(seed) - 1]


def a4_sat_resource(fraction, k, key):
    return int(A4_SAT_CFG["resources"][str(fraction)][str(k)][key])


def a4_sat_retry_mem(fraction, k, attempt):
    factor = A4_SAT_RETRY_FACTORS[min(int(attempt), len(A4_SAT_RETRY_FACTORS)) - 1]
    return min(int(round(a4_sat_resource(fraction, k, "mem_mb") * factor)), A4_SAT_RETRY_CAP)


def a4_sat_db(wildcards):
    return A4_CFG["ora"]["databases"][wildcards.database]


# Cells excluded because the data cannot support them; see the comment on
# saturation.exclude_cells in archs4.yaml for why.  Filtering here keeps every
# downstream target list (models, ORAs, validation) consistent automatically.
A4_SAT_EXCLUDED = {
    (int(f), int(k), int(s))
    for f, k, s in A4_SAT_CFG.get("exclude_cells", [])
}

A4_SAT_CELLS = [
    (f, k, s)
    for f in A4_SAT_FRACTIONS
    for k in A4_SAT_KS
    for s in A4_SAT_SEEDS
    if (f, k, s) not in A4_SAT_EXCLUDED
]


def a4_sat_base_present():
    """Cells whose CLAMPbase already exists, so no refit is needed.

    Checked rather than hardcoded so the split stays correct as the few
    missing cells get filled in.
    """
    return [c for c in A4_SAT_CELLS if os.path.exists(a4_sat_base_z(*c))]


def a4_sat_base_missing():
    present = set(a4_sat_base_present())
    return [c for c in A4_SAT_CELLS if c not in present]


A4_SAT_BASE_ORA_PRESENT = [
    a4_sat_ora_dir(f, k, s, "CLAMPbase", db)
    for (f, k, s) in a4_sat_base_present()
    for db in A4_SAT_DATABASES
]
A4_SAT_BASE_ORA_ALL = [
    a4_sat_ora_dir(f, k, s, "CLAMPbase", db)
    for (f, k, s) in A4_SAT_CELLS
    for db in A4_SAT_DATABASES
]
A4_SAT_FULL_ORA = [
    a4_sat_ora_dir(f, k, s, "CLAMPfull", db)
    for (f, k, s) in A4_SAT_CELLS
    for db in A4_SAT_DATABASES
]
A4_SAT_VALIDATED = [a4_sat_validated(*c) for c in A4_SAT_CELLS]


rule fit_bp_saturation_base:
    """CLAMPbase at a forced K, for the cells the earlier sweep never produced.

    Inputs are ancient() because this rule is "produce if absent", not
    "produce if stale": the SVD/FBM are fixed upstream data and the CLAMPbase
    fits that already exist predate this code by months.  Without it, editing
    this script would queue refits of models that are already correct.
    """
    input:
        svd=lambda wc: ancient(f"{a4_sat_input_cell(wc.fraction, wc.seed)}/svd.rds"),
        fbm=lambda wc: ancient(f"{a4_sat_input_cell(wc.fraction, wc.seed)}/fbm_subsampled.bk"),
        subsample=lambda wc: ancient(f"{a4_sat_input_cell(wc.fraction, wc.seed)}/subsample_info.rds"),
        metadata=ancient(f"{A4_PROD}/00_preprocess/metadata_filtered.rds"),
        script=ancient("scripts/saturation/fit_clampbase.R"),
        runner=ancient("scripts/archs4/coverage/run_with_metrics.py"),
    output:
        base_rds=f"{A4_SAT_BASE_ROOT}/{A4_SAT_CFG['base_leaf']}/CLAMPbase.rds",
        base_z=f"{A4_SAT_BASE_ROOT}/{A4_SAT_CFG['base_leaf']}/CLAMPbase/Z.csv",
    log:
        f"{A4_SAT_BASE_ROOT}/{A4_SAT_CFG['base_leaf']}/clampbase.log"
    params:
        out_dir=lambda wc: a4_sat_base_cell(wc.fraction, wc.k, wc.seed),
        metrics_dir=lambda wc: f"{a4_sat_base_cell(wc.fraction, wc.k, wc.seed)}/attempts",
        rng_seed=lambda wc: a4_sat_rng_seed(wc.seed),
        initial_mem=lambda wc: a4_sat_resource(wc.fraction, wc.k, "mem_mb"),
    threads: A4_SAT_THREADS
    resources:
        mem_mb=lambda wc, attempt: a4_sat_retry_mem(wc.fraction, wc.k, attempt),
        fit_mem_mb=lambda wc, attempt: a4_sat_retry_mem(wc.fraction, wc.k, attempt),
        runtime=lambda wc: a4_sat_resource(wc.fraction, wc.k, "runtime"),
        fit_slots=1,
    wildcard_constraints:
        fraction=A4_SAT_FRACTION_PATTERN,
        k=A4_SAT_K_PATTERN,
        seed=A4_SAT_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.runner} --metrics-dir {params.metrics_dir} "
        "--requested-mem-mb {resources.mem_mb} --initial-mem-mb {params.initial_mem} "
        "--requested-cpus {threads} --requested-runtime-min {resources.runtime} "
        "--seed {params.rng_seed} -- "
        "Rscript {input.script} --fraction {wildcards.fraction} --clamp-k {wildcards.k} "
        "--seed-index {wildcards.seed} --seed {params.rng_seed} "
        "--metadata {input.metadata} --subsample-info {input.subsample} "
        "--fbm-backing {input.fbm} --svd {input.svd} "
        "--requested-mem-mb {resources.mem_mb} --requested-cpus {threads} "
        "--requested-runtime-min {resources.runtime} "
        "--out-dir {params.out_dir} > {log} 2>&1"


rule fit_bp_saturation_full:
    """CLAMPfull at a forced K with the GO:BP prior, chained from CLAMPbase."""
    input:
        svd=lambda wc: ancient(f"{a4_sat_input_cell(wc.fraction, wc.seed)}/svd.rds"),
        fbm=lambda wc: ancient(f"{a4_sat_input_cell(wc.fraction, wc.seed)}/fbm_subsampled.bk"),
        subsample=lambda wc: ancient(f"{a4_sat_input_cell(wc.fraction, wc.seed)}/subsample_info.rds"),
        metadata=ancient(f"{A4_PROD}/00_preprocess/metadata_filtered.rds"),
        base=lambda wc: ancient(a4_sat_base_rds(wc.fraction, wc.k, wc.seed)),
        prior=A4_SAT_CFG["prior_gmt"],
        script="scripts/archs4/coverage/fit_clampfull.R",
        runner="scripts/archs4/coverage/run_with_metrics.py",
    output:
        model_dir=directory(
            f"{A4_SAT_MODEL_ROOT}/rs{{fraction}}/k{{k}}/seed{{seed}}/{A4_SAT_CFG['model_name']}"
        ),
    log:
        f"{A4_SAT_MODEL_ROOT}/rs{{fraction}}/k{{k}}/seed{{seed}}/fit.log"
    params:
        metrics_dir=lambda wc: f"{A4_SAT_MODEL_ROOT}/rs{wc.fraction}/k{wc.k}/seed{wc.seed}/attempts",
        rng_seed=lambda wc: a4_sat_rng_seed(wc.seed),
        initial_mem=lambda wc: a4_sat_resource(wc.fraction, wc.k, "mem_mb"),
        max_iter=A4_SAT_CFG["max_iter"],
        multiplier=A4_CFG["clamp"]["multiplier"],
    threads: A4_SAT_THREADS
    resources:
        mem_mb=lambda wc, attempt: a4_sat_retry_mem(wc.fraction, wc.k, attempt),
        fit_mem_mb=lambda wc, attempt: a4_sat_retry_mem(wc.fraction, wc.k, attempt),
        runtime=lambda wc: a4_sat_resource(wc.fraction, wc.k, "runtime"),
        fit_slots=1,
    wildcard_constraints:
        fraction=A4_SAT_FRACTION_PATTERN,
        k=A4_SAT_K_PATTERN,
        seed=A4_SAT_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.runner} --metrics-dir {params.metrics_dir} "
        "--requested-mem-mb {resources.mem_mb} --initial-mem-mb {params.initial_mem} "
        "--requested-cpus {threads} --requested-runtime-min {resources.runtime} "
        "--seed {params.rng_seed} --prior {input.prior} -- "
        "Rscript {input.script} --dataset archs4 --fraction {wildcards.fraction} "
        "--seed-index {wildcards.seed} --seed {params.rng_seed} --metadata {input.metadata} "
        "--subsample-info {input.subsample} --fbm-backing {input.fbm} --svd {input.svd} "
        "--clamp-k {wildcards.k} --base-model {input.base} --prior-gmt {input.prior} "
        "--requested-mem-mb {resources.mem_mb} --requested-cpus {threads} "
        "--requested-runtime-min {resources.runtime} "
        "--max-iter {params.max_iter} --multiplier {params.multiplier} "
        "--out-dir {output.model_dir} > {log} 2>&1"


rule validate_bp_saturation_model:
    input:
        model_dir=lambda wc: a4_sat_model_dir(wc.fraction, wc.k, wc.seed),
        script="scripts/archs4/coverage/validate_model.py",
    output:
        f"{A4_SAT_MODEL_ROOT}/rs{{fraction}}/k{{k}}/seed{{seed}}/validated.json"
    resources:
        mem_mb=2000,
        runtime=60,
    wildcard_constraints:
        fraction=A4_SAT_FRACTION_PATTERN,
        k=A4_SAT_K_PATTERN,
        seed=A4_SAT_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.script} --model-dir {input.model_dir} --output {output}"


rule ora_bp_saturation_base:
    """Score a CLAMPbase.  No prior involved, so no refit is ever needed."""
    input:
        z=lambda wc: ancient(a4_sat_base_z(wc.fraction, wc.k, wc.seed)),
        subsample=lambda wc: ancient(f"{a4_sat_input_cell(wc.fraction, wc.seed)}/subsample_info.rds"),
        database=lambda wc: a4_sat_db(wc)["path"],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(
            f"{A4_SAT_ORA_ROOT}/rs{{fraction}}/k{{k}}/seed{{seed}}/CLAMPbase/{{database}}"
        ),
    log:
        f"{A4_SAT_ORA_ROOT}/rs{{fraction}}/k{{k}}/seed{{seed}}/CLAMPbase/{{database}}.log"
    params:
        database_type=lambda wc: a4_sat_db(wc)["type"],
        database_label=lambda wc: a4_sat_db(wc)["label"],
        exclude=lambda wc: (
            f"--exclude-term-regex '{a4_sat_db(wc)['exclude_term_regex']}'"
            if "exclude_term_regex" in a4_sat_db(wc) else ""
        ),
        sheet=lambda wc: f"--sheet {a4_sat_db(wc)['sheet']}" if "sheet" in a4_sat_db(wc) else "",
        columns=lambda wc: (
            f"--term-column {a4_sat_db(wc)['term_column']} --gene-column {a4_sat_db(wc)['gene_column']}"
            if "term_column" in a4_sat_db(wc) else ""
        ),
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        ora_mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
        ora_slots=1,
    wildcard_constraints:
        fraction=A4_SAT_FRACTION_PATTERN,
        k=A4_SAT_K_PATTERN,
        seed=A4_SAT_SEED_PATTERN,
        database=A4_SAT_DATABASE_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --out-dir {output.ora_dir} "
        "--dataset archs4 --fraction {wildcards.fraction} --seed {wildcards.seed} "
        "--model CLAMPbase --subsample-info {input.subsample} --database {wildcards.database} "
        "--database-label '{params.database_label}' --database-type {params.database_type} "
        "--database-path {input.database} {params.exclude} {params.sheet} {params.columns} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {A4_CFG[ora][min_size]} "
        "--max-size {A4_CFG[ora][max_size]} --pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} "
        "--qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} --fdr {A4_CFG[ora][fdr][0]} "
        "> {log} 2>&1"


rule ora_bp_saturation_full:
    input:
        validated=lambda wc: a4_sat_validated(wc.fraction, wc.k, wc.seed),
        database=lambda wc: a4_sat_db(wc)["path"],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(
            f"{A4_SAT_ORA_ROOT}/rs{{fraction}}/k{{k}}/seed{{seed}}/CLAMPfull/{{database}}"
        ),
    log:
        f"{A4_SAT_ORA_ROOT}/rs{{fraction}}/k{{k}}/seed{{seed}}/CLAMPfull/{{database}}.log"
    params:
        z=lambda wc: f"{a4_sat_model_dir(wc.fraction, wc.k, wc.seed)}/Z.csv",
        manifest=lambda wc: f"{a4_sat_model_dir(wc.fraction, wc.k, wc.seed)}/manifest.json",
        database_type=lambda wc: a4_sat_db(wc)["type"],
        database_label=lambda wc: a4_sat_db(wc)["label"],
        exclude=lambda wc: (
            f"--exclude-term-regex '{a4_sat_db(wc)['exclude_term_regex']}'"
            if "exclude_term_regex" in a4_sat_db(wc) else ""
        ),
        sheet=lambda wc: f"--sheet {a4_sat_db(wc)['sheet']}" if "sheet" in a4_sat_db(wc) else "",
        columns=lambda wc: (
            f"--term-column {a4_sat_db(wc)['term_column']} --gene-column {a4_sat_db(wc)['gene_column']}"
            if "term_column" in a4_sat_db(wc) else ""
        ),
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        ora_mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
        ora_slots=1,
    wildcard_constraints:
        fraction=A4_SAT_FRACTION_PATTERN,
        k=A4_SAT_K_PATTERN,
        seed=A4_SAT_SEED_PATTERN,
        database=A4_SAT_DATABASE_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {params.z} --out-dir {output.ora_dir} "
        "--dataset archs4 --fraction {wildcards.fraction} --seed {wildcards.seed} "
        "--model CLAMPfull --model-manifest {params.manifest} --database {wildcards.database} "
        "--database-label '{params.database_label}' --database-type {params.database_type} "
        "--database-path {input.database} {params.exclude} {params.sheet} {params.columns} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {A4_CFG[ora][min_size]} "
        "--max-size {A4_CFG[ora][max_size]} --pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} "
        "--qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} --fdr {A4_CFG[ora][fdr][0]} "
        "> {log} 2>&1"


rule aggregate_bp_saturation:
    input:
        A4_SAT_FULL_ORA + A4_SAT_BASE_ORA_ALL,
        script="scripts/saturation/aggregate_saturation.R",
    output:
        saturation_long=f"{A4_SAT_BIO}/saturation_long.csv",
        panel_ready=f"{A4_SAT_BIO}/saturation_panel_ready.csv",
    log:
        f"{A4_SAT_BIO}/aggregate.log"
    resources:
        mem_mb=16000,
        runtime=60,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --ora-root {A4_SAT_ORA_ROOT} "
        "--saturation-out {output.saturation_long} --panel-out {output.panel_ready} "
        "> {log} 2>&1"


rule saturation_report_archs4:
    input:
        saturation_long=rules.aggregate_bp_saturation.output.saturation_long,
        panel_ready=rules.aggregate_bp_saturation.output.panel_ready,
        notebook=f"{A4_SAT_NB}/00_saturation.ipynb",
    output:
        complete=touch(f"{A4_SAT_BIO}/notebook.complete"),
    log:
        notebook=f"{A4_SAT_NB}/00_saturation.executed.ipynb",
    params:
        saturation_dir=A4_SAT_BIO,
        ora_root=A4_SAT_ORA_ROOT,
    resources:
        mem_mb=16000,
        runtime=60,
    conda: "clamp-analyses"
    notebook:
        f"{A4_SAT_NB}/00_saturation.ipynb"


# ---- aggregate targets ----

rule saturation_bp_base_ora:
    """Phase 1: score every CLAMPbase already on disk.  No fitting."""
    input:
        A4_SAT_BASE_ORA_PRESENT,


rule saturation_bp_base_models:
    """Phase 2: fill in the CLAMPbase cells the earlier sweep never produced."""
    input:
        [a4_sat_base_rds(*c) for c in a4_sat_base_missing()],


rule saturation_bp_models:
    input:
        A4_SAT_VALIDATED,


rule saturation_bp_ora:
    input:
        A4_SAT_FULL_ORA + A4_SAT_BASE_ORA_ALL,


rule saturation_bp:
    input:
        rules.saturation_report_archs4.output.complete,


rule archs4_saturation:
    input:
        rules.saturation_report_archs4.output.complete,
