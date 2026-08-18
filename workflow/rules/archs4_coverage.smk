A4_COV_CFG = A4_CFG["coverage"]
A4_COV_BIO = f"{A4_BIO}/00_coverage"
A4_COV_NB = f"{A4_BIO_NB}/00_coverage"
A4_COV_MODEL_ROOT = A4_COV_CFG["model_root"]
A4_COV_ORA_ROOT = A4_COV_CFG["ora_root"]
A4_COV_INPUT_ROOT = A4_COV_CFG["input_root"]
A4_COV_LEVELS = [int(level) for level in A4_COV_CFG["levels"]]
A4_COV_SEEDS = [int(seed) for seed in A4_COV_CFG["seeds"]]
# At 100% there is no subsampling, so every seed's SVD/CLAMPbase (and, for the
# comparators, CLAMPbase itself) is identical or not seed-varied to begin
# with -- only this one seed's CLAMPfull is ever used downstream. Fitting and
# scoring the other seeds there would just recompute discarded data.
A4_COV_REFERENCE_SEED = 1
A4_COV_RNG_SEEDS = [int(seed) for seed in A4_COV_CFG["rng_seeds"]]
A4_COV_DATABASES = list(A4_CFG["ora"]["databases"])
A4_COV_COMPARATORS = list(A4_COV_CFG["comparators"])
A4_COV_LEVEL_PATTERN = "|".join(map(str, A4_COV_LEVELS))
A4_COV_SEED_PATTERN = "|".join(map(str, A4_COV_SEEDS))
A4_COV_DATABASE_PATTERN = "|".join(A4_COV_DATABASES)
A4_COV_COMPARATOR_PATTERN = "|".join(A4_COV_COMPARATORS)
A4_COV_RETRY_FACTORS = [float(value) for value in A4_COV_CFG["routing"]["retry_mem_factors"]]
A4_COV_RETRY_CAP = int(A4_COV_CFG["routing"]["retry_mem_cap_mb"])


def a4_cov_cell(fraction, seed):
    fraction = str(fraction)
    level_dir = A4_COV_CFG["levels"][fraction]
    return (
        f"{A4_COV_INPUT_ROOT}/{level_dir}/"
        f"study_coverage_rs{fraction}_seed_{seed}"
    )


def a4_cov_model_dir(dataset, fraction, seed):
    return (
        f"{A4_COV_MODEL_ROOT}/{dataset}/rs{fraction}/seed{seed}/"
        f"{A4_COV_CFG['model_name']}"
    )


def a4_cov_validated(dataset, fraction, seed):
    return f"{A4_COV_MODEL_ROOT}/{dataset}/rs{fraction}/seed{seed}/validated.json"


def a4_cov_ora_dir(dataset, fraction, seed, model, database):
    return f"{A4_COV_ORA_ROOT}/{dataset}/rs{fraction}/seed{seed}/{model}/{database}"


def a4_cov_rng_seed(seed):
    return A4_COV_RNG_SEEDS[int(seed) - 1]


def a4_cov_initial_mem(fraction):
    return int(A4_COV_CFG["resources"][str(fraction)]["mem_mb"])


def a4_cov_retry_mem(fraction, attempt):
    factor = A4_COV_RETRY_FACTORS[min(int(attempt), len(A4_COV_RETRY_FACTORS)) - 1]
    return min(int(round(a4_cov_initial_mem(fraction) * factor)), A4_COV_RETRY_CAP)


def a4_cov_db(wildcards):
    return A4_CFG["ora"]["databases"][wildcards.database]


def a4_cov_base_z(wildcards):
    if wildcards.dataset == "archs4":
        return f"{a4_cov_cell(wildcards.fraction, wildcards.seed)}/CLAMPbase/Z.csv"
    return A4_COV_CFG["comparators"][wildcards.dataset]["base_z"]


def a4_cov_sample_arg(wildcards, full):
    if full:
        manifest = f"{a4_cov_model_dir(wildcards.dataset, wildcards.fraction, wildcards.seed)}/manifest.json"
        return f"--model-manifest {manifest}"
    if wildcards.dataset == "archs4":
        return f"--subsample-info {a4_cov_cell(wildcards.fraction, wildcards.seed)}/subsample_info.rds"
    samples = A4_COV_CFG["comparators"][wildcards.dataset]["samples_rds"]
    return f"--samples-rds {samples}"


def a4_cov_seeds_for(fraction):
    """Seeds to fit/score by default at this fraction.

    Below 100% each seed draws its own random subsample, so its SVD,
    CLAMPbase and CLAMPfull are all genuinely independent and worth scoring
    separately. At 100% there is nothing left to subsample -- every seed
    would fit from the same full compendium -- so only the reference seed is
    requested by default.
    """
    return [A4_COV_REFERENCE_SEED] if int(fraction) == 100 else A4_COV_SEEDS


A4_COV_ARCH_VALIDATED = [
    a4_cov_validated("archs4", fraction, seed)
    for fraction in A4_COV_LEVELS
    for seed in a4_cov_seeds_for(fraction)
]
A4_COV_COMPARATOR_VALIDATED = [
    a4_cov_validated(dataset, 100, seed)
    for dataset in A4_COV_COMPARATORS
    for seed in a4_cov_seeds_for(100)
]
A4_COV_ARCH_FULL_ORA = [
    a4_cov_ora_dir("archs4", fraction, seed, "CLAMPfull", database)
    for fraction in A4_COV_LEVELS
    for seed in a4_cov_seeds_for(fraction)
    for database in A4_COV_DATABASES
]
A4_COV_COMPARATOR_FULL_ORA = [
    a4_cov_ora_dir(dataset, 100, seed, "CLAMPfull", database)
    for dataset in A4_COV_COMPARATORS
    for seed in a4_cov_seeds_for(100)
    for database in A4_COV_DATABASES
]
A4_COV_ARCH_BASE_ORA = [
    a4_cov_ora_dir("archs4", fraction, seed, "CLAMPbase", database)
    for fraction in A4_COV_LEVELS
    for seed in a4_cov_seeds_for(fraction)
    for database in A4_COV_DATABASES
]
A4_COV_COMPARATOR_BASE_ORA = [
    a4_cov_ora_dir(dataset, 100, 0, "CLAMPbase", database)
    for dataset in A4_COV_COMPARATORS
    for database in A4_COV_DATABASES
]
A4_COV_FULL_ORA = A4_COV_ARCH_FULL_ORA + A4_COV_COMPARATOR_FULL_ORA
A4_COV_BASE_ORA = A4_COV_ARCH_BASE_ORA + A4_COV_COMPARATOR_BASE_ORA


rule fit_bp_archs4_coverage:
    input:
        metadata=f"{A4_PROD}/00_preprocess/metadata_filtered.rds",
        subsample=lambda wc: f"{a4_cov_cell(wc.fraction, wc.seed)}/subsample_info.rds",
        fbm=lambda wc: f"{a4_cov_cell(wc.fraction, wc.seed)}/fbm_subsampled.bk",
        svd=lambda wc: f"{a4_cov_cell(wc.fraction, wc.seed)}/svd.rds",
        k=lambda wc: f"{a4_cov_cell(wc.fraction, wc.seed)}/CLAMP_K.rds",
        base=lambda wc: f"{a4_cov_cell(wc.fraction, wc.seed)}/CLAMPbase.rds",
        prior=A4_COV_CFG["prior_gmt"],
        script="scripts/archs4/coverage/fit_clampfull.R",
        runner="scripts/archs4/coverage/run_with_metrics.py",
    output:
        model_dir=directory(
            f"{A4_COV_MODEL_ROOT}/archs4/rs{{fraction}}/seed{{seed}}/{A4_COV_CFG['model_name']}"
        ),
    log:
        f"{A4_COV_MODEL_ROOT}/archs4/rs{{fraction}}/seed{{seed}}/fit.log"
    params:
        rng_seed=lambda wc: a4_cov_rng_seed(wc.seed),
        initial_mem=lambda wc: a4_cov_initial_mem(wc.fraction),
        max_iter=A4_CFG["clamp"]["max_iter"],
        multiplier=A4_CFG["clamp"]["multiplier"],
    threads: 16
    resources:
        mem_mb=lambda wc, attempt: a4_cov_retry_mem(wc.fraction, attempt),
        fit_mem_mb=lambda wc, attempt: a4_cov_retry_mem(wc.fraction, attempt),
        runtime=lambda wc: int(A4_COV_CFG["resources"][str(wc.fraction)]["runtime"]),
        large_fit_slot=lambda wc: 1 if int(wc.fraction) >= 50 else 0,
        fit_slots=1,
    wildcard_constraints:
        fraction=A4_COV_LEVEL_PATTERN,
        seed=A4_COV_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.runner} --metrics-dir {A4_COV_MODEL_ROOT}/archs4/rs{wildcards.fraction}/seed{wildcards.seed}/attempts "
        "--requested-mem-mb {resources.mem_mb} --initial-mem-mb {params.initial_mem} "
        "--requested-cpus {threads} --requested-runtime-min {resources.runtime} --seed {params.rng_seed} --prior {input.prior} -- "
        "Rscript {input.script} --dataset archs4 --fraction {wildcards.fraction} "
        "--seed-index {wildcards.seed} --seed {params.rng_seed} --metadata {input.metadata} "
        "--subsample-info {input.subsample} --fbm-backing {input.fbm} --svd {input.svd} "
        "--k {input.k} --base-model {input.base} --prior-gmt {input.prior} "
        "--requested-mem-mb {resources.mem_mb} --requested-cpus {threads} "
        "--requested-runtime-min {resources.runtime} "
        "--max-iter {params.max_iter} --multiplier {params.multiplier} "
        "--out-dir {output.model_dir} > {log} 2>&1"


rule fit_bp_comparator_coverage:
    input:
        fbm=lambda wc: A4_COV_CFG["comparators"][wc.dataset]["fbm_rds"],
        genes=lambda wc: A4_COV_CFG["comparators"][wc.dataset]["genes_rds"],
        samples=lambda wc: A4_COV_CFG["comparators"][wc.dataset]["samples_rds"],
        svd=lambda wc: A4_COV_CFG["comparators"][wc.dataset]["svd"],
        k=lambda wc: A4_COV_CFG["comparators"][wc.dataset]["k"],
        base=lambda wc: A4_COV_CFG["comparators"][wc.dataset]["base"],
        prior=A4_COV_CFG["prior_gmt"],
        script="scripts/archs4/coverage/fit_clampfull.R",
        runner="scripts/archs4/coverage/run_with_metrics.py",
    output:
        model_dir=directory(
            f"{A4_COV_MODEL_ROOT}/{{dataset}}/rs100/seed{{seed}}/{A4_COV_CFG['model_name']}"
        ),
    log:
        f"{A4_COV_MODEL_ROOT}/{{dataset}}/rs100/seed{{seed}}/fit.log"
    params:
        rng_seed=lambda wc: a4_cov_rng_seed(wc.seed),
        max_iter=lambda wc: A4_COV_CFG["comparators"][wc.dataset]["max_iter"],
        multiplier=A4_CFG["clamp"]["multiplier"],
    threads: 8
    resources:
        mem_mb=lambda wc, attempt: min(
            int(round(32000 * A4_COV_RETRY_FACTORS[min(int(attempt), len(A4_COV_RETRY_FACTORS)) - 1])),
            A4_COV_RETRY_CAP,
        ),
        fit_mem_mb=lambda wc, attempt: min(
            int(round(32000 * A4_COV_RETRY_FACTORS[min(int(attempt), len(A4_COV_RETRY_FACTORS)) - 1])),
            A4_COV_RETRY_CAP,
        ),
        runtime=5760,
        large_fit_slot=0,
        fit_slots=1,
    wildcard_constraints:
        dataset=A4_COV_COMPARATOR_PATTERN,
        seed=A4_COV_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.runner} --metrics-dir {A4_COV_MODEL_ROOT}/{wildcards.dataset}/rs100/seed{wildcards.seed}/attempts "
        "--requested-mem-mb {resources.mem_mb} --initial-mem-mb 32000 "
        "--requested-cpus {threads} --requested-runtime-min {resources.runtime} --seed {params.rng_seed} --prior {input.prior} -- "
        "Rscript {input.script} --dataset {wildcards.dataset} --fraction 100 "
        "--seed-index {wildcards.seed} --seed {params.rng_seed} --fbm-rds {input.fbm} "
        "--genes-rds {input.genes} --samples-rds {input.samples} --svd {input.svd} "
        "--k {input.k} --base-model {input.base} --prior-gmt {input.prior} "
        "--requested-mem-mb {resources.mem_mb} --requested-cpus {threads} "
        "--requested-runtime-min {resources.runtime} "
        "--max-iter {params.max_iter} --multiplier {params.multiplier} "
        "--out-dir {output.model_dir} > {log} 2>&1"


rule validate_bp_coverage_model:
    input:
        model_dir=lambda wc: a4_cov_model_dir(wc.dataset, wc.fraction, wc.seed),
        script="scripts/archs4/coverage/validate_model.py",
    output:
        f"{A4_COV_MODEL_ROOT}/{{dataset}}/rs{{fraction}}/seed{{seed}}/validated.json"
    resources:
        mem_mb=2000,
        runtime=1440,
    wildcard_constraints:
        dataset="archs4|" + A4_COV_COMPARATOR_PATTERN,
        fraction=A4_COV_LEVEL_PATTERN,
        seed=A4_COV_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.script} --model-dir {input.model_dir} --output {output}"


rule ora_bp_full_coverage:
    input:
        validated=lambda wc: a4_cov_validated(wc.dataset, wc.fraction, wc.seed),
        database=lambda wc: a4_cov_db(wc)["path"],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(
            f"{A4_COV_ORA_ROOT}/{{dataset}}/rs{{fraction}}/seed{{seed}}/CLAMPfull/{{database}}"
        ),
    log:
        f"{A4_COV_ORA_ROOT}/{{dataset}}/rs{{fraction}}/seed{{seed}}/CLAMPfull/{{database}}.log"
    params:
        z=lambda wc: f"{a4_cov_model_dir(wc.dataset, wc.fraction, wc.seed)}/Z.csv",
        sample_arg=lambda wc: a4_cov_sample_arg(wc, True),
        database_type=lambda wc: a4_cov_db(wc)["type"],
        database_label=lambda wc: a4_cov_db(wc)["label"],
        exclude=lambda wc: (
            f"--exclude-term-regex '{a4_cov_db(wc)['exclude_term_regex']}'"
            if "exclude_term_regex" in a4_cov_db(wc) else ""
        ),
        sheet=lambda wc: f"--sheet {a4_cov_db(wc)['sheet']}" if "sheet" in a4_cov_db(wc) else "",
        columns=lambda wc: (
            f"--term-column {a4_cov_db(wc)['term_column']} --gene-column {a4_cov_db(wc)['gene_column']}"
            if "term_column" in a4_cov_db(wc) else ""
        ),
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        ora_mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
        ora_slots=1,
    wildcard_constraints:
        dataset="archs4|" + A4_COV_COMPARATOR_PATTERN,
        fraction=A4_COV_LEVEL_PATTERN,
        seed=A4_COV_SEED_PATTERN,
        database=A4_COV_DATABASE_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {params.z} --out-dir {output.ora_dir} "
        "--dataset {wildcards.dataset} --fraction {wildcards.fraction} --seed {wildcards.seed} "
        "--model CLAMPfull {params.sample_arg} --database {wildcards.database} "
        "--database-label '{params.database_label}' --database-type {params.database_type} "
        "--database-path {input.database} {params.exclude} {params.sheet} {params.columns} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {A4_CFG[ora][min_size]} "
        "--max-size {A4_CFG[ora][max_size]} --pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} "
        "--qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} --fdr {A4_CFG[ora][fdr][0]} "
        "> {log} 2>&1"


rule ora_bp_base_coverage:
    input:
        z=a4_cov_base_z,
        database=lambda wc: a4_cov_db(wc)["path"],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(
            f"{A4_COV_ORA_ROOT}/{{dataset}}/rs{{fraction}}/seed{{seed}}/CLAMPbase/{{database}}"
        ),
    log:
        f"{A4_COV_ORA_ROOT}/{{dataset}}/rs{{fraction}}/seed{{seed}}/CLAMPbase/{{database}}.log"
    params:
        sample_arg=lambda wc: a4_cov_sample_arg(wc, False),
        database_type=lambda wc: a4_cov_db(wc)["type"],
        database_label=lambda wc: a4_cov_db(wc)["label"],
        exclude=lambda wc: (
            f"--exclude-term-regex '{a4_cov_db(wc)['exclude_term_regex']}'"
            if "exclude_term_regex" in a4_cov_db(wc) else ""
        ),
        sheet=lambda wc: f"--sheet {a4_cov_db(wc)['sheet']}" if "sheet" in a4_cov_db(wc) else "",
        columns=lambda wc: (
            f"--term-column {a4_cov_db(wc)['term_column']} --gene-column {a4_cov_db(wc)['gene_column']}"
            if "term_column" in a4_cov_db(wc) else ""
        ),
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        ora_mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
        ora_slots=1,
    wildcard_constraints:
        dataset="archs4|" + A4_COV_COMPARATOR_PATTERN,
        fraction=A4_COV_LEVEL_PATTERN,
        seed="0|" + A4_COV_SEED_PATTERN,
        database=A4_COV_DATABASE_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --out-dir {output.ora_dir} "
        "--dataset {wildcards.dataset} --fraction {wildcards.fraction} --seed {wildcards.seed} "
        "--model CLAMPbase {params.sample_arg} --database {wildcards.database} "
        "--database-label '{params.database_label}' --database-type {params.database_type} "
        "--database-path {input.database} {params.exclude} {params.sheet} {params.columns} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {A4_CFG[ora][min_size]} "
        "--max-size {A4_CFG[ora][max_size]} --pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} "
        "--qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} --fdr {A4_CFG[ora][fdr][0]} "
        "> {log} 2>&1"


rule aggregate_bp_coverage:
    input:
        A4_COV_FULL_ORA + A4_COV_BASE_ORA,
        script="scripts/archs4/coverage/aggregate_coverage.R",
    output:
        coverage_long=f"{A4_COV_BIO}/coverage_long.csv",
        cross_dataset=f"{A4_COV_BIO}/cross_dataset_coverage.csv",
        panel_ready=f"{A4_COV_BIO}/coverage_panel_ready.csv",
    log:
        f"{A4_COV_BIO}/aggregate.log"
    resources:
        mem_mb=16000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --ora-root {A4_COV_ORA_ROOT} "
        "--coverage-out {output.coverage_long} --cross-out {output.cross_dataset} "
        "--panel-out {output.panel_ready} "
        "> {log} 2>&1"


rule coverage_report_archs4:
    input:
        coverage_long=rules.aggregate_bp_coverage.output.coverage_long,
        cross_dataset=rules.aggregate_bp_coverage.output.cross_dataset,
        panel_ready=rules.aggregate_bp_coverage.output.panel_ready,
        notebook=f"{A4_COV_NB}/00_coverage.ipynb",
    output:
        complete=touch(f"{A4_COV_BIO}/notebook.complete"),
    log:
        notebook=f"{A4_COV_NB}/00_coverage.executed.ipynb",
    params:
        coverage_dir=A4_COV_BIO,
    resources:
        mem_mb=16000,
        runtime=1440,
    conda: "clamp-analyses"
    notebook:
        f"{A4_COV_NB}/00_coverage.ipynb"


rule coverage_bp_archs4_models:
    input:
        A4_COV_ARCH_VALIDATED,


rule coverage_bp_comparator_models:
    input:
        A4_COV_COMPARATOR_VALIDATED,


rule coverage_bp_models:
    input:
        A4_COV_ARCH_VALIDATED + A4_COV_COMPARATOR_VALIDATED,


rule coverage_bp_ora:
    input:
        A4_COV_FULL_ORA + A4_COV_BASE_ORA,


rule coverage_bp_archs4_base_ora:
    input:
        A4_COV_ARCH_BASE_ORA,


rule coverage_bp_comparator_ora:
    input:
        A4_COV_COMPARATOR_FULL_ORA + A4_COV_COMPARATOR_BASE_ORA,


rule coverage_bp_comparator_base_ora:
    input:
        A4_COV_COMPARATOR_BASE_ORA,


rule coverage_bp:
    input:
        rules.coverage_report_archs4.output.complete,


rule archs4_coverage:
    input:
        rules.coverage_report_archs4.output.complete,
