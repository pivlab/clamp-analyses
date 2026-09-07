A4_TRAIT_CFG = A4_CFG["traits"]


# ============================================================
# Coverage: aggregate GLS trait-recovery summaries (CLAMPfull_bp and
# CLAMPbase, native full rank) and render the report notebook, alongside the
# pathway-coverage report in the same 00_coverage output/notebook folder.
# ============================================================

rule aggregate_archs4_coverage_traits:
    input:
        clampfull_dir=A4_TRAIT_CFG["coverage"]["clampfull_bp_dir"],
        clampbase_dir=A4_TRAIT_CFG["coverage"]["clampbase_dir"],
        finals_clampfull_dir=A4_TRAIT_CFG["finals"]["clampfull_bp_dir"],
        finals_clampbase_dir=A4_TRAIT_CFG["finals"]["clampbase_dir"],
        script="scripts/archs4/traits/aggregate_coverage_traits.R",
    output:
        traits_long=f"{A4_COV_BIO}/coverage_trait_recovery.csv",
        finals_long=f"{A4_COV_BIO}/final_trait_recovery.csv",
    log:
        f"{A4_COV_BIO}/traits_aggregate.log"
    params:
        fdr=A4_TRAIT_CFG["fdr"],
    resources:
        mem_mb=8000,
        runtime=30,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --clampfull-dir {input.clampfull_dir} "
        "--clampbase-dir {input.clampbase_dir} "
        "--finals-clampfull-dir {input.finals_clampfull_dir} "
        "--finals-clampbase-dir {input.finals_clampbase_dir} --fdr {params.fdr} "
        "--traits-out {output.traits_long} --finals-out {output.finals_long} > {log} 2>&1"


rule coverage_traits_report_archs4:
    input:
        traits_long=rules.aggregate_archs4_coverage_traits.output.traits_long,
        finals_long=rules.aggregate_archs4_coverage_traits.output.finals_long,
        notebook=f"{A4_COV_NB}/01_coverage_traits.ipynb",
    output:
        complete=touch(f"{A4_COV_BIO}/traits.complete"),
    log:
        notebook=f"{A4_COV_NB}/01_coverage_traits.executed.ipynb",
    resources:
        mem_mb=8000,
        runtime=30,
    conda: "clamp-analyses"
    notebook:
        f"{A4_COV_NB}/01_coverage_traits.ipynb"


# ============================================================
# Saturation: aggregate GLS trait-recovery summaries (CLAMPfull_bp and
# CLAMPbase, forced K=1728) and render the report notebook, alongside the
# pathway-saturation report in the same 01_saturation output/notebook folder.
# ============================================================

rule aggregate_archs4_saturation_traits:
    input:
        clampfull_dir=A4_TRAIT_CFG["saturation"]["clampfull_bp_dir"],
        clampbase_dir=A4_TRAIT_CFG["saturation"]["clampbase_dir"],
        script="scripts/archs4/traits/aggregate_saturation_traits.R",
    output:
        traits_long=f"{A4_SAT_BIO}/saturation_trait_recovery_k1728.csv",
    log:
        f"{A4_SAT_BIO}/traits_aggregate.log"
    params:
        fdr=A4_TRAIT_CFG["fdr"],
    resources:
        mem_mb=8000,
        runtime=30,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --clampfull-dir {input.clampfull_dir} "
        "--clampbase-dir {input.clampbase_dir} --fdr {params.fdr} "
        "--traits-out {output.traits_long} > {log} 2>&1"


rule saturation_traits_report_archs4:
    input:
        traits_long=rules.aggregate_archs4_saturation_traits.output.traits_long,
        notebook=f"{A4_SAT_NB}/01_saturation_traits.ipynb",
    output:
        complete=touch(f"{A4_SAT_BIO}/traits.complete"),
    log:
        notebook=f"{A4_SAT_NB}/01_saturation_traits.executed.ipynb",
    resources:
        mem_mb=8000,
        runtime=30,
    conda: "clamp-analyses"
    notebook:
        f"{A4_SAT_NB}/01_saturation_traits.ipynb"


rule archs4_traits:
    input:
        rules.coverage_traits_report_archs4.output.complete,
        rules.saturation_traits_report_archs4.output.complete,
