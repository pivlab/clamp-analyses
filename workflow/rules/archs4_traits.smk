A4_TRAIT_CFG = A4_CFG["traits"]


# ============================================================
# Coverage: aggregate GLS trait-recovery summaries (CLAMPfull_bp and
# CLAMPbase, native full rank) and render the report notebook, alongside the
# pathway-coverage report in the same 00_coverage output/notebook folder.
# ============================================================

rule aggregate_archs4_coverage_traits:
    input:
        # The per-model GLS summaries phenoplier.smk writes into the
        # `traits.*` directories below (one per coverage cell / final model),
        # so aggregating pulls the GLS runs, and behind them the model fits.
        clampfull=PHENOPLIER_COV_FULL,
        clampbase=PHENOPLIER_COV_BASE,
        finals_clampfull=PHENOPLIER_FIN_FULL,
        finals_clampbase=PHENOPLIER_FIN_BASE,
        script="scripts/archs4/traits/aggregate_coverage_traits.R",
    output:
        traits_long=f"{A4_COV_BIO}/coverage_trait_recovery.csv",
        finals_long=f"{A4_COV_BIO}/final_trait_recovery.csv",
    log:
        f"{A4_COV_BIO}/traits_aggregate.log"
    params:
        # The R scripts scan a directory; it is the one the summaries above
        # live in (archs4.yaml: traits.*).
        clampfull_dir=lambda wc, input: os.path.dirname(input.clampfull[0]),
        clampbase_dir=lambda wc, input: os.path.dirname(input.clampbase[0]),
        finals_clampfull_dir=lambda wc, input: os.path.dirname(input.finals_clampfull[0]),
        finals_clampbase_dir=lambda wc, input: os.path.dirname(input.finals_clampbase[0]),
        fdr=A4_TRAIT_CFG["fdr"],
    resources:
        mem_mb=8000,
        runtime=30,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --clampfull-dir {params.clampfull_dir} "
        "--clampbase-dir {params.clampbase_dir} "
        "--finals-clampfull-dir {params.finals_clampfull_dir} "
        "--finals-clampbase-dir {params.finals_clampbase_dir} --fdr {params.fdr} "
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
        clampfull=PHENOPLIER_SAT_FULL,
        clampbase=PHENOPLIER_SAT_BASE,
        script="scripts/archs4/traits/aggregate_saturation_traits.R",
    output:
        traits_long=f"{A4_SAT_BIO}/saturation_trait_recovery_k1728.csv",
    log:
        f"{A4_SAT_BIO}/traits_aggregate.log"
    params:
        clampfull_dir=lambda wc, input: os.path.dirname(input.clampfull[0]),
        clampbase_dir=lambda wc, input: os.path.dirname(input.clampbase[0]),
        fdr=A4_TRAIT_CFG["fdr"],
    resources:
        mem_mb=8000,
        runtime=30,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --clampfull-dir {params.clampfull_dir} "
        "--clampbase-dir {params.clampbase_dir} --fdr {params.fdr} "
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
