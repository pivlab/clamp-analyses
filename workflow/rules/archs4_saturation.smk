A4_SAT_BIO = f"{A4_BIO}/01_saturation"
A4_SAT_NB = f"{A4_BIO_NB}/01_saturation"

# ============================================================
# Pathway coverage against model rank K, at several study-subsampling levels.
#
# Coverage asks how much data CLAMP needs; saturation asks how much model
# capacity, and whether the two interact.  Same by-study subsampling and the same
# GPU ORA scorer as archs4_coverage.smk.
# ============================================================

rule saturation_report_archs4:
    input:
        saturation_dir=A4_SAT_BIO,
        notebook=f"{A4_SAT_NB}/00_saturation.ipynb",
    output:
        saturation_long=f"{A4_SAT_BIO}/saturation_long.csv",
        panel_ready=f"{A4_SAT_BIO}/saturation_panel_ready.csv",
        complete=touch(f"{A4_SAT_BIO}/notebook.complete"),
    log:
        notebook=f"{A4_SAT_NB}/00_saturation.executed.ipynb",
    params:
        saturation_dir=A4_SAT_BIO,
        fdr=A4_CFG["ora"]["fdr"],
    conda: "clamp-analyses"
    notebook:
        f"{A4_SAT_NB}/00_saturation.ipynb"


rule archs4_saturation:
    input:
        rules.saturation_report_archs4.output.complete,
