NB_ROOT   = "nbs/01_model_building/05_pseudobulk"
PROJ_ROOT = "nbs/03_model_biology/00_archs4/04_pseudobulk"
MODELS    = config["paths"]["models"]
BENCH_DIR = config["paths"]["benchmark"]
PROJ_DIR  = config["paths"]["projection"]

PROJ_NB_PREFIX = {
    "Brain_Mathys2023": "00",
    "Brain_Xiong2023":  "02",
    "PBMC_1k1k":        "04",
    "PBMC_Perez2022":   "06",
    "BRCA_Bassez2021":  "08",
    "CRC_Pelka2021":    "10",
}


rule preprocess:
    input:
        "data/pseudobulk/{dataset}/bulk_expr.csv"
    output:
        norm = f"{MODELS}/{{dataset}}/norm.csv",
        k    = f"{MODELS}/{{dataset}}/k.csv",
        nb   = f"{MODELS}/{{dataset}}/notebooks/00_preprocess.ipynb"
    params:
        mean_cutoff = config["preprocess"]["mean_cutoff"],
        var_cutoff  = config["preprocess"]["var_cutoff"]
    conda: "clamp-analyses"
    shell:
        """
        papermill {NB_ROOT}/00_preprocess.ipynb {output.nb} \
          -p DATASET {wildcards.dataset} \
          -p MEAN_CUTOFF {params.mean_cutoff} \
          -p VAR_CUTOFF {params.var_cutoff} \
          -k ir
        """


rule clamp:
    input:
        norm = rules.preprocess.output.norm,
        k    = rules.preprocess.output.k
    output:
        base_B   = f"{MODELS}/{{dataset}}/CLAMPbase/B.csv",
        full_B   = f"{MODELS}/{{dataset}}/CLAMPfull/B.csv",
        base_rds = f"{MODELS}/{{dataset}}/CLAMPbase/CLAMPbase.rds",
        full_rds = f"{MODELS}/{{dataset}}/CLAMPfull/CLAMPfull.rds",
        nb       = f"{MODELS}/{{dataset}}/notebooks/01_CLAMP.ipynb"
    params:
        max_iter = config["clamp"]["max_iter"]
    conda: "clamp-analyses"
    shell:
        """
        papermill {NB_ROOT}/01_CLAMP.ipynb {output.nb} \
          -p DATASET {wildcards.dataset} \
          -p MAX_ITER {params.max_iter} \
          -k ir
        """


rule plier:
    input:
        norm = rules.preprocess.output.norm,
        k    = rules.preprocess.output.k
    output:
        B   = f"{MODELS}/{{dataset}}/PLIER/B.csv",
        rds = f"{MODELS}/{{dataset}}/PLIER/PLIER.rds",
        nb  = f"{MODELS}/{{dataset}}/notebooks/02_PLIER.ipynb"
    conda: "clamp-analyses"
    shell:
        """
        papermill {NB_ROOT}/02_PLIER.ipynb {output.nb} \
          -p DATASET {wildcards.dataset} \
          -k ir
        """


rule pca_nmf_ica:
    input:
        norm = rules.preprocess.output.norm,
        k    = rules.preprocess.output.k
    output:
        pca_B   = f"{MODELS}/{{dataset}}/PCA/B.csv",
        nmf_B   = f"{MODELS}/{{dataset}}/NMF/B.csv",
        ica_B   = f"{MODELS}/{{dataset}}/ICA/B.csv",
        pca_pkl = f"{MODELS}/{{dataset}}/PCA/pca_model.pkl",
        nmf_pkl = f"{MODELS}/{{dataset}}/NMF/nmf_model.pkl",
        ica_pkl = f"{MODELS}/{{dataset}}/ICA/ica_model.pkl",
        nb      = f"{MODELS}/{{dataset}}/notebooks/03_PCA_NMF_ICA.ipynb"
    conda: "clamp-analyses"
    shell:
        """
        papermill {NB_ROOT}/03_PCA_NMF_ICA.ipynb {output.nb} \
          -p DATASET {wildcards.dataset} \
          -k clamp-analyses
        """


rule flashier:
    input:
        norm = rules.preprocess.output.norm,
        k    = rules.preprocess.output.k
    output:
        B   = f"{MODELS}/{{dataset}}/flashier/B.csv",
        rds = f"{MODELS}/{{dataset}}/flashier/flashier_model.rds",
        nb  = f"{MODELS}/{{dataset}}/notebooks/04_flashier.ipynb"
    params:
        backfit_maxiter = config["flashier"]["backfit_maxiter"]
    conda: "clamp-analyses"
    shell:
        """
        papermill {NB_ROOT}/04_flashier.ipynb {output.nb} \
          -p DATASET {wildcards.dataset} \
          -p BACKFIT_MAXITER {params.backfit_maxiter} \
          -k ir
        """


rule mofa_flex_prior:
    input:
        norm = rules.preprocess.output.norm,
        k    = rules.preprocess.output.k
    output:
        B   = f"{MODELS}/{{dataset}}/MOFA_FLEX_PRIOR/B_matrix.csv",
        pkl = f"{MODELS}/{{dataset}}/MOFA_FLEX_PRIOR/model.pkl",
        nb  = f"{MODELS}/{{dataset}}/notebooks/05_MOFA_FLEX_PRIOR.ipynb"
    params:
        max_epochs = config["mofa_flex"]["max_epochs"],
        seed       = config["mofa_flex"]["seed"]
    conda: "clamp-analyses"
    shell:
        """
        papermill {NB_ROOT}/05_MOFA_FLEX_PRIOR.ipynb {output.nb} \
          -p DATASET {wildcards.dataset} \
          -p MAX_EPOCHS {params.max_epochs} \
          -p SEED {params.seed} \
          -k clamp-analyses
        """



rule cogaps:
    input:
        norm = rules.preprocess.output.norm,
        k    = rules.preprocess.output.k
    output:
        B   = f"{MODELS}/{{dataset}}/CoGAPS/B.csv",
        rds = f"{MODELS}/{{dataset}}/CoGAPS/cogaps_model.rds",
        nb  = f"{MODELS}/{{dataset}}/notebooks/07_CoGAPS.ipynb"
    params:
        n_iterations = config["cogaps"]["n_iterations"]
    conda: "clamp-analyses"
    shell:
        """
        papermill {NB_ROOT}/07_CoGAPS.ipynb {output.nb} \
          -p DATASET {wildcards.dataset} \
          -p N_ITERATIONS {params.n_iterations} \
          -k ir
        """


rule gss:
    input:
        norm = rules.preprocess.output.norm,
        k    = rules.preprocess.output.k
    output:
        B   = f"{MODELS}/{{dataset}}/GSS/B.csv",
        rds = f"{MODELS}/{{dataset}}/GSS/RAVmodel.rds",
        nb  = f"{MODELS}/{{dataset}}/notebooks/08_GSS.ipynb"
    conda: "clamp-analyses"
    shell:
        """
        papermill {NB_ROOT}/08_GSS.ipynb {output.nb} \
          -p DATASET {wildcards.dataset} \
          -k ir
        """


rule project_archs4:
    input:
        norm  = rules.preprocess.output.norm,
        k     = rules.preprocess.output.k,
        model = config["paths"]["archs4_model"]
    output:
        proj = f"{PROJ_DIR}/{{dataset}}/projection.csv",
        nb   = f"{PROJ_DIR}/{{dataset}}/notebooks/projection.ipynb"
    params:
        prefix = lambda wc: PROJ_NB_PREFIX[wc.dataset]
    conda: "clamp-analyses"
    shell:
        """
        papermill {PROJ_ROOT}/{params.prefix}_{wildcards.dataset}_projection.ipynb {output.nb} \
          -p DATASET {wildcards.dataset} \
          -k ir
        """


rule check_outputs:
    input:
        expand(f"{MODELS}/{{ds}}/CLAMPbase/B.csv",              ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/CLAMPfull/B.csv",              ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/PLIER/B.csv",                  ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/PCA/B.csv",                    ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/NMF/B.csv",                    ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/ICA/B.csv",                    ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/flashier/B.csv",               ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/MOFA_FLEX_PRIOR/B_matrix.csv", ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/GSS/B.csv",                    ds=config["datasets"]),
    output:
        nb = f"{BENCH_DIR}/check_outputs.ipynb"
    conda: "clamp-analyses"
    shell:
        """
        papermill {NB_ROOT}/09_check_outputs.ipynb {output.nb} -k ir
        """


rule benchmark:
    input:
        expand(f"{MODELS}/{{ds}}/CLAMPbase/B.csv",              ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/CLAMPfull/B.csv",              ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/PLIER/B.csv",                  ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/PCA/B.csv",                    ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/NMF/B.csv",                    ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/ICA/B.csv",                    ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/flashier/B.csv",               ds=config["datasets"]),
        expand(f"{MODELS}/{{ds}}/MOFA_FLEX_PRIOR/B_matrix.csv", ds=config["datasets"]),
        # expand(f"{MODELS}/{{ds}}/CoGAPS/B.csv",                 ds=config["datasets"]),  # disabled: slow
        expand(f"{MODELS}/{{ds}}/GSS/B.csv",                    ds=config["datasets"]),
        f"{BENCH_DIR}/check_outputs.ipynb",
    output:
        long        = f"{BENCH_DIR}/benchmark_long.csv",
        summary     = f"{BENCH_DIR}/benchmark_summary.csv",
        assignments = f"{BENCH_DIR}/clampfull_lv_assignments.csv",
        nb          = f"{BENCH_DIR}/notebooks/00_method_benchmark.ipynb"
    conda: "clamp-analyses"
    shell:
        """
        papermill nbs/03_model_biology/02_pseudobulk/00_method_benchmark.ipynb \
          {output.nb} -k ir
        """


rule biology:
    input:
        full_rds    = expand(f"{MODELS}/{{ds}}/CLAMPfull/CLAMPfull.rds", ds=config["datasets"]),
        assignments = f"{BENCH_DIR}/clampfull_lv_assignments.csv",
    output:
        top_lvs  = f"{BENCH_DIR}/01_biology/top_lvs_per_celltype.csv",
        pathways = f"{BENCH_DIR}/01_biology/lv_pathways.csv",
        cm_ora   = f"{BENCH_DIR}/01_biology/lv_cellmarker_ora.csv",
        loadings = f"{BENCH_DIR}/01_biology/gene_loadings.csv",
        nb       = f"{BENCH_DIR}/notebooks/01_biology.ipynb"
    conda: "clamp-analyses"
    shell:
        """
        papermill nbs/03_model_biology/02_pseudobulk/01_biology.ipynb \
          {output.nb} -k ir
        """
