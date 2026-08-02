Deferred: CLAMPfull with the combined Hallmark+Reactome+GO:CC+MSigDB C8 pathway prior on
GTEx (formerly `scripts/gtex/clamp_hall.R`, rule `clamp_hall_gtex`). Not part of the
current GTEx model-building pipeline in `scripts/gtex/`/`workflow/rules/gtex.smk` -- this
belongs to a later section of the paper. Parked here with its `common.R` dependency so it
can be re-wired into a rule/script when that section is worked on.

Depends on `clamp_gtex`'s outputs (gtex_genes.rds, gtex_samples.rds, gtex_fbm_filt.rds,
gtex_svdRes.rds, CLAMPbase.rds, CLAMP_K_gtex.rds) and the local pathway GMT files under
`data/pathways/`.
