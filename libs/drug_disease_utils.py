"""
Utility functions for drug-disease association predictions using CLAMP.

Adapted from the PhenoPlier project (greenelab/phenoplier), replacing
MultiplierProjection / entity.Trait dependencies with direct pandas operations
and an explicit EFO → DOID mapping step.

Original: https://github.com/greenelab/phenoplier/blob/main/libs/drug_disease.py
"""

import pandas as pd


def map_traits_to_doid(data, preferred_doids, ukb_efo, efo_xrefs, do_xrefs):
    """
    Maps trait columns (UKB full codes) to Disease Ontology IDs (DOID).
    For traits mapping to multiple DOIDs, prefers those in `preferred_doids`.
    When a DOID appears from multiple traits, keeps the maximum score.
    """
    doid_efo = efo_xrefs[efo_xrefs['target_id_type'] == 'DOID']

    trait_to_doid = {}
    for trait in data.columns:
        if trait not in ukb_efo.index:
            continue
        rows = ukb_efo.loc[trait]
        if isinstance(rows, pd.Series):
            rows = rows.to_frame().T

        all_efo_codes = set()
        for term_codes in rows['term_codes'].dropna():
            for code in str(term_codes).split(','):
                code = code.strip()
                if code:
                    all_efo_codes.add(code)

        all_doids = set()
        for efo_code in all_efo_codes:
            mask = doid_efo['term_id'] == efo_code
            all_doids.update(doid_efo[mask]['target_id'].values)
            if efo_code.startswith('EFO:'):
                efo_num = efo_code[4:]
                mask2 = (do_xrefs['resource'] == 'EFO') & (do_xrefs['resource_id'] == efo_num)
                all_doids.update(do_xrefs[mask2]['doid_code'].values)

        if not all_doids:
            continue

        preferred = sorted(all_doids & preferred_doids)
        trait_to_doid[trait] = preferred[0] if preferred else sorted(all_doids)[0]

    data_mapped = data.loc[:, list(trait_to_doid.keys())].rename(columns=trait_to_doid)
    data_mapped = data_mapped.T.groupby(level=0).max().T
    return data_mapped


def _zero_nontop_genes(trait_vector, n_top, use_abs=True):
    """Zeros all but the top `n_top` values (genes or LVs) in a Series."""
    values = trait_vector.abs() if use_abs else trait_vector
    top_idx = values.sort_values(ascending=False).head(n_top).index
    result = trait_vector.copy()
    result[~result.index.isin(top_idx)] = 0.0
    return result


def predict_dotprod_neg(
    drug_gene_data,
    gene_trait_data_filename,
    gene_trait_data,
    output_dir,
    base_method_name,
    preferred_doid_list,
    ukb_efo,
    efo_xrefs,
    do_xrefs,
    n_top_conditions=None,
    use_abs=True,
):
    """
    Computes drug-disease predictions as: score = -1 * drug^T * disease

    Saves an HDF5 file with keys:
    - full_prediction: all traits
    - prediction: DOID-mapped traits (for gold-standard comparison)
    - metadata

    File naming: {stem}-{all_genes|top_N_genes}-prediction_scores.h5
    """
    output_subdir = output_dir / 'dotprod_neg'
    output_subdir.mkdir(exist_ok=True, parents=True)

    suffix = 'all_genes' if n_top_conditions is None else f'top_{n_top_conditions}_genes'
    stem = gene_trait_data_filename.stem
    output_file = output_subdir / f'{stem}-{suffix}-prediction_scores.h5'

    print(f'  predicting {suffix}...')

    disease_data = gene_trait_data.copy()
    if n_top_conditions is not None:
        disease_data = disease_data.apply(
            lambda x: _zero_nontop_genes(x, n_top_conditions, use_abs)
        )

    scores = -1.0 * drug_gene_data.T.dot(disease_data)
    print(f'    shape: {scores.shape}')

    with pd.HDFStore(output_file, mode='w', complevel=4) as store:
        scores.index.name = 'drug'
        scores.columns.name = 'trait'
        full_pred = (
            scores.unstack()
            .reset_index()
            .rename(columns={0: 'score'})
        )
        full_pred['trait'] = full_pred['trait'].astype('category')
        full_pred['drug'] = full_pred['drug'].astype('category')
        assert full_pred.shape == full_pred.dropna().shape
        store.put('full_prediction', full_pred, format='table')

        scores_doid = map_traits_to_doid(
            scores, preferred_doid_list, ukb_efo, efo_xrefs, do_xrefs
        )
        assert scores_doid.index.is_unique
        assert scores_doid.columns.is_unique

        scores_doid.index.name = 'drug'
        scores_doid.columns.name = 'trait'
        doid_pred = (
            scores_doid.unstack()
            .reset_index()
            .rename(columns={0: 'score'})
        )
        doid_pred['trait'] = doid_pred['trait'].astype('category')
        doid_pred['drug'] = doid_pred['drug'].astype('category')
        assert doid_pred.shape == doid_pred.dropna().shape
        store.put('prediction', doid_pred, format='table')

        meta = pd.DataFrame({
            'method': [base_method_name],
            'n_top_genes': [-1.0 if n_top_conditions is None else float(n_top_conditions)],
            'data': [stem],
        })
        store.put('metadata', meta, format='table')

    print(f'    saved to: {output_file}')
