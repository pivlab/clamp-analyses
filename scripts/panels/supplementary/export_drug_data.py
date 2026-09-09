"""Export the main-repo drug-disease scores and ROC coordinates for Supplement 6."""
from pathlib import Path
import pandas as pd
from sklearn.metrics import roc_curve, roc_auc_score
ROOT = Path(__file__).resolve().parents[3]
source = ROOT / 'output/03_model_biology/02_archs4/02_drug_diseases_canonical/aggregate/max'
out = ROOT / 'output/99_panels/supp6/source_data'
out.mkdir(parents=True, exist_ok=True)
scores = pd.read_pickle(source / 'predictions_results_aggregated.pkl')
names = {'module_based_archs4':'ARCHS4','module_based_recount2':'recount2','module_based_gtex':'GTEx','gene_based':'Gene-based'}
scores['model'] = scores.method.map(names)
assert set(scores.model) == set(names.values())
scores.to_csv(out / 'prediction_scores.csv', index=False)
curves = []
summary = []
for model in names.values():
    sub = scores[scores.model == model]
    fpr, tpr, thresholds = roc_curve(sub.true_class, sub.score)
    curves.append(pd.DataFrame({'model':model, 'fpr':fpr, 'tpr':tpr, 'threshold':thresholds}))
    summary.append({'model':model, 'AUROC':roc_auc_score(sub.true_class, sub.score), 'n_pairs':len(sub), 'n_positive':int(sub.true_class.sum())})
pd.concat(curves).to_csv(out / 'roc_coordinates.csv', index=False)
pd.DataFrame(summary).to_csv(out / 'performance_summary.csv', index=False)
print(f'Exported {len(scores)} scores for four methods.')
