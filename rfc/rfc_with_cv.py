'''
Random Forest Classifier with cross-validation and feature importance logging
using W&B.

This module trains a Random Forest classifier on encoded PTM sequences,
logs model performance metrics, visualizes top features, and saves
a heatmap of feature importances.
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score, matthews_corrcoef
import wandb
from wandb.sklearn import plot_precision_recall, plot_feature_importances
from wandb.sklearn import plot_class_proportions, plot_learning_curve, plot_roc


def get_splits(dataset_size):
    '''
    Calculate the number of splits for cross-validation based on dataset size.

    Args:
        dataset_size (int): Number of samples in the dataset.

    Returns:
        int: Number of splits for RepeatedStratifiedKFold.
    '''
    return int(np.round(dataset_size / (0.2 * dataset_size)))


def main(X_dict, y_path, class_imbalance, hydro):
    '''
    Train Random Forest on encoded PTM features, perform cross-validation,
    log metrics to W&B, and visualize feature importance.

    Args:
        X_csv (str): Path to CSV file containing encoded features.
        y_path (str): Path to class labels file (multi-label or binary).
        class_imbalance (float): Imbalance factor for the dataset.
        hydro (bool): Whether hydrogen encoding is included.
    '''
    X = pd.read_csv(X_dict)
    feature_names = list(X.columns)
    feature_array = np.array(feature_names)

    y = pd.read_csv(y_path, delimiter='\t', header=None)
    y = y.astype('category')
    y = y.to_numpy().ravel()
    label_names = ['no_ptm', 'has_ptm']

    ptm = y_path.split("/")[-2]

    cv = RepeatedStratifiedKFold(n_splits=get_splits(X.shape[0]), n_repeats=10, random_state=42)

    for i, (train_index, test_index) in enumerate(cv.split(X, y)):
        print(i)
        if i != 0:
            break
        X_train, y_train = X.iloc[train_index,:], y[train_index]
        X_test, y_test = X.iloc[test_index,:], y[test_index]

        rfc = RandomForestClassifier(n_jobs=-1, n_estimators=100, random_state=42)
        rfc.fit(X_train, y_train)

        y_pred = rfc.predict(X_test)
        y_probas = rfc.predict_proba(X_test)
        importances = rfc.feature_importances_
        model_params = rfc.get_params()

        # wandb.init(
        #     project='bachelor-ptm4',
        #     name=ptm,
        #     config=model_params,
        #     reinit=True
        # )

        # wandb.config.update({'test_size' : 0.2,
        #                 'train_len' : len(X_train),
        #                 'test_len' : len(X_test),
        #                 'ptm': ptm,
        #                 'class_imbalance': class_imbalance,
        #                 'with_hydrogen': hydro,
        #                 'criterion' : 'gini'})

        # top_n = 10
        # sorted_idx = np.argsort(importances)[::-1]
        # top_features = [(feature_array[i], importances[i]) for i in sorted_idx[:top_n]]
        # table = wandb.Table(data=top_features, columns=['Feature', 'Importance'])

        # wandb.log({
        #     'class_proportions': plot_class_proportions(y_train, y_test, label_names),
        #     'learning_curve': plot_learning_curve(rfc, X_train, y_train),
        #     'roc_curve': plot_roc(y_test, y_probas, label_names),
        #     'precision_recall': plot_precision_recall(y_test, y_probas, label_names),
        #     'feature_importance_builtin': plot_feature_importances(rfc),
        #     'accuracy': accuracy_score(y_test, y_pred),
        #     'confusion_matrix': wandb.plot.confusion_matrix(
        #         y_true=y_test,
        #         preds=y_pred,
        #         class_names=label_names
        #     ),
        #     'top_feature_importance': wandb.plot.bar(
        #         table, 'Feature', 'Importance', title=f'Top {top_n} Feature Importances'
        #     ),
        #     'classification_report': classification_report(y_test, y_pred, target_names=label_names, output_dict=True),
        #     'MCC': matthews_corrcoef(y_test, y_pred)
        # })
        # wandb.finish()

    n_atoms = 10 if hydro == 'with_hydrogen' else 8
    n_positions = importances.shape[0] // n_atoms

    reshaped_importances = importances.reshape(n_atoms, n_positions, order='F')
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.matshow(reshaped_importances, cmap='Greys', aspect='auto')

    if n_atoms == 10:
        atom_labels = ['H', 'C', 'O', 'N', 'S', 'H', 'C', 'O', 'N', 'S']
    else:
        atom_labels = ['C', 'O', 'N', 'S', 'C', 'O', 'N', 'S']

    ax.set_yticks(range(len(atom_labels)))
    ax.set_yticklabels(atom_labels)
    ax.tick_params(axis='both', labelsize=14)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Feature Importance', rotation=270, labelpad=15, fontsize=14)

    # ax.set_title(f'Feature Importance Heatmap ({ptm})', pad=20)
    dir_name = 'data/feature_importances'
    plt.savefig(f'{dir_name}/fi_{ptm}.pdf')
    plt.clf()
    print('saved')
