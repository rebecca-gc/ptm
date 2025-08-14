import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score, matthews_corrcoef
import wandb
from wandb.sklearn import plot_precision_recall, plot_feature_importances
from wandb.sklearn import plot_class_proportions, plot_learning_curve, plot_roc


def get_splits(dataset_size):
    return int(np.round(dataset_size / (0.2 * dataset_size)))


def main(X_dict, y_path, class_imbalance, hydro):
    flatten_dict = {}
    for k in X_dict.keys():
        flatten_dict[k] = X_dict[k].to_numpy().flatten(order='F')
    X = pd.DataFrame.from_dict(flatten_dict, orient='index')
    X = X.reset_index(drop=True)
    feature_names = list(X.columns)
    feature_array = np.array(feature_names)

    y = pd.read_csv(y_path, delimiter='\t', header=None)
    y = y.astype('category')
    y = y.to_numpy().ravel()
    label_names = ['no_ptm', 'has_ptm']

    ptm = y_path.split("/")[-2]

    cv = RepeatedStratifiedKFold(n_splits=get_splits(X.shape[0]), n_repeats=10, random_state=42)

    for i, (train_index, test_index) in enumerate(cv.split(X, y)):
        X_train, y_train = X.iloc[train_index,:], y[train_index]
        X_test, y_test = X.iloc[test_index,:], y[test_index]
        rfc = RandomForestClassifier(n_jobs=-1, n_estimators=100, random_state=42)
        rfc.fit(X_train, y_train)

        y_pred = rfc.predict(X_test)
        y_probas = rfc.predict_proba(X_test)
        importances = rfc.feature_importances_
        model_params = rfc.get_params()

        wandb.init(
            project='bachelor-ptm4',
            name=ptm,
            config=model_params,
            reinit=True
        )

        wandb.config.update({'test_size' : 0.2,
                        'train_len' : len(X_train),
                        'test_len' : len(X_test),
                        'ptm': ptm,
                        'class_imbalance': class_imbalance,
                        'with_hydrogen': hydro})
        
        top_n = 10
        sorted_idx = np.argsort(importances)[::-1]
        top_features = [(feature_array[i], importances[i]) for i in sorted_idx[:top_n]]

        table = wandb.Table(data=top_features, columns=['Feature', 'Importance'])

        wandb.log({
            'class_proportions': plot_class_proportions(y_train, y_test, label_names),
            'learning_curve': plot_learning_curve(rfc, X_train, y_train),
            'roc_curve': plot_roc(y_test, y_probas, label_names),
            'precision_recall': plot_precision_recall(y_test, y_probas, label_names),
            'feature_importance_builtin': plot_feature_importances(rfc),
            'accuracy': accuracy_score(y_test, y_pred),
            'confusion_matrix': wandb.plot.confusion_matrix(
                y_true=y_test,
                preds=y_pred,
                class_names=label_names
            ),
            'top_feature_importance': wandb.plot.bar(
                table, 'Feature', 'Importance', title=f'Top {top_n} Feature Importances'
            ),
            'classification_report': classification_report(y_test, y_pred, target_names=label_names, output_dict=True),
            'MCC': matthews_corrcoef(y_test, y_pred)
        })
        
        wandb.finish()

    reshaped_importances = importances.reshape(X_dict[0].shape, order='F')
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.matshow(reshaped_importances, cmap='Greys', aspect='auto')

    if reshaped_importances.shape[0] == 8:
        atom_labels = ['C', 'O', 'N', 'S', 'C', 'O', 'N', 'S']
    elif reshaped_importances.shape[0] == 10:
        atom_labels = ['H', 'C', 'O', 'N', 'S', 'H', 'C', 'O', 'N', 'S']

    ax.set_yticks(range(len(atom_labels)))
    ax.set_yticklabels(atom_labels)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Feature Importance', rotation=270, labelpad=15)

    ax.set_title(f'Feature Importance Heatmap ({ptm})', pad=20)

    #plt.tight_layout()
    dir_name = 'data/feature_importances'
    plt.savefig(f'{dir_name}/fi_{ptm}.jpg')
    plt.clf()


if __name__ == '__main__':
    main()
