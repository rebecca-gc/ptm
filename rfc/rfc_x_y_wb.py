import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
import wandb
from wandb.sklearn import plot_precision_recall, plot_feature_importances
from wandb.sklearn import plot_class_proportions, plot_learning_curve, plot_roc


def main(X_dict, y_path):
    flatten_dict = {}
    for k in X_dict.keys():
        flatten_dict[k] = X_dict[k].to_numpy().flatten(order='F')
    X = pd.DataFrame.from_dict(flatten_dict, orient='index')
    X = X.reset_index(drop=True)
    print(X.shape)
    feature_names = list(X.columns)
    feature_array = np.array(feature_names)

    y = pd.read_csv(y_path, delimiter='\t', header=None)
    y = y.astype('category')
    y = y.to_numpy().ravel()
    print(y.shape)
    label_names = ["no_ptm", "has_ptm"]

    test_size = 0.2
    n_estimators = 100

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=23, test_size=test_size)

    rf = RandomForestClassifier(n_jobs=-1, n_estimators=n_estimators, random_state=23)
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)
    y_probas = rf.predict_proba(X_test)
    importances = rf.feature_importances_
    model_params = rf.get_params()

    if 'db_seqs_classes' in y_path:
        db_name = y_path.split("_")[-2].split("/")[-1]
    else:
        db_name = f'{y_path.split("/")[-2]}_all-all'

    wandb.init(
        project="bachelor-ptm2",
        name=f"rf_{db_name}_testsize{test_size}_estimators{n_estimators}",
        config=model_params
    )

    wandb.config.update({"test_size" : test_size,
                    "train_len" : len(X_train),
                    "test_len" : len(X_test),
                    "database": db_name.split('-')[1],
                    "ptm": db_name.split('-')[0]})
    
    top_n = 10
    sorted_idx = np.argsort(importances)[::-1]
    top_features = [(feature_array[i], importances[i]) for i in sorted_idx[:top_n]]

    table = wandb.Table(data=top_features, columns=["Feature", "Importance"])

    wandb.log({
        "class_proportions": plot_class_proportions(y_train, y_test, label_names),
        "learning_curve": plot_learning_curve(rf, X_train, y_train),
        "roc_curve": plot_roc(y_test, y_probas, label_names),
        "precision_recall": plot_precision_recall(y_test, y_probas, label_names),
        "feature_importance_builtin": plot_feature_importances(rf),
        "accuracy": accuracy_score(y_test, y_pred),
        "confusion_matrix": wandb.plot.confusion_matrix(
            y_true=y_test,
            preds=y_pred,
            class_names=label_names
        ),
        "top_feature_importance": wandb.plot.bar(
            table, "Feature", "Importance", title=f"Top {top_n} Feature Importances"
        ),
        "classification_report": classification_report(y_test, y_pred, target_names=label_names, output_dict=True)
    })
    
    wandb.finish()

    reshaped_importances = importances.reshape(X_dict[0].shape, order='F')
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.matshow(reshaped_importances, cmap='viridis_r', aspect='auto')

    if reshaped_importances.shape[0] == 8:
        atom_labels = ["C", "O", "N", "S", "C", "O", "N", "S"]
    elif reshaped_importances.shape[0] == 10:
        atom_labels = ["H", "C", "O", "N", "S", "H", "C", "O", "N", "S"]

    ax.set_yticks(range(len(atom_labels)))
    ax.set_yticklabels(atom_labels)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Feature Importance", rotation=270, labelpad=15)

    ax.set_title(f"Feature Importance Heatmap ({db_name})", pad=20)

    plt.tight_layout()
    dir_name = 'data/feature_importances'
    try:
        os.makedirs(dir_name)
    except FileExistsError:
        print(f"Directory '{dir_name}' already exists.")
    plt.savefig(f'{dir_name}/fi_{db_name}.jpg')
    plt.clf()


if __name__ == '__main__':
    main()
