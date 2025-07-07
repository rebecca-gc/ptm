import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, f1_score


# dirs = ['glycosylation', 's_nitrosylation', 'acetylation', 'methylation']

# results = 'results.txt'
# if os.path.exists(results):
#     os.remove(results)

# for dir in dirs:
#     encodings = os.path.join('../data', dir, 'iCAN_level_2_without_hydrogen.csv')
#     X = pd.read_csv(encodings, delimiter=',')
#     print(X.head())
#     print(X.shape)

#     classes = os.path.join('../data', dir, 'classes.txt')
#     y = pd.read_csv(classes, delimiter='\t', header=None)
#     y = y.astype('category')
#     y = y.to_numpy().ravel()
#     print(y.shape)

#     X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=23, test_size=0.1)

#     rf = RandomForestClassifier(random_state=23)
#     rf.fit(X_train, y_train)
#     y_pred = rf.predict(X_test)

#     with open(results, 'a') as f:
#         f.write(dir)
#         f.write(rf.score(X_test, y_test))
#         f.write(classification_report(y_test, y_pred))
#         f.write('\n-----\n-----\n')
#         print(rf.score(X_test, y_test))
#         print(classification_report(y_test, y_pred))
#         features = pd.DataFrame(rf.feature_importances_, index=X.columns)
#         # print(features.describe())

dirs = ['data/s_nitrosylation', 'data/methylation']

mean_accuracy_df = pd.DataFrame(np.nan, index=dirs, columns=range(50))
f1_score_df = pd.DataFrame(np.nan, index=dirs, columns=range(50))
proba_class1_df = pd.DataFrame(np.nan, index=dirs, columns=range(50))

def get_splits(dataset_size):
    return int(np.round(dataset_size / (0.2 * dataset_size)))

for dir in dirs:
    print(f'rfc for {dir}')
    encodings = os.path.join(dir, 'iCAN_level_2_without_hydrogen.csv')
    X = pd.read_csv(encodings, delimiter=',')

    classes = os.path.join(dir, 'classes.txt')
    y = pd.read_csv(classes, delimiter='\t', header=None)
    y = y.astype('category')
    y = y.to_numpy().ravel()

    cv = RepeatedStratifiedKFold(n_splits=get_splits(X.shape[0]), n_repeats=10, random_state=42)
            
    for i, (train_index, test_index) in enumerate(cv.split(X, y)):
        X_train, y_train = X.iloc[train_index,:], y[train_index]
        X_test, y_test = X.iloc[test_index,:], y[test_index]
        rfc = RandomForestClassifier(n_jobs=-1, n_estimators=100, random_state=42)
        rfc.fit(X_train, y_train)
        mean_accuracy = rfc.score(X_test, y_test)
        y_pred = rfc.predict(X_test)
        f1 = f1_score(y_test, y_pred)

        mean_accuracy_df.iat[dir, i] = mean_accuracy
        f1_score_df.iat[dir, i] = f1

    results_path = os.path.join('..', 'Results', 'csv')
    if os.path.exists(results_path) == False:
        os.mkdir(results_path)

    mean_accuracy_path = os.path.join(results_path, 'mean_accuracy_level_' + '.csv')
    f1_score_path = os.path.join(results_path, 'f1_score_level_' + '.csv')

    mean_accuracy_df.to_csv(mean_accuracy_path, index=True, header=True)
    f1_score_df.to_csv(f1_score_path, index=True, header=True)
