import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

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

def main(output):
    results = os.path.join(output, 'results.txt')
    if os.path.exists(results):
        os.remove(results)

    print(f'rfc for {output}')
    encodings = os.path.join(output, 'iCAN_level_2_without_hydrogen.csv')
    X = pd.read_csv(encodings, delimiter=',')
    #print(X.head())
    print(X.shape)

    classes = os.path.join(output, 'classes.txt')
    y = pd.read_csv(classes, delimiter='\t', header=None)
    y = y.astype('category')
    y = y.to_numpy().ravel()
    print(y.shape)

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=23, test_size=0.1)

    rf = RandomForestClassifier(n_jobs=-1, n_estimators=100, random_state=23)
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)

    with open(results, 'a') as f:
        f.write(f'Mean accuracy: {rf.score(X_test, y_test)}\n\n')
        f.write(classification_report(y_test, y_pred))
        # print(rf.score(X_test, y_test))
        # print(classification_report(y_test, y_pred))
        f.write('\nTop 10 features with the most influence:\n')
        features = pd.DataFrame(rf.feature_importances_, index=X.columns)
        top10 = features.sort_values(by=features.columns[0], ascending=False).head(10)
        f.write(top10.to_string())
        # print(features.describe())


if __name__ == '__main__':
    main()
