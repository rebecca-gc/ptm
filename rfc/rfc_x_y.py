import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report


def main(X, y_path):
    results = 'rfc/database_results.txt'

    #print(f'rfc for {output}')
    #encodings = os.path.join(output, 'iCAN_level_2_without_hydrogen.csv') # !!!!
    #X = pd.read_csv(encodings, delimiter=',')
    #print(X.head())
    print(X.shape)

    #classes = os.path.join(output, 'classes.txt') # !!!!
    y = pd.read_csv(y_path, delimiter='\t', header=None)
    y = y.astype('category')
    y = y.to_numpy().ravel()
    print(y.shape)

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=23, test_size=0.1)

    rf = RandomForestClassifier(n_jobs=-1, n_estimators=100, random_state=23)
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)

    with open(results, 'a') as f:
        f.write(f'{y_path.split("_")[-2].split("/")[-1]}\nMean accuracy: {rf.score(X_test, y_test)}\n\n')
        f.write(classification_report(y_test, y_pred))
        # print(rf.score(X_test, y_test))
        # print(classification_report(y_test, y_pred))
        f.write('\nTop 10 features with the most influence:\n')
        features = pd.DataFrame(rf.feature_importances_, index=X.columns)
        top10 = features.sort_values(by=features.columns[0], ascending=False).head(10)
        f.write(top10.to_string())
        f.write('\n')
        # print(features.describe())


if __name__ == '__main__':
    main()
