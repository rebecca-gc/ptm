import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

encodings = 'data/iCAN_level_2_without_hydrogen.csv'
X = pd.read_csv(encodings, delimiter=',')
print(X.head())
print(X.shape)

classes = 'data/classes.txt'
y = pd.read_csv(classes, delimiter='\t', header=None)
y = y.astype('category')
y = y.to_numpy().ravel()
print(y.shape)

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=23, test_size=0.1)

rf = RandomForestClassifier(random_state=23)
rf.fit(X_train, y_train)
y_pred = rf.predict(X_test)

print(rf.score(X_test, y_test))
print(classification_report(y_test, y_pred))
features = pd.DataFrame(rf.feature_importances_, index=X.columns)
print(features.describe())