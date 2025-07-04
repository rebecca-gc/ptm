import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from Bio import SeqIO


def aa_freq(seq):
    seq = seq.upper()
    return [seq.count(aa) / len(seq) for aa in aa_list]

dirs = ['test', 'glycosylation', 's_nitrosylation', 'acetylation', 'methylation']

results = 'results_freq_encoding.txt'
if os.path.exists(results):
    os.remove(results)

for dir in dirs:
    seqs_path = os.path.join('../data', dir, 'seqs.fasta')
    seqs = [str(record.seq) for record in SeqIO.parse(seqs_path, 'fasta')]

    aa_list = list("ACDEFGHIKLMNPQRSTVWY")  # 20 Standard-Aminosäuren

    X = np.array([aa_freq(seq) for seq in seqs])

    classes = os.path.join('../data', dir, 'classes.txt')
    y = pd.read_csv(classes, delimiter='\t', header=None)
    y = y.astype('category')
    y = y.to_numpy().ravel()

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=23, test_size=0.2)

    rf = RandomForestClassifier(random_state=23)
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)

    with open(results, 'a') as f:
        f.write(f'{dir}\nMean accuracy: {rf.score(X_test, y_test)}\n')
        f.write(classification_report(y_test, y_pred))
        f.write('\nTop 10 features with the most influence:\n')
        X = pd.DataFrame(X, columns=aa_list) 
        features = pd.DataFrame(rf.feature_importances_, index=X.columns)
        top10 = features.sort_values(by=features.columns[0], ascending=False).head(10)
        f.write(top10.to_string())
        f.write('\n\n-----\n-----\n\n')
        