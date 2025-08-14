import os
import sys
import numpy as np
import time
import random
import threading
import multiprocessing
from joblib import Parallel, delayed
from Bio import SeqIO
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from collections import Counter

sys.path.append(os.path.abspath('data_preprocess'))
sys.path.append(os.path.abspath('Source'))

import disease
import ican


def draw_progress(progress_dict, total_steps):
    ESC = '\033'
    CLEAR = f'{ESC}[2J'
    HIDE_CURSOR = f'{ESC}[?25l'
    SHOW_CURSOR = f'{ESC}[?25h'

    sys.stdout.write(CLEAR + HIDE_CURSOR)
    sys.stdout.flush()

    try:
        while True:
            sys.stdout.write(f'{ESC}[H')
            for key in sorted(progress_dict.keys()):
                current = progress_dict[key]
                total = total_steps[key]
                percent = current / total
                bar_length = 40
                filled = int(bar_length * percent)
                bar = '█' * filled + '-' * (bar_length - filled)
                sys.stdout.write(f'{key:30} |{bar}| {percent*100:5.1f}%\n')
            sys.stdout.flush()
            time.sleep(0.1)
            if all(progress_dict[k] >= total_steps[k] for k in progress_dict):
                break
    finally:
        sys.stdout.write(SHOW_CURSOR)
        sys.stdout.flush()


def ican_parallel(seq_file, queue):
    output = os.path.dirname(seq_file)
    sys.argv = ['ican.py', '--output_path=data', '--alphabet_mode=with_hydrogen', seq_file]

    X = ican.main(
        queue=queue,
        smiles_key=f'smiles',
        encode_key=f'encode',
    )

    y = np.loadtxt('data/multi_label_test.txt', dtype=int)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=23)

    rfc = RandomForestClassifier(n_estimators=100, random_state=23)

    rfc.fit(X_train, y_train)

    y_pred = rfc.predict(X_test)

    print("Subset accuracy:", accuracy_score(y_test, y_pred))


def run_parallel_with_bars(seq_file):
    def count_fasta_entries(file_path):
        return sum(1 for _ in SeqIO.parse(file_path, 'fasta'))

    steps_total = {}

    x = count_fasta_entries(seq_file)
    steps_total['smiles'] = x
    steps_total['encode'] = x

    manager = multiprocessing.Manager()
    progress_dict = manager.dict({key: 0 for key in steps_total})
    queue = manager.Queue()

    threading.Thread(target=draw_progress, args=(progress_dict, steps_total), daemon=True).start()

    # Queue-Listener
    def progress_updater():
        while True:
            key, delta = queue.get()
            progress_dict[key] += delta
            if all(progress_dict[k] >= steps_total[k] for k in steps_total):
                break

    threading.Thread(target=progress_updater, daemon=True).start()

    ican_parallel(seq_file, queue)

    time.sleep(0.5)

def filter(seq_file):
    seqs = list(SeqIO.parse(seq_file, 'fasta'))
    lens = [len(record.seq) for record in seqs]
    x = np.percentile(lens,95)
    filtered_indices = [i for i, len in enumerate(lens) if len <= x]
    y = np.loadtxt('data/multi_label.txt', dtype=int)

    # Keep labels corresponding to filtered indices
    y_filtered = y[filtered_indices]

    np.savetxt('data/multi_label_filtered.txt', y_filtered, fmt='%d')

    combos = [tuple(row) for row in y_filtered]
    counts = Counter(combos)

    max_samples_per_combo = 25
    combo_counts = {combo: 0 for combo in counts.keys()}
    indices_to_keep = []

    random.seed(42)
    shuffled_indices = list(range(len(y_filtered)))
    random.shuffle(shuffled_indices)

    for idx in shuffled_indices:
        combo = combos[idx]
        if combo_counts[combo] < max_samples_per_combo:
            indices_to_keep.append(idx)
            combo_counts[combo] += 1
    
    filtered_sequences = [seqs[filtered_indices[i]] for i in indices_to_keep]
    SeqIO.write(filtered_sequences, 'data/sequences_filtered.fasta', 'fasta')

    y_final = y_filtered[indices_to_keep]
    np.savetxt('data/multi_label_final.txt', y_final, fmt='%d')


def main():
    #disease.file_generator()
    
    seq_file = 'data/disease_seqs.fasta'

    #filter(seq_file)

    # y = np.loadtxt('data/multi_label_final.txt', dtype=float)
    # combos = [tuple(row) for row in y]
    # counts = Counter(combos)

    # # Sorted by frequency
    # for combo, count in counts.most_common():
    #     print(combo, count)

    seq_file = 'data/sequences_filtered_test.fasta'

    run_parallel_with_bars(seq_file)


if __name__ == '__main__':
    main()