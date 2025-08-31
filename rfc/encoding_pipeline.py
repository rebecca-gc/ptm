'''
Module for encoding PTM and non-PTM sequences using iCAN and training
Random Forest classifiers with cross-validation. Supports parallel
processing with progress bars.
'''

import os
import sys
import time
import threading
import multiprocessing
from joblib import Parallel, delayed
from Bio import SeqIO

sys.path.append(os.path.abspath('data_preprocess'))
sys.path.append(os.path.abspath('Source'))

import rfc_with_cv
import data_pipeline
import ican


def draw_progress(progress_dict, total_steps):
    '''
    Display dynamic terminal progress bars for multiple PTM encoding steps.

    Args:
        progress_dict (multiprocessing.Manager.dict): Dictionary holding
            current progress for each key.
        total_steps (dict): Dictionary with total steps for each key.
    '''
    ESC = '\033'
    CLEAR = f'{ESC}[2J'
    HIDE_CURSOR = f'{ESC}[?25l'
    SHOW_CURSOR = f'{ESC}[?25h'

    sys.stdout.write(CLEAR + HIDE_CURSOR)
    sys.stdout.flush()

    try:
        while True:
            sys.stdout.write(f'{ESC}[H')
            for key, total in total_steps.items():
                current = progress_dict[key]
                percent = current / total
                bar_length = 40
                filled = int(bar_length * percent)
                progress_bar = '█' * filled + '-' * (bar_length - filled)
                sys.stdout.write(f'{key:30} |{progress_bar}| {percent*100:5.1f}%\n')
            sys.stdout.flush()
            time.sleep(0.1)
            if all(progress_dict[k] >= total_steps[k] for k in progress_dict):
                break
    finally:
        sys.stdout.write(SHOW_CURSOR)
        sys.stdout.flush()


def ican_parallel(seq_file, queue):
    '''
    Run iCAN encoding and Random Forest classification for a single
    multi-FASTA file, reporting progress to a queue.

    Args:
        seq_file (str): Path to the multi-FASTA file to encode.
        queue (multiprocessing.Queue): Queue to report progress increments.
    '''
    output_dir = os.path.dirname(seq_file)
    ptm = output_dir.split('/')[-1]
    csv_file = os.path.join(output_dir, f'iCAN_encoding_{ptm}')
    sys.argv = ['ican.py', f'--output_path={csv_file}', '--alphabet_mode=with_hydrogen', seq_file]

    ican.main(
        queue=queue,
        smiles_key=f'{ptm}/smiles',
        encode_key=f'{ptm}/encode',
    )

    y_path = seq_file.replace('seqs.fasta', 'classes.txt')
    # rfc_with_cv.main(csv_file, y_path, '1', 'with_hydrogen')


def run_parallel_with_bars(ptms_dir):
    '''
    Run iCAN encoding in parallel for all PTM directories, with live
    progress bars for each PTM's smiles and encoding steps.

    Args:
        ptms_dir (str): Path to the parent directory containing PTM subdirectories.
    '''
    def count_fasta_entries(file_path):
        return sum(1 for _ in SeqIO.parse(file_path, 'fasta'))

    steps_total = {}
    seqs = []

    for ptm_dir in os.listdir(ptms_dir):
        path = os.path.join(ptms_dir, ptm_dir)
        if os.path.isdir(path):
            fasta_file = os.path.join(ptms_dir, ptm_dir, 'seqs.fasta')
            seqs.append(fasta_file)
            x = count_fasta_entries(fasta_file)
            steps_total[f'{ptm_dir}/smiles'] = x
            steps_total[f'{ptm_dir}/encode'] = x

    manager = multiprocessing.Manager()
    progress_dict = manager.dict({key: 0 for key in steps_total})
    queue = manager.Queue()

    threading.Thread(target=draw_progress, args=(progress_dict, steps_total), daemon=True).start()

    # Queue-Listener updates the progress dictionary
    def progress_updater():
        while True:
            key, delta = queue.get()
            progress_dict[key] += delta
            if all(progress_dict[k] >= steps_total[k] for k in steps_total):
                break

    threading.Thread(target=progress_updater, daemon=True).start()

    Parallel(n_jobs=-1)(
        delayed(ican_parallel)(seq, queue) for seq in seqs
    )

    time.sleep(0.5)


def main():
    '''
    Main function to run the encoding pipeline.

    Steps:
        1. Optionally run the preprocessing pipeline (commented out).
        2. Run parallel iCAN encoding and Random Forest classification
           with live progress bars for all PTM multi-FASTA datasets.
    '''
    # data_pipeline.main()
    run_parallel_with_bars('data/ptms')


if __name__ == '__main__':
    main()
