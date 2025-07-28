import os
import sys
import time
import threading
import multiprocessing
from joblib import Parallel, delayed
from Bio import SeqIO

sys.path.append(os.path.abspath('data_preprocess'))
sys.path.append(os.path.abspath('Source'))

import data_pipeline
import rfc_x_y_wb
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
    seq = seq_file.split('/')[-1]
    sys.argv = ['ican.py', f'--output_path={output}', '--alphabet_mode=with_hydrogen', seq_file]

    X = ican.main(
        queue=queue,
        smiles_key=f'{seq}/smiles',
        encode_key=f'{seq}/encode',
    )

    y = seq_file.replace('seqs.fasta', 'classes.txt')
    rfc_x_y_wb.main(X, y)

def run_parallel_with_bars(ptms_dir):
    def count_fasta_entries(file_path):
        return sum(1 for _ in SeqIO.parse(file_path, 'fasta'))
        
    steps_total = {}
    seqs = []

    for ptm_dir in os.listdir(ptms_dir):
        path = os.path.join(ptms_dir, ptm_dir, 'db_seqs_classes')
        if os.path.isdir(path):
            #fasta_file = os.path.join(ptms_dir, ptm_dir, 'seqs.fasta')
            #seqs.append(fasta_file)
            #x = count_fasta_entries(fasta_file)
            #steps_total['seqs.fasta/smiles'] = x
            #steps_total['seqs.fasta/encode'] = x
            for seq in os.listdir(path):
                if seq.endswith('.fasta'):
                    fasta_file = os.path.join(path, seq)
                    seqs.append(fasta_file)
                    x = count_fasta_entries(fasta_file)
                    steps_total[f'{seq}/smiles'] = x
                    steps_total[f'{seq}/encode'] = x

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

    Parallel(n_jobs=-1)(
        delayed(ican_parallel)(seq, queue) for seq in seqs
    )

    time.sleep(0.5)


def main():
    data_pipeline.main()
    run_parallel_with_bars('data/ptms')


if __name__ == '__main__':
    main()
