import os
import sys
import time
import threading
import multiprocessing
from joblib import Parallel, delayed

sys.path.append(os.path.abspath("data_preprocess"))
sys.path.append(os.path.abspath("Source"))

import rfc
import data_pipeline
import ican


dirs = ['test', 'glycosylation', 's_nitrosylation', 'acetylation', 'methylation']
STEPS_PER_TASK = 100


def draw_progress(progress_dict, total_steps):
    ESC = "\033"
    CLEAR = f"{ESC}[2J"
    HIDE_CURSOR = f"{ESC}[?25l"
    SHOW_CURSOR = f"{ESC}[?25h"

    sys.stdout.write(CLEAR + HIDE_CURSOR)
    sys.stdout.flush()

    try:
        while True:
            sys.stdout.write(f"{ESC}[H")  # Cursor nach oben links
            for key in sorted(progress_dict.keys()):
                current = progress_dict[key]
                total = total_steps[key]
                percent = current / total
                bar_length = 40
                filled = int(bar_length * percent)
                bar = '█' * filled + '-' * (bar_length - filled)
                sys.stdout.write(f"{key:15} |{bar}| {percent * 100:5.1f}%\n")
            sys.stdout.flush()
            time.sleep(0.1)
            if all(progress_dict[k] >= total_steps[k] for k in progress_dict):
                break
    finally:
        sys.stdout.write(SHOW_CURSOR)
        sys.stdout.flush()


def ican_parallel(dir, queue, steps):
    seqs = os.path.join('data', dir, 'seqs.fasta')
    output = os.path.join('data', dir)
    sys.argv = ['ican.py', f'--output_path={output}', seqs]

    ican.main(queue=queue, task_name=dir)

    rfc.main(output)


def run_parallel_with_bars():
    steps_per_task = {dir: STEPS_PER_TASK for dir in dirs}

    manager = multiprocessing.Manager()
    progress_dict = manager.dict({dir: 0 for dir in dirs})
    queue = manager.Queue()

    threading.Thread(target=draw_progress, args=(progress_dict, steps_per_task), daemon=True).start()

    def progress_updater():
        while True:
            name, delta = queue.get()
            progress_dict[name] += delta
            if all(progress_dict[k] >= steps_per_task[k] for k in dirs):
                break

    threading.Thread(target=progress_updater, daemon=True).start()

    Parallel(n_jobs=-1)(
        delayed(ican_parallel)(dir, queue, steps_per_task[dir]) for dir in dirs
    )

    time.sleep(0.5)


def main():
    # data_pipeline.main()
    run_parallel_with_bars()


if __name__ == '__main__':
    main()
