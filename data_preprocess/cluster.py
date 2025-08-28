'''
Module to cluster multi-FASTA files using CD-HIT.

This script applies CD-HIT to cluster sequences in PTM and non-PTM datasets
to remove redundancy. Clustered sequences are saved in new multi-FASTA files.
'''

import subprocess
from pathlib import Path


def run_cd_hit(input_path, output_path, c=0.4, n=2, M=16000, T=8):
    '''
    Run CD-HIT to cluster sequences in a multi-FASTA file.

    Args:
        input_path (Path): Path to the input multi-FASTA file.
        output_path (Path): Path where the clustered multi-FASTA will be saved.
        c (float): Sequence identity threshold (default 0.4).
        n (int): Word length (default 2, recommended for protein sequences).
        M (int): Memory limit in MB (default 16000).
        T (int): Number of threads (default 8).
    '''
    cmd = [
        'cd-hit',
        '-i', str(input_path),
        '-o', str(output_path),
        '-c', str(c),
        '-n', str(n),
        '-M', str(M),
        '-T', str(T)
    ]
    subprocess.run(cmd, check=True)


def main(merged_path_str):
    '''
    Cluster PTM and corresponding non-PTM multi-FASTA files using CD-HIT.

    Args:
        merged_path_str (str): Path to the merged PTM multi-FASTA file.
    '''
    merged_path = Path(merged_path_str)
    clustered_path = merged_path.with_name(merged_path.stem.replace('merged', 'clustered') + merged_path.suffix)

    run_cd_hit(merged_path, clustered_path)

    ptm_name = merged_path.parent.name
    base_dir = merged_path.parents[1]
    no_ptm_dir = base_dir.with_name('no_ptm')
    no_ptm_path = no_ptm_dir / f'filtered_no_{ptm_name}.fasta'
    clustered_no_ptm_path = no_ptm_path.with_name(no_ptm_path.stem.replace('filtered', 'clustered') + no_ptm_path.suffix)

    run_cd_hit(no_ptm_path, clustered_no_ptm_path)
