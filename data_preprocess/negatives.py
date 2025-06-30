# download data from each Database into PTM folder

import os
import requests
from Bio import SeqIO


NO_PTM_DIR = 'data/no_ptm'

DIRS = ['data/glycosylation',
        'data/s_nitrosylation',
        'data/acetylation',
        'data/methylation']


def swiss_prot(files):
    urls = ['https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+NOT+%28ft_carbohyd%3AGlcNAc%29%29',
            'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+NOT+%28keyword%3AKW-0702%29%29',
            'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+NOT+%28keyword%3AKW-0007%29%29',
            'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+NOT+%28keyword%3AKW-0488%29%29']

    for i, url in enumerate(urls):
        filepath = os.path.join(NO_PTM_DIR, files[i])

        with requests.get(url, stream=True, timeout=10) as request:
            request.raise_for_status()
            with open(filepath, 'wb') as f:
                for chunk in request.iter_content(chunk_size=2**20):
                    f.write(chunk)


def filter_false_negatives(files):
    for i, file in enumerate(files):
        compare_path = os.path.join(DIRS[i], 'merged.fasta')
        target_path = os.path.join(NO_PTM_DIR, file)

        compare_seqs = {str(record.seq) for record in SeqIO.parse(compare_path, 'fasta')}
        target_seqs = {str(record.seq) for record in SeqIO.parse(target_path, 'fasta')}

        common_seqs = compare_seqs & target_seqs

        print(f'Count of common sequences: {len(common_seqs)}')

        filtered_path = os.path.join(NO_PTM_DIR, f'filtered_{file}')
        filtered3000_path = os.path.join(NO_PTM_DIR, f'filtered3000_{file}')

        if os.path.exists(filtered_path):
            os.remove(filtered_path)

        if os.path.exists(filtered3000_path):
            os.remove(filtered3000_path)

        with open(filtered_path, 'w') as filtered:
            for record in SeqIO.parse(target_path, 'fasta'):
                if record.seq not in common_seqs:
                    SeqIO.write(record, filtered, 'fasta')

        with open(filtered3000_path, 'w') as filtered:
            for record in SeqIO.parse(target_path, 'fasta'):
                if record.seq not in common_seqs and len(record.seq) <= 3000:
                    SeqIO.write(record, filtered, 'fasta')                            

        filtered_seqs = {str(record.seq) for record in SeqIO.parse(filtered_path, 'fasta')}
        print(f'PTM: {len(compare_seqs)}, NO_PTM: {len(filtered_seqs)}')


def main():
    files = ['no_glyco.fasta', 'no_s_nitro.fasta', 'no_acet.fasta', 'no_methyl.fasta']
    print('Starting negative downloads...')
    # SwissProt https://www.uniprot.org/
    swiss_prot(files)

    print('Downloads of negative sequences successful\n')

    filter_false_negatives(files)

    print('\nFiltered out all false negatives')


if __name__ == '__main__':
    main()
