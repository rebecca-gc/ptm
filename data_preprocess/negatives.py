# download data from each Database into PTM folder

import os
import requests
from Bio import SeqIO


NO_PTM_DIR = "/Users/rebeccagrevens/Documents/ptm/data/no_ptm"

DIRS = ['/Users/rebeccagrevens/Documents/ptm/data/glycosylation',
        '/Users/rebeccagrevens/Documents/ptm/data/s_nitrosylation',
        '/Users/rebeccagrevens/Documents/ptm/data/acetylation']


def swiss_prot(files):
    urls = ["https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+NOT+%28ft_carbohyd%3AGlcNAc%29%29",
            "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+NOT+%28keyword%3AKW-0702%29%29",
            "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+NOT+%28keyword%3AKW-0007%29%29"]

    for i in range(0,len(urls)):
        filepath = os.path.join(NO_PTM_DIR, files[i])

        with requests.get(urls[i], stream=True, timeout=10) as request:
            request.raise_for_status()
            with open(filepath, 'wb') as f:
                for chunk in request.iter_content(chunk_size=2**20):
                    f.write(chunk)
    print("SwissProt download successful")


def filter_false_negatives(files):
    for i in range(0,len(files)):
        compare_path = os.path.join(DIRS[i], 'merged.fasta')
        target_path = os.path.join(NO_PTM_DIR, files[i])

        compare_seqs = {str(record.seq) for record in SeqIO.parse(compare_path, "fasta")}
        target_seqs = {str(record.seq) for record in SeqIO.parse(target_path, "fasta")}

        common_seqs = compare_seqs & target_seqs

        print(f"Count of common sequences: {len(common_seqs)}")

        filtered_path = os.path.join(NO_PTM_DIR, f"filtered_{files[i]}")

        with open(filtered_path, "w") as filtered:
            for record in SeqIO.parse(target_path, "fasta"):
                if record.seq not in common_seqs:
                    SeqIO.write(record, filtered, "fasta")
        
        if os.path.exists(target_path):
            os.remove(target_path)
        
        filtered_seqs = {str(record.seq) for record in SeqIO.parse(filtered_path, "fasta")}
        print(f"PTM: {len(compare_seqs)}, NO_PTM: {len(filtered_seqs)}")


def main():
    files = ['no_glyco.fasta', 'no_s_nitro.fasta', 'no_acet.fasta']
    print("Starting downloads...")
    # SwissProt https://www.uniprot.org/
    # swiss_prot(files)

    print("All downloads successful")

    filter_false_negatives(files)

    print("Filtered out all false negatives")


if __name__ == '__main__':
    main()
