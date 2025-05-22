# download data from each Database into PTM folder

import os
import gzip
import shutil
import requests
from Bio import Entrez


GLYCO_DIR = '/Users/rebeccagrevens/Documents/ptm/data/glycosylation'


def swiss_prot():
    filepath = os.path.join(GLYCO_DIR, 'swissProt.fasta')

    url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28ft_carbohyd%3AGlcNAc%29%29"

    with requests.get(url, stream=True, timeout=10) as request:
        request.raise_for_status()
        with open(filepath, 'wb') as f:
            for chunk in request.iter_content(chunk_size=2**20):
                f.write(chunk)
        print("SwissProt download succesful")


def ncbi():
    filepath = os.path.join(GLYCO_DIR, 'ncbi.fasta')

    Entrez.email = "rebeg02@zedat.fu-berlin.de"
    query = '"Homo sapiens"[Organism] AND glycosylated[All Fields] '
    handle = Entrez.esearch(db="protein", term=query, retmax=4000)
    record = Entrez.read(handle)
    ids = record["IdList"]

    handle = Entrez.efetch(db="protein", id=",".join(ids), rettype="fasta", retmode="text")
    fasta_data = handle.read()
    with open(filepath, "w") as f:
        f.write(fasta_data)
    print("NCBI download succesful")


def db_ptm():
    zippath = os.path.join(GLYCO_DIR, 'dbPTM.gz')
    txtpath = os.path.join(GLYCO_DIR, 'dbPTM.txt')
    filepath = os.path.join(GLYCO_DIR, 'dbPTM.fasta')

    urls = ["https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/C-linked%20Glycosylation.gz",
            "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/N-linked%20Glycosylation.gz",
            "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/O-linked%20Glycosylation.gz",
            "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/S-linked%20Glycosylation.gz"]

    for url in urls:
        response = requests.get(url, timeout=30)
        if response.ok:
            with open(zippath, "wb") as f:
                f.write(response.content)
        else:
            print("Failed to download file:", response.status_code)

        with gzip.open(zippath, 'rb') as f_in:
            with open(txtpath, 'ab') as f_out:
                shutil.copyfileobj(f_in, f_out)

    uniprot_ids = []

    with open(txtpath, 'r') as f:
        for line in f:
            if "\x00" in line:
                continue
            idx = line.find("HUMAN")
            if idx != -1:
                after_human = line[idx + len("HUMAN"):].strip()
                words = after_human.split()
                if words:
                    uniprot_ids.append(words[0])

    uniprot_ids = list(set(uniprot_ids))

    url = "https://rest.uniprot.org/uniprotkb/stream"
    chunk_size = 500
    for i in range(0, len(uniprot_ids), chunk_size):
        chunk = uniprot_ids[i:(i + chunk_size)]
        params = {
            "format": "fasta",
            "query": " OR ".join([f"accession:{uid}" for uid in chunk])
        }
        response = requests.get(url, params=params, timeout=30)
        if response.ok:
            with open(filepath, "a") as f:
                f.write(response.text)
        else:
            print("Batch download failed:", response.status_code)

    if os.path.exists(zippath):
        os.remove(zippath)
    if os.path.exists(txtpath):
        os.remove(txtpath)

    print("dbPTM download succesful")


def main():
    # SwissProt
    # swiss_prot()

    # National Library of Medicine
    # ncbi()

    # dbPTM
    db_ptm()

    print("All downloads succesful")


if __name__ == '__main__':
    main()
