# download data from each Database into PTM folder

import os
import gzip
import shutil
import requests
from Bio import Entrez


S_NITRO_DIR = '../data/s_nitrosylation'


def get_uniprot_seqs(uniprot_ids, filepath):
    if os.path.exists(filepath):
        os.remove(filepath)
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


def get_uniprot_seqs_from_names(uniprot_names,filepath): # very slow
    if os.path.exists(filepath):
        os.remove(filepath)

    base_search_url = "https://rest.uniprot.org/uniprotkb/search"
    base_fasta_url = "https://rest.uniprot.org/uniprotkb/{}.fasta"

    not_found = 0
    downloaded = 0
    failed = 0
    d_failed = 0

    with open(filepath, "w") as out_f:
        for name in uniprot_names:
            params = {
                "query": f'(protein_name:"{name}" OR gene:"{name}") AND organism_id:9606', # leave out reviewed:True ?
                "fields": "accession",
                "format": "json",
                "size": 1 # what if I want all matches, not just the first?
            }

            search_resp = requests.get(base_search_url, params=params, timeout=30)
            if search_resp.ok:
                results = search_resp.json()
                if results.get("results"):
                    accession = results["results"][0]["primaryAccession"]

                    fasta_resp = requests.get(base_fasta_url.format(accession), timeout=30)
                    if fasta_resp.ok:
                        out_f.write(fasta_resp.text)
                        # print(f"Downloaded: {name} ({accession})")
                        downloaded += 1
                    else:
                        # print(f"FASTA download failed for {name}")
                        d_failed += 1
                else:
                    # print(f"No UniProt entry found for {name}")
                    not_found += 1
            else:
                # print(f"Search failed for {name}")
                failed += 1

    print(f"\nDidn't find {not_found} entries")
    print(f"{failed} entries failed")
    print(f"Downloaded {downloaded} entries and {d_failed} failed to download\n")


def swiss_prot():
    filepath = os.path.join(S_NITRO_DIR, 'swissProt.fasta')

    url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0702%29%29"
    with requests.get(url, stream=True, timeout=10) as request:
        request.raise_for_status()
        with open(filepath, 'wb') as f:
            for chunk in request.iter_content(chunk_size=2**20):
                f.write(chunk)
    print("SwissProt download successful")


def ncbi():
    filepath = os.path.join(S_NITRO_DIR, 'ncbi.fasta')

    Entrez.email = "rebeg02@zedat.fu-berlin.de"
    query = '"Homo sapiens"[Organism] AND S-nitrosylation[All Fields]'
    handle = Entrez.esearch(db="protein", term=query, retmax=4000)
    record = Entrez.read(handle)
    ids = record["IdList"]

    handle = Entrez.efetch(db="protein", id=",".join(ids), rettype="fasta", retmode="text")
    fasta_data = handle.read()
    with open(filepath, "w") as f:
        f.write(fasta_data)
    print("NCBI download successful")


def db_ptm():
    zippath = os.path.join(S_NITRO_DIR, 'dbPTM.gz')
    txtpath = os.path.join(S_NITRO_DIR, 'dbPTM.txt')
    filepath = os.path.join(S_NITRO_DIR, 'dbPTM.fasta')

    url = "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/S-nitrosylation.gz"

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

    get_uniprot_seqs(uniprot_ids,filepath)

    if os.path.exists(zippath):
        os.remove(zippath)
    if os.path.exists(txtpath):
        os.remove(txtpath)

    print("dbPTM download successful")


def ptmd():
    zippath = os.path.join(S_NITRO_DIR, 'ptmd.zip')
    txtpath = os.path.join(S_NITRO_DIR, 'S-Nitrosylation')
    filepath = os.path.join(S_NITRO_DIR, 'ptmd.fasta')

    url = "https://ptmd.biocuckoo.cn/Download/S-Nitrosylation.zip"

    response = requests.get(url, timeout=30)
    if response.ok:
        with open(zippath, "wb") as f:
            f.write(response.content)
    else:
        print("Failed to download file:", response.status_code)

    shutil.unpack_archive(zippath, S_NITRO_DIR, 'zip')

    uniprot_ids = []

    with open(txtpath, 'r') as f:
        for line in f:
            if "S-Nitrosylation" in line:
                uniprot_ids.append(line[0:6])

    uniprot_ids = list(set(uniprot_ids))

    get_uniprot_seqs(uniprot_ids,filepath)

    if os.path.exists(zippath):
        os.remove(zippath)
    if os.path.exists(txtpath):
        os.remove(txtpath)

    print("ptmd download successful")


def ptm_code2():
    zippath = os.path.join(S_NITRO_DIR, 'PTMcode2.zip')
    txtpath = os.path.join(S_NITRO_DIR, 'PTMcode2_associations_within_proteins.txt')
    filepath = os.path.join(S_NITRO_DIR, 'PTMcode2.fasta')

    url = "https://ptmcode.embl.de/data/PTMcode2_associations_within_proteins.txt.gz"

    response = requests.get(url, timeout=30)
    if response.ok:
        with open(zippath, "wb") as f:
            f.write(response.content)
    else:
        print("Failed to download file:", response.status_code)

    with gzip.open(zippath, 'rb') as f_in:
        with open(txtpath, 'ab') as f_out:
            shutil.copyfileobj(f_in, f_out)

    uniprot_names = []

    with open(txtpath, 'r') as f:
        for line in f:
            if "nitrosylation" in line:
                idx = line.find("Homo sapiens")
                if idx != -1:
                    before_human = line[:idx].strip()
                    words = before_human.split()
                    if words:
                        uniprot_names.append(words[0])

    uniprot_names = list(set(uniprot_names))
    get_uniprot_seqs_from_names(uniprot_names,filepath)

    if os.path.exists(zippath):
        os.remove(zippath)
    if os.path.exists(txtpath):
        os.remove(txtpath)

    print("PTMcode2 download successful")


def main():
    print("Starting downloads...")
    # SwissProt https://www.uniprot.org/
    # swiss_prot()

    # National Library of Medicine https://www.ncbi.nlm.nih.gov/
    # ncbi()

    # dbPTM https://biomics.lab.nycu.edu.tw/dbPTM/index.php
    # db_ptm()

    # ptmd https://ptmd.biocuckoo.cn/
    # ptmd()

    # PTMcode2 https://ptmcode.embl.de/index.cgi
    # ptm_code2() # only downloads associations within proteins (for now?)

    print("All downloads successful\n")


if __name__ == '__main__':
    main()
