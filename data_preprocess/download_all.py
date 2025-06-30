# download data from each Database into PTM folder

import os
import gzip
import shutil
from io import StringIO
import requests
import pandas as pd
from Bio import Entrez
from bs4 import BeautifulSoup


def get_uniprot_seqs(uniprot_ids, filepath):
    if os.path.exists(filepath):
        os.remove(filepath)
    url = 'https://rest.uniprot.org/uniprotkb/stream'
    chunk_size = 500
    for i in range(0, len(uniprot_ids), chunk_size):
        chunk = uniprot_ids[i:(i + chunk_size)]
        params = {
            'format': 'fasta',
            'query': ' OR '.join([f'accession:{uid}' for uid in chunk])
        }
        response = requests.get(url, params=params, timeout=30)
        if response.ok:
            with open(filepath, 'a') as f:
                f.write(response.text)
        else:
            print('Batch download failed:', response.status_code)


def get_uniprot_seqs_from_names(uniprot_names,filepath): # very slow
    if os.path.exists(filepath):
        os.remove(filepath)

    base_search_url = 'https://rest.uniprot.org/uniprotkb/search'
    base_fasta_url = 'https://rest.uniprot.org/uniprotkb/{}.fasta'

    not_found = 0
    downloaded = 0
    failed = 0
    d_failed = 0

    with open(filepath, 'w') as out_f:
        for name in uniprot_names:
            params = {
                'query': f'(protein_name:"{name}" OR gene:"{name}") AND organism_id:9606', # leave out reviewed:True ?
                'fields': 'accession',
                'format': 'json',
                'size': 1 # what if I want all matches, not just the first?
            }

            search_resp = requests.get(base_search_url, params=params, timeout=30)
            if search_resp.ok:
                results = search_resp.json()
                if results.get('results'):
                    accession = results['results'][0]['primaryAccession']

                    fasta_resp = requests.get(base_fasta_url.format(accession), timeout=30)
                    if fasta_resp.ok:
                        out_f.write(fasta_resp.text)
                        # print(f'Downloaded: {name} ({accession})')
                        downloaded += 1
                    else:
                        # print(f'FASTA download failed for {name}')
                        d_failed += 1
                else:
                    # print(f'No UniProt entry found for {name}')
                    not_found += 1
            else:
                # print(f'Search failed for {name}')
                failed += 1

    print(f'\nDid not find {not_found} entries')
    print(f'{failed} entries failed')
    print(f'Downloaded {downloaded} entries and {d_failed} failed to download\n')

#https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28protein_name%3Abla%29+OR+%28gene%3Abla%29%29
#https://rest.uniprot.org/uniparc/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28%28protein_name%3Aecl%29+OR+%28gene%3Abla%29%29%29
# def get_records(uniprot_names,filepath):
#     if os.path.exists(filepath):
#         os.remove(filepath)
#     chunk_size = 250
#     for i in range(0, len(uniprot_names), chunk_size):
#         chunk = uniprot_names[i:(i + chunk_size)]
#         query = '+OR+'.join([f'%28organism_id%3A9606%29+AND+%28%28protein_name%3A{name}%29+OR+%28gene%3A{name}%29%29' for name in chunk])
#         url = f'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28{query}%29'
#         #print(url)
#         response = requests.get(url, timeout=30)
#         if response.ok:
#             with open(filepath, 'a') as f:
#                 f.write(response.text)
#         else:
#             print('Batch download failed:', response.status_code)


def swiss_prot(ptm_dir, url):
    filepath = os.path.join(ptm_dir, 'swissProt.fasta')

    with requests.get(url, stream=True, timeout=10) as request:
        request.raise_for_status()
        with open(filepath, 'wb') as f:
            for chunk in request.iter_content(chunk_size=2**20):
                f.write(chunk)
    print('SwissProt download successful')


def ncbi(ptm_dir, query):
    filepath = os.path.join(ptm_dir, 'ncbi.fasta')

    Entrez.email = 'rebeg02@zedat.fu-berlin.de'
    handle = Entrez.esearch(db='protein', term=query, retmax=10000)
    record = Entrez.read(handle)
    ids = record['IdList']

    handle = Entrez.efetch(db='protein', id=','.join(ids), rettype='fasta', retmode='text')
    fasta_data = handle.read()
    with open(filepath, 'w') as f:
        f.write(fasta_data)
    print('NCBI download successful')


def db_ptm(ptm_dir, urls):
    zippath = os.path.join(ptm_dir, 'dbPTM.gz')
    txtpath = os.path.join(ptm_dir, 'dbPTM.txt')
    filepath = os.path.join(ptm_dir, 'dbPTM.fasta')

    for url in urls:
        response = requests.get(url, timeout=30)
        if response.ok:
            with open(zippath, 'wb') as f:
                f.write(response.content)
        else:
            print('Failed to download file:', response.status_code)

        with gzip.open(zippath, 'rb') as f_in:
            with open(txtpath, 'ab') as f_out:
                shutil.copyfileobj(f_in, f_out)

    uniprot_ids = []

    with open(txtpath, 'r') as f:
        for line in f:
            if '\x00' in line:
                continue
            idx = line.find('HUMAN')
            if idx != -1:
                after_human = line[idx + len('HUMAN'):].strip()
                words = after_human.split()
                if words:
                    uniprot_ids.append(words[0])

    uniprot_ids = list(set(uniprot_ids))

    get_uniprot_seqs(uniprot_ids,filepath)

    if os.path.exists(zippath):
        os.remove(zippath)
    if os.path.exists(txtpath):
        os.remove(txtpath)

    print('dbPTM download successful')


def ptmd(ptm_dir, ptm, url, ptmd_word):
    zippath = os.path.join(ptm_dir, 'ptmd.zip')
    txtpath = os.path.join(ptm_dir, ptm)
    filepath = os.path.join(ptm_dir, 'ptmd.fasta')

    response = requests.get(url, timeout=30)
    if response.ok:
        with open(zippath, 'wb') as f:
            f.write(response.content)
    else:
        print('Failed to download file:', response.status_code)

    shutil.unpack_archive(zippath, ptm_dir, 'zip')

    uniprot_ids = []

    with open(txtpath, 'r') as f:
        for line in f:
            if ptmd_word in line:
                uniprot_ids.append(line[0:6])

    uniprot_ids = list(set(uniprot_ids))

    get_uniprot_seqs(uniprot_ids,filepath)

    if os.path.exists(zippath):
        os.remove(zippath)
    if os.path.exists(txtpath):
        os.remove(txtpath)

    print('ptmd download successful')


def ptm_code2(ptm_dir, ptm_code2_word):
    zippath = os.path.join(ptm_dir, 'PTMcode2.zip')
    txtpath = os.path.join(ptm_dir, 'PTMcode2_associations_within_proteins.txt')
    filepath = os.path.join(ptm_dir, 'PTMcode2.fasta')

    url = 'https://ptmcode.embl.de/data/PTMcode2_associations_within_proteins.txt.gz'

    response = requests.get(url, timeout=30)
    if response.ok:
        with open(zippath, 'wb') as f:
            f.write(response.content)
    else:
        print('Failed to download file:', response.status_code)

    with gzip.open(zippath, 'rb') as f_in:
        with open(txtpath, 'ab') as f_out:
            shutil.copyfileobj(f_in, f_out)

    uniprot_names = []

    with open(txtpath, 'r') as f:
        for line in f:
            if ptm_code2_word in line:
                idx = line.find('Homo sapiens')
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

    print('PTMcode2 download successful')


def qptm(ptm_dir, ptm):
    filepath = os.path.join(ptm_dir, 'qPTM.fasta')
    txtpath = 'local_data/qPTM_all_data.txt'

    uniprot_ids = []

    with open(txtpath, 'r') as f:
        for line in f:
            words = line.strip().split('\t')
            if words[0] == 'Human' and words[5] == ptm:
                uniprot_ids.append(words[2])

    uniprot_ids = list(set(uniprot_ids))

    get_uniprot_seqs(uniprot_ids,filepath)

    print('qPTM download successful')


def unipep(ptm_dir):
    filepath = os.path.join(ptm_dir, 'unipep.fasta')

    url = 'https://db.systemsbiology.net/sbeams/cgi/Glycopeptide/browse_glycopeptides.cgi'

    response = requests.get(url, timeout=30)
    if response.ok:
        soup = BeautifulSoup(response.text, 'html.parser')
        table = soup.find('table')
        table_html = str(table)
        table_io = StringIO(table_html)
        df = pd.read_html(table_io)[0]
        df = df.dropna()
        protein_symbols = df[2].tolist()
        del protein_symbols[0]
        protein_symbols = list(set(protein_symbols))
        get_uniprot_seqs_from_names(protein_symbols,filepath)
    else:
        print('Failed to download file:', response.status_code)

    print('Unipep download successful')


def databases(ptm_dir, swiss_prot_url, query, db_ptm_urls, ptm, ptmd_url, ptmd_word, ptm_code2_word):
    swiss_prot(ptm_dir, swiss_prot_url)
    ncbi(ptm_dir, query)
    db_ptm(ptm_dir, db_ptm_urls)
    ptmd(ptm_dir, ptm, ptmd_url, ptmd_word)
    ptm_code2(ptm_dir, ptm_code2_word)
    if ptm != 'S-Nitrosylation':
        qptm(ptm_dir,ptm)
        if ptm == 'Glycosylation':
            unipep(ptm_dir)


def main():
    print('Starting downloads...')

    databases('data/glycosylation',
              'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28ft_carbohyd%3AGlcNAc%29%29',
              '"Homo sapiens"[Organism] AND glycosylated[All Fields]',
              ['https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/C-linked%20Glycosylation.gz',
              'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/N-linked%20Glycosylation.gz',
              'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/O-linked%20Glycosylation.gz',
              'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/S-linked%20Glycosylation.gz'],
              'Glycosylation',
              'https://ptmd.biocuckoo.cn/Download/Glycosylation.zip',
              'glycosylation',
              'glycosylation')

    databases('data/s_nitrosylation',
              'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0702%29%29',
              '"Homo sapiens"[Organism] AND S-nitrosylation[All Fields]',
              ['https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/S-nitrosylation.gz'],
              'S-Nitrosylation',
              'https://ptmd.biocuckoo.cn/Download/S-Nitrosylation.zip',
              'S-Nitrosylation',
              'nitrosylation')

    databases('data/acetylation',
              'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0007%29%29',
              '"Homo sapiens"[Organism] AND acetylated[All Fields]',
              ['https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Acetylation.gz'],
              'Acetylation',
              'https://ptmd.biocuckoo.cn/Download/Acetylation.zip',
              'Acetylation',
              'acetylation')

    databases('data/methylation',
              'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0488%29%29',
              '"Homo sapiens"[Organism] AND Methylation[All Fields]',
              ['https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Methylation.gz'],
              'Methylation',
              'https://ptmd.biocuckoo.cn/Download/Methylation.zip',
              'Methylation',
              'methylation')

    print('All downloads successful\n')


if __name__ == '__main__':
    main()
