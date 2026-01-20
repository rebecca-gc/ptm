"""
Module for downloading protein sequence datasets from various sources (UniProt, NCBI, dbPTM,
PTMD, PTMcode2, qPTM, UniPep). Handles extraction of UniProt IDs/names and retrieval of
FASTA sequences for specific post-translational modifications (PTMs).
"""

import gzip
import os
import shutil
from io import StringIO

import pandas as pd
import requests
from Bio import Entrez
from bs4 import BeautifulSoup


def get_uniprot_seqs(uniprot_ids, filepath):
    """
    Download UniProt sequences by accession IDs in batches and save to a multi-FASTA file.

    Parameters
    ----------
    uniprot_ids : list of str
        UniProt accession IDs.
    filepath : str
        Path to the output multi-FASTA file.
    """
    if os.path.exists(filepath):
        os.remove(filepath)
    url = "https://rest.uniprot.org/uniprotkb/stream"
    chunk_size = 500
    for i in range(0, len(uniprot_ids), chunk_size):
        chunk = uniprot_ids[i : (i + chunk_size)]
        params = {
            "format": "fasta",
            "query": " OR ".join([f"accession:{uid}" for uid in chunk]),
        }
        response = requests.get(url, params=params, timeout=30)
        if response.ok:
            with open(filepath, "a") as f:
                f.write(response.text)
        else:
            print("Batch download failed:", response.status_code)


def get_uniprot_seqs_from_names(uniprot_names, filepath):
    """
    Search UniProt by protein/gene names, download corresponding sequences, and save to multi-FASTA.

    Parameters
    ----------
    uniprot_names : list of str
        Protein or gene names to search in UniProt.
    filepath : str
        Path to the output multi-FASTA file.
    """
    if os.path.exists(filepath):
        os.remove(filepath)

    base_search_url = "https://rest.uniprot.org/uniprotkb/search"
    base_fasta_url = "https://rest.uniprot.org/uniprotkb/{}.fasta"

    not_found = downloaded = failed = d_failed = 0

    with open(filepath, "w") as out_f:
        for name in uniprot_names:
            params = {
                "query": f'(protein_name:"{name}" OR gene:"{name}") AND organism_id:9606',
                "fields": "accession",
                "format": "json",
                "size": 1,
            }

            search_resp = requests.get(
                base_search_url, params=params, timeout=30
            )
            if search_resp.ok:
                results = search_resp.json()
                if results.get("results"):
                    accession = results["results"][0]["primaryAccession"]
                    fasta_resp = requests.get(
                        base_fasta_url.format(accession), timeout=30
                    )
                    if fasta_resp.ok:
                        out_f.write(fasta_resp.text)
                        downloaded += 1
                    else:
                        d_failed += 1
                else:
                    not_found += 1
            else:
                failed += 1

    print(f"\nDid not find {not_found} entries")
    print(f"{failed} entries failed")
    print(
        f"Downloaded {downloaded} entries and {d_failed} failed to download\n"
    )


def swiss_prot(ptm_dir, url):
    """
    Download reviewed SwissProt sequences for a PTM and save as multi-FASTA.

    Parameters
    ----------
    ptm_dir : str
        Directory to store the output file.
    url : str
        Download URL for SwissProt FASTA.
    """
    filepath = os.path.join(ptm_dir, "swissProt.fasta")
    with requests.get(url, stream=True, timeout=10) as request:
        request.raise_for_status()
        with open(filepath, "wb") as f:
            for chunk in request.iter_content(chunk_size=2**20):
                f.write(chunk)
    print("SwissProt download successful")


def ncbi(ptm_dir, query):
    """
    Download protein sequences from NCBI for a given PTM query.

    Parameters
    ----------
    ptm_dir : str
        Directory to store the output file.
    query : str
        NCBI search query (e.g., PTM type).
    """
    filepath = os.path.join(ptm_dir, "ncbi.fasta")
    Entrez.email = "rebeg02@zedat.fu-berlin.de"
    handle = Entrez.esearch(db="protein", term=query, retmax=10000)
    record = Entrez.read(handle)
    ids = record["IdList"]
    handle = Entrez.efetch(
        db="protein", id=",".join(ids), rettype="fasta", retmode="text"
    )
    fasta_data = handle.read()
    with open(filepath, "w") as f:
        f.write(fasta_data)
    print("NCBI download successful")


def db_ptm(ptm_dir, urls):
    """
    Download and extract UniProt IDs from dbPTM, then fetch corresponding FASTA sequences.

    Parameters
    ----------
    ptm_dir : str
        Directory to store the output files.
    urls : list of str
        Download URLs for dbPTM .gz files.
    """
    zippath = os.path.join(ptm_dir, "dbPTM.gz")
    txtpath = os.path.join(ptm_dir, "dbPTM.txt")
    filepath = os.path.join(ptm_dir, "dbPTM.fasta")

    for url in urls:
        response = requests.get(url, timeout=30)
        if response.ok:
            with open(zippath, "wb") as f:
                f.write(response.content)
        else:
            print("Failed to download file:", response.status_code)
        with gzip.open(zippath, "rb") as f_in, open(txtpath, "ab") as f_out:
            shutil.copyfileobj(f_in, f_out)

    uniprot_ids = []
    with open(txtpath, "r") as f:
        for line in f:
            if "\x00" in line:
                continue
            idx = line.find("HUMAN")
            if idx != -1:
                words = line[idx + len("HUMAN") :].strip().split()
                if words:
                    uniprot_ids.append(words[0])

    get_uniprot_seqs(list(set(uniprot_ids)), filepath)

    if os.path.exists(zippath):
        os.remove(zippath)
    if os.path.exists(txtpath):
        os.remove(txtpath)

    print("dbPTM download successful")


def ptmd(ptm_dir, ptm, url, ptmd_word):
    """
    Download and extract UniProt IDs from PTMD, then fetch sequences.

    Parameters
    ----------
    ptm_dir : str
        Directory to store the output files.
    ptm : str
        PTM type.
    url : str
        Download URL for PTMD.
    ptmd_word : str
        Keyword to filter PTMD entries.
    """
    zippath = os.path.join(ptm_dir, "ptmd.zip")
    txtpath = os.path.join(ptm_dir, ptm)
    filepath = os.path.join(ptm_dir, "ptmd.fasta")

    response = requests.get(url, timeout=30)
    if response.ok:
        with open(zippath, "wb") as f:
            f.write(response.content)
    else:
        print("Failed to download file:", response.status_code)

    shutil.unpack_archive(zippath, ptm_dir, "zip")

    uniprot_ids = []
    with open(txtpath, "r") as f:
        for line in f:
            if ptmd_word in line:
                uniprot_ids.append(line[0:6])

    get_uniprot_seqs(list(set(uniprot_ids)), filepath)

    if os.path.exists(zippath):
        os.remove(zippath)
    if os.path.exists(txtpath):
        os.remove(txtpath)

    print("ptmd download successful")


def ptm_code2(ptm_dir, ptm_code2_word):
    """
    Download PTMcode2 data, extract human proteins, and fetch corresponding FASTA sequences.

    Parameters
    ----------
    ptm_dir : str
        Directory to store the output files.
    ptm_code2_word : str
        Keyword to filter PTMcode2 entries for the target PTM.
    """
    zippath = os.path.join(ptm_dir, "PTMcode2.zip")
    txtpath = os.path.join(ptm_dir, "PTMcode2_associations_within_proteins.txt")
    filepath = os.path.join(ptm_dir, "PTMcode2.fasta")
    url = "https://ptmcode.embl.de/data/PTMcode2_associations_within_proteins.txt.gz"

    response = requests.get(url, timeout=30)
    if response.ok:
        with open(zippath, "wb") as f:
            f.write(response.content)
    else:
        print("Failed to download file:", response.status_code)

    with gzip.open(zippath, "rb") as f_in, open(txtpath, "ab") as f_out:
        shutil.copyfileobj(f_in, f_out)

    uniprot_names = []
    with open(txtpath, "r") as f:
        for line in f:
            if ptm_code2_word in line:
                idx = line.find("Homo sapiens")
                if idx != -1:
                    words = line[:idx].strip().split()
                    if words:
                        uniprot_names.append(words[0])

    get_uniprot_seqs_from_names(list(set(uniprot_names)), filepath)

    if os.path.exists(zippath):
        os.remove(zippath)
    if os.path.exists(txtpath):
        os.remove(txtpath)

    print("PTMcode2 download successful")


def qptm(ptm_dir, ptm):
    """
    Download sequences from qPTM for a given PTM type.

    Parameters
    ----------
    ptm_dir : str
        Directory to store the output file.
    ptm : str
        Target PTM type (e.g., 'Glycosylation').
    """
    filepath = os.path.join(ptm_dir, "qPTM.fasta")
    txtpath = "local_data/qPTM_all_data.txt"

    uniprot_ids = []
    with open(txtpath, "r") as f:
        for line in f:
            words = line.strip().split("\t")
            if words[0] == "Human" and words[5] == ptm:
                uniprot_ids.append(words[2])

    get_uniprot_seqs(list(set(uniprot_ids)), filepath)

    print("qPTM download successful")


def unipep(ptm_dir):
    """
    Download UniPep glycopeptide data and fetch corresponding UniProt sequences.

    Parameters
    ----------
    ptm_dir : str
        Directory to store the output file.
    """
    filepath = os.path.join(ptm_dir, "unipep.fasta")
    url = "https://db.systemsbiology.net/sbeams/cgi/Glycopeptide/browse_glycopeptides.cgi"

    response = requests.get(url, timeout=30)
    if response.ok:
        soup = BeautifulSoup(response.text, "html.parser")
        table = soup.find("table")
        table_io = StringIO(str(table))
        df = pd.read_html(table_io)[0].dropna()
        protein_symbols = list(set(df[2].tolist()[1:]))
        get_uniprot_seqs_from_names(protein_symbols, filepath)
    else:
        print("Failed to download file:", response.status_code)

    print("Unipep download successful")


def databases(
    ptm_dir,
    swiss_prot_url,
    query,
    db_ptm_urls,
    ptm,
    ptmd_url,
    ptmd_word,
    ptm_code2_word,
):
    """
    Orchestrates downloading all databases for a specific PTM.

    Parameters
    ----------
    ptm_dir : str
        Directory to store database files.
    swiss_prot_url : str
        URL for SwissProt FASTA download.
    query : str
        NCBI query string for the PTM.
    db_ptm_urls : list of str
        List of dbPTM download URLs.
    ptm : str
        Target PTM type.
    ptmd_url : str
        URL for PTMD download.
    ptmd_word : str
        Keyword to filter PTMD entries.
    ptm_code2_word : str
        Keyword to filter PTMcode2 entries.
    """
    swiss_prot(ptm_dir, swiss_prot_url)
    ncbi(ptm_dir, query)
    db_ptm(ptm_dir, db_ptm_urls)
    ptmd(ptm_dir, ptm, ptmd_url, ptmd_word)
    ptm_code2(ptm_dir, ptm_code2_word)

    if ptm != "S-Nitrosylation":
        qptm(ptm_dir, ptm)
        if ptm == "Glycosylation":
            unipep(ptm_dir)


def main():
    """
    Main entry point: sets up directories and downloads all PTM datasets sequentially.
    """
    print("Starting downloads...")

    ptm_configs = [
        (
            "glycosylation",
            "KW-0325",
            '"Homo sapiens"[Organism] AND glycosylated[All Fields]',
            [
                "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/C-linked%20Glycosylation.gz",
                "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/N-linked%20Glycosylation.gz",
                "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/O-linked%20Glycosylation.gz",
                "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/S-linked%20Glycosylation.gz",
            ],
            "https://ptmd.biocuckoo.cn/Download/Glycosylation.zip",
            "Glycosylation",
            "glycosylation",
            "glycosylation",
        ),
        (
            "s_nitrosylation",
            "KW-0702",
            '"Homo sapiens"[Organism] AND S-nitrosylation[All Fields]',
            [
                "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/S-nitrosylation.gz"
            ],
            "https://ptmd.biocuckoo.cn/Download/S-Nitrosylation.zip",
            "S-Nitrosylation",
            "S-Nitrosylation",
            "nitrosylation",
        ),
        (
            "acetylation",
            "KW-0007",
            '"Homo sapiens"[Organism] AND acetylated[All Fields]',
            [
                "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Acetylation.gz"
            ],
            "https://ptmd.biocuckoo.cn/Download/Acetylation.zip",
            "Acetylation",
            "Acetylation",
            "acetylation",
        ),
        (
            "methylation",
            "KW-0488",
            '"Homo sapiens"[Organism] AND Methylation[All Fields]',
            [
                "https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Methylation.gz"
            ],
            "https://ptmd.biocuckoo.cn/Download/Methylation.zip",
            "Methylation",
            "Methylation",
            "methylation",
        ),
    ]

    for (
        ptm_name,
        kw,
        query,
        db_urls,
        ptmd_url,
        ptm,
        ptmd_word,
        ptm_code2_word,
    ) in ptm_configs:
        dir_name = f"test/{ptm_name}/databases"
        os.makedirs(dir_name, exist_ok=True)
        databases(
            dir_name,
            f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=((organism_id:9606) AND (reviewed:true) AND (keyword:{kw}))",
            query,
            db_urls,
            ptm,
            ptmd_url,
            ptmd_word,
            ptm_code2_word,
        )

    print("All downloads successful\n")


if __name__ == "__main__":
    main()
