'''
Module to retrieve disease information from UniProt and generate
disease-related multi-label datasets for PTMs.

Functions include:
- Extracting UniProt IDs from multi-FASTA files
- Downloading disease annotations via UniProt API
- Visualizing disease frequency
- Generating multi-label datasets combining PTM and disease info
'''

import os
import csv
from io import StringIO
from collections import Counter
import re
import requests
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO


def get_ids(filepath):
    '''
    Extract UniProt IDs from a multi-FASTA file.

    Args:
        filepath (str): Path to a multi-FASTA file.

    Returns:
        tuple: (list of UniProt IDs, list of other SeqIO records)
    '''
    uniprot_ids = []
    other_records = []
    for record in SeqIO.parse(filepath, 'fasta'):
        if record.id.startswith('sp|') or record.id.startswith('tr|'):
            name = record.id.split('|')[1]
            if name[-2] == '.':
                name = name[:-2]
            uniprot_ids.append(name)
        else:
            other_records.append(record)
    return uniprot_ids, other_records


def get_records(uniprot_ids, chunk_size=100):
    '''
    Download disease annotations for a list of UniProt IDs.

    Args:
        uniprot_ids (list): List of UniProt IDs.

    Returns:
        set: Set of tuples containing (accession, entry name, MIM string, sequence)
    '''
    records = set()
    
    base_url = 'https://rest.uniprot.org/uniprotkb/stream'
    fields = 'accession,id,cc_disease,sequence'
    
    for i in range(0, len(uniprot_ids), chunk_size):
        chunk = uniprot_ids[i:i + chunk_size]
        query = ' OR '.join([f'accession:{uid}' for uid in chunk])
        
        params = {
            'format': 'tsv',
            'fields': fields,
            'query': query
        }
        
        response = requests.get(base_url, params=params, timeout=30)
        
        if response.ok:
            tsv_data = StringIO(response.text)
            reader = csv.DictReader(tsv_data, delimiter='\t')
            for row in reader:
                acc = row.get('Entry', '')
                name = row.get('Entry Name', '')
                disease = row.get('Involvement in disease', '')
                mim_ids = re.findall(r'\[MIM:(\d+)', disease)
                mim_string = ' '.join([f'MIM:{m}' for m in mim_ids])
                seq = row.get('Sequence', '')
                records.add((acc, name, mim_string, seq))
        else:
            print(f'(Disease) Batch download failed: {response.status_code} - {response.text[:200]}')

    return records


def vis_disease(dir_path):
    '''
    Visualize the frequency of top 10 diseases in a PTM multi-FASTA file.

    Args:
        dir_path (str): Directory containing 'merged.fasta' for a PTM.
    '''
    diseases = []
    for record in SeqIO.parse(f'{dir_path}/merged.fasta', 'fasta'):
        if 'MIM:' in record.description:
            parts = record.description.split(' MIM:')
            for disease_part in parts[1:]:
                diseases.append(disease_part)

    counts = Counter(diseases)
    top10 = counts.most_common(10)
    labels = [d for d, _ in top10]
    values = [c for _, c in top10]

    if 'glycosylation' in dir_path:
        c = '#66c2a5'
        e = '#1b9e77'
    elif 's_nitrosylation' in dir_path:
        c = '#e78ac3'
        e = '#e7298a'
    elif 'acetylation' in dir_path:
        c = '#fc8d62'
        e = '#d95f02'
    elif 'methylation' in dir_path:
        c = '#8da0cb'
        e = '#7570b3'
    else:
        c = 'black'
        e = 'black'

    plt.figure(figsize=(5, 4))
    plt.bar(labels, values, color=c, edgecolor=e)
    plt.xlabel('MIM ID of Diseases', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    # plt.title(f'Frequency of Top 10 Diseases (PTM: {dir_path.split("/")[-1]})')
    plt.tight_layout()
    plt.savefig(f'{dir_path}/top10diseases_{dir_path.split("/")[-1]}.pdf')
    plt.clf()


def disease_stacked(ptms_dir):
    '''
    Plot the top 10 most frequent MIM IDs for each PTM as a stacked bar chart.

    Args:
        ptms_dir (str): Directory containing PTM subdirectories with 'merged.fasta' files.
    '''

    diseases_per_ptm = dict()
    for ptm in os.listdir(ptms_dir):
        ptm_path = os.path.join(ptms_dir, ptm)
        if os.path.isdir(ptm_path):
            diseases = []
            for record in SeqIO.parse(f'{ptm_path}/merged.fasta', 'fasta'):
                if 'MIM:' in record.description:
                    parts = record.description.split(' MIM:')
                    for disease_part in parts[1:]:
                        diseases.append(disease_part)
            diseases_per_ptm[ptm] = diseases

    all_mims = []
    for ptm, mim_list in diseases_per_ptm.items():
        all_mims.extend(mim_list)
    
    n = 10
    mim_counts = Counter(all_mims)
    top_mims = [mim for mim, _ in mim_counts.most_common(n)]

    df = pd.DataFrame(0, index=top_mims, columns=diseases_per_ptm.keys())
    for ptm, mim_list in diseases_per_ptm.items():
        counts = Counter(mim_list)
        for mim in top_mims:
            df.loc[mim, ptm] = counts.get(mim, 0)
    
    face_colors = {
        'glycosylation': '#66c2a5',
        's_nitrosylation': '#e78ac3',
        'acetylation': '#fc8d62',
        'methylation': '#8da0cb'
    }
    edge_colors = {
        'glycosylation': '#1b9e77',
        's_nitrosylation': '#e7298a',
        'acetylation': '#d95f02',
        'methylation': '#7570b3'
    }

    fig, ax = plt.subplots(figsize=(10,5))
    bottom = [0]*len(top_mims)
    for ptm in df.columns:
        ax.bar(
            df.index,
            df[ptm],
            bottom=bottom,
            color=face_colors.get(ptm, 'grey'),
            edgecolor=edge_colors.get(ptm, 'black'),
            label=ptm
        )
        bottom = [i+j for i,j in zip(bottom, df[ptm])]
    
    ax.set_xlabel('MIM ID of Diseases')
    ax.set_ylabel('Frequency')
    # ax.set_title(f'Top {n} MIM IDs by PTM')
    ax.legend()
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{ptms_dir}/top10diseases.pdf')
    plt.clf()


def main(filepath):
    '''
    Wrapper function to extract disease records from a PTM multi-FASTA file.

    Args:
        filepath (str): Path to a PTM multi-FASTA file.

    Returns:
        tuple: (set of disease records, list of non-UniProt records)
    '''
    ids, other_records = get_ids(filepath)
    records = get_records(ids)
    return records, other_records
