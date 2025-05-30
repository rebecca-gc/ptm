# get disease information

import os
import requests
from Bio import SeqIO
import csv
from io import StringIO
from collections import Counter
import matplotlib.pyplot as plt


GLYCO_DIR = '../data/glycosylation'
S_NITRO_DIR = '../data/s_nitrosylation'
ACET_DIR = '../data/acetylation'


def get_ids(filepath):
    uniprot_ids = []
    count = 0
    for record in SeqIO.parse(filepath, "fasta"):
        if record.id.startswith("sp|") or record.id.startswith("tr|"):
            uniprot_ids.append(record.id.split("|")[1])
        else:
            count += 1
    print(f"{count} sequences not in Uniprot")
    return uniprot_ids


def get_records(uniprot_ids):
    filepath = "./disease.fasta"
    if os.path.exists(filepath):
        os.remove(filepath)

    url = "https://rest.uniprot.org/uniprotkb/stream"
    chunk_size = 500
    records = set()
    for i in range(0, len(uniprot_ids), chunk_size):
        chunk = uniprot_ids[i:(i + chunk_size)]
        query = " OR ".join([f"accession:{uid}" for uid in chunk])
        fields = "accession,id,cc_disease,sequence"
        url = f"https://rest.uniprot.org/uniprotkb/stream?format=tsv&query={query}&fields={fields}"
        
        response = requests.get(url, timeout=30)
        if response.ok:
            tsv_data = StringIO(response.text)
            reader = csv.DictReader(tsv_data, delimiter='\t')
            for row in reader:
                acc = row.get("Entry", "")
                id = row.get("Entry Name", "")
                disease = row.get("Involvement in disease", "")
                seq = row.get("Sequence", "")
                stuff = [acc,id,disease,seq]
                records.add(tuple(stuff))
                """with open(filepath, "a") as f:
                    f.write(f">sp|{acc}|{id} {disease}\n{seq}\n")"""
        else:
            print("(Disease) Batch download failed:", response.status_code)

    return records


def vis_disease(filepath):
    diseases = []
    for record in SeqIO.parse(f"{filepath}/merged_d.fasta", "fasta"):
        if "DISEASE" in record.description:
            parts = record.description.split("DISEASE: ")
            for disease_part in parts[1:]:
                diseases.append(disease_part)

    counts = Counter(diseases)
    top10 = counts.most_common(10)
    labels = [d[:20] for d, c in top10]
    values = [c for d, c in top10]
    plt.bar(labels, values, color='skyblue', edgecolor='black')
    plt.xlabel('Diseases')
    plt.ylabel('Absolute Frequency')
    plt.xticks(rotation=45, ha='right')
    plt.title(f'Absolute Frequency of Top 10 Diseases (PTM: {filepath.split("/")[-1]})')
    plt.tight_layout()
    plt.savefig(f"{filepath}/top10diseases.jpg")


def main(filepath):
    ids =  get_ids(filepath)
    records = get_records(ids)
    return records


if __name__ == '__main__':
    main()
