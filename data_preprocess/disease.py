# get disease information

import csv
from io import StringIO
from collections import Counter
import re
import requests
from Bio import SeqIO
import matplotlib.pyplot as plt


def get_ids(filepath):
    uniprot_ids = []
    other_records = []
    count = 0
    for record in SeqIO.parse(filepath, "fasta"):
        if record.id.startswith("sp|") or record.id.startswith("tr|"):
            name = record.id.split("|")[1]
            if name[-2] == '.':
                name = name[:-2]
            uniprot_ids.append(name)
        else:
            count += 1
            other_records.append(record)
    print(f"{count} sequences not in Uniprot")
    return uniprot_ids, other_records


def get_records(uniprot_ids):
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
                name = row.get("Entry Name", "")
                disease = row.get("Involvement in disease", "")
                mim_ids = re.findall(r'\[MIM:(\d+)', disease)
                mim_string = ""
                for m in mim_ids:
                    mim_string += f" MIM:{m}"
                seq = row.get("Sequence", "")
                record = [acc,name,mim_string,seq]
                records.add(tuple(record))
        else:
            print("(Disease) Batch download failed:", response.status_code)

    return records


def vis_disease(dir_path):
    diseases = []
    for record in SeqIO.parse(f"{dir_path}/merged.fasta", "fasta"):
        if "MIM:" in record.description:
            parts = record.description.split(" MIM:")
            for disease_part in parts[1:]:
                diseases.append(disease_part)

    counts = Counter(diseases)
    top10 = counts.most_common(10)
    labels = [d for d, c in top10]
    values = [c for d, c in top10]
    plt.bar(labels, values, color='skyblue', edgecolor='black')
    plt.xlabel('MIM ID of Diseases')
    plt.ylabel('Absolute Frequency')
    plt.xticks(rotation=45, ha='right')
    plt.title(f'Absolute Frequency of Top 10 Diseases (PTM: {dir_path.split("/")[-1]})')
    plt.tight_layout()
    plt.savefig(f"{dir_path}/top10diseases.jpg")
    plt.clf()


def main(filepath):
    ids, other_records =  get_ids(filepath)
    records = get_records(ids)
    return records, other_records


if __name__ == '__main__':
    main()
