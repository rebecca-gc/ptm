from Bio import SeqIO

file1 = "path"
file_target = "path"

seqs1 = {str(record.seq) for record in SeqIO.parse(file1, "fasta")}
seqs_target = {str(record.seq) for record in SeqIO.parse(file_target, "fasta")}

common_seqs = seqs1 & seqs_target

print(f"Count of common sequences: {len(common_seqs)}")

with open("filtered.fasta", "w") as filtered:
    for record in SeqIO.parse(file_target, "fasta"):
        if record.seq not in common_seqs:
            SeqIO.write(record, filtered, "fasta")

print("Removed common sequences and saved as filtered.fasta")