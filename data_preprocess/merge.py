# merge existing FASTA files (from databases) and remove duplicate sequences

import os
import disease


GLYCO_DIR = '../data/glycosylation'
S_NITRO_DIR = '../data/s_nitrosylation'
ACET_DIR = '../data/acetylation'


def main():
    dirs = [GLYCO_DIR, S_NITRO_DIR, ACET_DIR]
    for d in dirs:
        merged_file = os.path.join(d, "merged.fasta")
        if os.path.exists(merged_file):
            os.remove(merged_file)

        records = []
        new_records = []
        rec_without_disease = []

        for filename in os.listdir(d):
            if filename.endswith(".fasta") and filename != "merged.fasta" and filename != "seqs.fasta":
                filepath = os.path.join(d, filename)
                with_d, without_d = disease.main(filepath)
                rec_without_disease.append(without_d)
                for record in with_d:
                    new_records.append(record)
                all_records = records + new_records
                all_records = list(set(all_records))
                removed = (len(records) + len(new_records)) - len(all_records)
                records = all_records
                print(f"{filename}: {len(new_records)} sequences. {removed} were already there")
                new_records = []

        with open(merged_file, "w") as merged:
            for rec in records:
                merged.write(f">{rec[0]}|{rec[1]} {rec[2]}\n{rec[3]}\n")
            for recs in rec_without_disease:
                for r in recs:
                    merged.write(f">{r.description}\n{r.seq}\n")

        print("Merged succesfully\n")
        disease.vis_disease(d)


if __name__ == '__main__':
    main()
