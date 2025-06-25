import os
import requests
import pandas as pd
import matplotlib.pyplot as plt


d = '../data/all_ptm'

ptms = ['phosphorylation', 'acetylation', 'methylation', 'hydroxylation', 'nitrosylation', 'glycosylation', 'ubl_conjugation', 'lipoprotein']

urls = ['https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0597%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0007%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0488%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0379%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0702%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0325%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0832%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0449%29%29']

# for i, url in enumerate(urls):
#     filepath = os.path.join(d, f'{ptms[i]}.tsv')
#     with requests.get(url, stream=True, timeout=10) as request:
#         request.raise_for_status()
#         with open(filepath, 'wb') as f:
#             for chunk in request.iter_content(chunk_size=2**20):
#                 f.write(chunk)

dfs = []
for filename in os.listdir(d):
    filepath = os.path.join(d, filename)
    temp_df = pd.read_csv(filepath, sep='\t')
    temp_df['ptm'] = filename.split('.')[0]
    dfs.append(temp_df)

df_all = pd.concat(dfs, ignore_index=True)

df_all['date_created'] = pd.to_datetime(df_all['Date of creation'])
df_all['year'] = df_all['date_created'].dt.year

entries_per_year_ptm = df_all.groupby(['year', 'ptm']).size().unstack(fill_value=0)

cumulative_counts = entries_per_year_ptm.cumsum()

cumulative_counts.plot(kind='bar', stacked=True, figsize=(10, 6))
plt.ylabel('Cumulative number of entries')
plt.title('Cumulative Entries per PTM by Year')
plt.legend(title='PTM')
plt.tight_layout()
plt.savefig('../data/all_ptm/cum_counts.jpg')
plt.clf()

cumulative_percentages = cumulative_counts.div(cumulative_counts.sum(axis=1), axis=0) * 100

cumulative_percentages.plot(kind='bar', stacked=True, figsize=(10, 6))
plt.ylabel('Cumulative % of entries')
plt.title('Cumulative Distribution of PTMs per Year (Percent Stacked)')
plt.legend(title='PTM', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('../data/all_ptm/cum_counts_percentage.jpg')
plt.clf()
