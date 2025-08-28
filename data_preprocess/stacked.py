'''
Module to visualize cumulative counts of different PTMs over the years
using UniProt data.

Downloads TSV files for selected PTMs, aggregates them by year, and
plots both cumulative counts and percentage distributions.
'''

import os
import requests
import pandas as pd
import matplotlib.pyplot as plt


def download_ptm_data(ptms, urls, output_dir):
    '''
    Download PTM data from UniProt as TSV files.

    Args:
        ptms (list): List of PTM names.
        urls (list): Corresponding UniProt query URLs.
        output_dir (str): Directory where TSV files will be saved.
    '''
    os.makedirs(output_dir, exist_ok=True)

    for i, url in enumerate(urls):
        filepath = os.path.join(output_dir, f'{ptms[i]}.tsv')
        with requests.get(url, stream=True, timeout=10) as request:
            request.raise_for_status()
            with open(filepath, 'wb') as f:
                for chunk in request.iter_content(chunk_size=2**20):
                    f.write(chunk)


def aggregate_ptm_data(tsv_dir):
    '''
    Aggregate PTM entries per year from downloaded TSV files.

    Args:
        tsv_dir (str): Directory containing PTM TSV files.

    Returns:
        pd.DataFrame: DataFrame of cumulative counts per PTM per year.
        pd.DataFrame: DataFrame of cumulative percentages per PTM per year.
    '''
    dfs = []
    for filename in os.listdir(tsv_dir):
        if filename.endswith('.tsv'):
            filepath = os.path.join(tsv_dir, filename)
            temp_df = pd.read_csv(filepath, sep='\t')
            temp_df['ptm'] = filename.split('.')[0]
            dfs.append(temp_df)

    df_all = pd.concat(dfs, ignore_index=True)
    df_all['date_created'] = pd.to_datetime(df_all['Date of creation'])
    df_all['year'] = df_all['date_created'].dt.year

    entries_per_year_ptm = df_all.groupby(['year', 'ptm']).size().unstack(fill_value=0)
    cumulative_counts = entries_per_year_ptm.cumsum()
    cumulative_percentages = cumulative_counts.div(cumulative_counts.sum(axis=1), axis=0) * 100

    return cumulative_counts, cumulative_percentages


def plot_cumulative_counts(cumulative_counts, output_dir):
    '''
    Plot stacked bar chart of cumulative PTM counts.

    Args:
        cumulative_counts (pd.DataFrame): Cumulative counts per PTM per year.
        output_dir (str): Directory to save plots.
    '''
    os.makedirs(output_dir, exist_ok=True)

    cumulative_counts.plot(kind='bar', stacked=True, figsize=(10, 6))
    plt.ylabel('Cumulative number of entries')
    plt.title('Cumulative Entries per PTM by Year')
    plt.legend(title='PTM')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cum_counts.jpg'))
    plt.clf()


def plot_cumulative_percentages(cumulative_percentages, output_dir):
    '''
    Plot stacked bar chart of cumulative PTM percentages.

    Args:
        cumulative_percentages (pd.DataFrame): Cumulative percentages per PTM per year.
        output_dir (str): Directory to save plots.
    '''
    cumulative_percentages.plot(kind='bar', stacked=True, figsize=(10, 6))
    plt.ylabel('Cumulative % of entries')
    plt.title('Cumulative Distribution of PTMs per Year (Percent Stacked)')
    plt.legend(title='PTM', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cum_counts_percentage.jpg'))
    plt.clf()


def main():
    '''
    Download UniProt data for multiple PTMs, aggregate entries by year, 
    and generate cumulative count and percentage plots.

    Steps:
        1. Download TSV files for selected PTMs from UniProt.
        2. Aggregate entries per PTM per year.
        3. Plot stacked bar charts for cumulative counts.
        4. Plot stacked bar charts for cumulative percentages.

    Output:
        - '../data/all_ptm/cum_counts.jpg'
        - '../data/all_ptm/cum_counts_percentage.jpg'
    '''
    output_dir = '../data/all_ptm'

    ptms = ['phosphorylation', 'acetylation', 'methylation', 'hydroxylation', 'nitrosylation', 'glycosylation', 'ubl_conjugation', 'lipoprotein']

    urls = [
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0597%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0007%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0488%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0379%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0702%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0325%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0832%29%29',
        'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cdate_created&format=tsv&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29+AND+%28keyword%3AKW-0449%29%29'
    ]

    download_ptm_data(ptms, urls, output_dir)
    cumulative_counts, cumulative_percentages = aggregate_ptm_data(output_dir)
    plot_cumulative_counts(cumulative_counts, output_dir)
    plot_cumulative_percentages(cumulative_percentages, output_dir)


if __name__ == '__main__':
    main()
