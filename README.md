# Project Title
### Classfication of Human Protein Post-Translational Modifications encoded with iCAN using a Random Forest Classifier

## Table of Contents

- [Authors](https://github.com/rebecca-gc/ptm#authors)
- [Dependencies](https://github.com/rebecca-gc/ptm#dependencies)
- [Data](https://github.com/rebecca-gc/ptm#data)
- [Code](https://github.com/rebecca-gc/ptm#code)

## Authors

- [Georges Hattab](https://github.com/ghattab)
- [Aleksandar Anžel](https://github.com/AAnzel)
- [Rebecca Grevens Carpi](https://github.com/rebecca-gc)

Created by the [Visualization group](https://visualization.group/), which is part of the Centre for Artifical Intelligence in Public Health Research (ZKI-PH) of the Robert Koch-Institute. This repository was developed as part of a bachelor thesis.

## Dependencies
The code is written in Python 3.11.13 with the following libraries installed:

|Library|Version|
|---|---|
|biopython|1.85|
|bs4|4.13.4|
|cd-hit|4.8.1|
|joblib|1.5.1|
|NumPy7|2.3.0|
|matplotlib|3.10.2|
|pandas|2.3.0|
|requests|2.32.4|
|scikit-learn|1.7.0|
|venn|0.1.3|
|wandb|0.21.0|

### Conda Environment
It is recommended to use Conda to build the environment with all dependencies. This process requires users to install Anaconda or Miniconda package managers. Users can install Conda using the instructions on the following link [https://docs.anaconda.com/miniconda/#quick-command-line-install](https://docs.anaconda.com/miniconda/#quick-command-line-install). If Conda is installed, the following sequence of instructions (for Linux-based systems) installs all dependencies and activates the environment:

```bash
conda env create -f environment.yml
conda activate rebecca_ptm_env
```

## Data
The protein sequences used to predict the Post-Translational Modifications are collected from seven different databases.
- SwissProt
- NCBI
- dbPTM
- ptmd
- PTMCode2
- qPTM (Note: requires data access inquiry and has to be downloaded manually)
- Unipep

These databases were integrated, filtered, and processed into positive and negative datasets for machine learning classification. Preprocessing steps included merging, deduplication, false negative removal, clustering, and sequence-length filtering. Full details are described in the thesis.

## Code
The repository is organized into two main directories:

|Script|Description|
|---|---|
|[data_preprocess/](./data_preprocess/)|Scripts for retrieving, cleaning, and preparing the datasets.|
|[data_preprocess/data_pipeline.py](./data_preprocess/data_pipeline.py)|Orchestrates the entire preprocessing workflow: merging, deduplication, filtering, clustering, and length cutoff. Produces final positive/negative multi-FASTA files and labels.|
|[data_preprocess/download_all.py](./data_preprocess/download_all.py)|Downloads protein sequences in FASTA format from all databases.|
|[data_preprocess/merge.py](./data_preprocess/merge.py)|Merges multi-FASTA files from multiple databases into unified PTM-specific files, ensuring consistency and removing duplicates.|
|[data_preprocess/disease.py](./data_preprocess/class_generator.py)|Extracts disease associations from UniProt annotations and links them to protein sequences.|
|[data_preprocess/negatives.py](./data_preprocess/negatives.py)|Generates negative datasets by filtering out sequences that contain the PTM of interest while allowing other PTMs.|
|[data_preprocess/cluster.py](./data_preprocess/cluster.py)|Runs CD-HIT to cluster sequences at 40% similarity, retaining one representative per cluster to reduce redundancy and avoid data leakage.|
|[data_preprocess/class_generator.py](./data_preprocess/class_generator.py)|Generates the classes.txt label file and a consolidated seqs.fasta containing both positive and negative sequences for machine learning input.|
|[data_preprocess/venn_diagrams.py](./data_preprocess/venn_diagrams.py)|General visualization script; generates Venn diagrams of sequence overlaps across PTMs and diseases.|
|[data_preprocess/stacked.py](./data_preprocess/stacked.py)|Creates stacked bar plots showing the distribution of 8 main PTMs in UniProt over time.|

|Script|Description|
|---|---|
|[rfc/](./rfc/)|Scripts related to feature encoding and classification with Random Forests.|
|[rfc/encoding_pipeline.py](./rfc/encoding_pipeline.py)|Converts sequences into machine learning features via iCAN for Random Forest training.|
|[rfc/rfc_with_cv.py](./rfc/rfc_with_cv.py)|Trains and evaluates the Random Forest classifier with cross-validation. Reports performance metrics.|

## Running

The code should be run from the project's root directory. For processing the data users should run:

```bash
python data_preprocess/data_pipeline.py
```
