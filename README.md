# Project Title
### Machine Learning Data for Classification using iCAN for Post Translational Modifications in Human Proteins

## Table of Contents

- [Authors](https://github.com/rebecca-gc/ptm#authors)
- [Dependencies](https://github.com/rebecca-gc/ptm#dependencies)
- [Data](https://github.com/rebecca-gc/ptm#data)
- [Code](https://github.com/rebecca-gc/ptm#code)

## Authors

- [Aleksandar Anžel](https://github.com/AAnzel)
- [Georges Hattab](https://github.com/ghattab)
- [Rebecca Grevens](https://github.com/rebecca-gc)

Created by the [Visualisation group](https://visualization.group/), which is part of the Centre for Artifical Intelligence in Public Health Research (ZKI-PH) of the Robert Koch-Institute.

## Dependencies
The code is written in Python 3.11.12 with the following libraries installed:

|Library|Version|
|---|---|
|biopython|1.78|

## Data
The example datasets which are used to predict the Post Translational Modifications are collected from the [UniProt Database](https://www.uniprot.org/) with the filter (reviewed:true) to obtain only experimentally verified proteins.

## Code
|Script|Description|
|---|---|
|[data_preprocess/](./data_preprocess/)|contains all scripts used for dataset retrieval and modification.|
|[data_preprocess/download.py](./data_preprocess/download.py)|contains the code for downloading protein sequences as FASTA files from different databases|
|[data_preprocess/remove_common_seqs.py](./data_preprocess/remove_common_seqs.py)|contains the code for filtering out false negatives.|
|[data_preprocess/generator_x_balanced.py](./data_preprocess/generator_x_balanced.py)|contains the code to generate a balanced classes.txt and seqs.fasta for two given FASTA files.|
