# RUST - Ribo-seq Unit Seq Transform 

This program was developed by [Patrick O'Connor](https://pubmed.ncbi.nlm.nih.gov/?term=%22O%E2%80%99Connor+PBF%22%5BAuthor%5D) as part of his work in [LAPTI](http://lapti.ucc.ie) in University College Cork. 

This work resulted in the publication of [Comparative survey of the relative impact of mRNA features on local ribosome profiling read density](https://doi.org/10.1038/ncomms12915) in Nat Comms. 

## Abstract 
Ribosome profiling (Ribo-seq), a promising technology for exploring ribosome decoding rates, is characterized by the presence of infrequent high peaks in ribosome footprint density and by long alignment gaps. Here, to reduce the impact of data heterogeneity we introduce a simple normalization method, Ribo-seq Unit Step Transformation (RUST). RUST is robust and outperforms other normalization techniques in the presence of heterogeneous noise. We illustrate how RUST can be used for identifying mRNA sequence features that affect ribosome footprint densities globally. We show that a few parameters extracted with RUST are sufficient for predicting experimental densities with high accuracy. Importantly the application of RUST to 30 publicly available Ribo-seq data sets revealed a substantial variation in sequence determinants of ribosome footprint frequencies, questioning the reliability of Ribo-seq as an accurate representation of local ribosome densities without prior quality control. This emphasizes our incomplete understanding of how protocol parameters affect ribosome footprint densities.

# RUST v1.3

RUST version 1.3 only differs from v1.2 by version of python and some minor error handling 
The work in this repo primarily aims to bring RUST to python 3.10 and to improve the installation experience. 

## Installation

Clone this repository 
```
git clone https://github.com/JackCurragh/RUST.git 
```

Create a conda environment with python 3.10
```
conda create -n RUST python=3.10
conda activate RUST
```

Install required packages using pip
```
python -m pip install git+https://github.com/JackCurragh/RUST.git
```
