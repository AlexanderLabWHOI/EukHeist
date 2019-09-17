# EukHeist: Capturing environmental eukaryotic genomes
## Snakemake-based workflow for the retrieval of Eukaryotic (and other) genomes from metagenomes


# Introduction
Molecular and genomic approaches, particularly those applied to whole, mixed communities (e.g. metagenomics, metatranscriptomics), have shed light on the ecological roles, evolutionary histories, and physiological capabilities of these organisms. We developed a scalable and reproducible pipeline to facilitate the retrieval, taxonomic assignment, and annotation of eukaryotic metagenome assembled genomes (MAGs) from mixed community metagenomes. The below pipeline uses **EukHeist** to retrieve eukaryotic (and other) metagenome assembled genomes from the *Tara Oceans* global dataset.

_General workflow_
**Enter flow chart here**


# Setup
***
### Dependencies
Create a conda environment to run the EukHeist pipeline:

```conda env create --name EukHeist --file environment.yaml```

This conda environment runs the snakemake pipeline manager. Snakemake is then dependent on the environments available in ```EukHeist/envs/``` or through snakemake wrappers to run modules within the pipeline.


### Data
Download Tara Expedition metagenomic and metatranscriptomic data. [Use this pipeline](https://github.com/AlexanderLabWHOI/tara-download-snakemake/blob/master/Snakefile)
**This needs to get sorted**

