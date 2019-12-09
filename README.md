# EukHeist: Capturing environmental eukaryotic genomes
## Snakemake-based workflow for the retrieval of Eukaryotic (and other) genomes from metagenomes


# Introduction
Molecular and genomic approaches, particularly those applied to whole, mixed communities (e.g. metagenomics, metatranscriptomics), have shed light on the ecological roles, evolutionary histories, and physiological capabilities of these organisms. We developed a scalable and reproducible pipeline to facilitate the retrieval, taxonomic assignment, and annotation of eukaryotic metagenome assembled genomes (MAGs) from mixed community metagenomes. The below pipeline uses **EukHeist** to retrieve eukaryotic (and other) metagenome assembled genomes from the *Tara Oceans* global dataset.

_General workflow_
**Enter flow chart here**


# Setup
Inputs from user:   
* 2 input directories storing raw metatranscriptome and metagenome fastq files. If read pairs are nested by sample ID, that is OK.
```
EukHeist/raw_dir/
├── metagenome
└── metatranscriptome

```
* Provide 2 *sample data tables* (one for metaT and one for metaG) that list all sample IDs for input data. See example data table ```NAME```. _(optional)_
```
example sample data table
```
* Generate an *Assembly group table* which lists a unique name for each group of samples you wish to assemble together. Based on the provided example data set and data table(```NAME```), we've included scripts you can modify to generate your assembly group table. 
```
example assembly group table
HEADERS:
Assembly_group	Sample_list

*Sample-list must be comma separated with a space

```
* Modify ```config.yaml``` to tell EukHeist the location of raw read directories, assembly group table, and where you want results to be stored.
As an example, my input fastq files are located somewhere else, while the assembly group table is located in the input directory in this repo. 
```
copy config file
```

##outstanding questions/to-do:
- 


### Dependencies
Create a conda environment to run the EukHeist pipeline:

```conda env create --name EukHeist --file environment.yaml```

This conda environment runs the snakemake pipeline manager. Snakemake is then dependent on the environments available in ```EukHeist/envs/``` or through snakemake wrappers to run modules within the pipeline.


### Data
Download Tara Expedition metagenomic and metatranscriptomic data. [Use this pipeline](https://github.com/AlexanderLabWHOI/tara-download-snakemake/blob/master/Snakefile)
**This needs to get sorted**

### Configure working directories


## TO DO
* update so that input files can be R1_001.fastq.gz or 1.fastq.gz - or else? 
* Nested config.yaml file type?
* User needs to supply a SAMPLE LIST and an ASSEMBLY list
* Need to include a check at the beginning that (a) ensures all fastq files in SAMPLE LIST are present and not duplicated file names, (b) all samples listed in the ASSMEBLY GROUPS are also present in the SAMPLE LIST and not included more than once.
* Example input file with test data can have an accompanying scripts (maybe even in R and python!) that make a SAMPLE LIST from what is present in a given directory (see make-manifest-current.R) and a second one that sorts into assembly lists by parsing part of the file name???

