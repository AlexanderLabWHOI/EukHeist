# EukHeist: Capturing environmental eukaryotic genomes
## Snakemake-based workflow for the retrieval of Eukaryotic (and other) genomes from metagenomes


# Introduction
Molecular and genomic approaches, particularly those applied to whole, mixed communities (e.g. metagenomics, metatranscriptomics), have shed light on the ecological roles, evolutionary histories, and physiological capabilities of these organisms. We developed a scalable and reproducible pipeline to facilitate the retrieval, taxonomic assignment, and annotation of eukaryotic metagenome assembled genomes (MAGs) from mixed community metagenomes. The below pipeline uses **EukHeist** to retrieve eukaryotic (and other) metagenome assembled genomes from the *Tara Oceans* global dataset.

_General workflow_
**Enter flow chart here**


# Setup
Inputs from user:   
1.  2 input directories storing raw metatranscriptome and metagenome fastq files. If read pairs are nested by sample ID, that is OK.
```
 EukHeist/raw_dir/
 ├── metagenome
 └── metatranscriptome
```
Using the test data:
```
 EukHeist/raw_dir/
 │── metagenome
   ├── ERR1726828  #This is an example of nested
   │   ├── ERR1726828_1.fastq.gz
   │   ├── ERR1726828_2.fastq.gz
   │   └── md5sum.tab
   ├── ERR599214
   │   ├── ERR599214_1.fastq.gz
   │   ├── ERR599214_2.fastq.gz
   │   └── md5sum.tab
   └── ERR868421
       ├── ERR868421_1.fastq.gz
       ├── ERR868421_2.fastq.gz
       └── md5sum.tab
```

2. Based on the above data sets, provide 2 *sample data tables* (one for metaT and one for metaG) that list all sample IDs for input data. See example data table ```NAME```. _(optional)_
```
SAMPLEID        FULLPATH        R1      R2      OMIC    ASSEMBLY_GROUPING
ERR1711907      /vortexfs1/omics/alexander/shu/EukHeist/raw_dir/metatranscriptome/ERR1711907    ERR1711907_1.fastq.gz   ERR1711907_2.fastq.gz   METATRANSCRIPTOME       
ERR1719262      /vortexfs1/omics/alexander/shu/EukHeist/raw_dir/metatranscriptome/ERR1719262    ERR1719262_1.fastq.gz   ERR1719262_2.fastq.gz   METATRANSCRIPTOME       
ERR1740135      /vortexfs1/omics/alexander/shu/EukHeist/raw_dir/metatranscriptome/ERR1740135    ERR1740135_1.fastq.gz   ERR1740135_2.fastq.gz   METATRANSCRIPTOME  
```
To automatically generate this using the contents in ```metagenome``` and ```metatranscriptome```, see Rscripts that find the directory's fastq files and create the above table.
```
# Example R script:
cd ./EukHeist/raw_dir/metatranscriptome/

# Enable an R environment to run R
conda activate r_3.5.1
Rscript generate-metaT-samplelist.r

# Output sample list: 'samplelist-metaT.txt'
## Repeat for metagenome

```
*The last column* of the sample list table

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

