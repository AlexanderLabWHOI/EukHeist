# EukHeist: Capturing environmental eukaryotic genomes
## Snakemake-based workflow for the retrieval of Eukaryotic (and other) genomes from metagenomes


# Introduction
Molecular and genomic approaches, particularly those applied to whole, mixed communities (e.g. metagenomics, metatranscriptomics), have shed light on the ecological roles, evolutionary histories, and physiological capabilities of these organisms. We developed a scalable and reproducible pipeline to facilitate the retrieval, taxonomic assignment, and annotation of eukaryotic metagenome assembled genomes (MAGs) from mixed community metagenomes. The below pipeline uses **EukHeist** to retrieve eukaryotic (and other) metagenome assembled genomes from the *Tara Oceans* global dataset.

_General workflow_

![flowchart](input/flowchart1.png)

# Setup
Inputs from user:   
## 1.  2 input directories storing raw metatranscriptome and metagenome fastq files.
```
 EukHeist/raw_dir/
 ├── metagenome
 └── metatranscriptome
```
To use the test data, run:
```
bash download-test-data.sh
```
This is download raw reads to the metagenome and metatranscriptome directories.

## 2. Provide 2 *sample data tables* (one for metaT and one for metaG) that list all sample IDs for input data. Using the test data, enable an R environment and run the provided R script in the metagenome and metatranscriptome directories (provided).
The final product will look like this:
```
SAMPLEID        FULLPATH        R1      R2      OMIC    ASSEMBLY_GROUPING
ERR1711907      /vortexfs1/omics/alexander/shu/EukHeist/raw_dir/metatranscriptome/ERR1711907    ERR1711907_1.fastq.gz   ERR1711907_2.fastq.gz   METATRANSCRIPTOME       
ERR1719262      /vortexfs1/omics/alexander/shu/EukHeist/raw_dir/metatranscriptome/ERR1719262    ERR1719262_1.fastq.gz   ERR1719262_2.fastq.gz   METATRANSCRIPTOME       
ERR1740135      /vortexfs1/omics/alexander/shu/EukHeist/raw_dir/metatranscriptome/ERR1740135    ERR1740135_1.fastq.gz   ERR1740135_2.fastq.gz   METATRANSCRIPTOME  
```
To automatically generate this using the contents in ```metagenome``` and ```metatranscriptome```, see Rscripts that find the directory's fastq files and create the above table.
```
cd ./EukHeist/raw_dir/metatranscriptome/

# Enable an R environment to run R
conda activate r_3.5.1
Rscript generate-metaT-samplelist.r

# Output sample list: 'samplelist-metaT.txt'
## Repeat for metagenome in /EukHeist/raw_dir/metagenome/
### Output from metagenome script is 'samplelist-metaG.txt'
```
*The last column* of the sample list table _ASSEMBLY_GROUP_ lists how you want the samples to be assembled for the metagenomic and metatranscriptomic pipelines. See final versions of this in ```EukHeist/input/samplelist-meta*-wgroups.txt```

## 3. Generate an *Assembly group table* which lists a unique name for each group of samples you wish to assemble together. Based on the provided example data set and data table(```NAME```), we've included scripts you can modify to generate your assembly group table. You can also use the provided script to generate this file for you, based on the _ASSEMBLY_GROUP_ column from the ```samplelist-metaG-wgroups.txt``` and ```samplelist-metaT-wgroups.txt``` files.
```
ASSEMBLY_GROUPING       SAMPLE_LIST
Group1  ERR1726828,ERR599214
Group2  ERR868421
```

To generate automatically, run Rscript ```/EukHeist/input/generate-assembly-group-tables.r```.
```
cd /EukHeist/input/
# Ensure metaG and metaT files are present, named:
## samplelist-metaG-wgroups.txt
## samplelist-metaT-wgroups.txt

Rscript generate-assembly-group-tables.r
# Output file: assembly-list-metaG.txt and assembly-list-metaT.txt

```

## 4. Modify ```config.yaml``` to tell EukHeist the location of raw read directories, assembly group table, and where you want results to be stored.
As an example, my input fastq files are located somewhere else, while the assembly group table is located in the input directory in this repo. 
```
inputDIR: /vortexfs1/omics/alexander/data/TARA	#Location and full path to raw sequences

scratch: 	#Location and full path to output scratch directory

metaG_ena_table: input/ENA_tables/PRJEB4352_metaG_wenv_PE-TEST.txt	## Update with new name for these files?
metaT_ena_table: input/ENA_tables/PRJEB6609_metaT_wenv_PE-TEST.txt
metaG_sample_list: input/SampleList_ForAssembly_metaG_python-TEST.txt	#Location of assembly file
metaT_sample_list: input/SampleList_ForAssembly_metaT_python-TEST.txt 

```

## 5. Create a conda environment to run the EukHeist pipeline:

```conda env create --name EukHeist --file environment.yaml```

This conda environment runs the snakemake pipeline manager. Snakemake is then dependent on the environments available in ```EukHeist/envs/``` or through snakemake wrappers to run modules within the pipeline.   


Explanation of working directory:

```
EukHeist
├── cluster.yaml         # Specifications to submit Snakemake jobs through SLURM on HPC
├── config-test.yaml     # Test configuration file
├── config.yaml          # Configuration file, modified above. Edit to customize Snakemake pipeline
├── environmentv0.2.yaml # Starting environment to load Snakemake
├── environment.yaml     # Old environment file - to delete
├── envs		 # All conda environments required for snakemake rules
│   
├── input		 # Required input files for snakefile to run and scripts
│   ├── adapters
│   ├── assembly-list-metaG.txt
│   ├── assembly-list-metaT.txt
│   ├── generate-assembly-group-tables.r
│   ├── samplelist-metaG-wgroups.txt
│   └── samplelist-metaT-wgroups.txt
├── raw_dir			
│   ├── metagenome
│   └── metatranscriptome
├── README.md
├── rules
│   └── normalization.smk
├── Snakefile		 # Complete snakemake file to run full pipeline
├── start-up
└── submit_script	 # Submit scripts to enable running on HPC / with slurm.
```

## 6. Data
Download Tara Expedition metagenomic and metatranscriptomic data. [Use this pipeline](https://github.com/AlexanderLabWHOI/tara-download-snakemake/blob/master/Snakefile)
**This needs to get sorted**  

Also instructions for test data?


## 7. Other considerations - program specifications


## 8. Test run snakemake

```
snakemake -np --use-conda

# Or if running with HPC
bash submit_script/dry_submit_snakemake.sh
```

## 9. Execute full snakemake pipeline

```
snakemake --use-conda

# Or if running with HPC 
bash submit_script/submit_snakemake.sh
```


### Troubleshooting snakemake
When throwing an error, snakemake will list log files. Each time snakemake is executed, a log file is created in ```CURRENT_SNAKEMAKE_DIR/.snakemake/log/```. These are dated and provide the printed output. Some common errors and steps to diagnose.   

*Compatibility with snakemake and conda* Note the version of snakemake and anaconda you are using. Upon conda's switch from _source activate_ to _conda activate_, snakemake was not calling on the conda environments properly. Errors associted with these were ```returned non-zero exit status 127``` and an error about *line 56* in *thread.py* like this: ```LOCATION_OF_YOUR_CONDA_ENVS/.conda/envs/snake-18S/lib/python3.6/concurrent/futures/thread.py", line 56, in run```
Update your version of snakemake. Versions listed above are compatible. This error will also be generated when there is an incomptaible conda environment called, such as an outdated [snakemake wrapper](https://snakemake-wrappers.readthedocs.io/en/stable/). Try updating that - but ensure you have first updated the snakemake version.   

Check all log files, found in ```./snakemake/log``` and the log files generated as output. Most helpful are the output files from each slurm job, ```slurm-XXX.out```. Look at the most recent slurm out files after a job fails to find more informative errors.    

See ```snakemake -h``` for additional commands to clean up working snakemake environment, list steps, or restarting attempts, etc. When running, snakemake "locks" the working directory, so only one snakemake command can be run at a time. If you have to cancel, make sure to run ```snakemake --unlock``` to clear it. See other flags to clean up your working environment after updated conda, snakemake, or other environment versions (```--cleanup-shadow```, ```--cleanup-conda```).
To look for additional code error that may result in Syntax errors, try adding this to snakemake execution:
* ```--summary``` or ```--detailed-summary```
* ```--printshellcmds```
* ```--debug```

## TO DO
* update so that input files can be R1_001.fastq.gz or 1.fastq.gz - or else? 
* Nested config.yaml file type?
* User needs to supply a SAMPLE LIST and an ASSEMBLY list - rather than the assessions that assume metaG vs metaT
* Need to include a check at the beginning that (a) ensures all fastq files in SAMPLE LIST are present and not duplicated file names, (b) all samples listed in the ASSMEBLY GROUPS are also present in the SAMPLE LIST and not included more than once.
* snakemake checkpoints?
