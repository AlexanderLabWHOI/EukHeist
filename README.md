# EukHeist: Capturing environmental eukaryotic genomes
## Snakemake-based workflow for the retrieval of Eukaryotic (and other) genomes from metagenomes


## 1.0 Introduction
Molecular and genomic approaches, particularly those applied to whole, mixed communities (e.g. metagenomics, metatranscriptomics), have shed light on the ecological roles, evolutionary histories, and physiological capabilities of these organisms. We developed a scalable and reproducible pipeline to facilitate the retrieval, taxonomic assignment, and annotation of eukaryotic metagenome assembled genomes (MAGs) from mixed community metagenomes. The below pipeline uses **EukHeist** to retrieve eukaryotic (and other) metagenome assembled genomes from the *Tara Oceans* global dataset.

_General workflow_

![flowchart](input/flowchart1.png)


_User workflow_
* get set up on HPC with snakemake
* input sample list, assembly group
* Run tests and confirm rules in pipeline you plan to use


## 2.0 Setup

Locate and organize metagenomic and metatranscriptomic fastq files. You will need to know the full path for all files, individual sample IDs, and an idea of how the assemblies should be grouped. In this step, we are creating input file lists that will tell EukHeist where to look for our input fastq reads, what to name them, and how to group the assemblies or mapping.


### 2.1 Structure input fastq files in individual directories

Example file structure:
```
 EukHeist/raw_dir/
 ├── metagenome
 └── metatranscriptome
```

To use the test data, run:
```bash download-test-data.sh```

This will download raw reads, for a test run, to the metagenome and metatranscriptome directories.

### 2.2 Create _sample list_ files for metagenome and metatranscriptome data

The final sample list files should look like this, with complete paths under "FULLPATH" and the assembly grouping lists how you want the samples to be assembled for the metagenomic and metatranscriptomic pipelines.

**Example sample list files**

```
SAMPLEID        SAMPLENAME      OMIC    FULLPATH        R1      R2      ASSEMBLY_GROUPING
ERR1163068      CTD1200_DNA     METAGENOMIC     /../../../ERR1163068 ERR1163068_1.fastq.gz   ERR1163068_2.fastq.gz   CTD1200_DNA
ERR1163069      CTD1200_DNA     METAGENOMIC     /../../../ERR1163069 ERR1163069_1.fastq.gz   ERR1163069_2.fastq.gz   CTD1200_DNA
ERR1163070      CTD1200_DNA     METAGENOMIC     /../../../ERR1163070 ERR1163070_1.fastq.gz   ERR1163070_2.fastq.gz   CTD1200_DNA
```

See example file "input-samplelist-example-metagenomic.txt". These should be _.txt_ files and separated with a tab.



### 2.3 Generate an *Assembly group* file

The first column, 'ASSEMBLY_GROUPING' lists the unique names for each group specified from the sample list file. The second column lists the sample IDs (in this example, the accession numbers), that are associated with the ASSEMBLY_GROUPING. These sample IDs are collapsed with commas.

```
ASSEMBLY_GROUPING       SAMPLE_LIST
Group1  ERR1726828,ERR599214
Group2  ERR868421
```


Create sample list files that end in ```metatranscriptomic.txt``` and metagenomic.txt```, and run this R script to automatically generate the Assembly group files. 

Run Rscript ```/../generate-assembly-group-tables.r```.



## 3.0 Set up EukHeist

Start by cloning this repo.
```
git clone https://github.com/AlexanderLabWHOI/EukHeist.git
```


Then, create a conda environment to run the EukHeist pipeline:

```conda env create --name EukHeist --file environment.yaml```

This conda environment runs the snakemake pipeline manager. Snakemake is then dependent on the environments available in ```EukHeist/envs/``` or through snakemake wrappers to run modules within the pipeline.   


### 3.1

Modify ```config.yaml``` to tell EukHeist the location of raw read directories, assembly group table, and where you want results to be stored.

As an example, my input fastq files are located somewhere else, while the assembly group table is located in the input directory in this repo. 
```
inputDIR: /vortexfs1/omics/alexander/data/TARA	#Location and full path to raw sequences

scratch: 	#Location and full path to output scratch directory

metaG_ena_table: input/ENA_tables/PRJEB4352_metaG_wenv_PE-TEST.txt	## Update with new name for these files?
metaT_ena_table: input/ENA_tables/PRJEB6609_metaT_wenv_PE-TEST.txt
metaG_sample_list: input/SampleList_ForAssembly_metaG_python-TEST.txt	#Location of assembly file
metaT_sample_list: input/SampleList_ForAssembly_metaT_python-TEST.txt 
```

Explanation of ```config.yaml``` file:
_need to figure out if the file structure can be changed upfront?_
```
metaG_accession: PRJEB4352 
metaT_accession: PRJEB6609

metaG_ena_table: input/ENA_tables/PRJEB4352_metaG_wenv_PE-TEST.txt
metaT_ena_table: input/ENA_tables/PRJEB6609_metaT_wenv_PE-TEST.txt

#metaG_ena_table: input/ENA_tables/PRJEB4352_metaG_wenv_PE.txt
#metaT_ena_table: input/ENA_tables/PRJEB6609_metaT_wenv_PE.txt

inputDIR: /vortexfs1/omics/alexander/data/TARA
outputDIR: /vortexfs1/omics/alexander/akrinos/output-data-TARA

#outputDIR: /vortexfs1/omics/alexander/data/TARA/PRJEB4352-snakmake-output 

scratch:  /vortexfs1/scratch/akrinos/tara
adapters: input/adapters/illumina-adapters.fa
metaG_sample_list: input/SampleList_ForAssembly_metaG_python-TEST.txt
metaT_sample_list: input/SampleList_ForAssembly_metaT_python-TEST.txt 

#metaG_sample_list: input/SampleList_ForAssembly_metaG_python.txt
#metaT_sample_list: input/SampleList_ForAssembly_metaT_python.txt 

megahit_other: --continue --k-list 29,39,59,79,99,119 
megahit_cpu: 80
megahit_min_contig: 1000
megahit_mem: .95
restart-times: 0
max-jobs-per-second: 1
max-status-checks-per-second: 10
local-cores: 1
rerun-incomplete: true
keep-going: true
```


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
