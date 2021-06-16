import io 
import os
import pandas as pd
import numpy as np
import pathlib
import yaml
import subprocess
from snakemake.exceptions import print_exception, WorkflowError

with open("hierarchy_cluster.yaml") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    
def createSampleTable(tag, category, input_dir, folder):
    """
    Creates a sample table and assembly group data for the input files if the relevant flag is
    specified. Otherwise, the appropriate entries from the configuration file are returned, if 
    they exist.
    """
    if config[tag]["create_sample_table"]:
        create_call = ["python",os.path.join("scripts", "createinputfiles.py"),
                       "-d", os.path.join(INPUTDIR, folder),
                      "-c", category, "-o", os.path.join(input_dir, "tab" + "_" + tag)]
        if "samples" in config[tag]:
            create_call.append("-s")
            create_call.extend(config[tag]["samples"])
        if "assembly_groups" in config[tag]:
            create_call.append("-a")
            create_call.extend(config[tag]["assembly_groups"])
        subprocess.Popen(create_call).wait()

        return os.path.join(input_dir, "tab" + "_" + tag + "_sample_file.tsv"), \
               os.path.join(input_dir, "tab" + "_" + tag + "_assembly_file.tsv")
    else:
        try:
            sample_table = config[tag]["sample_data_table"]
            assembly_group_table = config[tag]["assembly_group_table"]
        except:
            print("You did not specify that the sample table should be created, nor did you " +
                  "include sample data tables and assembly group tables.")
            sys.exit(1)
        return sample_table, assembly_group_table

#----SET VARIABLES----#

#----STATIC REQUIRED VARIABLES----#
INPUTDIR = config["directories"]["input"]
ADAPTERS = config["adapters"]
SCRATCHDIR = config["directories"]["scratch"]
OUTPUTDIR = config["directories"]["output"]

METAG_FOLDER = config["metaG"]["folder_name"]
METAT_FOLDER = config["metaT"]["folder_name"]
METAG_SAMPLE_TABLE, METAG_ASSEMBLY_TABLE = createSampleTable("metaG", "METAGENOME", INPUTDIR, METAG_FOLDER)
METAT_SAMPLE_TABLE, METAT_ASSEMBLY_TABLE = createSampleTable("metaT", "METATRANSCRIPTOME", INPUTDIR, METAT_FOLDER)
    
METAG_SAMPLES = pd.read_table(METAG_SAMPLE_TABLE)
METAT_SAMPLES = pd.read_table(METAT_SAMPLE_TABLE)
METAG_STUDY = list(set(METAG_SAMPLES["SAMPLEID"].tolist()))
METAT_STUDY = list(set(METAT_SAMPLES["SAMPLEID"].tolist()))
METAG_SAMPLELIST = pd.read_table(METAG_ASSEMBLY_TABLE, index_col="ASSEMBLY_GROUPING")
METAG_ASSEMBLYGROUP = list(METAG_SAMPLELIST.index)
METAT_SAMPLELIST = pd.read_table(METAT_ASSEMBLY_TABLE, index_col="ASSEMBLY_GROUPING")
METAT_ASSEMBLYGROUP= list(METAT_SAMPLELIST.index)
ASSEMBLYGROUP = METAG_ASSEMBLYGROUP
USEFILE = False # whether to use the paths in the sample file to find the files if True, or use the metaT/metaG folder if False

#----COMPUTE VAR----#
MEGAHIT_CPU = config["megahit"]["cpu"]
MEGAHIT_MIN_CONTIG = config["megahit"]["min_contig"]
MEGAHIT_MEM = config["megahit"]["mem"]
MEGAHIT_OTHER = config["megahit"]["other"]

#----ACCESSION NUMS----#
metaG_run_accession = list(METAG_SAMPLES.SAMPLEID)
metaT_run_accession = list(METAT_SAMPLES.SAMPLEID)

#----FUNCTIONS FOR PIPELINE----#

def identify_read_groups(assembly_group_name, STUDY, FORWARD=True):
    ERR_list = str(METAG_SAMPLELIST.loc[assembly_group_name, 'SAMPLE_LIST']).split(',')
    ident_col = "R2"
    if FORWARD: 
        ident_col = "R1"
        
    if USEFILE:
        outlist = [os.path.join((METAG_SAMPLES.loc[(METAG_SAMPLES["SAMPLEID"] == curr_sample), 
                                                   "FULLPATH"]).to_string(index=False).strip(),
                   (METAG_SAMPLES.loc[METAG_SAMPLES["SAMPLEID"] == curr_sample, 
                                      ident_col]).to_string(index=False).strip()) for curr_sample in ERR_list]
    else:
        outlist = [os.path.join(INPUTDIR, METAG_FOLDER,
                   (METAG_SAMPLES.loc[METAG_SAMPLES["SAMPLEID"] == curr_sample, 
                                      ident_col]).to_string(index=False).strip()) for curr_sample in ERR_list]
        
    return(outlist)


def get_sample_list(assembly_group, SAMPLELIST, STUDY):
    outlist = []
    for assembly_group_name in assembly_group: 
        ERR_list = SAMPLELIST.loc[assembly_group_name, 'SAMPLE_LIST'].split(',')
        for E in ERR_list:
            out = SCRATCHDIR + "/mapping/{}/{}/{}.bam".format(assembly_group_name,STUDY, E) 
            outlist.append(out)
    return(outlist)


def get_sample_list_onegroup(assembly_group_name, SAMPLELIST, STUDY):
    outlist = []
    ERR_list = SAMPLELIST.loc[assembly_group_name, 'SAMPLE_LIST'].split(',')
    for E in ERR_list:
        out = SCRATCHDIR + "/mapping/{}/{}/{}.bam".format(assembly_group_name, STUDY, E)
        outlist.append(out)
    return(outlist)        

# Create output directory
pathlib.Path(OUTPUTDIR).mkdir(parents=True, exist_ok=True)