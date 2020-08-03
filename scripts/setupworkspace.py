import io 
import os
import pandas as pd
import numpy as np
import pathlib
import yaml
from snakemake.exceptions import print_exception, WorkflowError

with open("../hierarchy_cluster.yaml") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

#----SET VARIABLES----#

METAG_FOLDER = config["metaG"]["folder_name"]
METAT_FOLDER = config["metaT"]["folder_name"]
METAG_SAMPLES = pd.read_table(config["metaG"]["sample_data_table"])
METAT_SAMPLES = pd.read_table(config["metaT"]["sample_data_table"])
INPUTDIR = config["directories"]["input"]
ADAPTERS = config["adapters"]
METAG_STUDY = list(set(METAG_SAMPLES["SAMPLEID"].tolist()))
METAT_STUDY = list(set(METAT_SAMPLES["SAMPLEID"].tolist()))
SCRATCHDIR = config["directories"]["scratch"]
OUTPUTDIR = config["directories"]["output"]
METAG_SAMPLELIST = pd.read_table(config["metaG"]["assembly_group_table"], index_col="ASSEMBLY_GROUPING")
METAG_ASSEMBLYGROUP = list(METAG_SAMPLELIST.index)
METAT_SAMPLELIST = pd.read_table(config["metaT"]["assembly_group_table"], index_col="ASSEMBLY_GROUPING")
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

#----FUNCTIONS----#

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