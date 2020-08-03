import io 
import os
import pandas as pd
import numpy as np
import pathlib
import yaml
import argparse
import math 

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sample_names', nargs = "+", default = [], dest = "sample_names",
                    help = "The samples to investigate.")
parser.add_argument('-a', '--assembly_groups', nargs = "+", default = [], dest = "assembly_groups",
                    help = "The assembly group corresponding to each sample.")
parser.add_argument('-d', '--sample_dir', type = str, default = "samples", dest = "sample_dir",
                    help = "The directory in which the sample files can be found.")
parser.add_argument('-c', '--category', type = str, default = "METAGENOME", dest = "category",
                    help = "Whether these samples are metagenome or metatranscriptome.")
parser.add_argument('-o', '--output', type = str, default = os.path.join("samples", "sample_tab"), dest = "output",
                    help = "The prefix directory for the output files.")

args = parser.parse_args()

sample_file = pd.DataFrame(columns = ["SAMPLEID", "FULLPATH", "R1", "R2", "OMIC", "ASSEMBLY_GROUPING"])
assembly_dict = dict()
sample_names = args.sample_names
assembly_groups = args.assembly_groups

available_files = os.listdir(args.sample_dir)
if len(sample_names) < 1: # if no sample names are provided
    shortest_name = [s.split("_")[0] for s in available_files]
    longer_name = ["_".join(s.split("_")[0:-1]) for s in available_files]
    if (len(set(shortest_name)) >= math.floor(len(available_files) / 2)):
        sample_names = shortest_name
    elif (len(set(shortest_name)) >= math.floor(len(available_files) / 2)):
        sample_names = longer_name
    else:
        sample_names = available_files
    assembly_groups = ["Group1"] * len(sample_names)
        
for ind in range(len(sample_names)):
    s = sample_names[ind]
    a = assembly_groups[ind] 
    
    matches = [curr for curr in available_files if s in curr]
    if len(matches) == 2:
        sorted_matches = sorted(matches, reverse = False)
        curr_frame = pd.DataFrame({"SAMPLEID": s, "FULLPATH": args.sample_dir, 
                      "R1": sorted_matches[0], "R2": sorted_matches[1],
                      "OMIC": args.category, "ASSEMBLY_GROUPING": a}, index = [0])
        if a in assembly_dict:
            assembly_dict[a].append(s)
        else:
            assembly_dict[a] = [s]
            
        sample_file = pd.concat([sample_file, curr_frame], ignore_index = True)
    

assembly_file = pd.DataFrame({"ASSEMBLY_GROUPING": list(assembly_dict.keys()), 
                              "SAMPLE_LIST": list(assembly_dict.values())})
   
os.system("mkdir -p " + os.path.dirname(args.output))
assembly_file.to_csv(args.output + "_assembly_file.tsv", sep = "\t")
sample_file.to_csv(args.output + "_sample_file.tsv", sep = "\t")