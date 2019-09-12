#! /usr/bin/env python
"""
Split tables to identify PE vs SE reads in PRJEB files supplied from ENA for download through the snakemake pipeline. 
"""
import sys
import argparse
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument('projectfile')
    p.add_argument('--column',default='fastq_ftp', required=False)
    p.add_argument('--delimiter', default=';', required=False)
    args = p.parse_args() 
    #read in file and split the ftp column
    df = pd.read_table(args.projectfile) 
    df[['c1', 'c2']]=df[args.column].str.split(args.delimiter, expand=True)
    #test if PE read by looking for fastq2 file
    pe_ind = df.c2.notnull()    
    se_ind = df.c2.isnull()
    basename = args.projectfile.split('.')[0]
    PE_out = basename + '_PE.txt'
    SE_out = basename + '_SE.txt'
    PE_READS = df.loc[pe_ind].drop(['c1', 'c2'], axis=1)
    SE_READS = df.loc[se_ind].drop(['c1', 'c2'], axis=1)  
    PE_READS.to_csv(PE_out, sep='\t') 
    SE_READS.to_csv(SE_out, sep='\t')
    
if __name__ == '__main__':
        sys.exit(main())
