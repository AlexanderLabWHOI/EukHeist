configfile: "config.yaml"  

import io 
import os
import pandas as pd
import numpy as np
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

#----SET VARIABLES----#
METAG_ACCESSION = config["metaG_accession"]
METAT_ACCESSION = config["metaT_accession"]
METAG_SAMPLES = pd.read_table(config["metaG_ena_table"])
METAT_SAMPLES = pd.read_table(config["metaT_ena_table"])
INPUTDIR = config["inputDIR"]
ADAPTERS = config["adapters"]
METAG_STUDY = list(set(METAG_SAMPLES["study_accession"].tolist()))
METAT_STUDY = list(set(METAT_SAMPLES["study_accession"].tolist()))
SCRATCHDIR = config["scratch"]
OUTPUTDIR = config["outputDIR"]
METAG_SAMPLELIST = pd.read_table(config["metaG_sample_list"], index_col="Assembly_group")
METAG_ASSEMBLYGROUP = list(METAG_SAMPLELIST.index)
METAT_SAMPLELIST = pd.read_table(config["metaT_sample_list"], index_col="Assembly_group")
METAT_ASSEMBLYGROUP= list(METAT_SAMPLELIST.index)
ASSEMBLYGROUP = METAG_ASSEMBLYGROUP

#----COMPUTE VAR----#
MEGAHIT_CPU = config["megahit_cpu"]
MEGAHIT_MIN_CONTIG = config["megahit_min_contig"]
MEGAHIT_MEM = config["megahit_mem"]
MEGAHIT_OTHER = config["megahit_other"]

#----ACCESSION NUMS----#
metaG_run_accession = list(METAG_SAMPLES.run_accession)
metaT_run_accession = list(METAT_SAMPLES.run_accession)


#----FUNCTIONS----#

def identify_read_groups(assembly_group_name, STUDY, FORWARD=True):
    outlist=[] 
    ERR_list = METAG_SAMPLELIST.loc[assembly_group_name, 'ERR_list'].split(', ')
    if FORWARD: 
        num = 1
    else: 
        num = 2
    for E in ERR_list: 
        outlist.append(SCRATCHDIR + "/trimmed/{}/{}_{}.trimmed.fastq.gz".format(STUDY,E, num)) 
    return(outlist)

def get_sample_list(assembly_group, SAMPLELIST, STUDY):
    outlist = []
    for assembly_group_name in assembly_group: 
        ERR_list = SAMPLELIST.loc[assembly_group_name, 'ERR_list'].split(', ')
        for E in ERR_list:
            out = SCRATCHDIR + "/mapping/{}/{}/{}.bam".format(assembly_group_name,STUDY, E) 
            outlist.append(out)
    return(outlist)


def get_sample_list_onegroup(assembly_group_name, SAMPLELIST, STUDY):
    outlist = []
    ERR_list = SAMPLELIST.loc[assembly_group_name, 'ERR_list'].split(', ')
    for E in ERR_list:
        out = SCRATCHDIR + "/mapping/{}/{}/{}.bam".format(assembly_group_name, STUDY, E)
        outlist.append(out)
    return(outlist)        

#----QC DATA FILE----#

assert(len(METAT_STUDY)==1), 'This metaT ena table contains more than one study accession'
assert(len(METAG_STUDY)==1), 'This metaG ena table contains more than one study accession' 
assert(METAG_STUDY[0]==METAG_ACCESSION), 'The study accession provided in the config file does not match the study accession provided in the metaG ena table.'
assert(METAT_STUDY[0]==METAT_ACCESSION), 'The study accession provided in the config file does not match the study accession provided in the ena table.'
#assert(np.all(METAG_SAMPLELIST.index == METAT_SAMPLELIST.index)), 'The sample list provided for the MetaG does not match the MetaT'

pathlib.Path(OUTPUTDIR).mkdir(parents=True, exist_ok=True)

#----DEFINE RULES----#

localrules: multiqc, copy_bwa_index 

rule all: 
    input:
        # QC DATA
        fastqcZIP_rawG = expand("{base}/qc/fastqc/{study}/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, study = METAG_ACCESSION, sample=metaG_run_accession, num = [1,2]),
        fastqcZIP_rawT = expand("{base}/qc/fastqc/{study}/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, study = METAT_ACCESSION, sample=metaT_run_accession, num = [1,2]),  
        fastqcZIP_trimmedG = expand("{base}/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.zip", base = OUTPUTDIR, study = METAG_ACCESSION, sample=metaG_run_accession, num = [1,2]),  
        fastqcZIP_trimmedT = expand("{base}/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.zip", base = OUTPUTDIR, study = METAT_ACCESSION, sample=metaT_run_accession, num = [1,2]),  
        #MULTIQC
        html_rawG = OUTPUTDIR + "/qc/rawG_multiqc.html",
        stats_rawG = OUTPUTDIR + "/qc/rawG_multiqc_general_stats.txt",
        html_trimmedG = OUTPUTDIR + "/qc/trimmedG_multiqc.html",
        stats_trimmedG = OUTPUTDIR + "/qc/trimmedG_multiqc_general_stats.txt",
        html_rawT = OUTPUTDIR + "/qc/rawT_multiqc.html",
        stats_rawT = OUTPUTDIR + "/qc/rawT_multiqc_general_stats.txt",
        html_trimmedT = OUTPUTDIR + "/qc/trimmedT_multiqc.html",
        stats_trimmedT = OUTPUTDIR + "/qc/trimmedT_multiqc_general_stats.txt",
        #TRIM DATA
        trimmedDataG = expand("{base}/trimmed/{study}/{sample}_{num}.trimmed.fastq.gz", base = SCRATCHDIR, study = METAG_ACCESSION, sample=metaG_run_accession, num = [1,2]), 
        trimmedDataT = expand("{base}/trimmed/{study}/{sample}_{num}.trimmed.fastq.gz", base = SCRATCHDIR, study = METAT_ACCESSION, sample=metaT_run_accession, num = [1,2]), 

        #CALCULATE SOURMASH
        signatureG = expand("{base}/sourmash/{study}/{sample}.10k.sig", base = SCRATCHDIR, study = METAG_ACCESSION, sample = metaG_run_accession),
        signatureT = expand("{base}/sourmash/{study}/{sample}.10k.sig", base = SCRATCHDIR, study = METAT_ACCESSION, sample = metaT_run_accession),

        #ASSEMBLE
        assembly = expand("{base}/megahit/{assembly_group}/final.contigs.fa", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP),  
        assembly_copy = expand("{base}/bwa_index/{assembly_group}.fa", base = SCRATCHDIR, assembly_group = METAG_ASSEMBLYGROUP),  

        #BWA INDEX
        bwa_index = expand("{base}/bwa_index/{assembly_group}.fa.{bwa_tail}", base = SCRATCHDIR, assembly_group = METAG_ASSEMBLYGROUP, bwa_tail = ["amb", "ann", "bwt", "pac", "sa"]), 

        #BWA MAPPING:
        bwa_memG = get_sample_list(METAG_ASSEMBLYGROUP, METAG_SAMPLELIST, METAG_ACCESSION), 
        bwa_memT = get_sample_list(METAT_ASSEMBLYGROUP, METAT_SAMPLELIST, METAT_ACCESSION), 

        #BINNING 

        #METABAT2 
        jgi_abund = expand("{base}/metabat2/{assembly_group}/jgi_abund.txt", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP),
        metabat2_bins = expand("{base}/metabat2/{assembly_group}/{assembly_group}_bin", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP),
        #EUKREP
        eukrep =  expand("{base}/eukrep/{assembly_group}/euk.final.contigs.fa", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP), 
        metabat2_bins_euk = expand("{base}/metabat2_euk/{assembly_group}/{assembly_group}_eukbin", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP),
        
        #PROTEIN PREDICITION
        #PRODIGAL
        proteins = expand("{base}/prodigal/{assembly_group}/proteins.faa", base = OUTPUTDIR, assembly_group = METAG_ASSEMBLYGROUP), 

rule fastqc:
    input:
        INPUTDIR + "/{study}/{sample}/{sample}_{num}.fastq.gz"     
    output:
        html = OUTPUTDIR + '/qc/fastqc/{study}/{sample}_{num}_fastqc.html', 
        zip = OUTPUTDIR + '/qc/fastqc/{study}/{sample}_{num}_fastqc.zip'
    params: ""
    log: 
        OUTPUTDIR + '/logs/fastqc/{study}/{sample}_{num}.log'
    wrapper:
        "0.27.1/bio/fastqc"

rule trimmomatic: 
    input:
        r1 = INPUTDIR + "/{study}/{sample}/{sample}_1.fastq.gz", 
        r2 = INPUTDIR + "/{study}/{sample}/{sample}_2.fastq.gz" 
    output:
        r1 = SCRATCHDIR + "/trimmed/{study}/{sample}_1.trimmed.fastq.gz",
        r2 = SCRATCHDIR + "/trimmed/{study}/{sample}_2.trimmed.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired = SCRATCHDIR + "/trimmed/{study}/{sample}_1.unpaired.fastq.gz",
        r2_unpaired = SCRATCHDIR + "/trimmed/{study}/{sample}_2.unpaired.fastq.gz"
    log:
        OUTPUTDIR +  "/logs/trimmomatic/{study}/{sample}.log"
    params:
        trimmer=["ILLUMINACLIP:{}:2:30:7".format(ADAPTERS), "LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:2", "MINLEN:50"],
        extra=""
    wrapper:
        "0.27.1/bio/trimmomatic/pe"

rule fastqc_trimmed:
    input:
        SCRATCHDIR + "/trimmed/{study}/{sample}_{num}.trimmed.fastq.gz" 
    output:
        html = OUTPUTDIR + '/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.html', 
        zip = OUTPUTDIR + '/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.zip'
    params: ""
    log: 
        OUTPUTDIR + '/logs/fastqc/{study}/{sample}_{num}.trimmed.log'
    wrapper:
        "0.27.1/bio/fastqc"

rule multiqc:
    input:
        rawG = expand("{base}/qc/fastqc/{study}/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, study = METAG_ACCESSION, sample = metaG_run_accession, num = [1,2]), 
        trimmedG = expand("{base}/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.zip", base = OUTPUTDIR, study = METAG_ACCESSION, sample = metaG_run_accession, num = [1,2]), 
        rawT = expand("{base}/qc/fastqc/{study}/{sample}_{num}_fastqc.zip", base = OUTPUTDIR, study = METAT_ACCESSION, sample = metaT_run_accession, num = [1,2]), 
        trimmedT = expand("{base}/qc/fastqc/{study}/{sample}_{num}.trimmed_fastqc.zip", base = OUTPUTDIR, study = METAT_ACCESSION, sample = metaT_run_accession, num = [1,2]) 
    output:
        html_rawG = OUTPUTDIR + "/qc/rawG_multiqc.html", 
        stats_rawG = OUTPUTDIR + "/qc/rawG_multiqc_general_stats.txt",
        html_trimmedG = OUTPUTDIR + "/qc/trimmedG_multiqc.html", 
        stats_trimmedG = OUTPUTDIR + "/qc/trimmedG_multiqc_general_stats.txt",
        html_rawT = OUTPUTDIR + "/qc/rawT_multiqc.html", 
        stats_rawT = OUTPUTDIR + "/qc/rawT_multiqc_general_stats.txt",
        html_trimmedT = OUTPUTDIR + "/qc/trimmedT_multiqc.html", 
        stats_trimmedT = OUTPUTDIR + "/qc/trimmedT_multiqc_general_stats.txt"
    conda: 
        "envs/multiqc-env.yaml"
    shell: 
        """
        multiqc -n multiqc.html {input.rawG}
        mv multiqc.html {output.html_rawG}
        mv multiqc_data/multiqc_general_stats.txt {output.stats_rawG} 
        rm -rf multiqc_data

        multiqc -n multiqc.html {input.trimmedG}
        mv multiqc.html {output.html_trimmedG}
        mv multiqc_data/multiqc_general_stats.txt {output.stats_trimmedG} 
        rm -rf multiqc_data
        
        multiqc -n multiqc.html {input.trimmedT}
        mv multiqc.html {output.html_trimmedT}
        mv multiqc_data/multiqc_general_stats.txt {output.stats_trimmedT} 
        rm -rf multiqc_data
        
        multiqc -n multiqc.html {input.trimmedT}
        mv multiqc.html {output.html_trimmedT}
        mv multiqc_data/multiqc_general_stats.txt {output.stats_trimmedT} 
        rm -rf multiqc_data
        """ 

rule compute_sigs:
    input:
        r1 = SCRATCHDIR + "/trimmed/{study}/{sample}_1.trimmed.fastq.gz",
        r2 = SCRATCHDIR + "/trimmed/{study}/{sample}_2.trimmed.fastq.gz" 
    output: 
        SCRATCHDIR + "/sourmash/{study}/{sample}.10k.sig"
    conda: 
        "envs/sourmash.yaml"
    log:
         OUTPUTDIR +  "/logs/sourmash/{study}/sourmash_{sample}.log"
    shell: 
        """
        zcat {input.r1} {input.r2} | sourmash compute -k 21,31,51\
            --scaled 10000  --track-abundance \
            -o {output} - 2> {log}
        """

rule megahit_assembly: 
    input: r1 = lambda wildcards: identify_read_groups("{assembly_group}".format(assembly_group=wildcards.assembly_group), METAG_ACCESSION), 
           r2 = lambda wildcards: identify_read_groups("{assembly_group}".format(assembly_group=wildcards.assembly_group), METAG_ACCESSION, FORWARD=False) 
    output: 
       OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa"  
    conda: 
        "envs/megahit.yaml"
    log: 
        OUTPUTDIR + "/logs/megahit/{assembly_group}.log" 
    params: 
        inputr1 = lambda wildcards, input: ','.join(input.r1),
        inputr2 = lambda wildcards, input: ','.join(input.r2),
        min_contig_len = MEGAHIT_MIN_CONTIG,
        cpu_threads = MEGAHIT_CPU, 
        memory = MEGAHIT_MEM, 
        other_options = MEGAHIT_OTHER,
        megahit_output_name = lambda wildcards: "{}/megahit/{}".format(OUTPUTDIR, wildcards.assembly_group),
        megahit_output_prefix = lambda wildcards: "{}".format(wildcards.assembly_group),
    shell: 
        """
        megahit -1 {params.inputr1} -2 {params.inputr2} --min-contig-len {params.min_contig_len} --memory {params.memory} -t {params.cpu_threads} --out-dir {params.megahit_output_name} {params.other_options}  >> {log} 2>&1
        """

rule bwa_index:
    input:
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa"
    output:
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa.amb",
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa.ann",
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa.bwt",
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa.pac",
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa.sa"
    log:
        OUTPUTDIR + "/logs/bwa_index/{assembly_group}.log"
    params:
        algorithm="bwtsw"
    conda:
        "envs/metabat-env.yaml"
 
    shell:
        """
        bwa index {input} 2> {log}
        """

rule copy_bwa_index:
    input: 
        OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa"
    output:
        SCRATCHDIR + "/bwa_index/{assembly_group}.fa"
    shell:
        """
        cp {input} {output}
        """

rule bwa_mem:
    input:
        amb = SCRATCHDIR + "/bwa_index/{assembly_group}.fa.amb", 
        ann = SCRATCHDIR + "/bwa_index/{assembly_group}.fa.ann", 
        bwt = SCRATCHDIR + "/bwa_index/{assembly_group}.fa.bwt", 
        pac = SCRATCHDIR + "/bwa_index/{assembly_group}.fa.pac",  
        sa = SCRATCHDIR + "/bwa_index/{assembly_group}.fa.sa", 
        reference = SCRATCHDIR + "/bwa_index/{assembly_group}.fa", 
        r1 = SCRATCHDIR + "/trimmed/{study}/{sample}_1.trimmed.fastq.gz",
        r2 = SCRATCHDIR + "/trimmed/{study}/{sample}_2.trimmed.fastq.gz", 
    output:
        SCRATCHDIR + "/mapping/{assembly_group}/{study}/{sample}.bam"
    log:
        OUTPUTDIR + "/logs/bwa_mem/{assembly_group}/{study}/{sample}.log"
    params: 
        extra="", 
        #pipe_cmd = "samtools sort -o {output} -",  
        threads = 8
    conda: 
        "envs/metabat-env.yaml"
    shell:
        """ 
        bwa mem -t {params.threads} {params.extra} {input.reference} {input.r1} {input.r2} | samtools sort -o {output} - >> {log} 2>&1
        """ 

rule metabat_abundance:
    input:
        lambda wildcards: get_sample_list_onegroup("{assembly_group}".format(assembly_group=wildcards.assembly_group), METAG_SAMPLELIST, METAG_ACCESSION)
    output:
        OUTPUTDIR + "/metabat2/{assembly_group}/jgi_abund.txt"
    conda:
         "envs/metabat-env.yaml"
    log:
        OUTPUTDIR + "/logs/metabat2/{assembly_group}.abun.log"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input} > {log} 2>&1
        """

rule metabat_binning:
    input:
        assembly = OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa",
        depth = OUTPUTDIR + "/metabat2/{assembly_group}/jgi_abund.txt"
    output:
        OUTPUTDIR + "/metabat2/{assembly_group}/{assembly_group}_bin"
    conda:
         "envs/metabat-env.yaml"
    params: 
        other = "--saveCls",
        threads = 8
    log:
        OUTPUTDIR + "/logs/metabat2/{assembly_group}.bin.log"
    shell:
        """
        metabat2 {params.other} --numThreads {params.threads} -i {input.assembly} -a {input.depth} -o {output} > {log} 2>&1
        """

rule eukrep:
    input: 
        assembly = OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa",
    output:
        OUTPUTDIR + "/eukrep/{assembly_group}/euk.final.contigs.fa"
    conda: 
        "envs/EukRep.yaml"
    log:
        OUTPUTDIR + "/logs/eukrep/{assembly_group}.eukrep.log"
    params: 
        prok = OUTPUTDIR + "/eukrep/{assembly_group}/prok.final.contigs.fa",
        min_contig = 1000
    shell: 
        """
        EukRep -i {input} -o {output} --prokarya {params.prok} --min {params.min_contig} > {log} 2>&1
        """

rule metabat_binning_euk:
    input:
        assembly = OUTPUTDIR + "/eukrep/{assembly_group}/euk.final.contigs.fa",
        depth = OUTPUTDIR + "/metabat2/{assembly_group}/jgi_abund.txt"
    output:
        OUTPUTDIR + "/metabat2_euk/{assembly_group}/{assembly_group}_eukbin"
    conda:
         "envs/metabat-env.yaml"
    params:
        other = "--saveCls",
        threads = 8
    log:
        OUTPUTDIR + "/logs/metabat2/{assembly_group}.eukbin.log"
    shell:
        """
        metabat2 {params.other} --numThreads {params.threads} -i {input.assembly} -a {input.depth} -o {output} > {log} 2>&1
        """

rule prodigal:
    input:
        assembly = OUTPUTDIR + "/megahit/{assembly_group}/final.contigs.fa",
    output: 
        proteins = OUTPUTDIR + "/prodigal/{assembly_group}/proteins.faa",
        genes = OUTPUTDIR + "/prodigal/{assembly_group}/genes.gff"
    conda:
        "envs/prodigal.yaml"
    shell:
        """
        prodigal -i {input.assembly} -f gff -o {output.genes} -a {output.proteins} -p meta
        """
