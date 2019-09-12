
rule normalized_samples: 
    input: 
        r1 = OUTPUTDIR + "/trimmed/{sample}_1.trimmed.fastq.gz",   
        r2 = OUTPUTDIR + "/trimmed/{sample}_2.trimmed.fastq.gz"
    output: 
        r1 = OUTPUTDIR + "/normalized/{sample}_1.trimmed.normalized.fastq.gz",
        r2 = OUTPUTDIR + "/normalized/{sample}_2.trimmed.normalized.fastq.gz", 
        inhist = OUTPUTDIR + "/normalized/{sample}.inhist",
        outhist = OUTPUTDIR + "/normalized/{sample}.outhist",
    params: 
        bb_targetDepth = 30, 
        bb_minDepth = 2,  
        bb_otherparams = "", 
        bb_threads = "", 
        r1 = OUTPUTDIR + "/normalized/{sample}_1.trimmed.normalized.fastq",
        r2 = OUTPUTDIR + "/normalized/{sample}_2.trimmed.normalized.fastq",
    conda: 
        'envs/bbmap-env.yaml'
    log: 
        OUTPUTDIR + "/logs/normalized/bbnorm_{sample}.log" 
    shell: 
        """ 
        bbnorm.sh in1={input.r1} in2={input.r2} out={params.r1} out2={params.r2} target={params.bb_targetDepth} min={params.bb_minDepth} hist={output.inhist} histout={output.outhist} {params.bb_otherparams} 2> {log} 
        pigz {params.r1}
        pigz {params.r2}
        """

rule interleave_reads:
    input: 
        r1 =  OUTPUTDIR + "/trimmed/{sample}_1.trimmed.fastq.gz", 
        r2 = OUTPUTDIR + "/trimmed/{sample}_2.trimmed.fastq.gz"
    output: 
        temp(SCRATCHDIR + "/{sample}.trimmed.interleaved.fastq.gz") 
    conda: 
        "envs/trim_low_abund.yaml"
    log: 
         OUTPUTDIR +  "/logs/errtrim/interleave_PE_{sample}.log"
    shell: 
        '''  
            interleave-reads.py --gzip {input.r1} {input.r2} -o {output} 2> {log} 
        '''

rule split_reads: 
    input: 
        SCRATCHDIR + "/{sample}.trimmed.interleaved.errtrim.fastq.gz",
    output:  
        r1 = OUTPUTDIR + "/errtrim/{sample}_1.trimmed.errtrim.fastq.gz",
        r2 = OUTPUTDIR + "/errtrim/{sample}_2.trimmed.errtrim.fastq.gz",        
        r0= OUTPUTDIR + "/errtrim/{sample}_SE.trimmed.errtrim.fastq.gz",
    params:   
        tmp_r1 = SCRATCHDIR + "/{sample}_1.trimmed.errtrim.fastq.gz", 
        tmp_r2 = SCRATCHDIR + "/{sample}_2.trimmed.errtrim.fastq.gz",
        tmp_r0 = SCRATCHDIR + "/{sample}_SE.trimmed.errtrim.fastq.gz",
    conda: 
        'envs/trim_low_abund.yaml'
    log: 
        OUTPUTDIR +  "/logs/errtrim/split_PE_{sample}.log"
    shell:
        ''' 
        split-paired-reads.py {input} --gzip -0 {params.tmp_r0} -1 {params.tmp_r1} -2 {params.tmp_r2} 2> {log} 
        mv {params.tmp_r1} {output.r1}
        mv {params.tmp_r2} {output.r2}
        mv {params.tmp_r0} {output.r0} 
        '''

rule trim_low_abund:
    input:
         SCRATCHDIR + "/{sample}.trimmed.interleaved.fastq.gz"  
    output: 
         temp(SCRATCHDIR +  "/{sample}.trimmed.interleaved.errtrim.fastq.gz")
    params: 
        memory = "20e9",  
        other = "-C 2 -Z 18 -V -k31 --gzip"
    conda: 
        'envs/trim_low_abund.yaml'
    log:
         OUTPUTDIR +  "/logs/errtrim/trim_low_abund_{sample}.log" 
    shell:
        """  
        trim-low-abund.py -M {params.memory} {params.other} -o {output} {input} 2> {log} 
        """


