rule DownloadSRA:
    output:
        "SRA_Fastq/{SRR_accension}.fastq.gz"
    log:
        "logs/DownloadSRA/{SRR_accension}.log"
    shell:
        """
        fastq-dump --gzip {wildcards.SRR_accension} -O SRA_Fastq/ > {log}
        """

# rule GatherSRA:
#     input:
#         SRA_fastq_for_download

def GetInput_ConcatFastqPerSample(wildcards):
    SRA_fastq_list = expand("SRA_Fastq/{SRA_Accension}.fastq.gz", SRA_Accension=[x for x in FastQList.loc[FastQList['Sample'] == wildcards.sample ]["SRA"].tolist() if str(x) != 'nan'] )
    Other_fastq_list = expand("{filename}", filename=[x for x in FastQList.loc[FastQList['Sample'] == wildcards.sample ]["R1"].tolist() if str(x) != 'nan'] )
    return SRA_fastq_list + Other_fastq_list
    

def GetInput_ConcatFastqPerSample_R2(wildcards):
    SRA_fastq_list = expand("SRA_Fastq/{SRA_Accension}.fastq.gz", SRA_Accension=[x for x in FastQList.loc[FastQList['Sample'] == wildcards.sample ]["SRA"].tolist() if str(x) != 'nan'] )
    Other_fastq_list = expand("{filename}", filename=[x for x in FastQList.loc[FastQList['Sample'] == wildcards.sample ]["R2"].tolist() if str(x) != 'nan'] )
    return SRA_fastq_list + Other_fastq_list

rule ConcatFastqPerSample:
    input:
        GetInput_ConcatFastqPerSample
    output:
        "Sample_Fastq/{sample}.fastq.gz"
    shell:
        "cat {input} > {output}"

rule ConcatFastqPerSample_PE:
    input:
        R1 = GetInput_ConcatFastqPerSample,
        R2 = GetInput_ConcatFastqPerSample_R2
    output:
        R1 = "Sample_Fastq_PE/{sample}.R1.fastq.gz",
        R2 = "Sample_Fastq_PE/{sample}.R2.fastq.gz"
    shell:
        """
        cat {input.R1} > {output.R1}
        cat {input.R2} > {output.R2}
        """


