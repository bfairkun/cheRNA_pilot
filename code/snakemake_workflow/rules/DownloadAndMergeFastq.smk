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
    

rule ConcatFastqPerSample:
    input:
        GetInput_ConcatFastqPerSample
    output:
        "Sample_Fastq/{sample}.fastq.gz"
    shell:
        "cat {input} > {output}"
