Human_genomeDir = config["Human_ref"]["genome_dir"][:-1]

rule STAR_make_index_human:
    input:
        fasta=config["Human_ref"]["genome_fasta"],
        gtf=config["Human_ref"]["genome_gtf"]
    output:
        index = config["Human_ref"]["genome_dir"] + "chrLength.txt",
    log:
        "logs/STAR_make_index_human.log"
    params:
        genomeDir = Human_genomeDir
    shell:
        """
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN 4 --genomeDir {params.genomeDir} --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} &> {log}
        """

rule STAR_alignment_FirstPass:
    input:
        index_human = config["Human_ref"]["genome_dir"] + "chrLength.txt",
        R1 = "Sample_Fastq/{sample}.fastq.gz"
    log:
        "logs/STAR_FirstPass/{sample}.log"
    threads: 8
    params:
        index = Human_genomeDir
    output:
        SJout = "Alignments/FirstPass/{sample}/SJ.out.tab"
    shell:
        """
        ulimit -v 31538428657
        STAR --runThreadN {threads} --genomeDir {params.index} --readFilesIn {input.R1} --readFilesCommand zcat --outSAMmultNmax 1 --outSAMtype None --outFileNamePrefix Alignments/FirstPass/{wildcards.sample}/ &> {log}
        """

rule STAR_make_index_human_AfterFirstPass:
    input:
        fasta=config["Human_ref"]["genome_fasta"],
        gtf=config["Human_ref"]["genome_gtf"],
        FirstPass_sjdb=expand( "Alignments/FirstPass/{sample}/SJ.out.tab", sample=SampleList),
    output:
        index = "Alignments/FirstPass_sjdb_index/chrLength.txt"
    log:
        "logs/STAR_make_index_human_AfterFirstPass.log"
    shell:
        """
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN 4 --genomeDir Alignments/FirstPass_sjdb_index/ --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} --sjdbFileChrStartEnd {input.FirstPass_sjdb} &> {log}
        """

rule STAR_alignment_SecondPass:
    input:
        index = "Alignments/FirstPass_sjdb_index/chrLength.txt",
        R1 = "Sample_Fastq/{sample}.fastq.gz"
    log:
        "logs/STAR_SecondPass/{sample}.log"
    threads: 8
    output:
        SJout = "Alignments/SecondPass/{sample}/SJ.out.tab",
        bam = "Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam",
        bai =  "Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        """
        ulimit -v 31538428657
        STAR --runThreadN {threads} --genomeDir Alignments/FirstPass_sjdb_index --readFilesIn <(zcat {input.R1} | fastx_trimmer -l 46) --outSAMtype BAM SortedByCoordinate --outWigType wiggle --outFileNamePrefix Alignments/SecondPass/{wildcards.sample}/ &> {log}
        samtools index {output.bam}
        """

rule unzip_AnnotatedSpliceBed:
    input:
        "../../data/GencodeHg38_all_introns.corrected.bed.gz"
    output:
        "Misc/GencodeHg38_all_introns.corrected.bed"
    shell:
        "zcat {input} > {output}"

rule AnnotateSplicingTypes:
    input:
        SJout = "Alignments/SecondPass/{sample}/SJ.out.tab",
        AnnotatedIntrons = "Misc/GencodeHg38_all_introns.corrected.bed"
    output:
        SJout_AS_annotated = "Alignments/SecondPass/{sample}/SJ.out.annotated.tab"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$4=="2" {{print $1,$2-1,$3,".",$7,"-"}} $4=="1" {{print $1,$2-1,$3,".",$7,"+"}}' {input.SJout} | scripts/AnnotateSplicingType.py -I - -A {input.AnnotatedIntrons} -O {output}
        """

rule IntersectIntronsWithGenes:
    input:
        SJout_AS_annotated = "Alignments/SecondPass/{sample}/SJ.out.annotated.tab",
        Genes = "../../data/GtexGenes.bed"
    output:
        config["gitinclude_output"] + "SJoutAnnotatedAndIntersected/{sample}.tab.gz"
    shell:
        """
        bedtools intersect -a {input.SJout_AS_annotated} -b {input.Genes} -wao -s | gzip - > {output}
        """

# SedReplace = ("Alignments/SecondPass/").replace("/", 
# "\/")
# rule PowerAnalysisSubreadCount:
#     input:
#         expand( "Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam", sample=SampleList),
#     output:
#         Output = "PowerAnalysis/Subread/{species}.subread.txt.gz",
#         OutputGit = config["gitinclude_output"] + "PowerAnalysisFullCountTable.{species}.subread.txt.gz",
#         OutputSummary = "PowerAnalysis/Subread/{species}.subread.txt.summary"
#     log:
#         "logs/SubreadCount.txt"
#     params:
#         Gtf=config["Human_ref"]["genome_gtf"]
#     shell:
#         """
#         ~/software/subread-1.6.3-Linux-x86_64/bin/featureCounts -a {params.Gtf} -o {wildcards.species}.temp.txt {input} &> {log}
#         # Format the header for the count table to just contain sample name
#         # instead of bam filepath
#         cat {wildcards.species}.temp.txt | sed -e '2s/{SedReplace}//g' | sed -e '2s/\/Aligned.sortedByCoord.out.bam//g' | gzip - > {output.Output}
#         cat {wildcards.species}.temp.txt.summary | sed -e '1s/{SedReplace}//g' | sed -e '1s/\/Aligned.sortedByCoord.out.bam//g' > {output.OutputSummary}
#         rm {wildcards.species}.temp.txt {wildcards.species}.temp.txt.summary
#         cp {output.Output} {output.OutputGit}
#         """
