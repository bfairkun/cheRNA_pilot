rule STAR_alignment_FirstPassQuick:
    input:
        index_human = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/STARIndex/chrLength.txt",
        R1 = "Sample_Fastq/{sample}.fastq.gz"
    log:
        "logs/STAR_FirstPass/{sample}.log"
    threads: 8
    params:
        index = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/STARIndex/",
        clippingParam = GetClippingParameterForLexogenDataset
    output:
        "Alignments/FirstPassQuick/{sample}/Aligned.sortedByCoord.out.bam"
    shell:
        """
        ulimit -v 31538428657
        STAR --readMapNumber 10000000 --runThreadN {threads} --limitOutSJcollapsed 2000000 --genomeDir {params.index} --readFilesIn <(zcat {input.R1} {params.clippingParam})   --outSAMmultNmax 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix Alignments/FirstPassQuick/{wildcards.sample}/ &> {log}
        """


rule QualimapRnaseqQuick:
    input:
        gtf="/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
        bam="Alignments/FirstPassQuick/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "qualimap_quick/{sample}/rnaseq_qc_results.txt"
    log:
        "logs/QualimapQucikRnaseq/{sample}.log"
    params:
        libtype = GetQualimapLibtype
    shell:
        """
        qualimap rnaseq -bam {input.bam} -gtf {input.gtf} -p {params.libtype} --java-mem-size=8G -outdir qualimap_quick/{wildcards.sample}/ &> {log}
        """

rule MultiQCQuick:
    input:
        expand("qualimap_quick/{sample}/rnaseq_qc_results.txt", sample=SampleList),
        expand("Alignments/FirstPassQuick/{sample}/Aligned.sortedByCoord.out.bam", sample=SampleList),
    log: "logs/MultiqcQuick.log"
    output:
        "MultiqcQuick/multiqc_report.html"
    shell:
        """
        multiqc -f -o MultiqcQuick/ qualimap_quick/ Alignments/FirstPassQuick/ &> {log}
        """

rule idxStat:
    input:
        "Alignments/FirstPassQuick/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        idx = "QuickAnalyses/idxstats/{sample}.tab",
        bai = "Alignments/FirstPassQuick/{sample}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        """
        samtools index {input}
        samtools idxstats {input} > {output.idx}
        """

rule ConcatIdxStat:
    input:
        expand(  "QuickAnalyses/idxstats/{sample}.tab" , sample=SampleList),
    output:
        "../../output/QuickAnalyses/IdxStats.tab"
    shell:
        """
        awk -F'\\t' '{{ print $1, $2, $3, FILENAME }}' {input} > {output}
        """

rule QuickFeatureCounts:
    input:
        bam = expand("Alignments/FirstPassQuick/{sample}/Aligned.sortedByCoord.out.bam", sample=SampleList),
        gtf="/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Annotations/gencode.v34.primary_assembly.annotation.gtf",
    output:
        "QuickAnalyses/FeatureCounts.txt"
    shell:
        """
        featureCounts -M -a {input.gtf} -o {output} {input.bam}
        """

SedReplace = ("Alignments/FirstPassQuick/").replace("/", 
"\/")
rule reformatQuickFeatureCounts:
    input:
        "QuickAnalyses/FeatureCounts.txt"
    output:
        "../../output/QuickAnalyses/FeatureCounts.txt.gz"
    shell:
        """
        cat {input} | sed -e '2s/{SedReplace}//g' | sed -e '2s/\/Aligned.sortedByCoord.out.bam//g' | gzip - > {output}
        """

rule QuickBigiwgs:
    input:
        bam = "Alignments/FirstPassQuick/{sample}/Aligned.sortedByCoord.out.bam",
        fai = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/Fasta/GRCh38.primary_assembly.genome.fa.fai"
    output:
        "QuickAnalyses/Bigwigs/{sample}.bw"
    shell:
        """
        scripts/BamToBigwig.sh {input.fai} {input.bam} {output}
        """
