# An example collection of Snakemake rules imported in the main Snakefile.
rule DownloadChromHMM:
    output:
        # temp("Misc/Gm12878.chromHMM.hg19.bed.gz")
        "Misc/Gm12878.chromHMM.hg19.bed.gz"
    shell:
        "wget -O {output} http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationChromhmmGm12878.bed.gz"

rule LiftOverChromHMM:
    input:
        InputBedGz = "Misc/Gm12878.chromHMM.hg19.bed.gz",
        chain = config["Misc"]["hg19Tohg38Liftover"]
    output:
        OutputBed = "../../data/Gm12878.chromHMM.hg38.bed.gz"
    shell:
        "liftOver {input.InputBedGz} {input.chain} /dev/stdout /dev/null | gzip - > {output.OutputBed}"

rule DownloadGencodeAnnotations:
    output:
        "Misc/GencodeOriginalFiles/gencode.v31.annotation.gtf.gz"
    shell:
        """
        wget -O {output} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
        """

# rule FilterForHighSupportTranscripts:
#     input:
#         "Misc/GencodeOriginalFiles/gencode.v31.annotation.gtf.gz"
#     output:
#         ""

rule GetUniqIntrons:
    input:
        bed = "Misc/GencodeHg38_all_introns.corrected.bed",
        faidx = config["Human_ref"]["genome_fasta"] + ".fai"
    output:
        "Misc/GencodeHg38_all_introns.corrected.uniq.bed"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '{{ print $1, $2, $3,$5"_"$1"_"$2"_"$3"_"$6,".", $6 }}' {input.bed} | sort | uniq | bedtools sort -i - -faidx {input.faidx} > {output}
        """

rule MakeWindowsForLongIntrons:
    """
    For calculating downward slope over introns. Like in Jonkers & Lis. Split
    introns into equal number of bins
    """
    input:
        bed = "Misc/GencodeHg38_all_introns.corrected.uniq.bed",
        faidx = config["Human_ref"]["genome_fasta"] + ".fai"
    output:
        "Misc/GencodeHg38_all_introns.corrected.uniq.bed.LongIntronWindows.bed"
    shell:
        """
        set +o pipefail;
        bedtools intersect -a {input.bed} -b {input.bed} -c -s -sorted | awk -F '\\t' '$NF==1 && ($3-$2)>10000' | bedtools makewindows -b - -n 100 -i srcwinnum | bedtools sort -i - -faidx {input.faidx} | awk -F'\\t' -v OFS='\\t' '{{ split($4,a,"_") }} a[5]=="+" {{print $1,$2,$3,$4,".",a[5]}} a[5]=="-" {{print $1,$2,$3, a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]"_"101-a[6],".", a[5] }}' > {output}
        """

rule CountReadsInLongIntronWindows:
    """
    Unstranded counting. Because of different lib types.
    I dunno why this rule is sometimes throwing non-zero exit codes since the output files seem to look fine. Need to set more permissive shell options to let snakemake run successfully.. I suspect samtools view is throwing a non-zero exit code for reasons I don't understand
    """
    input:
        bed = "Misc/GencodeHg38_all_introns.corrected.uniq.bed.LongIntronWindows.bed",
        faidx = config["Human_ref"]["genome_fasta"] + ".fai",
        bam = "Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam",
        bai =  "Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        bed = "LongIntronWindowCounts/{sample}.bed.gz"
    log:
        "logs/CountReadsInLongIntronWindows/{sample}.log"
    shell:
        """
        set +e
        set +o pipefail
        (samtools view -bh -F256 {input.bam} | bedtools intersect -sorted -g {input.faidx} -a {input.bed} -b - -c -split -F 0.5 | gzip - > {output}) &> {log}
        """

rule MoveLongIntronCountsToOutput:
    input:
        bed = "LongIntronWindowCounts/{sample}.bed.gz"
    output:
        bed = "../output/LongIntronWindowCounts/{sample}.bed.gz"
    log:
        "logs/MoveLongIntronCountsToOutput/{sample}.log"
    shell:
        """
        cp {input} {output}
        """

# rule CountReadsInLongIntronWindows_gzip:
#     input:
#         bed = "LongIntronWindowCounts/{sample}.bed"
#     output:
#         bed = "LongIntronWindowCounts/{sample}.bed.gz"
#     shell:
#         """
#         cat {input} | gzip - > {output}
#         """



rule MakeWindowsForLongIntrons_equalSized:
    """
    For calculating downward slope over introns. Like in Jonkers & Lis. Split
    introns into equal length bins
    """
    input:
        bed = "Misc/GencodeHg38_all_introns.corrected.uniq.bed",
        faidx = config["Human_ref"]["genome_fasta"] + ".fai"
    output:
        "Misc/GencodeHg38_all_introns.corrected.uniq.bed.LongIntronWindows_equalLength.bed"
    params:
        WinLen = 100
    shell:
        """
        set +o pipefail;
        bedtools intersect -a {input.bed} -b {input.bed} -c -s -sorted | awk -F '\\t' '$NF==1 && ($3-$2)>10000' | bedtools makewindows -b - -n 100 -i srcwinnum | bedtools sort -i - -faidx {input.faidx} | awk -F'\\t' -v OFS='\\t' '{{ split($4,a,"_") }} a[5]=="+" {{print $1,$2,$3,$4,".",a[5]}} a[5]=="-" {{print $1,$2,$3, a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]"_"101-a[6],".", a[5] }}' > {output}
        """

rule DownloadVeloso_Txn_Rates:
    output:
        "../../output/VelosoTxnRates.xlsx"
    shell:
        """
        wget -O {output} https://genome.cshlp.org/content/suppl/2014/04/16/gr.171405.113.DC1/Supplemental_Table1.xlsx
        """
