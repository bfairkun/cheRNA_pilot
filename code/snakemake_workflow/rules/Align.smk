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

rule hisat2_make_index:
    input:
        fasta=config["Human_ref"]["genome_fasta"]
    output:
        index = config["Human_ref"]["genome_dir_hisat"] + "genome.1.ht2",
    log:
        "logs/hisat2_make_index.log"
    params:
        index_prefix = config["Human_ref"]["genome_dir_hisat"] + "genome",
    threads: 4
    shell:
        """
        hisat2-build -p {threads} {input} {params.index_prefix} &> {log}
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

rule STAR_alignment_SecondPass_PE:
    input:
        index = "Alignments/FirstPass_sjdb_index/chrLength.txt",
        R1 = "Sample_Fastq_PE/{sample}.R1.fastq.gz",
        R2 = "Sample_Fastq_PE/{sample}.R2.fastq.gz"
    log:
        "logs/STAR_SecondPass_PE/{sample}.log"
    threads: 8
    output:
        SJout = "Alignments/SecondPass_PE/{sample}/SJ.out.tab",
        bam = "Alignments/SecondPass_PE/{sample}/Aligned.out.bam",
        # R1Unmapped = "Alignments/SecondPass_PE/{sample}/Unmapped.out.mate1.fastq.gz",
        # R2Unmapped = "Alignments/SecondPass_PE/{sample}/Unmapped.out.mate2.fastq.gz"
    shell:
        """
        ulimit -v 31538428657
        STAR --runThreadN {threads} --genomeDir Alignments/FirstPass_sjdb_index --readFilesIn {input.R1} {input.R2} --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix Alignments/SecondPass_PE/{wildcards.sample}/ --outReadsUnmapped Fastx --alignEndsType EndToEnd &> {log}
        samtools view -bh Alignments/SecondPass_PE/{wildcards.sample}/Aligned.out.bam > {output.bam}
        samtools index {output.bam}
        """

rule sortBam:
    input:
        "{Filepath}.bam"
    output:
        sortedbam = "{Filepath}.sorted.bam",
        index = "{Filepath}.sorted.bam.bai"
    shell:
        """
        samtools sort {input} > {output.sortedbam}
        samtools index {output.sortedbam}
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

rule AlignLariatJunctions_PE:
    input:
        R1Unmapped = "Alignments/SecondPass_PE/{sample}/Unmapped.out.mate1",
        R2Unmapped = "Alignments/SecondPass_PE/{sample}/Unmapped.out.mate2"
    log:
        "logs/AlignLariatJunctions_PE/{sample}.log"
    output:
    params:
        index_prefix = config["Human_ref"]["genome_dir_hisat"] + "genome",
    shell:
        """
        bash scripts/AlignLariatJunctions_SingleSample.sh {params.index_prefix} {output.}
        """

rule unzip_AnnotatedSpliceBed:
    input:
        "../../data/GencodeHg38_all_introns.corrected.bed.gz"
    output:
        "Misc/GencodeHg38_all_introns.corrected.bed"
    params:
        faidx = config["Human_ref"]["genome_fasta"] + ".fai"
    shell:
        "zcat {input} | bedtools sort -i - -faidx {params.faidx} > {output}"

rule Make3ssUpstreamAndDownstreamWindows:
    input:
        "Misc/GencodeHg38_all_introns.corrected.bed"
    output:
        upstream = "Misc/GencodeHg38_all_introns.upstream3ssWindows.bed",
        downstream = "Misc/GencodeHg38_all_introns.downstream3ssWindows.bed",
    params:
        faidx = config["Human_ref"]["genome_fasta"] + ".fai"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$6=="+" {{print $1,$3-25,$3,$1"."$3".+",".",$6}} $6=="-" {{print $1,$2,$2+25,$1"."$2".-",".",$6}}' Misc/GencodeHg38_all_introns.corrected.bed | sort | uniq | bedtools sort -i - -faidx {params.faidx} > {output.upstream}
        bedtools shift -i {output.upstream} -p 25 -m -25 -g {params.faidx} | bedtools sort -i - -faidx {params.faidx} > {output.downstream}
        """

rule CountUpstreamAndDownstreamCoverage:
    """
    Output a bedfile with 7th column for coverage in window upstream of 3ss and
    8th column as downstream coverage
    """
    input:
        upstream = "Misc/GencodeHg38_all_introns.upstream3ssWindows.bed",
        downstream = "Misc/GencodeHg38_all_introns.downstream3ssWindows.bed",
        bam = "Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam",
        bai =  "Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        "3ssCoverage/{sample}.bed.gz"
    params:
        faidx = config["Human_ref"]["genome_fasta"] + ".fai"
    shell:
        """
        paste <(samtools view -bh -F256 {input.bam} | bedtools intersect -sorted -g {params.faidx} -a {input.upstream} -b - -c -split) <(samtools view -bh -F256 {input.bam} | bedtools intersect -sorted -g {params.faidx} -a {input.downstream} -b - -c -split) | awk -F'\\t' -v OFS='\\t' '{{ print $1,$2,$3,$4,$5,$6,$7,$14 }}' | gzip - > {output}
        """

rule MoveUpstreamAndDownstreamCoverageToGitRepo:
    input:
        "3ssCoverage/{sample}.bed.gz"
    output:
        config["gitinclude_output"] + "3ssCoverageBeds/{sample}.bed.gz"
    shell: "cp {input} {output}"

rule AnnotateSplicingTypes:
    input:
        SJout = "Alignments/SecondPass/{sample}/SJ.out.tab",
        AnnotatedIntrons = "Misc/GencodeHg38_all_introns.corrected.bed"
    output:
        SJOutASBed = "Alignments/JunctionBeds/{sample}.junc",
        SJout_AS_annotated = "Alignments/SecondPass/{sample}/SJ.out.annotated.tab"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$4=="2" {{print $1,$2-1,$3,".",$7,"-"}} $4=="1" {{print $1,$2-1,$3,".",$7,"+"}}' {input.SJout} | tee {output.SJOutASBed} | scripts/AnnotateSplicingType.py -I - -A {input.AnnotatedIntrons} -O {output.SJout_AS_annotated}
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

rule faidx:
    input: config["Human_ref"]["genome_fasta"]
    output: config["Human_ref"]["genome_fasta"] + ".fai"
    shell: "samtools faidx {input}"

rule BedgraphAndBigwigs:
    input:
        bam = "Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam",
        fai =config["Human_ref"]["genome_fasta"] + ".fai"
    output:
        PlusBg = "Alignments/SecondPass/{sample}/Plus.bg",
        MinusBg = "Alignments/SecondPass/{sample}/Minus.bg",
        UnstrandBg = "Alignments/SecondPass/{sample}/Unstranded.bg",
        PlusBw = "Alignments/SecondPass/{sample}/Plus.bw",
        MinusBw = "Alignments/SecondPass/{sample}/Minus.bw",
        UnstrandBw = "Alignments/SecondPass/{sample}/Unstranded.bw",
    shell:
        """
        samtools view -bh -F256 {input.bam} | bedtools genomecov -max 10000 -ibam - -bg -strand + -split -scale $(bc <<< "scale=6;1000/$(samtools idxstats {input.bam} | awk '{{sum+=$2}} END{{printf sum}}')") | sort -k1,1 -k2,2n > {output.PlusBg}
        samtools view -bh -F256 {input.bam} | bedtools genomecov -max 10000 -ibam - -bg -strand - -split -scale $(bc <<< "scale=6;1000/$(samtools idxstats {input.bam} | awk '{{sum+=$2}} END{{printf sum}}')") | sort -k1,1 -k2,2n > {output.MinusBg}
        samtools view -bh -F256 {input.bam} | bedtools genomecov -max 10000 -ibam - -bg -split -scale $(bc <<< "scale=6;1000/$(samtools idxstats {input.bam} | awk '{{sum+=$2}} END{{printf sum}}')") | sort -k1,1 -k2,2n > {output.UnstrandBg}
        bedGraphToBigWig {output.PlusBg} {input.fai} {output.PlusBw}
        bedGraphToBigWig {output.MinusBg} {input.fai} {output.MinusBw}
        bedGraphToBigWig {output.UnstrandBg} {input.fai} {output.UnstrandBw}
        """

rule leacutter_cluster:
    input:
        expand("Alignments/JunctionBeds/{sample}.junc", sample=SampleList)
    output:
        juncfilelist = "leafcutter/junctionfiles.txt",
        counts = "leafcutter/leafcutter_perind.counts.gz",
        numers = "leafcutter/leafcutter_perind_numers.counts.gz"
    log:
    shell:
        """
        echo {input} | tr " " '\\n' > {output.juncfilelist}
        leafcutter_cluster.py -s True -j {output.juncfilelist} -r leafcutter/
        """

samplesForMetaPlots=["Sultan_rRNADeplete_Total", "Sultan_polyA_Total", "Sultan_rRNADepelete_nuclear", "Sultan_rRNADepelete_cytoplasmic", "19201_cheRNA_1", "NA19201_argonne"]
rule computeMatrix_metatranscript:
    input:
        UnstrandBw = expand("UnstrandedBigwigs/{sample}.bw", sample=samplesForMetaPlots),
        gtf=config["Human_ref"]["genome_gtf"],
    output:
        "Misc/computeMatrix.gz"
    threads: 2
    log:
        "logs/computeMatrix.log"
    shell:
        """
        computeMatrix scale-regions -S {input.UnstrandBw} -R {input.gtf} -a 1000 -b 1000 --metagene --missingDataAsZero -o {output} -p {threads} &> {log}
        """

rule computeMatrix_introns:
    input:
        UnstrandBw = expand("Alignments/SecondPass/{sample}/Unstranded.bw", sample=samplesForMetaPlots),
        regions="Misc/GencodeHg38_all_introns.corrected.bed"
    output:
        "Misc/computeMatrix.introns.gz"
    threads: 2
    log:
        "logs/computeMatrix_introns.log"
    shell:
        """
        computeMatrix scale-regions -S {input.UnstrandBw} -R {input.regions} --missingDataAsZero -o {output} -p {threads} &> {log}
        """

rule copyBigwigs:
    input:
        UnstrandBw = "Alignments/SecondPass/{sample}/Unstranded.bw",
    output:
        UnstrandedBigwigs="UnstrandedBigwigs/{sample}.bw"
    shell:
        "cp {input} {output}"

rule MakeTopHatJuncBedFiles:
    input:
        SJ = "Alignments/SecondPass/{sample}/SJ.out.tab",
    output:
        bed = "Junctions/{sample}.junctions.bed.gz",
    log:
        "logs/MakeTopHatJuncBedFiles/{sample}.log"
    shell:
        """
        bash scripts/SJ_to_junctions.sh {input.SJ} Junctions/{wildcards.sample}.junctions.bed
        gzip Junctions/{wildcards.sample}.junctions.bed
        """

# rule TopHatJuncBedFilesToBigBed:
#     input:
#         bed = "Alignments/SecondPass/{sample}/SJ.out.junctions.bed",
#         fai =config["Human_ref"]["genome_fasta"] + ".fai"
#     output:
#         Bigbed = "JunctionBigbeds/{sample}.bb"
#     shell:
#         """
#         ucsc-bedtobigbed {input.bed} {input.fai} {output.Bigbed}
#         """

        

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
