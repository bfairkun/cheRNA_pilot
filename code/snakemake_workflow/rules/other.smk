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

rule FilterForHighSupportTranscripts:
    input:
        "Misc/GencodeOriginalFiles/gencode.v31.annotation.gtf.gz"
    output:
        ""


