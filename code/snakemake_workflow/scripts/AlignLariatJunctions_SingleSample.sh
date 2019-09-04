#!/usr/bin/env bash
##run from results directory
##Output BED file is sorted for indexing and loading into IGV
##awk for converting SJ.out.tab to bed12 format
##based on code originally published by frymor at http://seqanswers.com/forums/showthread.php?t=62896
if [ "$1" == "-h" ]; then
  echo "USAGE: `basename $0`  <Unmapped.fastq.gz> <hisat2-index> <Output.bam> <Output.junctions.bed>"
  echo ""
  echo "Takes a single fastq file of antisense reads (R1 in most stranded RNA-seq protocols that use Illumina paired end sequencing) and attempts to align lariat junctions (inverted split alignments that cross the branchpoint and 5'ss) without using any filtering or priors based on splice site annotations. The general strategy is the same as implemented in Stepankiw et al (Nucleic Acids Res. 2015 Sep 30; 43(17): 8488–8501). Briefly, reads are split at each instance of GT (actually they are split at each instance of AC since this script assumes sequenced fastq reads are antisense sequence from the actual lariat) with a minimal overhang of 15 on each side of the split. Multiple split read pairs may derive from each original read. Then split read pairs are aligned as to the genome as paired end reads with hisat2 aligner with default parameters. Only the best alignment for each original read is kept. If the original fastq file is sense orientation (R2 in most stranded RNA-seq protocols), you should reverse complement it (fastx_reverse_complement) prior to running this script
  
  ARGUMENTS:

  <Unmapped.fastq.gz> Input fastq file
  <hisat2-index> Path to hisat2 index prefix (minus railing .X.ht2)
  <Output.bam> Output bam file of paired end alignments. Not sorted
  <Output.junctions.bed> Output bed file for easy visualization in IGV
  
  DEPENDENCIES:
  The following should be in \$PATH:
  
  hisat2
  fastx_reverse_complement"
  exit 0
fi


set -xe #Debug mode

FqIn=$1
hisatIndex=$2
bamOut=$3
junctionsOut=$4

temp=$(mktemp -dt "$(basename $0).XXXXXXXXXX")

# Split reads into inward facing read pairs
zcat $FqIn | paste -d"\t" - - - - | perl -F"\t" -lane 'while ($F[1] =~ /(?=[NAGCT]{13}(AC)[NAGCT]{15})/g) {if(defined($1)){printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $F[0],substr($F[1],$+[1]),$F[2],substr($F[3],$+[1]),$F[0],substr($F[1],0,$+[1]),$F[2],substr($F[3],0,$+[1])}}' | awk  -v OFS='\n' -F '\t' '{print $1,$2,$3,$4}' > $temp/myR1.fastq
zcat $FqIn | paste -d"\t" - - - - | perl -F"\t" -lane 'while ($F[1] =~ /(?=[NAGCT]{13}(AC)[NAGCT]{15})/g) {if(defined($1)){printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $F[0],substr($F[1],$+[1]),$F[2],substr($F[3],$+[1]),$F[0],substr($F[1],0,$+[1]),$F[2],substr($F[3],0,$+[1])}}' | awk  -v OFS='\n' -F '\t' '{print $5,$6,$7,$8}' | fastx_reverse_complement > $temp/myR2.fastq

# Align read pairs
hisat2 --no-spliced-alignment --no-softclip -I 20 -X 500000 --no-mixed --no-discordant -x $hisatIndex -1 $temp/myR1.fastq -2 $temp/myR2.fastq | samtools view -bh - > $temp/MyOutput.bam

# Keep only best alignments for each original read
cat <(samtools view -H $temp/MyOutput.bam) <(samtools view -F 4 $temp/MyOutput.bam | paste -d"~" - - | perl -F"\t" -lne '$_ =~ m/NM:i:(\d+)\t.+NM:i:(\d+)\t/; print $1+$2 . "\t$_"' | sort -k2,2 -k1n,1 | sort -u -k2,2 | perl -lne '$_ =~ m/\d+\t(.+)~(.+)/; print "$1\n$2"') | samtools view -bh > $bamOut


#Convert lariat junction PE alignments to junctions.bed similar to TopHat output
bedtools bamtobed -bedpe -mate1 -i $bamOut | awk -F'\t' -v OFS='\t' '$10=="-" {print $1,$2,$6,".",".",$10} $10=="+" {print $1,$5,$3,".",".",$10}' | sort | uniq -c | awk -v OFS='\t' '{print $2,$3,$4,".",$1,$7}' > $temp/BedPE.bed

awk \
{'if($6=="-") print ""$1"\t"$2"\t"$3"\tJUNC_"$1"_"$2"_"$3"_-\t"$5"\t-\t"$2"\t"$3"\t255,0,0\t2\t0,0\t0,"$3-$2; \
else \
if($6=="+") print ""$1"\t"$2"\t"$3"\tJUNC_"$1"_"$2"_"$3"_+\t"$5"\t+\t"$2"\t"$3"\t0,0,255\t2\t0,0\t0,"$3-$2'} \
$temp/BedPE.bed | sort -k 1,1 -k2,2n | awk 'BEGIN{print "track name=junctions description=\47TopHat junctions\47"} {print}' > $junctionsOut
