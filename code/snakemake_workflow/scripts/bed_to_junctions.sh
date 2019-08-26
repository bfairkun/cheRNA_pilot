#!/bin/bash
##run from results directory
##Output BED file is sorted for indexing and loading into IGV
##awk for converting SJ.out.tab to bed12 format
##based on code originally published by frymor at http://seqanswers.com/forums/showthread.php?t=62896

if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0`  <FileIn.bed> <FileOut.junctions.bed>"
  echo ""
  echo "Takes a stranded bed file of introns, with the number of junction reads in column5, and outputs a bed track like the one output by TopHat for viewing in IGV"
  exit 0
fi

bed=$1
junctions=$2

echo ${bed}
echo "Converting..."
awk \
{'if($6=="-") print ""$1"\t"$2"\t"$3"\tJUNC_"$1"_"$2"_"$3"_-\t"$5"\t-\t"$2"\t"$3"\t255,0,0\t2\t0,0\t0,"$3-$2; \
else \
if($6=="+") print ""$1"\t"$2"\t"$3"\tJUNC_"$1"_"$2"_"$3"_+\t"$5"\t+\t"$2"\t"$3"\t0,0,255\t2\t0,0\t0,"$3-$2'} \
${bed} | sort -k 1,1 -k2,2n | awk 'BEGIN{print "track name=junctions description=\47TopHat junctions\47"} {print}' > $junctions
echo "Complete"
