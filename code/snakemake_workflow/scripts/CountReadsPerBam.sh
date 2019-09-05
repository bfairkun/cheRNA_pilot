#!/usr/bin/env bash

#debug to stderr
set -xe

CountFile=$1
printf "Filename\tReadCount\n" > $CountFile
shift;
for filename in "$@" ; do
    #process item
    printf "$filename\t" >> $CountFile
    printf "$(samtools view -c -F256 $filename)\n" >> $CountFile
done
