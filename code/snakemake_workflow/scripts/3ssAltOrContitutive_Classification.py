#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 19 15:09:19 2018
@author: benfair
"""

import sys
import re
from signal import signal, SIGPIPE, SIG_DFL
import argparse
import sys
import collections

parser = argparse.ArgumentParser(description="""
        This command line script will read in a stranded 6column bed file of annotated introns (5' and 3' splice site pairings) and regard those as annotated introns. From there it will annotate an 6 column stranded bed file of input as to whether the introns are annotated, alt5'ss, alt3'ss, new intron (neither 5' or 3' splice sites are annotated), or new splice site pairing (both 5' and 3' splice sites are annotated but the pairing has not been annotated, as is the case in a novel exon skipping event).
        """)
parser.add_argument('-I','--FileIn',  default="-", help="REQUIRED. Stranded 6 column bed file of introns that need alternative splice site classification. Use 'stdin' or '-' to read from stdin", required=True)
parser.add_argument('-A','--AlternativeBed',  help="REQUIRED. Stranded 6 column bed file of annotated introns", required=True)
parser.add_argument('-C', '--ConstitutiveBed',  default="stdout", help="REQUIRED. Output file.", required=True)
parser.add_argument("--quiet", help="OPTIONAL. quiet the output verbosity", action="store_true")
args = parser.parse_args()


if args.FileIn in ("stdin", "-"):
    fhIn = sys.stdin
else:
    fhIn = open(args.FileIn, 'rU')


# fh_AnnotatedSpliceJunctionsBedFile= open(args.AnnotationsBed, 'rU')
signal(SIGPIPE, SIG_DFL)

def GetSpliceDonorAndAcceptorCoordinates(chrom, start, stop, strand):
    '''
    function returns a tuple pair of (SpliceDonor, SpliceAcceptor) where SpliceDonor and SpliceAcceptor are themselves a tuple of (chromosome, coordinate, strand).
    '''
    if strand=='+':
        SpliceDonorCoord=(chrom, start, strand)
        SpliceAcceptorCoord=(chrom, stop, strand)
    elif strand=='-':
        SpliceDonorCoord=(chrom, stop, strand)
        SpliceAcceptorCoord=(chrom, start, strand)
    else:
        print("Failed... Need strand information coded as + or -. This information should be stored in field6 of bed file")
    return (SpliceDonorCoord, SpliceAcceptorCoord)


# AnnotatedSpliceJunctionsBedFile='/Users/benfair/Downloads/TargetIntrons.sorted.bed'
# InputSpliceJunctionsBedFile='/Users/benfair/Downloads/SI_Tables_AlternativeCounts.bed.txt'


#Initialize sets of annotated splice sites
AnnotatedSpliceDonors=collections.Counter()
AnnotatedSpliceAcceptors=collections.Counter()
AnnotatedSplicePairs=set()

for line in fhIn:
    l=line.strip('\n').split('\t')
    chrom, start, stop, name, score, strand = l[0:6]
    Donor,Acceptor = GetSpliceDonorAndAcceptorCoordinates(chrom, start, stop, strand)

    #Add to sets and counters
    AnnotatedSpliceDonors[Donor] += 1
    AnnotatedSpliceAcceptors[Acceptor] += 1
    AnnotatedSplicePairs.add((Donor,Acceptor))
fhIn.close()

for i, sspair in enumerate(AnnotatedSplicePairs):
    if i<10:
        # print(sspair)
        pass

for x in list(AnnotatedSpliceDonors)[0:100]:
    print(x, AnnotatedSpliceDonors[x])
