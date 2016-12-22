#!/usr/bin/env python

"""

Date : May 9th 2014

Author : Julia di Iulio

Script used to filter (from the alignment obtained with the reads that contained the molecular barcode in the beginning of
their sequence) reads that contain mismatches in the 6 first nucleotides (nts). The rationale is that it is extremely
unlikely (probability = 1/4e6) that the molecular barcode (which is a random hexamer sequence) corresponds exactly to the
genomic sequence next to where the rest of the read aligned to. So among the reads that aligned (even though they still
contained the molecular barcode in the beginning of their sequence), we will extract the ones that did not have any
mismatches in the 6 first nts. Indeed, those reads are most likely due to mispriming event from the reverse transcriptase
primer. We will then remove those reads from the alignment obtained with the reads that did not contain the molecular
barcode in the beginning of their sequence.

use : python extractReadsWithMismatchesIn6FirstNct.py stdin stdout

"""

import sys, pysam, re

iBAM = pysam.Samfile("-", 'r') # reads from the standard input
oBAM = pysam.Samfile("-", 'w', template=iBAM) # output to the standard output

for line in iBAM:
    md=re.findall(r'\d+', [tag[1] for tag in line.tags if tag[0]=='MD'][0])
    if len(md) == 1 :       # if there are no mismatches
        oBAM.write(line)    # write the alignment into the output file
    else:
        if (not line.is_reverse) and (int(md[0]) >= 6): # if the first mismatch occurs after the 6th nt (from the 5' end)
            oBAM.write(line)                           # write the alignment into the output file
        elif (line.is_reverse) and (int(md[-1]) >= 6):  # same as above but for reads that align to the reverse strand
            oBAM.write(line)

# close files
iBAM.close()
oBAM.close()

