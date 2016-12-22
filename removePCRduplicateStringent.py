#!/usr/bin/env python

"""

Date : May 2nd 2014

Author : Julia di Iulio

script used to extract all reads starting at each position from a coverage file (only position with score > 0), to then
keep only 1 read per molecular barcode (to avoid getting signal from PCR duplicates). Ideally we should check the
distribution of the MB, but as it is unlikely to have very high coverage at a single position (due to the large size of the
human (un-)stable transcriptome), we are gonna chose the easy way and keep only 1 read per molecular barcode.
In the future, we could also take into consideration the insert size of the read, and the mismatches if any.

use : python removePCRduplicateStringent.py inBedGraph1 (with start coverage)                                                  [1]
                                            inTxt1  (file containing splicing intermediates (SI) positions with the folowing
                                                    nomenclature : chrX_strand_Start (where Start position are 1 based)        [2]
                                            inBAM1  (input BAM file with only unique alignment)                                [3]
                                            outBAM1 (containing no SI, and only non duplicated reads)                          [4]
                                            outBAM2 (containing only non duplicated reads, but still containing the SI if any) [5]
                                            outBedGraph1 (containing no SI, and only non duplicated reads)                          [6]
                                            outBedGraph2 (containing only non duplicated reads, but still containing the SI if any) [7]

"""

import sys, pysam, os, numpy, re


iBG1  = open(sys.argv[1], 'r')
iTXT1 = set(open(sys.argv[2], 'r').readlines())
iBAM1 = pysam.Samfile(sys.argv[3], 'rb')
oBAM1 = pysam.Samfile(sys.argv[4], 'wb', template=iBAM1)
oBAM2 = pysam.Samfile(sys.argv[5], 'wb', template=iBAM1)
oBG1 = open(sys.argv[6], 'w')
oBG2 = open(sys.argv[7], 'w')

if '_pos.bedGraph' in sys.argv[1]:
    std = 'pos'
elif '_neg.bedGraph' in sys.argv[1]:
    std = 'neg'

#iBG1.readline()

for bgLine in iBG1:
    chrom, st, ed , cov = bgLine.rstrip().split('\t')
    if float(cov) > 0:
        for j in range(0,(int(ed)-int(st))):
            MB = set()
            SI = 0
            woSI = 0
            wiSI = 0
            if chrom+"_"+std+"_"+str(int(st)+j+1)+"\n" in iTXT1:
                SI = 1
            for line in iBAM1.fetch(chrom, int(st)+j, int(st)+j+1):
                if (line.is_reverse and std=='pos' and line.aend-1 == int(st)) or \
                        (not line.is_reverse and std=='neg' and line.pos == int(st)):
                    mb = line.qname.split('_MolecularBarcode:')[1]
                    if mb not in MB:
                        MB.add(mb)
                        if SI == 0:
                            oBAM1.write(line)
                            woSI += 1/float([tag[1] for tag in line.tags if tag[0]=='NH'][0])
                        oBAM2.write(line)
                        wiSI += 1/float([tag[1] for tag in line.tags if tag[0]=='NH'][0])

            if SI == 0:
                oBG1.write("%s\t%s\t%s\t%.6f\n" % (chrom, str(int(st)+j), str(int(st)+j+1), float(woSI)))
            oBG2.write("%s\t%s\t%s\t%.6f\n" % (chrom, str(int(st)+j), str(int(st)+j+1), float(wiSI)))


iBG1.close()
iBAM1.close()
oBAM1.close()
oBAM2.close()
oBG1.close()
oBG2.close()
