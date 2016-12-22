#!/usr/bin/env python

"""
date : 17 april 2013

Author : Julia di Iulio

Script used to monitor coverage (either the whole read, the start of the read (which represents PolII position),
or the end of the read) taking into account the number of times the read mapped (given by the NH:i: tag in the
alignment file): by giving a weight inversely proportional to the number of times the read mapped.
This script is highly inspired from http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html

use: python customCoverage.py stdin (alignment file in sam format) Experiment Name [1]

"""


import sys, HTSeq

#d stands for float as we will use the NH tag to weight the coverage
cvg = HTSeq.GenomicArray( "auto", stranded=True, typecode='d' ) #cvg stands for whole read coverage
stt = HTSeq.GenomicArray( "auto", stranded=True, typecode='d' ) #stt stands for start of reads
end = HTSeq.GenomicArray( "auto", stranded=True, typecode='d' ) #end stands for end of reads

alignment_file = HTSeq.SAM_Reader(sys.stdin)
for algt in alignment_file:
    idxM = [] # creates an empty list which is gonna contain the indexes of the "matched" segment(s) of the read

    for i in range(0,len(algt.cigar)): # builds the index positions of each matched segment of the read
        if algt.cigar[i].type == 'M' :
            idxM.append(i)

    if algt.aligned:
        for m in idxM:
            # adds weighted coverage throughout the matched segement(s) of the read
            cvg[ algt.cigar[m].ref_iv ] += float(1)/float( algt.optional_field( "NH" ) )

        # extract the start position of the read as a genomic interval
        stt_iv = HTSeq.GenomicInterval( algt.iv.chrom , algt.iv.start_d , algt.iv.start_d + 1, algt.iv.strand )
        if algt.iv.strand == '+':
            # extract the end position of the read as a genomic interval (if the read aligned on the "+" strand)
            end_iv = HTSeq.GenomicInterval( algt.iv.chrom , algt.iv.end - 1, algt.iv.end, algt.iv.strand )
        elif algt.iv.strand == '-':
            # extract the end position of the read as a genomic interval (if the read aligned on the "-" strand)
            end_iv = HTSeq.GenomicInterval( algt.iv.chrom , algt.iv.start, algt.iv.start + 1, algt.iv.strand )

        # adds weighted coverage at the position corresponding to the start of the read
        stt[ stt_iv ] += float(1)/float( algt.optional_field( "NH" ) )
        # adds weighted coverage at the position corresponding to the end of the read
        end[ end_iv ] += float(1)/float( algt.optional_field( "NH" ) )


cvg.write_bedgraph_file( sys.argv[1] + "_cov_neg.bedGraph", "+" )
cvg.write_bedgraph_file( sys.argv[1] + "_cov_pos.bedGraph", "-" )
stt.write_bedgraph_file( sys.argv[1] + "_stt_neg.bedGraph", "+" )
stt.write_bedgraph_file( sys.argv[1] + "_stt_pos.bedGraph", "-" )
end.write_bedgraph_file( sys.argv[1] + "_end_neg.bedGraph", "+" )
end.write_bedgraph_file( sys.argv[1] + "_end_pos.bedGraph", "-" )

