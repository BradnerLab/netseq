#!/usr/bin/env bash

FILES=`echo -e "X\nY\nZ"` # replace X,Y,Z with the names of the samples.
                          # (the names have to match the starting fastq files. f.ex: X.fastq, Y.fastq, Z.fastq)


inDir=/define/the/analyses/directory 
mkdir -p ${inDir}

# the starting fastq files have to be put in this directory, as well as all the python scripts and additional files provided

##############################################################################################################################
##############################################################################################################################
######################################################## molecular Barcode Extraction #########################################
# Description:
#                       removal of 6 first nucleotides (molecular barcode) from FastQ files
#                       creates 2 additional files containing the molecular barcode/ligation hexamer counts
#
# Input Name(s) :
#                       <Sample Name>.fastq
#
# Output Name(s) :
#                       <Sample Name>_noBarcode.fastq # contains read's sequence without the 6 first sequenced nucleotides
#                       <Sample Name>_barcodeDistribution.txt   # file containing barcode counts
#                       <Sample Name>_ligationDistribution.txt  # file containing ligation counts (3 last nucleotides of barcode
#                                                                 and 3 first nucleotides of the read without molecular barcode (MB))
#
# Program/commands:
#                       extractMolecularBarcode.py

cd ${inDir}
echo "Extracting Molecular Barcode... "
for f in $FILES
do
    echo "Doing file "$f
    python extractMolecularBarcode.py ${f}.fastq ${f}_noBarcode.fastq ${f}_barcodeDistribution.txt ${f}_split_ligationDistribution.txt
done

##############################################################################################################################
##############################################################################################################################
######################################################## Alignment  ##########################################################
# Description:
#                       STAR run on reads not containing/containing molecular barcodes ...
#                       Required for further Reverse Transcriptase (RT) bias removal
#
# Inputs:
#                       <Sample Name>_noBarcode.fastq
#                       <Sample Name>.fastq
#
# Outputs:
#                       <Sample Name>_Aligned.sortedByCoord.out.bam # alignment file obtained from reads without MB
#                       <Sample Name>_withBC_Aligned.sortedByCoord.out.bam # alignment file obtained from reads still containing MB
#
# Program:
#                       STAR
#


echo "Aligning reads containing or not Molecular Barcode ... "
program=/path/to/STAR # provide the path to STAR (https://github.com/alexdobin/STAR)
pFile=Parameters.in
idxDir=/path/to/STAR/index/ # provide the path to the STAR index directory

for f in $FILES
do
    echo "Doing file "$f
    reads=${f}_noBarcode.fastq
    mkdir -p STAR/${f}/
    pushd STAR/${f}/
    $program --genomeDir ${idxDir} --readFilesIn ${reads} --runThreadN 4 --parametersFiles ${pFile} --outFileNamePrefix ${f}_
    popd

    reads=${f}.fastq
    mkdir -p STAR/${f}_withBC/
    pushd STAR/${f}_withBC/
    $program --genomeDir ${idxDir} --readFilesIn ${reads} --runThreadN 4 --parametersFiles ${pFile} --outFileNamePrefix ${f}_withBC_
    popd

done



##############################################################################################################################
##############################################################################################################################
########################################## Reverse transdcriptase mispriming event removal ###################################
# Description:
#               to remove artefactual reads (arising from RT !!)
#               that are mapping even when we leave the barcode sequence !
#
# Inputs:
#               <Sample Name>_Aligned.sortedByCoord.out.bam
#               <Sample Name>_withBC_Aligned.sortedByCoord.out.bam
#
#
# Outputs:
#               <Sample Name>_noRTbias.sorted.bam  # alignment not containing reads arising from reverse transcriptase mispriming
#
#
# Program:
#               samtools
#               extractReadsWithMismatchesIn6FirstNct.py : a python script that extract reads arising from RT bias
#

echo "removing reads with RT bias ... "

program=/path/to/samtools # provide the path to samtools (http://samtools.sourceforge.net/)
script=extractReadsWithMismatchesIn6FirstNct.py

for f in $FILES
do
    echo "Doing file "$f
    pushd STAR/
    Bam1=${f}/${f}_Aligned.sortedByCoord.out.bam
    Bam2=${f}_withBC/${f}_withBC_Aligned.sortedByCoord.out.bam
    folder=${f}
    $program view -H ${Bam1} > ${f}/headers.sam ;
    $program view ${Bam1} | sed 's/_MolecularBarcode/\t_MB/' | sort -k1,1 > ${f}/temp.sam ;
    $program view -h ${Bam2} | python ${script} | cut -f 1 | sort -k1,1 | uniq > ${f}/rtBiasnames.txt ;
    join -v 1 ${f}/temp.sam ${f}/rtBiasnames.txt | sed -e 's/ _MB/_MolecularBarcode/' -e 's/ /\t/g' | \
        cat ${f}/headers.sam - | $program view -Sb -o ${f}/temp_noRTbias.bam - ;
    $program sort ${f}/temp_noRTbias.bam ${f}/${f}_noRTbias.sorted ;
    $program index ${f}/${f}_noRTbias.sorted ; rm ${f}/temp.sam ; rm ${f}/temp_noRTbias.bam;
    $program view -b -h -q 50 ${f}/${f}_noRTbias.sorted.bam > ${f}/${f}_noRTbias.sorted_uniq.bam ;
    $program index ${f}/${f}_noRTbias.sorted_uniq.bam

    popd
done


##############################################################################################################################
##############################################################################################################################
########################################## Extracting coverage tract #########################################################
# Description:
#                       Monitors whole read coverage and pol II position
#                       Generates bedgraph file
# Inputs:
#                       <Sample Name>_noRTbias.sorted.bam
#
# Outputs:
#                       <SampleName>_noRTbias_cov_???.bedGraph (contains the coverage information of the whole read)
#                       <SampleName>_noRTbias_stt_???.bedGraph (monitors the start (stt) of the read (= PolII position))
#                       <SampleName>_noRTbias_end_???.bedGraph (monitors the end of the read)
#                       The ??? stands either for "pos" or "neg";
#                       The *_pos.bedGraph files contain info from polII transcribing in the pos strand,
#                       and inversely for the *_neg.bedGraph files
#
## Program:
#                       samtools
#                       customCoverage.py : a python script that extract reads arising from RT bias. This script uses the HTSeq
#                                           package written by Simon Anders http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html
#

echo "monitoring Coverage and PolII position ... "

program=/path/to/samtools # provide the path to samtools (http://samtools.sourceforge.net/)
script=customCoverage.py

for f in $FILES
do
    echo "Doing file "$f
    outDir=CoverageSTAR/${f}
    mkdir -p $outDir
    Bam=STAR/${f}/${f}_noRTbias.sorted.bam
    $program view -q 50 ${Bam} | sed 's/\tjM.*$//' | python $script ${outDir}/${f}_uniq_noRTbias_
    $program view ${Bam} | sed 's/\tjM.*$//' | python $script ${outDir}/${f}_noRTbias_

done


##############################################################################################################################
##############################################################################################################################
########################################## Remove PCR duplicates #############################################################
# Description:
#               to remove artefactual reads (arising from PCR duplicates)
#               that have the same genomic position and molecular barcode than another read
#
# Inputs:
#               <SampleName>_noRTbias.sorted.bam
#               <SampleName>_noRTbias_stt_???.bedGraph
#               SI_coordinates_1based.txt : file containing the splicing intermediates (SI) position (chromosome_strand_position)
#
#
# Outputs:
#               <SampleName>_noRTbias_noPCRdup_noSI_???.bam (containing no SI, and only non duplicated reads)
#               <SampleName>_noRTbias_noPCRdup_???.bam (containing only non duplicated reads, but still contains the SI if any)
#               <SampleName>_noRTbias_noPCRdup_noSI_stt_???.bedGraph (containing no SI, and only non duplicated reads)
#               <SampleName>_noRTbias_noPCRdup_stt_???.bedGraph (containing only non dup. reads, but still contains the SI if any)
#                       The ??? stands either for "pos" or "neg";
#                       The *_pos.bedGraph files contain info from polII transcribing in the pos strand,
#                       and inversely for the *_neg.bedGraph files
#
# Program:
#               samtools
#               removePCRduplicateStringent.py : a python script that extract reads arising from PCR duplication events
#

echo "removing PCR duplicates and Splicing intermediates from bam and bedGraph files ... "
program=/path/to/samtools # provide the path to samtools (http://samtools.sourceforge.net/)
script=removePCRduplicateStringent.py
SIfile=SI_coordinates_1based.txt

for f in $FILES
do
    echo "Doing file " $f
    inCovDir=CoverageSTAR/${f}
    inBamDir=STAR/${f}/

    iBamUniq=${inBamDir}/${f}_noRTbias.sorted_uniq.bam
    iCovUniqP=${inCovDir}/${f}_uniq_noRTbias_stt_pos.bedGraph
    iCovUniqN=${inCovDir}/${f}_uniq_noRTbias_stt_neg.bedGraph
    oCovUniqP=${inCovDir}/${f}_uniq_noRTbias_noPCRdup_stt_pos.bedGraph
    oCovUniqN=${inCovDir}/${f}_uniq_noRTbias_noPCRdup_stt_neg.bedGraph
    oCovUniqPnoSI=${inCovDir}/${f}_uniq_noRTbias_noPCRdup_noSI_stt_pos.bedGraph
    oCovUniqNnoSI=${inCovDir}/${f}_uniq_noRTbias_noPCRdup_noSI_stt_neg.bedGraph
    oBamUniqP=${inBamDir}/${f}_noRTbias_noPCRdup_uniq_pos.bam
    oBamUniqN=${inBamDir}/${f}_noRTbias_noPCRdup_uniq_neg.bam
    oBamUniqPnoSI=${inBamDir}/${f}_noRTbias_noPCRdup_noSI_uniq_pos.bam
    oBamUniqNnoSI=${inBamDir}/${f}_noRTbias_noPCRdup_noSI_uniq_neg.bam

    python $script $iCovUniqP $SIfile $iBamUniq $oBamUniqPnoSI $oBamUniqP $oCovUniqPnoSI $oCovUniqP

    python $script $iCovUniqN $SIfile $iBamUniq $oBamUniqNnoSI $oBamUniqN $oCovUniqNnoSI $oCovUniqN

    sed 's/$/\t+/' $oCovUniqPnoSI > ${oCovUniqPnoSI}.temp ; sed 's/$/\t-/' $oCovUniqNnoSI > ${oCovUniqNnoSI}.temp ;                          cat ${oCovUniqPnoSI}.temp ${oCovUniqNnoSI}.temp | sort-bed - > ${oCovUniqPnoSI%pos.bedGraph}all.bed ;                                    rm ${oCovUniqPnoSI}.temp ; rm ${oCovUniqNnoSI}.temp

    sed 's/$/\t+/' $oCovUniqP > ${oCovUniqP}.temp ; sed 's/$/\t-/' $oCovUniqN > ${oCovUniqN}.temp ;                                          cat ${oCovUniqP}.temp ${oCovUniqN}.temp | sort-bed - > ${oCovUniqP%pos.bedGraph}all.bed ;                                                rm ${oCovUniqP}.temp ; rm ${oCovUniqN}.temp

    $program sort ${oBamUniqNnoSI} ${oBamUniqNnoSI%.bam}.sorted ; $program sort ${oBamUniqPnoSI} ${oBamUniqPnoSI%.bam}.sorted ;              $program merge -f -r ${oBamUniqPnoSI%_pos.bam}.bam ${oBamUniqNnoSI%.bam}.sorted.bam ${oBamUniqPnoSI%.bam}.sorted.bam;                    $program index ${oBamUniqPnoSI%_pos.bam}.bam ;                                                                                           $program index ${oBamUniqNnoSI%.bam}.sorted.bam ; $program index ${oBamUniqPnoSI%.bam}.sorted.bam

    $program sort ${oBamUniqN} ${oBamUniqN%.bam}.sorted ; $program sort ${oBamUniqP} ${oBamUniqP%.bam}.sorted ;                              $program merge -f -r ${oBamUniqP%_pos.bam}.bam ${oBamUniqN%.bam}.sorted.bam ${oBamUniqP%.bam}.sorted.bam;                                $program index ${oBamUniqP%_pos.bam}.bam ;$program index ${oBamUniqN%.bam}.sorted.bam ; $program index ${oBamUniqP%.bam}.sorted.bam

done


##############################################################################################################################
##############################################################################################################################
########################################## Extract read counts per gene ######################################################
# Description:
#                       Extract counts (number of pol II) per feature
#                       Generates a table with the feature ID (column 1) and the counts (column 2).
#                       The script uses a modified version of HTSeq-count (written by Simon Anders), that enables to monitor
#                       reads mapping multiple times (giving a weight inversely proportional to the number of times the read
#                       mapped), reads mapping to overlapping features (by randomly attributin the read to 1 of the feature)
#                       and modified to monitor only the start of the read (which correpsonds to PolII position)
#
# Inputs:
#                       <Sample Name>_noRTbias.sorted.bam
#                       gencode.v16_intronExonSIannot.gtf : gtf file containing the feature annotation; in this case it contains
#                                                           specific gene_id for exons, introns and splicing intermediates position
#                                                           of the same gene (the suffix for introns is ".Intron" and the suffix for
#                                                           splicing intermediate positions is ".SpliceInterm"). Depending on the need
#                                                           for the analysis you may want to remove all the suffixes so that all reads
#                                                           mapping to a given gene (no matter whether it is in the exons, intron, etc)
#                                                           will be counted together.
#
# Outputs:
#                       <SampleName>_noRTbias_stt.txt          (contains weighted counts from all the reads,
#                                                                  for each feature provided in the gtf file)
#                       <SampleName>_uniq_noRTbias_stt.txt     (contains counts only from the uniquely aligned reads,
#                                                                  for each feature provided in the gtf file)
#                       <SampleName>_IE_uniq_noRTbias_stt.txt  (contains counts only from the uniquely aligned reads,
#                                                                  for total intronic (column 2) and exonic (column 3) regions
#                                                                  of the intron containing genes provided in the gtf file)
#
# Program:
#                       samtools
#                       htseq-count_custom.py : a python script that extract read counts per feature
#                                               (written by Simon Anders, and customized as explained above)

echo "extract PolII counts in every feature provided by the gtf file ... "
program=/path/to/samtools # provide the path to the samtools (http://samtools.sourceforge.net/)
script=htseq-count_custom.py
GTFfile=gencode.v16_intronExonSIannot.gtf

for f in $FILES
do
    echo "Doing file "$f
    cd ${inDir}
    outDir=HTSeqSTAR/${f}
    mkdir -p $outDir
    Bam=STAR/${f}/${f}_noRTbias.sorted.bam

    $program view ${Bam} | sed 's/\tjM.*$//' | python $script --mode=intersection-nonempty --stranded=reverse                                --type=exon --idattr=gene_id - ${GTFfile} > ${f}_noRTbias_stt.txt

    $program view -q 50 ${Bam} | sed 's/\tjM.*$//' | python $script --mode=intersection-nonempty --stranded=reverse                          --type=exon --idattr=gene_id - ${GTFfile} > ${f}_uniq_noRTbias_stt.txt

## this is used for analysis of intron versus exon, you then would need to normalize the read by the total respective size of
## the features in the gtf file, to be able to compare them.
#    grep Intron ${f}_uniq_noRTbias_stt.txt | sed 's/\.Intron\t/\t/g' | sort -k1,1 > ${f}_uniq_noRTbias_stt.txt.Intemp
#    grep -v Intron ${f}_uniq_noRTbias_stt.txt | sort -k1,1 > ${f}_uniq_noRTbias_stt.txt.Extemp
#    join -1 1 -2 1 ${f}_uniq_noRTbias_stt.txt.Intemp ${f}_uniq${Osuf1}stt.${nb}.txt.Extemp | sed 's/ /\t/g' \
#        > ${f}_IE_uniq_uniq_noRTbias_stt.txt ; rm ${f}_uniq_noRTbias_stt.txt.??temp
#
done

# Note: to plot the total number of polII per gene, you will need to add the exons and introns counts.


##############################################################################################################################
##############################################################################################################################
###################################### Extract traveling ratio per gene ######################################################
# Description:
#                       Extract polII per kb in the promoter proximal (PP) and gene body (GB) regions, as well as their ratio
#                       Generates a table with the gene ID (column 1), polII/kb in PP (column 2),  polII/kb in GB (column 3)
#                       and travel ratio, which is the ratio of column 2 and column 3 (column 4).
#
# Inputs:
#                       <Sample Name>_uniq_noRTbias_noPCRdup_noSI_stt_ (prefix of coverage files to be used)
#                       PP_coord.txt : coordinates of promoter proximal regions, with the following structure:
#                                      column 1: geneID; column 2: chromosome; column 3: start; column 4: end; column 6: strand
#                                      In the list provided, only the pol II genes where we were able to identify a signal for promoter
#                                      proximal pausing (PPP) were used. And the PP corresponds to the PPP position -80bp to PPP +250bp
#                                      which corresponds approximately to the used TSS -30bp to TSS +300bp.
#                       GB_coord.txt : coordinates of gene body region, with the following structure:
#                                      column 1: geneID; column 2: chromosome; column 3: start; column 4: end; column 6: strand
#                                      In the list provided, only the pol II genes where we were able to identify a signal for promoter
#                                      proximal pausing (PPP) were used. And the GB corresponds to the PPP +250bp up to the pA site.
#                                      which corresponds approximately to the used TSS +300bp up to the pA site.
#                       Note : different coordinates for PP and GB can be provided by the user, depending on which set of genes the
#                              user is interested (f.ex. only protein coding, non-overlaping, expressed, etc)
#
# Outputs:
#                       <SampleName>_TravelRatio.txt   (contains polII/kb in both PP and GB regions, and their ratio)
#                                                      structure of the file: gene ID (column 1), polII/kb in PP (column 2),
#                                                      polII/kb in GB (column 3) and travel ratio (column 4)
#
# Program:
#                       bedops (https://bedops.readthedocs.io/en/latest/, which will be called from within the r script)
#                       travelRatio.r : a R script that generates table with travel ratios.

echo "extract travel ratio in every gene provided ... "
program=/path/to/Rscript # provide the path to R (https://cran.r-project.org/)
script=travelRatio.r
PPcoord=PP_coord.txt
GBcoord=GB_coord.txt

for f in $FILES
do
    echo "Doing file "$f
    cd ${inDir}
    outDir=travelRatio/${f}/
    mkdir -p $outDir
    CoveragePrefix=CoverageSTAR/${f}/${f}_uniq_noRTbias_noPCRdup_noSI_stt_

    $program ${script} ${f} ${PPcoord} ${GBcoord} ${CoveragePrefix} ${outDir}
done

