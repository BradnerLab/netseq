#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
s = args[1]
PPPcoord = args[2]
GBcoord  = args[3]
inDirCov = args[4]
outDir   = args[5]

bashCommand='bedextract $1 $2 '

PPPcoordForTR = read.table(PPPcoord, sep='\t', row.names=1, col.names=c('id', 'chrom', 'st', 'ed', 'std'))
GBcoordForTR  = read.table(GBcoord, sep='\t', row.names=1, col.names=c('id', 'chrom', 'st', 'ed', 'std'))

sizeGBs = GBcoordForTR$ed - GBcoordForTR$st
names = rownames(GBcoordForTR)[which(sizeGBs>0)]

print(s)
mt = data.frame(row.names=names, PPP=rep(NA, length(names)), GB=rep(NA, length(names)))
for (j in 1:length(names))
{
    if (j%%100 == 1){print(paste(j))}
    # j is printed in order to have an idea of how slow/fast the algorythm is.
    geneID = names[j]
    tmp = PPPcoordForTR[geneID,c(1:3)]
    tmp2 = GBcoordForTR[geneID,c(1:3)]
    write.table(tmp, paste('tmpCoord_',s,'.txt', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
    write.table(tmp2, paste('tmp2Coord_',s,'.txt', sep=''), sep="\t", quote=F, row.names=F, col.names=F)

    if (PPPcoordForTR[geneID,4]=='+'){std="pos"}else{std="neg"}

    covFile = paste(inDirCov,std,'.bedGraph', sep='')
    temp = read.table(pipe(paste(bashCommand, covFile, ' tmpCoord_',s,'.txt', sep='')), sep="\t", header=FALSE,
                          col.names=c('chr', 'start', 'end', 'Cov'), colClasses=c('character', rep('numeric',3)))
    temp2 = read.table(pipe(paste(bashCommand, covFile, ' tmp2Coord_',s,'.txt', sep='')), sep="\t", header=FALSE,
                          col.names=c('chr', 'start', 'end', 'Cov'), colClasses=c('character', rep('numeric',3)))
    if (nrow(temp)>0){
        if (temp[1,'start'] < tmp[,'st']){temp[1,'start'] = tmp[,'st']}
        if (temp[nrow(temp),'end'] > tmp[,'ed']){temp[nrow(temp),'end'] = tmp[,'ed']}
        sizeTemp = tmp[,'ed'] -  tmp[,'st']
        temp = subset(temp, Cov>0)
        mt[geneID, 'PPP'] = sum(rep(temp[,'Cov'], temp[,'end'] - temp[,'start']))/sizeTemp
    }else{
        mt[geneID, 'PPP'] = 0
    }

    if (nrow(temp2)>0){
        if (temp2[1,'start'] < tmp2[,'st']){temp2[1,'start'] = tmp2[,'st']}
        if (temp2[nrow(temp2),'end'] > tmp2[,'ed']){temp2[nrow(temp2),'end'] = tmp2[,'ed']}
        sizeTemp2 = tmp2[,'ed'] -  tmp2[,'st']
        temp2 = subset(temp2, Cov>0)
        mt[geneID, 'GB'] = sum(rep(temp2[,'Cov'], temp2[,'end'] - temp2[,'start']))/sizeTemp2
    }else{
        mt[geneID, 'GB'] = 0
    }
}

mt$geneID = rownames(mt)
mt$TR = mt$PPP / mt$GB
mt = mt[,c('geneID', 'PPP', 'GB', 'TR')]
colnames(mt)[2] = 'PP'
write.table(mt, file=paste(outDir,s,'_TravelRatio.txt', sep=''), sep='\t', col.names=T, quote=F, row.names=F)



