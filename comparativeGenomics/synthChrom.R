library(GenomicRanges)
library(dplyr)

rm(list = ls())

load("Desktop/RTN_domains/R_objects/mappedGaps/mm10.hg19.netData.RData")



switchGenome <- function(queGenome.gr){
  GRanges(mcols(queGenome.gr)$queRanges, 
          queRanges = granges(queGenome.gr, use.mcols = FALSE), sData = mcols(queGenome.gr)$sData)
}


## get the loss bases 

ol <- findOverlaps(queGap.gr,queAncDna.gr)
queGapAnc.gr <- pintersect(queGap.gr[queryHits(ol)], queAncDna.gr[subjectHits(ol)], drop.nohit.ranges=TRUE)

ol <- findOverlaps(queGap.gr,
                   GenomicRanges::setdiff(queGap.gr,queAncDna.gr, ignore.strand = TRUE))
queGapNonAnc.gr <- pintersect(queGap.gr[queryHits(ol)], 
                              GenomicRanges::setdiff(queGap.gr,queAncDna.gr,ignore.strand = TRUE)[subjectHits(ol)], 
                              drop.nohit.ranges=TRUE)

# correct seperation of ancestral and non ancestral
all(!overlapsAny(queGapAnc.gr, queGapNonAnc.gr))



ol <- findOverlaps(refGap.gr,refAncDna.gr)
refGapAnc.gr <- pintersect(refGap.gr[queryHits(ol)], refAncDna.gr[subjectHits(ol)], drop.nohit.ranges=TRUE)

ol <- findOverlaps(refGap.gr,
                   GenomicRanges::setdiff(refGap.gr,refAncDna.gr, ignore.strand = TRUE))
refGapNonAnc.gr <- pintersect(refGap.gr[queryHits(ol)], 
                              GenomicRanges::setdiff(refGap.gr,refAncDna.gr, ignore.strand = TRUE)[subjectHits(ol)], 
                              drop.nohit.ranges=TRUE)

# correct seperation of ancestral and non ancestral
all(!overlapsAny(refGapAnc.gr, refGapNonAnc.gr))
    

refMissingGaps.gr <- GenomicRanges::setdiff(refFillGaps.gr, refGap.gr, ignore.strand = TRUE)
all(!overlapsAny(refMissingGaps.gr, refGap.gr))



# set up genomes

queIns.gr <- switchGenome(queGapNonAnc.gr)
queIns.gr <- sort(sortSeqlevels(queIns.gr))
queIns.gr$type = "queIns"

queDel.gr <- sort(sortSeqlevels(refGapAnc.gr))
queDel.gr$type = "queDel"

refDel.gr <- switchGenome(queGapAnc.gr)
refDel.gr <- sort(sortSeqlevels(refDel.gr))
refDel.gr$type = "refDel"

refIns.gr <- sort(sortSeqlevels(refGapNonAnc.gr))
refIns.gr$type = "refIns"


# annotate 
mcols(refMissingGaps.gr)$queRanges <- GRanges(seqnames = "empty", ranges = IRanges(start = 1, end = 1))
mcols(refMissingGaps.gr)$sData = "*"
mcols(refMissingGaps.gr)$type = "missingGap"


mcols(refSeqGaps.gr)$queRanges <- GRanges(seqnames = "empty", ranges = IRanges(start = 1, end = 1))
mcols(refSeqGaps.gr)$sData = "*"
mcols(refSeqGaps.gr)$type = "seqGap"


refAll <- c(refDel.gr, refIns.gr,
            #queDel.gr, queIns.gr, 
            refMissingGaps.gr, refSeqGaps.gr)



#refAll <- c(refDel.gr, refIns.gr,queDel.gr, queIns.gr, refSeqGaps.gr)


refAll <- sort(refAll, by = ~ seqnames + start)



refSynth <- refAll

seqlengths(refSynth) = NA 

end(refSynth[mcols(refSynth)$type == "refDel"]) = start(refSynth[mcols(refSynth)$type == "refDel"]) + width(refSynth[mcols(refSynth)$type == "refDel"]$queRanges)  -1

#end(refSynth[mcols(refSynth)$type == "queIns"]) = start(refSynth[mcols(refSynth)$type == "queIns"]) + width(refSynth[mcols(refSynth)$type == "queIns"]$queRanges)  -1


pDiff <- end(refSynth[1:(length(refSynth)-1)]) - start(refSynth[2:length(refSynth)]) +1
pDiff[pDiff < 0] = 0

synthChrInfo <- refChrInfo
for(chrChoice in as.character(unique(seqnames(refSynth)))){
  cSum <- c(0,cumsum(pDiff[as.character(seqnames(refSynth)) == chrChoice]))
  refSynth[as.character(seqnames(refSynth)) == chrChoice] <- shift(refSynth[as.character(seqnames(refSynth)) == chrChoice], cSum[1:(length(cSum) - 1)])
  addedBases <- pDiff[as.character(seqnames(refSynth)) == chrChoice]
  addedBases <- addedBases[1:(length(addedBases)-1)]
  synthChrInfo[synthChrInfo$chrom == chrChoice,]$size <- synthChrInfo[synthChrInfo$chrom == chrChoice,]$size + sum(addedBases, na.rm = TRUE)
}



cov <- coverage(refSynth)
sl <- IRanges::slice(cov, lower = 2)
all(unlist(lapply(sl, length)) == 0)



synthGenome <- GRanges(seqnames = Rle(synthChrInfo$chrom), 
                       ranges = IRanges(width = synthChrInfo$size, end = synthChrInfo$size))

synthBin.gr <- unlist(slidingWindows(synthGenome, width = 100000, step = 100000))

#synthDel <- refSynth[mcols(refSynth)$type == "del"]

ol <- findOverlaps(refSynth, synthBin.gr)
pInt <- pintersect(refSynth[queryHits(ol)], synthBin.gr[subjectHits(ol)])



df <- data.frame(width = width(pInt), type = pInt$type, binNo = subjectHits(ol))

sumDF <- summarise(group_by(df, type, binNo), total = sum(width))

mcols(synthBin.gr)$refIns <- 0
mcols(synthBin.gr)$refDel <- 0
mcols(synthBin.gr)$queIns <- 0
mcols(synthBin.gr)$queDel <- 0
mcols(synthBin.gr)$missingGap <- 0
mcols(synthBin.gr)$seqGap <- 0

mcols(synthBin.gr[sumDF$binNo[sumDF$type == "refIns"]])$refIns = sumDF$total[sumDF$type == "refIns"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "refDel"]])$refDel = sumDF$total[sumDF$type == "refDel"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "queIns"]])$queIns = sumDF$total[sumDF$type == "queIns"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "queDel"]])$queDel = sumDF$total[sumDF$type == "queDel"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "missingGap"]])$missingGap = sumDF$total[sumDF$type == "missingGap"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "seqGap"]])$seqGap = sumDF$total[sumDF$type == "seqGap"]





synthBin.gr[100:110]


df <- data.frame(synthBin.gr)
df[rowSums(df[,c("missingGap", "seqGap")]) > (.05 * 1e5),] <- NA


layout(c(1,2,3,4))
par(mar=c(1,1,1,1))
chrChoice = "chr1"

plot(df$start[df$seqnames == chrChoice],
     df$refDel[df$seqnames == chrChoice], 
     pch = 16, cex = .5, ylim = c(0,1e5))
legend("topright", "refDel", bty = "n", cex = 1.5)

plot(df$start[df$seqnames == chrChoice],
     df$refIns[df$seqnames == chrChoice], 
     pch = 16, cex = .5, ylim = c(0,1e5))
legend("topright", "refIns", bty = "n", cex = 1.5)

plot(df$start[df$seqnames == chrChoice],
     df$queDel[df$seqnames == chrChoice], 
     pch = 16, cex = .5, ylim = c(0,1e5))
legend("topright", "queDel", bty = "n", cex = 1.5)

plot(df$start[df$seqnames == chrChoice],
     df$queIns[df$seqnames == chrChoice], 
     pch = 16, cex = .5, ylim = c(0,1e5))
legend("topright", "queIns", bty = "n", cex = 1.5)

# how to get the original coordinates back

# maybe do a bootstrapping approach or something 
# I'll just shuffle all of the stuff around and look at the expected neighbor values

hist(df$refIns, breaks= 100)
hist(df$refDel, breaks= 100)
hist(df$queIns, breaks= 100)
hist(df$queDel, breaks= 100)


hist(rowSums(data.frame(mcols(synthBin.gr))), breaks = 100)

# so lets get other things involved

hist(log10(width(refMissingGaps.gr)), breaks = 100)
hist(log10(width(refSeqGaps.gr)), breaks = 100)

hist(synthBin.gr$missingGap, breaks = 100)

sum(synthBin.gr$missingGap)
sum(as.numeric(width(refMissingGaps.gr)), na.rm = TRUE)

refSeqGaps.gr

gMeanRM = NULL
for(chrChoice in synthChrInfo$chrom){
  chrMean <- rollapply(data= c(rep(NA, 15), df[df$seqnames == chrChoice,]$refIns,rep(NA, 15)), 
                     width = 31, na.rm = TRUE, FUN = mean)
  gMeanRM <- c(gMeanRM, chrMean)
}
gMean <- median(gMean, na.rm = TRUE)
hist(gMeanRM, breaks = 100)
abline(v = gMean)

gSdRM = NULL
for(chrChoice in synthChrInfo$chrom){
  chrSd <- rollapply(data= c(rep(NA, 15), df[df$seqnames == chrChoice,]$refIns,rep(NA, 15)), 
                       width = 31, na.rm = TRUE, FUN = sd)
  gSdRM <- c(gSdRM, chrSd)
}
gSd <- median(gSdRM, na.rm = TRUE)
hist(gSdRM, breaks = 100)
abline(v = gSd)

# maybe use bootstrap to see how often we get certain values

# If I do bootstrapping on our rolling mean, maybe then I can get some 

# the senstivity to small values 

Msamp <- replicate(sample(gMeanRM[!is.nan(gMeanRM)], replace = TRUE, size = length(gMeanRM[!is.nan(gMeanRM)])),n = 100)


hist(Msamp, breaks = 100)




pdf(file = "~/Desktop/mm10.pdf",height = 6,width = 12 ,onefile = TRUE)
for(chrChoice in synthChrInfo$chrom[-grep("_", synthChrInfo$chrom)]){
  layout(1:4)
  par(mar = c(2,3,0,1), oma = c(5,5,5,5))
  plot(df[df$seqnames == chrChoice,]$start,
       df[df$seqnames == chrChoice,]$refIns, 
       type = "p", pch = 16, cex = .7,ylim = c(0,100000))
  rm <- rollapply(data= c(rep(NA, 15), df[df$seqnames == chrChoice,]$refIns,rep(NA, 15)), 
                  width = 31, na.rm = TRUE, FUN = mean)
  lines(df[df$seqnames == chrChoice,]$start, rm, col = 2)
  mtext(paste(genomes[1],"gain (bp)"),side = 2,line = 3, cex = .8)
  
  plot(df[df$seqnames == chrChoice,]$start,
       df[df$seqnames == chrChoice,]$refDel, 
       type = "p", pch = 16, cex = .7,ylim = c(0,100000))
  rm <- rollapply(data= c(rep(NA, 15), df[df$seqnames == chrChoice,]$refDel,rep(NA, 15)), 
                  width = 31, na.rm = TRUE, FUN = mean)
  lines(df[df$seqnames == chrChoice,]$start, rm, col = 2)
  mtext(paste(genomes[1],"loss (bp)"),side = 2,line = 3, cex = .8)
  
  plot(df[df$seqnames == chrChoice,]$start,
       df[df$seqnames == chrChoice,]$queIns, 
       type = "p", pch = 16, cex = .7,ylim = c(0,100000))
  rm <- rollapply(data= c(rep(NA, 15), df[df$seqnames == chrChoice,]$queIns,rep(NA, 15)), 
                  width = 31, na.rm = TRUE, FUN = mean)
  lines(df[df$seqnames == chrChoice,]$start, rm, col = 2)
  mtext(paste(genomes[2],"gain (bp)"),side = 2,line = 3, cex = .8)
  
  plot(df[df$seqnames == chrChoice,]$start,
       df[df$seqnames == chrChoice,]$queDel, 
       type = "p", pch = 16, cex = .7,ylim = c(0,100000))
  rm <- rollapply(data= c(rep(NA, 15), df[df$seqnames == chrChoice,]$queDel,rep(NA, 15)), 
                  width = 31, na.rm = TRUE, FUN = mean)
  lines(df[df$seqnames == chrChoice,]$start, rm, col = 2)
  mtext(paste(genomes[2],"loss (bp)"),side = 2,line = 3, cex = .8)
  
  
  title(main = paste("mm10", chrChoice), outer = TRUE)
  
}
dev.off()



# bootstrap mean 

# bootstrap sd?


#important stats

resGenomeStat <- matrix(data = NA, ncol = 9, nrow = 2)
rownames(resGenomeStat) <- genomes
colnames(resGenomeStat) <- c("genome_size", "fills","ancestral" ,
                             "fills_percent_ancestral", "all_gaps", "placed_gaps", 
                             "placed_gaps_percentage","gain", "loss")
resGenomeStat <- as.data.frame(resGenomeStat)

# genome size
## mm10
resGenomeStat$genome_size[1] <- round( (sum(as.numeric(seqlengths(refSeqGaps.gr))) - sum(width(refSeqGaps.gr)))/1e6 )
## hg19
resGenomeStat$genome_size[2] <-round((sum(as.numeric(seqlengths(queSeqGaps.gr))) - sum(width(queSeqGaps.gr)))/1e6)

# all non gaps
## mm10
resGenomeStat$fills[1] <- round(sum(width(refFill.gr))/1e6)
## hg19
resGenomeStat$fills[2] <- round(sum(width(queFill.gr))/1e6)

#ancestral
resGenomeStat$ancestral[1] <- round(sum(width(refAncDna.gr))/1e6)
resGenomeStat$ancestral[2] <- round(sum(width(queAncDna.gr))/1e6)


# ancestral %
# mm10
resGenomeStat$fills_percent_ancestral[1] <- round(sum(width(GenomicRanges::intersect(refFill.gr, refAncDna.gr))) / sum(width(refFill.gr)) * 100, digits = 1)
# hg19
resGenomeStat$fills_percent_ancestral[2] <- round(sum(width(GenomicRanges::intersect(queFill.gr, queAncDna.gr))) / sum(width(queFill.gr)) * 100, digits = 1)


# all gaps
#mm10
resGenomeStat$all_gaps[1] <- round(sum(width(refFillGaps.gr))/1e6)
# hg19
resGenomeStat$all_gaps[2] <- round(sum(width(queFillGaps.gr))/1e6)


# placeable gaps 
#mm10
resGenomeStat$placed_gaps[1] <- round(sum(width(refGap.gr))/1e6)
#hg19
resGenomeStat$placed_gaps[2] <- round(sum(width(queGap.gr))/1e6)



# placeable gaps percentage
#mm10
resGenomeStat$placed_gaps_percentage[1] <- round(sum(width(refGap.gr))/sum(width(refFillGaps.gr))*100,digits = 1)
#hg19
resGenomeStat$placed_gaps_percentage[2] <- round(sum(width(queGap.gr))/sum(width(queFillGaps.gr))*100, digits = 1)


# gain
## mm10
resGenomeStat$gain[1] <- round(sum(width(refGapNonAnc.gr))/1e6)
# hg19
resGenomeStat$gain[2] <- round(sum(width(queGapNonAnc.gr))/1e6)

# loss
# mm10 
resGenomeStat$loss[1] <- round(sum(width(queGapAnc.gr))/1e6)
# hg19
resGenomeStat$loss[2] <- round(sum(width(refGapAnc.gr))/1e6)

t(resGenomeStat)


