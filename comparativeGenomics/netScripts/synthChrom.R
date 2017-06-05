pkgs = names(sessionInfo()$otherPkgs)
pkgs = paste('package:', pkgs, sep = "")
lapply(pkg, detach, character.only = TRUE, unload = TRUE, force = TRUE)

library(spdep)
library(igraph)
library(dplyr)
library(zoo)
library(GenomicRanges)

rm(list = ls())

genoExpandBreak <- function(x.gr, synthGenome, expandedSeqlengths){
  seqlengths(x.gr) <- seqlengths(synthGenome)
  ol <- findOverlaps(x.gr, synthGenome)
  pInt <- pintersect(x.gr[queryHits(ol)], synthGenome[subjectHits(ol)], drop.nohit.ranges=TRUE)
  seqlengths(pInt) <- expandedSeqlengths
  expanded.gr <- shift(pInt, shift = synthGenome$shift[subjectHits(ol)])
  return(expanded.gr)
}

genoExpandStretch <- function(x.gr, synthGenome, expandedSeqlengths){
  seqlengths(x.gr) <- seqlengths(synthGenome)
  olStart <- findOverlaps(x.gr, synthGenome, select = "first")
  olEnd <- findOverlaps(x.gr, synthGenome, select = "last")
  expanded.gr <- GRanges(seqnames = seqnames(x.gr), 
                         ranges = IRanges(start = start(x.gr) + synthGenome[olStart]$shift,
                                          end = end(x.gr) + synthGenome[olEnd]$shift))
  mcols(expanded.gr) <- mcols(x.gr)
  seqlengths(expanded.gr) <- expandedSeqlengths
  genome(expanded.gr) <- genome(x.gr)
  return(expanded.gr)
}





specRef = "hg19"
specQue = "mm10"

load(paste("Desktop/RTN_domains/R_objects/mappedGaps/",specRef,".",specQue,".netData.RData",sep = ""))
load(paste("Desktop/RTN_domains/R_objects/mappedGaps/",specRef,".stretch.RData", sep = ""))

load(paste("Desktop/RTN_domains/R_objects/mappedGaps/",specRef,".expand.breaks.RData", sep = ""))


# annotate extra files

refMissingGaps.gr <- GenomicRanges::setdiff(refFillGaps.gr, refGap.gr, ignore.strand = TRUE)
all(!overlapsAny(refMissingGaps.gr, refGap.gr))

mcols(refMissingGaps.gr)$queRanges <- GRanges(seqnames = "empty", ranges = IRanges(start = 1, end = 1))
mcols(refMissingGaps.gr)$sData = "*"
mcols(refMissingGaps.gr)$chainID = NA
mcols(refMissingGaps.gr)$type = "missingGap"

refSeqGaps.gr <- reduce(refSeqGaps.gr)
mcols(refSeqGaps.gr)$queRanges <- GRanges(seqnames = "empty", ranges = IRanges(start = 1, end = 1))
mcols(refSeqGaps.gr)$sData = "*"
mcols(refSeqGaps.gr)$chainID = NA
mcols(refSeqGaps.gr)$type = "seqGap"



missingGenome.gr <- sort(c(refMissingGaps.gr, refSeqGaps.gr))

missingGenome.gr <- genoExpandBreak(missingGenome.gr, newSynthRefShift, seqlengths(stretchedRef.gr))



refSynth <- sort(c(stretchedRef.gr,missingGenome.gr))


# check for overlappign ranges
cov <- coverage(refSynth)
sl <- IRanges::slice(cov, lower = 2)
all(unlist(lapply(sl, length)) == 0)


# create synthetic genome
synthGenome <- GRanges(seqnames = seqlevels(stretchedRef.gr), 
                       ranges = IRanges(width = seqlengths(stretchedRef.gr), end = seqlengths(stretchedRef.gr)))

# select bin size
binSize <- 2e5

# bin synthetic genome
synthBin.gr <- unlist(slidingWindows(synthGenome, width = binSize, step = binSize))

# get overlapping ranges
ol <- findOverlaps(refSynth, synthBin.gr)
pInt <- pintersect(refSynth[queryHits(ol)], synthBin.gr[subjectHits(ol)])


# sort into data frame
df <- data.frame(width = width(pInt), type = pInt$type, binNo = subjectHits(ol))
sumDF <- summarise(group_by(df, type, binNo), total = sum(width))

# setup bin columns
mcols(synthBin.gr)$refIns <- 0
mcols(synthBin.gr)$refDel <- 0
mcols(synthBin.gr)$queIns <- 0
mcols(synthBin.gr)$queDel <- 0
mcols(synthBin.gr)$missingGap <- 0
mcols(synthBin.gr)$seqGap <- 0

# sort into bin columns
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "refIns"]])$refIns = sumDF$total[sumDF$type == "refIns"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "refDel"]])$refDel = sumDF$total[sumDF$type == "refDel"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "queIns"]])$queIns = sumDF$total[sumDF$type == "queIns"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "queDel"]])$queDel = sumDF$total[sumDF$type == "queDel"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "missingGap"]])$missingGap = sumDF$total[sumDF$type == "missingGap"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "seqGap"]])$seqGap = sumDF$total[sumDF$type == "seqGap"]


# create df for plotting
df <- data.frame(synthBin.gr)
# remove bins with missing data
df[rowSums(df[,c("missingGap", "seqGap")]) > (.5 * binSize),c("refIns", "refDel", "queIns", "queDel")] <- NA

df[,c("refIns", "refDel", "queIns", "queDel")] <- df[,c("refIns", "refDel", "queIns", "queDel")]/(binSize - rowSums(df[,c("missingGap", "seqGap")]))
df[,c("refIns", "refDel", "queIns", "queDel")] <- df[,c("refIns", "refDel", "queIns", "queDel")]*binSize

# set ylims
ylim = c(0, .65*binSize)



# plot
chrPlot <- seqlevels(refSynth)[-grep("_", seqlevels(refSynth))]
chrPlot <- chrPlot[-grep("Y", chrPlot)]
chrPlot <- chrPlot[-grep("M", chrPlot)]




gapChoice = c("refIns", "refDel", "queIns", "queDel")
direction = c("gain", "loss", "gain", "loss")
ylimMax = c(120000, 40000,120000,120000)

pdf(file = paste("~/Desktop/RTN_domains/plots/netGainLoss/", genomes["ref"],"SlidingWindow.pdf", sep = ""),height = 6,width = 12 ,onefile = TRUE)
for(chrChoice in chrPlot){
  layout(1:4)
  par(mar = c(2,3,0,1), oma = c(5,5,5,5))
  
  for(i in 1:4){
  plot(df[df$seqnames == chrChoice,]$start,
       df[df$seqnames == chrChoice,gapChoice[i]], 
       type = "p", pch = 16, cex = .7,ylim = c(0,ylimMax[i]), col = "grey70")
  rm <- rollapply(data= c(rep(NA, 7), df[df$seqnames == chrChoice,gapChoice[i]],rep(NA, 7)), 
                  width = 15, na.rm = TRUE, FUN = mean)
  rsd <- rollapply(data= c(rep(NA, 7), df[df$seqnames == chrChoice,gapChoice[i]],rep(NA, 7)), 
                   width = 15, na.rm = TRUE, FUN = sd)
  lines(df[df$seqnames == chrChoice,]$start, rm, col = 2, lwd = 2)
  lines(df[df$seqnames == chrChoice,]$start, rm + (2*rsd), col = 4, lwd = 1)
  lines(df[df$seqnames == chrChoice,]$start,  rm - (2*rsd), col = 4, lwd = 1)
  
  mtext(paste(c(genomes["ref"], genomes["ref"], genomes["que"], genomes["que"])[i],direction[i] ,"(bp)"),side = 2,line = 3, cex = .8)
  abline(h = mean(df[,gapChoice[i]], na.rm = TRUE), lty = 2, col= 1, lwd = 1)
  
  }
  
  title(main = paste(genomes["ref"], chrChoice), outer = TRUE)
  
}
dev.off()


## GstarI approach


ol<-findOverlaps(synthBin.gr, maxgap = 3*binSize)
ol <- ol[!(isRedundantHit(ol))]
# remove NA hits
ol <- ol[!is.na(df$refIns[queryHits(ol)])]
ol <- ol[!is.na(df$refIns[subjectHits(ol)])]

G <- graph.data.frame(as.matrix(ol),directed=FALSE)
A <- as_adjacency_matrix(G,type="both",names=FALSE,sparse=TRUE)
wMat <- mat2listw(A)



df3 <- df

for(i in 1:4){
  score <- df[!is.na(df$refDel),gapChoice[i]]
  G <- localG(x = score,wMat)
  df3[!is.na(df3$refIns),gapChoice[i]] <- G
}



refGenome.gr <- GRanges(seqnames = Rle(refChrInfo$chrom),
                        ranges = IRanges(start = 1,end= refChrInfo$size+1))
refGenome.gr <- sort(sortSeqlevels(refGenome.gr))
seqlengths(refGenome.gr) <- seqlengths(refFill.gr)

refGenomeBin.gr <- unlist(slidingWindows(refGenome.gr, width = 10e6, 10e6))
refGenomeBin.gr <- resize(refGenomeBin.gr, width = 1, fix = "start")
#refGenomeBin.gr <- refGenomeBin.gr[start(refGenomeBin.gr) != 1]

refGenomeMarks.df <- data.frame(genoExpandStretch(x.gr = refGenomeBin.gr,synthGenome = newSynthRefShift, expandedSeqlengths = seqlengths(stretchedRef.gr)))




pdf(file = paste("~/Desktop/RTN_domains/plots/netGainLoss/",genomes["ref"] ,"Gstar.pdf", sep = ""),height = 6,width = 12 ,onefile = TRUE)
for(chrChoice in chrPlot){
  layout(1:4)
  par(mar = c(2,3,0,1), oma = c(5,5,5,5))
  
  for(i in 1:4){
    
    plot(df[df$seqnames == chrChoice,]$start,
         scale(df[,gapChoice[i]])[df$seqnames == chrChoice], 
         type = "p", pch = 16, cex = .7,ylim = c(-6,+6), col = "grey70", xaxt = "n")
    
    axis(side = 1, at = refGenomeMarks.df$start[refGenomeMarks.df$seqnames == chrChoice],
        labels = seq(from = 0, by = 10e6, length.out = sum(refGenomeMarks.df$seqnames == chrChoice)))
    
    lines(df3[df3$seqnames == chrChoice,]$start,
          df3[df3$seqnames == chrChoice,gapChoice[i]])
    
    mtext(paste(c(genomes["ref"], genomes["ref"], genomes["que"], genomes["que"])[i],direction[i] ,"(bp)"),side = 2,line = 3, cex = .8)
    abline(h = 0, lty = 2, col= 4, lwd = 1)
    abline(h = 2, lty = 2, col= 2, lwd = 1)
    abline(h = -2, lty = 2, col= 2, lwd = 1)
  }
  
  title(main = paste(genomes["ref"], chrChoice), outer = TRUE)
  
}
dev.off()


pdf(file = "Desktop/diffsMM10.pdf", width = 10, height  = 3)
for(chrChoice in chrPlot){

rmRef <- rollapply(data= c(rep(NA, 7), (df$refIns - df$refDel)[df$seqnames == chrChoice],rep(NA, 7)), 
                   width = 15, na.rm = TRUE, FUN = mean)
rsdRef <- rollapply(data= c(rep(NA, 7), (df$refIns - df$refDel)[df$seqnames == chrChoice],rep(NA, 7)), 
                    width = 15, na.rm = TRUE, FUN = sd)

rmQue <- rollapply(data= c(rep(NA, 7), (df$queIns - df$queDel)[df$seqnames == chrChoice],rep(NA, 7)), 
                   width = 15, na.rm = TRUE, FUN = mean)
rsdQue <- rollapply(data= c(rep(NA, 7), (df$queIns - df$queDel)[df$seqnames == chrChoice],rep(NA, 7)), 
                    width = 15, na.rm = TRUE, FUN = sd)


plot(df[df$seqnames == chrChoice, "start"],rmRef,
     ylim = c(-100e3,100e3), type = "n", xaxt = "n",
     main = chrChoice)
#lines(df[df$seqnames == chrChoice, "start"],rmQue, col = 2)
abline(h=0)
axis(side = 1, at = refGenomeMarks.df$start[refGenomeMarks.df$seqnames == chrChoice],
     labels = seq(from = 0, by = 10e6, length.out = sum(refGenomeMarks.df$seqnames == chrChoice)))


pX <- c(df[df$seqnames == chrChoice, "start"][complete.cases(df[df$seqnames == chrChoice,])],
        rev(df[df$seqnames == chrChoice, "start"][complete.cases(df[df$seqnames == chrChoice,])]))

pYref <- c( (rmRef + rsdRef)[complete.cases(df[df$seqnames == chrChoice,])],
          rev((rmRef - rsdRef)[complete.cases(df[df$seqnames == chrChoice,])]))
pYque <- c( (rmQue + rsdQue)[complete.cases(df[df$seqnames == chrChoice,])],
            rev((rmQue - rsdQue)[complete.cases(df[df$seqnames == chrChoice,])]))
polygon(pX, pYref, 
        density = -1 , col = alpha("grey",.5), border = NA)
polygon(pX, pYque, 
        density = -1 , col = alpha("grey",.5), border = NA)
pos <- df[df$seqnames == chrChoice, c("start","end")]
pos <- pos[!complete.cases(df[df$seqnames == chrChoice,]),]

lines(df[df$seqnames == chrChoice, "start"],rmQue, col = 2)
lines(df[df$seqnames == chrChoice, "start"],rmRef, col = 1)
rect(xleft = pos$start -1,xright = pos$end + 1, ytop = 100000, ybottom = -100000, col = "white", border = NA)

}

dev.off()



## how to integrate this data?
## What questions are worth asking?

## what is hapening on average?



pdf(file = "Desktop/turnOverMM10.pdf", width = 10, height  = 3)
for(chrChoice in chrPlot){
  
  rmRef <- rollapply(data= c(rep(NA, 7), (df$refIns + df$refDel)[df$seqnames == chrChoice],rep(NA, 7)), 
                     width = 15, na.rm = TRUE, FUN = mean)
  rsdRef <- rollapply(data= c(rep(NA, 7), (df$refIns + df$refDel)[df$seqnames == chrChoice],rep(NA, 7)), 
                      width = 15, na.rm = TRUE, FUN = sd)
  
  rmQue <- rollapply(data= c(rep(NA, 7), (df$queIns + df$queDel)[df$seqnames == chrChoice],rep(NA, 7)), 
                     width = 15, na.rm = TRUE, FUN = mean)
  rsdQue <- rollapply(data= c(rep(NA, 7), (df$queIns + df$queDel)[df$seqnames == chrChoice],rep(NA, 7)), 
                      width = 15, na.rm = TRUE, FUN = sd)
  
  
  plot(df[df$seqnames == chrChoice, "start"],rmRef,
       ylim = c(0,200e3), type = "n", xaxt = "n",
       main = chrChoice)
  #lines(df[df$seqnames == chrChoice, "start"],rmQue, col = 2)
  #abline(h=0)
  axis(side = 1, at = refGenomeMarks.df$start[refGenomeMarks.df$seqnames == chrChoice],
       labels = seq(from = 0, by = 10e6, length.out = sum(refGenomeMarks.df$seqnames == chrChoice)))
  
  
  pX <- c(df[df$seqnames == chrChoice, "start"][complete.cases(df[df$seqnames == chrChoice,])],
          rev(df[df$seqnames == chrChoice, "start"][complete.cases(df[df$seqnames == chrChoice,])]))
  
  pYref <- c( (rmRef + rsdRef)[complete.cases(df[df$seqnames == chrChoice,])],
              rev((rmRef - rsdRef)[complete.cases(df[df$seqnames == chrChoice,])]))
  pYque <- c( (rmQue + rsdQue)[complete.cases(df[df$seqnames == chrChoice,])],
              rev((rmQue - rsdQue)[complete.cases(df[df$seqnames == chrChoice,])]))
  polygon(pX, pYref, 
          density = -1 , col = alpha("grey",.5), border = NA)
  polygon(pX, pYque, 
          density = -1 , col = alpha("grey",.5), border = NA)
  pos <- df[df$seqnames == chrChoice, c("start","end")]
  pos <- pos[!complete.cases(df[df$seqnames == chrChoice,]),]
  
  lines(df[df$seqnames == chrChoice, "start"],rmQue, col = 2)
  lines(df[df$seqnames == chrChoice, "start"],rmRef, col = 1)
  rect(xleft = pos$start -1,xright = pos$end + 1, ytop = 300000, ybottom = -100000, col = "white", border = NA)
  
}

dev.off()















lines(df[df$seqnames == chrChoice, "start"],rmQue, col = 2)
lines(df[df$seqnames == chrChoice, "start"],rmQue + rsdQue, lty = 2, col=2)
lines(df[df$seqnames == chrChoice, "start"],rmQue - rsdQue, lty = 2, col = 2)

# it is almost like the exchange level is similar to human when mouse insertion goes up
# rather than this overall pattern of deletion

# regions of conserved size


plot(rep(1,10))
polygon(c(1:10,10:1), 
        c( 1 + c(rep(.2,5),0, rep(.2,4)), rev(1-c(rep(.2,5),0, rep(.2,4)) )))






library(ggplot2)

chrChoice = "chr4"
dfG <- filter(df, seqnames == chrChoice)
dfG <- data.frame(start = c(dfG$start,dfG$start), 
                  netGain = c(dfG$refIns - dfG$refDel, dfG$queIns - dfG$queDel),
                  genome = c( rep("ref", nrow(dfG)) , rep("que", nrow(dfG)) ))

ggplot(data = dfG, aes(x = start, y = netGain, col = genome))  + 
  geom_point(cex = .5) + 
  geom_smooth(span = .1, method = "loess", fullrange = TRUE)+
  scale_x_continuous(breaks = refGenomeMarks.df$start[refGenomeMarks.df$seqnames == chrChoice], 
                   labels=seq(from = 0, by = 10e6, length.out = sum(refGenomeMarks.df$seqnames == chrChoice)))






plot(df[df$seqnames == "chr9", "start"],
     (df$queIns)[df$seqnames == "chr9"], 
     ylim = c(0,100e3))
abline(h = 0)


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




sum(width(queFill.gr))


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
