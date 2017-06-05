

library(spdep)
library(igraph)
library(zoo)
library(dplyr)
library(GenomicRanges)

rm(list = ls())


devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")




specRef = "hg19"
specQue = "mm10"

load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/",specRef,".",specQue,".netData.RData",sep = ""))
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/",specRef,".stretch.RData", sep = ""))
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/shiftData/",specRef,".expand.breaks.RData", sep = ""))




# use our stretch function to breate the hnew marks
# could maybe use synthBin as expanded seqlengths
refGenomeMarks.df <- data.frame(genoExpandStretch(x.gr = refGenomeBin.gr,synthGenome = newSynthRefShift, 
                                                  expandedSeqlengths = seqlengths(stretchedRef.gr)))








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



