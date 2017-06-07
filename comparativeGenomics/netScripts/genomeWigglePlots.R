




library(spdep)
library(igraph)
library(zoo)
library(dplyr)
library(GenomicRanges)

rm(list = ls())


devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")




specRef = "mm10"
specQue = "hg19"

load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/",specRef,".",specQue,".netData.RData",sep = ""))
#load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/",specRef,".stretch.RData", sep = ""))
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/shiftData/",specRef,".expand.breaks.RData", sep = ""))
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/",specRef,".synthBin.RData", sep = ""))


# use our stretch function to breate the hnew marks
# could maybe use synthBin as expanded seqlengths

refGenome.gr <- GRanges(seqnames = Rle(refChrInfo$chrom),
                        ranges = IRanges(start = 1,end= refChrInfo$size))
seqlevels(refGenome.gr) <- refChrInfo$chrom
seqlengths(refGenome.gr) <- refChrInfo$size
genome(refGenome.gr) <- specRef

refGenome.gr <- sort(sortSeqlevels(refGenome.gr))
seqlengths(refGenome.gr) <- seqlengths(refFill.gr)

refGenomeBin.gr <- unlist(slidingWindows(refGenome.gr, width = 10e6, 10e6))
refGenomeBin.gr <- resize(refGenomeBin.gr, width = 1, fix = "start")

refGenomeMarks.df <- data.frame(genoExpandStretch(x.gr = refGenomeBin.gr,synthGenome = newSynthRefShift, 
                                                  expandedSeqlengths = seqlengths(synthBin.gr)))







binSize <- width(synthBin.gr)[1]
# create df for plotting
df <- data.frame(synthBin.gr)
# remove bins with missing data
df[rowSums(df[,c("missingGap", "seqGap")]) > (.5 * binSize),c("refIns", "refDel", "queIns", "queDel","gcContent","fill")] <- NA

df[,c("refIns", "refDel", "queIns", "queDel")] <- df[,c("refIns", "refDel", "queIns", "queDel")]/(binSize - rowSums(df[,c("missingGap", "seqGap")]))
df[,c("refIns", "refDel", "queIns", "queDel")] <- df[,c("refIns", "refDel", "queIns", "queDel")]*binSize

# if the fill content is too low we won't use gc content
df[df$fill < 5e3 & !is.na(df$gcContent),"gcContent"] <- NA


df$refTurnover <- df$refIns+df$refDel
df$refNetGrowth <- df$refIns-df$refDel

df$queTurnover <- df$queIns+df$queDel
df$queNetGrowth <- df$queIns-df$queDel


colChoice <- c("refIns", "refDel", "queIns", "queDel", 
               "gcContent", "fill", 
               "refTurnover", "refNetGrowth", "queTurnover", "queNetGrowth")

# next we get the df rolling means
dfRollMean <- NULL
dfSplit <- split(df,f=df$seqnames)
for(i in 1:length(dfSplit)){
  for(j in 1:length(colChoice)){
    rMean <- rollapply(data= c(rep(NA, 7), dfSplit[[i]][,colChoice[j]],rep(NA, 7)), 
                  width = 15, na.rm = TRUE, FUN = mean)
    rMean[is.na(dfSplit[[i]][,colChoice[j]])] <- NA
    dfSplit[[i]][,colChoice[j]] <- rMean
  }
  dfRollMean <- rbind(dfRollMean, dfSplit[[i]])
}


dfRollSd <- NULL
dfSplit <- split(df,f=df$seqnames)
for(i in 1:length(dfSplit)){
  for(j in 1:length(colChoice)){
    rSd <- rollapply(data= c(rep(NA, 7), dfSplit[[i]][,colChoice[j]],rep(NA, 7)), 
                       width = 15, na.rm = TRUE, FUN = sd)
    rSd[is.na(dfSplit[[i]][,colChoice[j]])] <- NA
    dfSplit[[i]][,colChoice[j]] <- rSd 
  }
  dfRollSd <- rbind(dfRollSd, dfSplit[[i]])
}




blueRed <- colorRampPalette(c("blue", "white", "red"))
heatmap(cor(as.matrix(dfRollMean[,colChoice]), use = "complete.obs"), col= blueRed(400), zlim=c(-1,1), scale="none")

newDF <- df[complete.cases(df),]
newDF$chrType = NA
newDF$chrType[newDF$seqnames == "chrX"] <- "sexChr"
newDF$chrType[newDF$seqnames != "chrX"] <- "autoChr"
lMod <- lm(data = newDF, refNetGrowth ~ gcContent + chrType + refTurnover)
summary(lMod)
layout(1)
plot(lMod)
plot(lMod$residuals)


plot(lMod)
summary(lMod)
plot(newDF$gcContent,lMod$residuals)

# set ylims
ylim = c(0, .65*binSize)



# plot
chrPlot <- seqlevels(synthBin.gr)[-grep("_", seqlevels(synthBin.gr))]
chrPlot <- chrPlot[-grep("Y", chrPlot)]
chrPlot <- chrPlot[-grep("M", chrPlot)]



# maybe remove this plot
if(specRef == "hg19"){
  laminB1 <- read.table(paste("~/Desktop/RTN_domains/data/UCSCtracks/", specRef,"/laminB1Lads", sep = ""), 
                        col.names = c("bin", "seqnames", "start", "end"))  
}

if(specRef == "mm10"){
  laminB1 <- read.table(paste("~/Desktop/RTN_domains/data/UCSCtracks/", specRef,"/laminB1Lads", sep = ""),
                        col.names = c("seqnames", "start", "end", "name", "score", "name2"))
}

laminB1.gr <- GRanges(laminB1)
start(laminB1.gr) <- start(laminB1.gr + 1)
seqlevels(laminB1.gr) <- refChrInfo$chrom
seqlengths(laminB1.gr) <- refChrInfo$size
laminB1.gr <- sort(sortSeqlevels(laminB1.gr))
laminB1.gr <- genoExpandStretch(x.gr = laminB1.gr,synthGenome = newSynthRefShift,expandedSeqlengths = seqlengths(synthBin.gr))


pdf(file = paste("~/Desktop/", genomes["ref"] ,"stats.pdf", sep = ""), width = 15, height = 6, onefile = TRUE)

for(chrChoice in chrPlot){
  
  dfChrMean <- dfRollMean[dfRollMean$seqnames == chrChoice,]
  dfChrSd <- dfRollSd[dfRollSd$seqnames == chrChoice,]
  
  pX <- dfChrMean[complete.cases(dfChrMean), "start"]
  pX <- c(pX, rev(pX))
  
  plotChoices <- c("gcContent", "refNetGrowth", "queNetGrowth", "refTurnover", "queTurnover")
  pY <- matrix(NA, nrow = 2*nrow(dfChrMean[complete.cases(dfChrMean),]), ncol = length(plotChoices))
  colnames(pY) <- plotChoices
  for(plotChoice in plotChoices){
    pYtop <- dfChrMean[,plotChoice] + (qnorm(.975) * dfChrSd[,plotChoice] / sqrt(15))
    pYbottom <- dfChrMean[,plotChoice] - (qnorm(.975) * dfChrSd[,plotChoice] / sqrt(15))
    pYtop <- pYtop[ complete.cases(dfChrMean)]
    pYbottom <- pYbottom[ complete.cases(dfChrMean)]
    pY[,plotChoice] <- c(pYtop, rev(pYbottom))
  }
  pY <- data.frame(pY)
  
  
  layout(matrix(1:6, byrow = T, ncol = 2), widths =c(2,10) )
  par(mar = c(3,4,1,2), oma = c(5,5,5,5))
  
  
  ## GC content
  
  d.plot <- density(df$gcContent ,na.rm = T)
  plot(-d.plot$y, d.plot$x, type = "l", axes = FALSE, ylim = c(.25,.7) , 
        yaxs="i",main = "", ylab = "")
  box()
  
  
  plot(dfChrMean[,"start"],dfChrMean[, "gcContent"] * 100 , type = "n", ylim = c(25,70), 
       xlim = c(0, max(seqlengths(synthBin.gr))),
       frame.plot = FALSE, xaxs = "i", yaxs = "i",
       xaxt = "n", ylab = "GC content (%)")
  
  polygon(pX, pY$gcContent*100, density = -1 , col = scales::alpha("black",.2), border = NA)
  lines(dfChrMean[,"start"],dfChrMean[, "gcContent"] * 100)
  
  
  
  rect(xleft = dfChrMean$start[!complete.cases(dfChrMean)], 
       xright = dfChrMean$end[!complete.cases(dfChrMean)],
       ybottom = 25, ytop = 70, col = "white", border = "white")
  
  ### laminB1
  rect(xleft = start(laminB1.gr[seqnames(laminB1.gr) == chrChoice]),
       xright = end(laminB1.gr[seqnames(laminB1.gr) == chrChoice]),
       ybottom =  25,ytop = 70, col = scales::alpha(colour = "darkblue", .3), border = FALSE)
  
  rect(xleft = 0, xright = max(seqlengths(synthBin.gr)[chrChoice]), ybottom = 25, ytop = 70)
  
  
  
  axis(side = 1, at = refGenomeMarks.df$start[refGenomeMarks.df$seqnames == chrChoice], 
       seq(from = 0, by = 10, length.out = sum(refGenomeMarks.df$seqnames == chrChoice)))
  
  
  ### growth
  
  d.plot <- density(df$refNetGrowth ,na.rm = T)
  plot(-d.plot$y, d.plot$x, type = "l", axes = FALSE, ylim = c(-7e4,7e4) , 
       yaxs="i",main = "", ylab = "")

  d.plot <- density(df$queNetGrowth ,na.rm = T)
  lines(-d.plot$y, d.plot$x, col = 2)
  box()

  
  plot(dfChrMean[,c("start", "refNetGrowth")], type = "n", ylim = c(-7e4,7e4),
       xlim = c(0, max(seqlengths(synthBin.gr))),
       frame.plot = FALSE, xaxs = "i", yaxs = "i", ylab = "net genome growth (bp)", xaxt = "n")
  
  polygon(pX, pY$refNetGrowth, density = -1 , col = scales::alpha("black",.2), border = NA)
  lines(dfChrMean[,c("start", "refNetGrowth")])
  polygon(pX, pY$queNetGrowth, density = -1 , col = scales::alpha("black",.2), border = NA)
  lines(dfChrMean[,c("start", "queNetGrowth")],  col = 2)
  
  rect(xleft = dfChrMean$start[!complete.cases(dfChrMean)], 
       xright = dfChrMean$end[!complete.cases(dfChrMean)],
       ybottom = -7e4, ytop = 7e4, col = "white", border = "white")
  rect(xleft = 0, xright = max(seqlengths(synthBin.gr)[chrChoice]), ybottom = -7e4, ytop = 7e4)
  
  lines(c(0,max(seqlengths(synthBin.gr)[chrChoice])), c(0,0), lty = 2)
  
  axis(side = 1, at = refGenomeMarks.df$start[refGenomeMarks.df$seqnames == chrChoice], 
       seq(from = 0, by = 10, length.out = sum(refGenomeMarks.df$seqnames == chrChoice)))
  
  
  ### turnover
  
  d.plot <- density(df$refTurnover ,na.rm = T)
  plot(-d.plot$y, d.plot$x, type = "l", axes = FALSE, ylim = c(0, 1.5e5) , 
       yaxs="i",main = "", ylab = "")
  
  d.plot <- density(df$queTurnover ,na.rm = T)
  lines(-d.plot$y, d.plot$x, col = 2)
  box()
  
  
  plot(dfChrMean[,c("start", "refTurnover")], type = "n", ylim = c(0, 1.5e5),
       xlim = c(0, max(seqlengths(synthBin.gr))),
       frame.plot = FALSE, xaxs = "i", yaxs = "i",
       xaxt = "n", ylab = "genome turnover (bp)")
  
  polygon(pX, pY$refTurnover, density = -1 , col = scales::alpha("black",.2), border = NA)
  lines(dfChrMean[,c("start", "refTurnover")])
  polygon(pX, pY$queTurnover, density = -1 , col = scales::alpha("black",.2), border = NA)
  lines(dfChrMean[,c("start", "queTurnover")],  col = 2)
  
  rect(xleft = dfChrMean$start[!complete.cases(dfChrMean)], 
       xright = dfChrMean$end[!complete.cases(dfChrMean)],
       ybottom = 0, ytop = 15e4, col = "white", border = "white")
  rect(xleft = 0, xright = max(seqlengths(synthBin.gr)[chrChoice]), ybottom = 0, ytop = 1.5e5)
  
  axis(side = 1, at = refGenomeMarks.df$start[refGenomeMarks.df$seqnames == chrChoice], 
       seq(from = 0, by = 10, length.out = sum(refGenomeMarks.df$seqnames == chrChoice)))
  
  lines(c(0,max(seqlengths(synthBin.gr)[chrChoice])), 
        c(mean(df$refTurnover, na.rm = TRUE),
          mean(df$refTurnover, na.rm = TRUE)), 
        lty = 2)
  lines(c(0,max(seqlengths(synthBin.gr)[chrChoice])), 
        c(mean(df$queTurnover, na.rm = TRUE),
          mean(df$queTurnover, na.rm = TRUE)), 
        lty = 2, col = 2)
  
  title(paste(genomes["ref"], chrChoice), outer = TRUE, cex.main = 2)
  mtext(text = "position (Mb)", side = 1, 
        at = max(seqlengths(synthBin.gr))/2, line = 3 , cex = 1)
  
}

dev.off()

















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





# plot turnOver, gain and GC


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



pdf(file = "Desktop/turnoverHuman.pdf", width = 10, height  = 3)
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
          density = -1 , col = scales::alpha("grey",.5), border = NA)
  polygon(pX, pYque, 
          density = -1 , col = scales::alpha("grey",.5), border = NA)
  pos <- df[df$seqnames == chrChoice, c("start","end")]
  pos <- pos[!complete.cases(df[df$seqnames == chrChoice,]),]
  
  lines(df[df$seqnames == chrChoice, "start"],rmQue, col = 2)
  lines(df[df$seqnames == chrChoice, "start"],rmRef, col = 1)
  rect(xleft = pos$start -1,xright = pos$end + 1, ytop = 300000, ybottom = -100000, col = "white", border = NA)
  
}

dev.off()


pdf(file = "Desktop/gcHuman.pdf", width = 10, height  = 3)
for(chrChoice in chrPlot){
  
  rmRef <- rollapply(data= c(rep(NA, 7), (df$gcContent)[df$seqnames == chrChoice],rep(NA, 7)), 
                     width = 15, na.rm = TRUE, FUN = mean)
  rsdRef <- rollapply(data= c(rep(NA, 7), (df$gcContent)[df$seqnames == chrChoice],rep(NA, 7)), 
                      width = 15, na.rm = TRUE, FUN = sd)
  
  rmQue <- rollapply(data= c(rep(NA, 7), (df$gcContent)[df$seqnames == chrChoice],rep(NA, 7)), 
                     width = 15, na.rm = TRUE, FUN = mean)
  rsdQue <- rollapply(data= c(rep(NA, 7), (df$gcContent)[df$seqnames == chrChoice],rep(NA, 7)), 
                      width = 15, na.rm = TRUE, FUN = sd)
  
  
  plot(df[df$seqnames == chrChoice, "start"],rmRef,
       ylim = c(.3,.6), type = "n", xaxt = "n",
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
          density = -1 , col = scales::alpha("grey",.5), border = NA)
  polygon(pX, pYque, 
          density = -1 , col = scales::alpha("grey",.5), border = NA)
  pos <- df[df$seqnames == chrChoice, c("start","end")]
  pos <- pos[!complete.cases(df[df$seqnames == chrChoice,]),]
  
  lines(df[df$seqnames == chrChoice, "start"],rmRef, col = 1)
  rect(xleft = pos$start -1,xright = pos$end + 1, ytop = 300000, ybottom = -100000, col = "white", border = NA)
  
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




