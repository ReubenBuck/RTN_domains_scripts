library(grid)
library(gridExtra)
library(gtable)
library(GenomicRanges)

rm(list = ls())

specRef = "hg19"
specQue = "mm10"

# get dta
load(paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/", specRef, "repNoRep.RData", sep = ""))
load(paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/", specRef, "synthBinNorm.RData", sep = ""))
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/",specRef,".",specQue,".netData.RData",sep = ""))
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/shiftData/",specRef,".expand.breaks.RData", sep = ""))

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")

# 50Mb marks
refGenome.gr <- GRanges(seqnames = Rle(refChrInfo$chrom),
                        ranges = IRanges(start = 1,end= refChrInfo$size))
seqlevels(refGenome.gr) <- refChrInfo$chrom
seqlengths(refGenome.gr) <- refChrInfo$size
genome(refGenome.gr) <- specRef

refGenome.gr <- sort(sortSeqlevels(refGenome.gr))
seqlengths(refGenome.gr) <- seqlengths(refFill.gr)

refGenomeBin.gr <- unlist(slidingWindows(refGenome.gr, width = 50e6, 50e6))
refGenomeBin.gr <- resize(refGenomeBin.gr, width = 1, fix = "start")

refGenomeMarks.df <- data.frame(genoExpandStretch(x.gr = refGenomeBin.gr,synthGenome = newSynthRefShift, 
                                                  expandedSeqlengths = seqlengths(synthBinNorm.gr)))


# plot missing Ranges
missingRange <- gaps(synthBinNorm.gr)
missingRange <- missingRange[strand(missingRange) == "*"]

# gap Names

if(specRef == "hg19"){
  gapType = c("refIns", "refDel", "queIns", "queDel")
}
if(specRef == "mm10"){
  gapType = c("queIns", "queDel", "refIns", "refDel")
}
gapNames = c("hg19 gain", "hg19 loss", "mm10 gain", "mm10 loss")
gapCols = c("aquamarine3", "darkblue", "red", "purple")
gapCols = RColorBrewer::brewer.pal(4, "Dark2")


# select chromosomes
chrAll <- seqlevels(synthBinNorm.gr)
chrAll <- chrAll[-grep("_", chrAll)]

chrAll <- chrAll[chrAll!="chrM"]
chrAll <- chrAll[chrAll!="chrY"]
#chrAll <- chrAll[chrAll!="chrX"]


# get important ranges
intSigRanges = NULL
for(i in 1:length(gapType)){
  if(gapType[i] == "refDel" & specRef == "hg19"){
    intSigRange <- repSigRanges[[gapType[i]]]
    intSigRange <- intSigRange[seqnames(intSigRange) %in% chrAll]
    intSigRanges <- c(intSigRanges, list(intSigRange))
  }else if(gapType[i] == "queDel" & specRef == "mm10"){
    intSigRange <- repSigRanges[[gapType[i]]]
    intSigRange <- intSigRange[seqnames(intSigRange) %in% chrAll]
    intSigRanges <- c(intSigRanges, list(intSigRange))
  }else{
    intSigRange <- intersect(noRepSigRanges[[gapType[i]]], repSigRanges[[gapType[i]]])
    intSigRange <- intSigRange[seqnames(intSigRange) %in% chrAll]
    intSigRanges <- c(intSigRanges, list(intSigRange))
  }
}
names(intSigRanges) <- gapType




refGenomeMarks.df <- refGenomeMarks.df[refGenomeMarks.df$seqnames %in% chrAll,]
missingRange <- missingRange[seqnames(missingRange) %in% chrAll]

##### Generate table

gapOL <- matrix(NA,nrow = 4, ncol = 4,
                dimnames = list(gapType,gapType))
for(i in 1:length(gapType)){
  for(j in 1:length(gapType)){
    gapOL[i,j] <- sum(width(intersect(intSigRanges[[gapType[i]]], intSigRanges[[gapType[j]]])))
  }
}
gapOL <- gapOL/200000


sigOL <- matrix(NA,nrow = 4, ncol = 4,
                dimnames = list(gapType,gapType))
EOL <- matrix(NA,nrow = 4, ncol = 4,
              dimnames = list(gapType,gapType))
for(i in gapType){
  for(j in gapType){
    if(i == j)next
    x =  as.table(
      rbind(
        c(gapOL[i,j], 
          gapOL[i,i] - gapOL[i,j]),
        c(gapOL[j,j] - gapOL[i,j], 
          length(synthBinNorm.gr[seqnames(synthBinNorm.gr) %in% chrAll]) - gapOL[i,i] - gapOL[j,j] + gapOL[i,j])
      )
    )
    dNames <- list(c("sig", "norm"), c("sig", "norm"))
    names(dNames) <- c(i,j)
    dimnames(x)  <- dNames
    sigOL[i,j] <- fisher.test(x)$p.value
    EOL[i,j] <- round(fisher.test(x)$estimate,2)
  }
}

# expected overlap

olMat <-round(gapOL/diag(gapOL) * 100, digits = 2)
olMat[sigOL < .05 & sigOL > .01 & !is.na(sigOL)] <- paste(olMat[sigOL < .05 & sigOL > .01 & !is.na(sigOL)], "*")
olMat[sigOL < .01 & !is.na(sigOL)] <- paste(olMat[sigOL < .01 & !is.na(sigOL)], "**")
olMat[olMat == "100"] <- ""
dimnames(olMat) <- list(gapNames, gsub(" ", "\n",gapNames))

OoE <- matrix(as.character(EOL), nrow = 4)
OoE[!is.na(OoE)] <- paste("\n(",OoE[!is.na(OoE)],")", sep = "")
# if i can colour it
olMat2 <- olMat

for(i in 1:4){
  for(j in 1:4){
    if(i == j){next}
    olMat2[i,j] <- paste(olMat[i,j], OoE[i,j], sep = "")
  }
}

pdf(file = paste("~/Desktop/RTN_domains/RTN_domain_plots/netGainLoss/hotspotPlot/genomeDistribution",specRef,".pdf", sep = ""), width = 8,height = 11)

plot(1,xlim = c(0,max(seqlengths(synthBinNorm.gr))), ylim = c(0,(length(chrAll))*5),
     axes = FALSE, ylab = "", xlab = "", type = "n", main = specRef)

mtext(text = chrAll, side = 2, at = ((length(chrAll):1)*5 )-2, las = 2)

# plot missing areas
rect(xleft = start(missingRange),
     xright = end(missingRange),
     ytop = (-(as.integer(as.factor(seqnames(missingRange)))) + length(chrAll) + 1) * 5   ,
     ybottom = (-(as.integer(as.factor(seqnames(missingRange)))) + length(chrAll) + 1) * 5 -4,
     border = FALSE, col = scales::alpha(alpha = 1, "grey85"), density = -1)

#plot chr
rect(xleft = 0, 
     xright = seqlengths(synthBinNorm.gr)[chrAll], 
     ytop = (length(chrAll):1)*5,
     ybottom = ((length(chrAll):1)*5 )-4)

# place hotspots
for(i in 1:4){
  rect(xleft = start(reduce(intSigRanges[[gapType[i]]])),
       xright = end(reduce(intSigRanges[[gapType[i]]])),
       ytop = (-(as.integer(as.factor(seqnames(reduce(intSigRanges[[gapType[i]]]))))) + length(chrAll) +1) * 5 - i +1  ,
       ybottom = (-(as.integer(as.factor(seqnames(reduce(intSigRanges[[gapType[i]]]))))) + length(chrAll) ) * 5 + (5-i),
       border = scales::alpha(gapCols[i], 1), col = scales::alpha(gapCols[i], 1), density = -1)
}

#place Marks
segments(x0 = refGenomeMarks.df$start,
         x1 = refGenomeMarks.df$start,
         y0= (-as.integer(as.factor(refGenomeMarks.df$seqnames)) + length(chrAll) + 1) * 5,
         y1= ((-as.integer(as.factor(refGenomeMarks.df$seqnames)) + length(chrAll) + 1)  * 5) - 4,
         lty = 1)  


legend("right", fill = gapCols, legend = gapNames, bty = "n", title = "Hotspots")


# grid themes
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.8)),
  colhead = list(fg_params=list(cex = 0.8)),
  rowhead = list(fg_params=list(cex = 0.8)))

t1 <- gridExtra::tableGrob(olMat2, theme = mytheme)
title <- textGrob("Hotspot overlap (%)",gp=gpar(fontsize=12))
padding <- unit(5,"mm")

table <- gtable_add_rows(
  t1, 
  heights = grobHeight(title) + padding,
  pos = 0)
table <- gtable_add_grob(
  table, 
  title, 
  1, 1, 1, ncol(table))


pushViewport(viewport(y= .285, x = .75 ))
grid.draw(table)


dev.off()


