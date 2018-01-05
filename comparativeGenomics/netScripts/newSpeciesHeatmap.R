
library(dplyr)
library(GenomicRanges)
library(igraph)
library(spdep)

rm(list = ls())




ZtoP <- function(Zs){
  2*pnorm(-abs(Zs))
}

ZsigAdj <- function(Zs, cutoff){
  score = p.adjust(p = ZtoP(Zs), method = "BH")
  return(score <= cutoff)
}





options(stringsAsFactors = FALSE)
genomes = c(ref = "hg19", que = "mm10")


# load genome information
loadPath <- paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/ingroupSpecies/",
                  genomes["ref"],".ingroup.RData", sep = "")

load(loadPath)

newInGroups$species <- factor(newInGroups$species, 
                              levels = c("micMur", "tarSyr", "papHam", "panTro", "hg19",
                                         "ochPri", "dipOrd", "rn", "mm10"))
MyaDivergence <- c(74, 67.1, 29.44, 6.65, 0, 82, 69.9, 20.9, 0)
names(MyaDivergence) <- levels(newInGroups$species)


newInGroups$MyaDivergence <- MyaDivergence[as.character(newInGroups$species)]

fullSpeciesNames <- c("pika", "kangaroo rat", "rat", "mouse", "mouse lemur", "tarsier", "baboon", "chimpanzee", "human")
names(fullSpeciesNames) <- c("ochPri", "dipOrd", "rn", "mm10", "micMur", "tarSyr", "papHam", "panTro", "hg19")

df <- data.frame(newInGroups)

a <- group_by(df, type, species, lineage, MyaDivergence) %>%
  summarize(sum(width))

a$timePeriod <- c(a$MyaDivergence[nrow(a)], a$MyaDivergence[1:(nrow(a) - 1)])
a$timePeriod[a$timePeriod == 0] = 90
a$timePeriod <- a$timePeriod - a$MyaDivergence
a$rate <- (a$`sum(width)`/1e6)/a$timePeriod
a$fullName <- fullSpeciesNames[a$species]

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")


load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/",genomes["ref"],".",genomes["que"],".netData.RData",sep = ""))
load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/",genomes["ref"],".stretch.RData", sep = ""))
load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/shiftData/",genomes["ref"],".expand.breaks.RData", sep = ""))



# annotate extra files

refMissingGaps.gr <- GenomicRanges::setdiff(refFillGaps.gr, refGap.gr, ignore.strand = TRUE)
all(!overlapsAny(refMissingGaps.gr, refGap.gr))
mcols(refMissingGaps.gr)$lineage <- NA
mcols(refMissingGaps.gr)$species <- NA
mcols(refMissingGaps.gr)$type = "missingGap"
mcols(refMissingGaps.gr)$MyaDivergence = NA

refSeqGaps.gr <- reduce(refSeqGaps.gr)
mcols(refSeqGaps.gr)$lineage <- NA
mcols(refSeqGaps.gr)$species <- NA
mcols(refSeqGaps.gr)$type = "seqGap"
mcols(refSeqGaps.gr)$MyaDivergence = NA

missingGenome.gr <- sort(c(refMissingGaps.gr, refSeqGaps.gr))
missingGenome.gr <- genoExpandBreak(missingGenome.gr, newSynthRefShift, seqlengths(stretchedRef.gr))

refSynth <- sort(c(newInGroups,missingGenome.gr))



## create genomic bins
binSize <- 2e5

bins <- unlist(slidingWindows(GRanges(seqinfo(refSynth)), width = binSize, step = binSize))
ol <- findOverlaps(bins, refSynth)
pInt <- pintersect(refSynth[subjectHits(ol)], bins[queryHits(ol)])

df <- data.frame(binID = queryHits(ol), width = width(pInt), mcols(pInt))

df0 <- group_by(df, binID, type, species, lineage) %>% 
  summarise(width = sum(width))
df0$speciesType <- paste(df0$species, df0$type, sep = ".")


library(reshape)
df1 <- cast(df0[,c("binID", "speciesType", "width")], binID ~ speciesType)
df1[is.na(df1)] <- 0


# attache missing data to bin info and clean
mcols(bins) <- df1
mcols(bins)$nonNbases <- width(bins) - rowSums(df1[,c("NA.missingGap", "NA.seqGap")])
bins <- bins[mcols(bins)$nonNbases > .75*binSize]
mcols(bins) <- (data.frame(mcols(bins)) / mcols(bins)$nonNbases) * binSize
mcols(bins) <- data.frame(mcols(bins))[,!(names(mcols(bins)) %in% c("binID", "NA.seqGap", "NA.missingGap", "nonNbases"))]

# here we have the bin data with all the important information

# we can start doign hotspot detection



# hotspot identification


#for(m in mehtods){
synthBinNorm.gr <- bins
gapChoice <- names(mcols(bins))
df <- data.frame(bins)

j = 3
ol<-findOverlaps(synthBinNorm.gr, maxgap = j*width(synthBinNorm.gr)[1])
ol <- ol[!(isRedundantHit(ol))]
# remove NA hits
#ol <- ol[!is.na(df$refIns[queryHits(ol)])]
#ol <- ol[!is.na(df$refIns[subjectHits(ol)])]
olMat <- data.frame(ol)
G <- graph.data.frame(d = olMat,directed=FALSE)
#weight <- olMat$subjectHits - olMat$queryHits
#weight <- -(weight - (max(weight) + 1)) / max(weight)
weight = rep(1, length(ol))
weight[isSelfHit(ol)] <- 0
E(G)$weight <- weight
#A <- as_adjacency_matrix(G,type="both",names=FALSE,sparse=TRUE,edges = FALSE)
A <- as_adjacency_matrix(G,type="both",names=FALSE,sparse=TRUE,edges = FALSE, attr = "weight")
wMat <- mat2listw(A)
wMat = nb2listw(include.self(wMat$neighbours))

dfGscore <- data.frame(synthBinNorm.gr)

for(i in 1:length(gapChoice)){
  score <- df[,gapChoice[i]]
  G <- localG(x = score,wMat)
  dfGscore[,gapChoice[i]] <- G
}
#nDistList <- c(nDistList, list(dfGscore))


sigRanges = NULL
sigRangeDF <- NULL
for(i in 1:length(gapChoice)){
  z <- ZsigAdj(dfGscore[,gapChoice[i]],cutoff = .05) & dfGscore[,gapChoice[i]] > 0
  sigRanges <- c(sigRanges , list((synthBinNorm.gr[z])))
  sigRangeDF <- cbind(sigRangeDF, z)
}
names(sigRanges) <- gapChoice
colnames(sigRangeDF) <- gapChoice

#sigRangeList <- c(sigRangeList, list(sigRanges))




# now we have significant ranges
# need to do the all vall overlap

oMat <- pMat <- matrix(NA, ncol = ncol(sigRangeDF), nrow = ncol(sigRangeDF),
                       dimnames = list(gapChoice,gapChoice))
for(i in 1:ncol(sigRangeDF)){
  for(j in 1:ncol(sigRangeDF)){
    tab <- table(sigRangeDF[,i], sigRangeDF[,j])
    tab <- tab[c("TRUE", "FALSE"),c("TRUE", "FALSE")]
    fish = fisher.test(tab)
    oMat[i,j] <- fish$estimate
    pMat[i,j] <- fish$p.value
  }
}

pMatAdj <- matrix(p.adjust(pMat, method = "bonferroni") < .05, ncol = ncol(pMat))



oMatLog <- log2(oMat)
oMatLog[is.infinite(oMatLog)] <- NA
oMatLog[!pMatAdj] <- NA
oMatLog[upper.tri(oMatLog)] <- NA


lim <- max(abs(oMatLog), na.rm = TRUE)



# here is where the plotting happens

df1 <- oMatLog

blueRed <- colorRampPalette(c("blue", "white", "red"))

colnames(df1) <- gsub(pattern = "ref", "", x = colnames(df1))
colnames(df1) <- gsub(pattern = "que", "", x = colnames(df1))


newNames <- c(paste(levels(newInGroups$species), "Ins", sep = "."), 
              paste(levels(newInGroups$species), "Del", sep = "."))

df2 <- df1[,as.character(newNames)]


corMat <- df2
colnames(corMat) <- colnames(df2)
rownames(corMat) <- colnames(df2)



refIns <- paste(c("micMur", "tarSyr", "papHam", "panTro", "hg19"), "Ins", sep = ".")
queIns <- paste(c("ochPri", "dipOrd", "rn", "mm10"), "Ins", sep = ".")
refDel <- paste(c("micMur", "tarSyr", "papHam", "panTro", "hg19"), "Del", sep = ".")
queDel <- paste(c("ochPri", "dipOrd", "rn", "mm10"), "Del", sep = ".")

gridLayout <- data.frame(rows = c("refIns", "refDel", "queIns", "queDel",
                                  "newPlot", "refDel", "queIns", "queDel",
                                  "newPlot", "newPlot", "queIns", "queDel",
                                  "newPlot", "newPlot", "newPlot", "queDel"),
                         cols = c("refIns", "refIns", "refIns", "refIns",
                                  "newPlot", "refDel", "refDel", "refDel",
                                  "newPlot", "newPlot", "queIns", "queIns",
                                  "newPlot", "newPlot", "newPlot", "queDel"))


timing <- a$MyaDivergence
names(timing) <- paste(a$species, a$type, sep = ".")
names(timing) <- gsub(pattern = "ref", "", names(timing))
names(timing) <- gsub(pattern = "que", "", names(timing))


pdf(file = paste("~/Documents/RTN_domain_writing/manuscriptRound2/bioArxiv/main/figs/newSpecies/newSpeciesHeatmap", genomes["ref"], ".pdf", sep = ""),
    height = 10, width = 10)
layout(matrix((4*4):1, ncol = 4, byrow = TRUE)[,4:1])
par(oma = c(5,5,5,5), mar = c(3,3,3,3))

for(i in 1:nrow(gridLayout)){
  
  if(i == 9){
    par(mar = c(12,3,3,3))
     image(x = -6:6, y = 1,
           matrix(-lim:lim), 
           col = blueRed(20), 
           main = "OR (log2)", 
           yaxt = "n", xlab = "", ylab = "")
    par(mar = c(3,3,3,3))
  }else if(gridLayout[i,1] == "newPlot"){
    plot.new()
  }else if(gridLayout$rows[i] == gridLayout$cols[i]){
    cols <- get(as.character(gridLayout$cols[i]))
    rows <- get(as.character(gridLayout$rows[i]))
    corMat1 <- corMat[rows, cols]
    corMat1[upper.tri(corMat1, diag = TRUE)] <- NA
    greyMat1 <- matrix(1,nrow = length(rows),ncol=length(cols))
    greyMat1[upper.tri(greyMat1, diag = TRUE)] <- NA
    image(1:nrow(corMat1), 1:ncol(corMat1), 
          greyMat1, 
          col = "grey90", zlim = c(-1,1),
          xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
    image(1:nrow(corMat1), 1:ncol(corMat1), add = TRUE,
          corMat1, 
          col = blueRed(20), zlim = c(-lim,lim),
          xaxt = "n", yaxt = "n", bty = "n",xlab = "", ylab = "")
    axis(side = 1, at = -.5:(nrow(corMat1) + 1), 
         labels = c("",90, round(timing[rownames(corMat1)])))
    axis(side = 2, at = -.5:(ncol(corMat1) + 1), 
         labels = c("",90, round(timing[colnames(corMat1)])), las = 2)
  }else{
    cols <- get(as.character(gridLayout$cols[i]))
    rows <- get(as.character(gridLayout$rows[i]))
    corMat1 <- corMat[rows, cols]
    greyMat1 <- matrix(1,nrow = length(rows),ncol=length(cols))
    image(1:nrow(corMat1), 1:ncol(corMat1), 
          greyMat1, 
          col = "grey90", zlim = c(-1,1),
          xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
    image(1:nrow(corMat1), 1:ncol(corMat1), 
          corMat1, 
          col = blueRed(20), zlim = c(-lim,lim),
          xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "", add = TRUE)
    axis(side = 1, at = -.5:(nrow(corMat1) + 1), 
         labels = c("",90, round(timing[rownames(corMat1)])))
    axis(side = 2, at = -.5:(ncol(corMat1) + 1), 
         labels = c("",90, round(timing[colnames(corMat1)])), las = 2)
  }
}

if(genomes["ref"] == "hg19"){
mtext(side = 1, outer = TRUE, at = seq(1/(4 *2),1 - (1/(4 *2)) ,length.out = 4),
      text = paste(c(genomes["ref"], genomes["ref"], genomes["que"], genomes["que"]),
                   c("gain", "loss", "gain", "loss")))

mtext(side = 4, outer = TRUE, at = seq(1/(4 *2),1 - (1/(4 *2)) ,length.out = 4),
      text = paste(c(genomes["ref"], genomes["ref"], genomes["que"], genomes["que"]),
                   c("gain", "loss", "gain", "loss")))
}else{
mtext(side = 1, outer = TRUE, at = seq(1/(4 *2),1 - (1/(4 *2)) ,length.out = 4),
      text = paste(c(genomes["que"], genomes["que"], genomes["ref"], genomes["ref"]),
                   c("gain", "loss", "gain", "loss")))

mtext(side = 4, outer = TRUE, at = seq(1/(4 *2),1 - (1/(4 *2)) ,length.out = 4),
      text = paste(c(genomes["que"], genomes["que"], genomes["ref"], genomes["ref"]),
                   c("gain", "loss", "gain", "loss")))
}
mtext(side = 1, outer = TRUE,line = 2, text = "MYA" )
mtext(side = 4, outer = TRUE,line = 2, text = "MYA" )

dev.off()





