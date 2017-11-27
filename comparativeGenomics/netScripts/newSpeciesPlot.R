
library(dplyr)
library(GenomicRanges)

rm(list = ls())


genomes = c(ref = "hg19", que = "mm10")

loadPath <- paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/ingroupSpecies/",
                  genomes["ref"],".ingroup.RData", sep = "")

load(loadPath)

newInGroups$species <- factor(newInGroups$species, 
                              levels = c("micMur", "tarSyr", "papHam", "panTro", "hg19",
                                         "ochPri", "dipOrd", "rn", "mm10"))
MyaDivergence <- c(74, 67.1, 29.44, 6.65, 0, 82, 69.9, 20.9, 0)
names(MyaDivergence) <- levels(newInGroups$species)


newInGroups$MyaDivergence <- MyaDivergence[as.character(newInGroups$species)]



df <- data.frame(newInGroups)

a <- group_by(df, type, species, lineage, MyaDivergence) %>%
  summarize(sum(width))

a$timePeriod <- c(a$MyaDivergence[nrow(a)], a$MyaDivergence[1:(nrow(a) - 1)])
a$timePeriod[a$timePeriod == 0] = 90
a$timePeriod <- a$timePeriod - a$MyaDivergence
a$rate <- (a$`sum(width)`/1e6)/a$timePeriod


genomeChoice <- data.frame(genome = c("hg19", "mm10"), 
                           lineage = c("primate", "rodent"),
                           insType = c("refIns", "queIns"),
                           delType = c("refDel", "queDel"),
                           que = c("mm10", "hg19"))


pdf(file = paste("~/Desktop/","TemporalDANTurnover.pdf", sep = ""), onefile = TRUE)

for(i in 1:2){
  
  dataPointsGain <- data.frame(x = a[a$type == genomeChoice$insType[i],"MyaDivergence"] + a$timePeriod[a$type == genomeChoice$insType[i]],
                               y = a[a$type == genomeChoice$insType[i] ,"rate"])
  dataPointsGain[nrow(dataPointsGain) +1,] <- c(0, dataPointsGain[nrow(dataPointsGain),"rate"])
  
  dataPointsLoss <- data.frame(x = a[a$type == genomeChoice$delType[i],"MyaDivergence"] + a$timePeriod[a$type == genomeChoice$delType[i]],
                               y = a[a$type == genomeChoice$delType[i] ,"rate"])
  dataPointsLoss[nrow(dataPointsLoss) +1,] <- c(0, dataPoints[nrow(dataPointsLoss),"rate"])
  
  primateMYA <- c(90,a$MyaDivergence[a$lineage == genomeChoice$lineage[i] & a$type == genomeChoice$insType[i]])
  primateNames <- c(genomeChoice$que[i], 
                    as.character(a$species[a$lineage == genomeChoice$lineage[i] & a$type == genomeChoice$insType[i]]))
  distMat <- matrix(data = primateMYA, nrow = length(primateNames), ncol = length(primateNames), 
                    dimnames = list(primateNames, primateNames))
  distMat <- as.dist(t(distMat))
  
  layout(c(1,2), heights = c(1.5,1))
  par(mar = c(0,5,5,5))
  plot(dataPointsGain, ylim = c(0,60), xlim = c(90,0), type = "s", 
       xaxs = "i", xaxt = "n", ylab = "DNA turnover (Mb / MY)",
       main = genomeChoice$genome[i], lwd = 2)
  lines(dataPointsLoss, col = 2, type = "s", lwd = 2)
  grid(ny = 0, nx = 9)
  legend("topright", lty = 1, col= c(1,2), legend = c("Gain", "Loss"), bty = "n")
  
  par(mar = c(5,5,0,5))
  plot(as.dendrogram(hclust(distMat)), horiz = TRUE,xlim = c(90,0), xaxs = "i", xlab = "MYA", xaxt = "n")
  axis(side = 1, at = seq(0,90,10))
  grid(ny = 0, nx = 9)
  
}

dev.off()





pdf(file = paste("~/Desktop/","AmountDANTurnover.pdf", sep = ""), onefile = TRUE)

for(i in 1:2){
  
  dataPointsGain <- data.frame(x = a[a$type == genomeChoice$insType[i],"MyaDivergence"] + a$timePeriod[a$type == genomeChoice$insType[i]],
                               y = a$`sum(width)`[a$type == genomeChoice$insType[i]] / 1e6)
  dataPointsGain[nrow(dataPointsGain) +1,] <- c(0, dataPointsGain[nrow(dataPointsGain),"y"])
  
  dataPointsLoss <- data.frame(x = a[a$type == genomeChoice$delType[i],"MyaDivergence"] + a$timePeriod[a$type == genomeChoice$delType[i]],
                               y = a$`sum(width)`[a$type == genomeChoice$delType[i]] / 1e6)
  dataPointsLoss[nrow(dataPointsLoss) +1,] <- c(0, dataPoints[nrow(dataPointsLoss),"y"])
  
  primateMYA <- c(90,a$MyaDivergence[a$lineage == genomeChoice$lineage[i] & a$type == genomeChoice$insType[i]])
  primateNames <- c(genomeChoice$que[i], 
                    as.character(a$species[a$lineage == genomeChoice$lineage[i] & a$type == genomeChoice$insType[i]]))
  distMat <- matrix(data = primateMYA, nrow = length(primateNames), ncol = length(primateNames), 
                    dimnames = list(primateNames, primateNames))
  distMat <- as.dist(t(distMat))
  
  layout(c(1,2), heights = c(1.5,1))
  par(mar = c(0,5,5,5))
  plot(dataPointsGain, ylim = c(0,460), xlim = c(90,0), type = "s", 
       xaxs = "i", xaxt = "n", ylab = "DNA turnover (Mb)",
       main = genomeChoice$genome[i], lwd = 2)
  lines(dataPointsLoss, col = 2, type = "s", lwd = 2)
  grid(ny = 0, nx = 9)
  legend("topright", lty = 1, col= c(1,2), legend = c("Gain", "Loss"), bty = "n")
  
  par(mar = c(5,5,0,5))
  plot(as.dendrogram(hclust(distMat)), horiz = TRUE,xlim = c(90,0), xaxs = "i", xlab = "MYA", xaxt = "n")
  axis(side = 1, at = seq(0,90,10))
  grid(ny = 0, nx = 9)
  
}

dev.off()




### so now its time to bin genomes and do correlation analysis 



bins <- unlist(slidingWindows(GRanges(seqinfo(newInGroups)), width = 1e7, step = 1e7))
ol <- findOverlaps(bins, newInGroups)
pInt <- pintersect(newInGroups[subjectHits(ol)], bins[queryHits(ol)])

df <- data.frame(binID = queryHits(ol), width = width(pInt), mcols(pInt))

df0 <- group_by(df, binID, type, species, lineage) %>% 
  summarise(width = sum(width))
df0$speciesType <- paste(df0$species, df0$type, sep = ".")


library(reshape)
df1 <- cast(df0[,c("binID", "speciesType", "width")], binID ~ speciesType)
df1[is.na(df1)] <- 0

blueRed <- colorRampPalette(c("blue", "white", "red"))

heatmap(cor(df1), scale = "none", col = blueRed(20), zlim = c(-1,1))

layout(1:5)
par(mar = c(0,1,0,1))
whiteRed <- colorRampPalette(c("white", "red"))
image(matrix(df1$hg19.refDel), col = whiteRed(20))
image(matrix(df1$panTro.refDel), col = whiteRed(20))
image(matrix(df1$papHam.refDel), col = whiteRed(20))
image(matrix(df1$tarSyr.refDel), col = whiteRed(20))
image(matrix(df1$micMur.refDel), col = whiteRed(20))


plot(df1$hg19.refIns)
plot(df1$panTro.refIns)
plot(df1$papHam.refIns)
plot(df1$tarSyr.refIns)
plot(df1$micMur.refIns)



# need to adjust bins for gaps and bases 






df <- data.frame(newInGroups)

a <- group_by(df, type, species, lineage) %>%
  summarize(sum(width))


hg19years <- c(15, 7, 40, 22)
mm10years <- c(8, 10 , 54)

# rn is rel3
# dipOrd is rel2
# ochPri is rel1
# 
# 
# micMur is rel1
# tarSyr is rel2
# papHam is rel3
# panTro is rel4