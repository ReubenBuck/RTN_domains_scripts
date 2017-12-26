
library(dplyr)
library(GenomicRanges)

rm(list = ls())

options(stringsAsFactors = FALSE)
genomes = c(ref = "mm10", que = "hg19")

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

if(genomes["ref"] == "hg19"){
  genomeChoice <- data.frame(genome = c("hg19", "mm10"), 
                           lineage = c("primate", "rodent"),
                           insType = c("refIns", "queIns"),
                           delType = c("refDel", "queDel"),
                           que = c("mm10", "hg19"))
}else if(genomes["ref"]=="mm10"){
  genomeChoice <- data.frame(genome = c("mm10", "hg19"), 
                             lineage = c("rodent", "primate"),
                             insType = c("refIns", "queIns"),
                             delType = c("refDel", "queDel"),
                             que = c("hg19", "mm10")) 
}



pdf(file = paste("~/Desktop/","TemporalDANTurnover",genomes["ref"] ,".pdf", sep = ""), onefile = TRUE)
for(i in 1:2){
  
  dataPointsGain <- data.frame(x = a[a$type == genomeChoice$insType[i],"MyaDivergence"] + a$timePeriod[a$type == genomeChoice$insType[i]],
                               y = a[a$type == genomeChoice$insType[i] ,"rate"])
  dataPointsGain[nrow(dataPointsGain) +1,] <- c(0, dataPointsGain[nrow(dataPointsGain),"rate"])
  
  dataPointsLoss <- data.frame(x = a[a$type == genomeChoice$delType[i],"MyaDivergence"] + a$timePeriod[a$type == genomeChoice$delType[i]],
                               y = a[a$type == genomeChoice$delType[i] ,"rate"])
  dataPointsLoss[nrow(dataPointsLoss) +1,] <- c(0, dataPointsLoss[nrow(dataPointsLoss),"rate"])
  
  primateMYA <- c(90,a$MyaDivergence[a$lineage == genomeChoice$lineage[i] & a$type == genomeChoice$insType[i]])
  primateNames <- c(genomeChoice$que[i], 
                    as.character(a$species[a$lineage == genomeChoice$lineage[i] & a$type == genomeChoice$insType[i]]))
  distMat <- matrix(data = primateMYA, nrow = length(primateNames), ncol = length(primateNames), 
                    dimnames = list(fullSpeciesNames[primateNames], fullSpeciesNames[primateNames]))
  distMat <- as.dist(t(distMat))
  
  layout(c(1,2), heights = c(1.5,1))
  par(mar = c(0,7,5,7))
  plot(dataPointsGain, ylim = c(0,60), xlim = c(90,0), type = "s", 
       xaxs = "i", xaxt = "n", ylab = "DNA turnover (Mb / MY)",
       main = genomeChoice$genome[i], lwd = 2)
  lines(dataPointsLoss, col = 2, type = "s", lwd = 2)
  grid(ny = 0, nx = 9)
  legend("topright", lty = 1, col= c(1,2), legend = c("Gain", "Loss"), bty = "n")
  
  par(mar = c(5,7,0,7))
  plot(as.dendrogram(hclust(distMat)), horiz = TRUE,xlim = c(90,0), xaxs = "i", xlab = "MYA", xaxt = "n")
  axis(side = 1, at = seq(0,90,10))
  grid(ny = 0, nx = 9)
  
}

dev.off()





pdf(file = paste("~/Desktop/","AmountDANTurnover",genomes["ref"] ,".pdf", sep = ""), onefile = TRUE)

for(i in 1:2){
  
  dataPointsGain <- data.frame(x = a[a$type == genomeChoice$insType[i],"MyaDivergence"] + a$timePeriod[a$type == genomeChoice$insType[i]],
                               y = a$`sum(width)`[a$type == genomeChoice$insType[i]] / 1e6)
  dataPointsGain[nrow(dataPointsGain) +1,] <- c(0, dataPointsGain[nrow(dataPointsGain),"y"])
  
  dataPointsLoss <- data.frame(x = a[a$type == genomeChoice$delType[i],"MyaDivergence"] + a$timePeriod[a$type == genomeChoice$delType[i]],
                               y = a$`sum(width)`[a$type == genomeChoice$delType[i]] / 1e6)
  dataPointsLoss[nrow(dataPointsLoss) +1,] <- c(0, dataPointsLoss[nrow(dataPointsLoss),"y"])
  
  primateMYA <- c(90,a$MyaDivergence[a$lineage == genomeChoice$lineage[i] & a$type == genomeChoice$insType[i]])
  primateNames <- c(genomeChoice$que[i], 
                    as.character(a$species[a$lineage == genomeChoice$lineage[i] & a$type == genomeChoice$insType[i]]))
  distMat <- matrix(data = primateMYA, nrow = length(primateNames), ncol = length(primateNames), 
                    dimnames = list(fullSpeciesNames[primateNames], fullSpeciesNames[primateNames]))
  distMat <- as.dist(t(distMat))
  
  layout(c(1,2), heights = c(1.5,1))
  par(mar = c(0,7,5,7))
  plot(dataPointsGain, ylim = c(0,460), xlim = c(90,0), type = "s", 
       xaxs = "i", xaxt = "n", ylab = "DNA turnover (Mb)",
       main = genomeChoice$genome[i], lwd = 2)
  lines(dataPointsLoss, col = 2, type = "s", lwd = 2)
  grid(ny = 0, nx = 9)
  legend("topright", lty = 1, col= c(1,2), legend = c("Gain", "Loss"), bty = "n")
  
  par(mar = c(5,7,0,7))
  plot(as.dendrogram(hclust(distMat)), horiz = TRUE,xlim = c(90,0), xaxs = "i", xlab = "MYA", xaxt = "n")
  axis(side = 1, at = seq(0,90,10))
  grid(ny = 0, nx = 9)
  
}

dev.off()









### so now its time to bin genomes and do correlation analysis 



bins <- unlist(slidingWindows(GRanges(seqinfo(newInGroups)), width = 1e6, step = 1e6))
ol <- findOverlaps(bins, newInGroups)
pInt <- pintersect(newInGroups[subjectHits(ol)], bins[queryHits(ol)])

df <- data.frame(binID = queryHits(ol), width = width(pInt), mcols(pInt))

df0 <- group_by(df, binID, type, species, lineage) %>% 
  summarise(width = sum(width))
df0$speciesType <- paste(df0$species, df0$type, sep = ".")


library(reshape)
df1 <- cast(df0[,c("binID", "speciesType", "width")], binID ~ speciesType)
df1[is.na(df1)] <- 0

# here we have bp per bin


# below I'm experimenting with ways to plot the relationships between this data 

# need to also chreat an aditive model

# also need to repeat this in mouse 

blueRed <- colorRampPalette(c("blue", "white", "red"))

colnames(df1) <- gsub(pattern = "ref", "", x = colnames(df1))
colnames(df1) <- gsub(pattern = "que", "", x = colnames(df1))


newNames <- c(paste(levels(newInGroups$species), "Ins", sep = "."), 
              paste(levels(newInGroups$species), "Del", sep = "."))

df2 <- df1[,as.character(newNames)]


df3 <- log10(df2)
df3 <- as.matrix(df3)
df3[is.infinite(df3)] <- NA
df3 <- df3[complete.cases(df3),]
corMat <- cor(df3)
colnames(corMat) <- colnames(df2)
rownames(corMat) <- colnames(df2)



refIns <- 1:5
queIns <- 6:9
refDel <- 10:14
queDel <- 15:18

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

layout(matrix((4*4):1, ncol = 4, byrow = TRUE)[,4:1])
par(oma = c(5,5,5,5), mar = c(3,3,3,3))

for(i in 1:nrow(gridLayout)){

  if(gridLayout[i,1] == "newPlot"){
    plot.new()
  }else if(gridLayout$rows[i] == gridLayout$cols[i]){
    cols <- get(as.character(gridLayout$cols[i]))
    rows <- get(as.character(gridLayout$rows[i]))
    corMat1 <- corMat[rows, cols]
    corMat1[upper.tri(corMat1, diag = TRUE)] <- NA
    image(1:nrow(corMat1), 1:ncol(corMat1),
          corMat1, 
          col = blueRed(20), zlim = c(-1,1),
          xaxt = "n", yaxt = "n", bty = "n",xlab = "", ylab = "")
    axis(side = 1, at = -.5:(nrow(corMat1) + 1), 
         labels = c("",90, timing[rownames(corMat1)]))
    axis(side = 2, at = -.5:(ncol(corMat1) + 1), 
         labels = c("",90, timing[colnames(corMat1)]), las = 2)
  }else{
    cols <- get(as.character(gridLayout$cols[i]))
    rows <- get(as.character(gridLayout$rows[i]))
    corMat1 <- corMat[rows, cols]
    image(1:nrow(corMat1), 1:ncol(corMat1), 
          corMat1, 
          col = blueRed(20), zlim = c(-1,1),
          xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
    axis(side = 1, at = -.5:(nrow(corMat1) + 1), 
         labels = c("",90, timing[rownames(corMat1)]))
    axis(side = 2, at = -.5:(ncol(corMat1) + 1), 
         labels = c("",90, timing[colnames(corMat1)]), las = 2)
  }
}

mtext(side = 1, outer = TRUE, at = seq(1/(4 *2),1 - (1/(4 *2)) ,length.out = 4),
      text = c("hg19 gain", "hg19 loss", "mm10 gain", "mm10 loss"))

mtext(side = 4, outer = TRUE, at = seq(1/(4 *2),1 - (1/(4 *2)) ,length.out = 4),
      text = c("hg19 gain", "hg19 loss", "mm10 gain", "mm10 loss"))


# now with a dat matrix we can get the interesting stuff out and plot it
# lets thing about the analysis required




# program this in terms of groups
# refIns 
# refDels
#Que Ins
# what does the image function to do our matirx


