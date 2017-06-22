
# use regioner
rm(list = ls())

ZtoP <- function(Zs){
  2*pnorm(-abs(Zs))
}

ZsigAdj <- function(Zs, nb, cutoff){
  score = p.adjustSP(p = ZtoP(Zs), nb = nb, method = "BH")
  return(score <= cutoff)
}



library(igraph)
library(spdep)
library(regioneR)



for(spec in c("hg19", "mm10")){
  load(file = paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/", 
                    spec,".synthBin.genome.variable.RData", sep = ""))
  
  
  df<- data.frame(synthBin.gr)
  
  colChoice <- c("dnasePeaks","ladCov","exon", "intron","ctcfMotif", 
                 "L1Motif", "prdm9Motif", "ancient", "new_L1", "new_SINE", 
                 "old_L1", "recHotSpot","cpgCov")
  
  remainRefBases <- width(synthBin.gr)[1] - (df$queIns + df$refDel + df$seqGap)
  df[remainRefBases > 0,colChoice] <- (df[remainRefBases > 0,colChoice]/remainRefBases[remainRefBases > 0]) * width(synthBin.gr)[1]
  
  df <- df[df$seqGap + df$missingGap < 100001,]
  df <- df[df$seqnames != "chrX",]
  df <- df[complete.cases(df),]
  
  # havn't completly normalised yet
  
  
  
  synthBinNorm.gr <- GRanges(df)
  
  seqinfo(synthBinNorm.gr) <- seqinfo(synthBin.gr)
  
  corDF <- data.frame(mcols(synthBinNorm.gr))
  
  
  corMap <- data.frame(cor(as.matrix(corDF), method = "spearman"))
  corMap <- corMap[,c("refDel","refIns","queDel", "queIns")]
  corMap <- corMap[!rownames(corMap) == "fill",]
  corMap <- corMap[!rownames(corMap) == "missingGap",]
  corMap <- corMap[!rownames(corMap) %in% c("refDel","refIns","queDel","queIns"),]
  
  corMap<- corMap[hclust(dist(corMap))$order,]
  
  corTestDF <- corMap
  corTestDF[,colnames(corMap)] <- NA
  for(i in rownames(corMap)){
    for(j in colnames(corMap)){
      corRes <- cor.test(corDF[,i], corDF[,j], method = "spearman")
      corTestDF[i,j] <- corRes$p.value
    }
    corTestDF[i,] <- p.adjust(corTestDF[i,], method = "fdr", n = length(corMat))
  }
  
  corTestMat <- as.matrix(corTestDF)
  corMat <- as.matrix(corMap)
  
  corMat[corTestMat>.05] <- NA
  
  assign(corMat,x = paste(spec, ".corMat", sep = ""))
  assign(corMap,x = paste(spec, ".corMap", sep = ""))
  
  
}

mm10.corMap <- mm10.corMap[rownames(hg19.corMap),]
mm10.corMat <- mm10.corMat[rownames(hg19.corMap),]


layout(1)
corMapOrder <- hclust(dist(cbind(hg19.corMap, mm10.corMap)))$order

hg19.corMat <- hg19.corMat[corMapOrder,]
mm10.corMat <- mm10.corMat[corMapOrder,]



layout(matrix(1:3,nrow = 1), widths = c(4,4,2))

par(mar = c(10,10,5,1))
# hg19
image(t(as.matrix(hg19.corMap)),
      col="grey", xaxt = "n", yaxt = "n", main = "hg19")
colFun = colorRampPalette(c("blue","lightblue", "white", "orange", "red"))
image(t(hg19.corMat),
      col=colFun(40), add = TRUE)
mtext(c("hg19 loss", "hg19 gain", "mm10 loss", "mm10 gqin"), at = seq(0,1,length.out = ncol(hg19.corMat)), las = 2, line = 1, side = 1)
mtext((rownames(hg19.corMat)), at = seq(0,1,length.out = nrow(hg19.corMat)), las = 2, line = 1, side = 2)

# mm10
par(mar = c(10,5,5,6))

mm10.corMat <- mm10.corMat[,c("queDel","queIns", "refDel", "refIns")]
image(t(as.matrix(mm10.corMap)),
      col="grey", xaxt = "n", yaxt = "n", main = "mm10")
colFun = colorRampPalette(c("blue","lightblue", "white", "orange", "red"))
image(t(mm10.corMat),
      col=colFun(40), add = TRUE)
mtext(c("hg19 loss", "hg19 gain", "mm10 loss", "mm10 gain"), at = seq(0,1,length.out = ncol(mm10.corMat)), las = 2, line = 1, side = 1)
#mtext((rownames(mm10.corMat)), at = seq(0,1,length.out = nrow(mm10.corMat)), las = 2, line = 1, side = 2)



par(mar = c(30,5,5,1))
image(t(matrix(seq(-1,1,length.out = 7), ncol = 1)), col = colFun(40),
      xaxt = "n", yaxt = "n")
axis(side = 2, at = seq(0,1,length.out = 5), labels = seq(-1,1,length.out = 5), las = 2)
# should include both 

# characterizing the DNA gain/loss landscape of the human and mouse genomes

# what factors are they correlated with 


#synthBinNorm.gr

# convert to 1MB bins and get the means


# so we are able to explain a certain amount of the variance

# maybe the method isn't working because we are misclassifing our important regions


# Next we could look at each enriched region?


# time to start writing






