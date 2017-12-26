# Test robustness of approaches 



library(igraph)
library(spdep)
library(GenomicRanges)

rm(list = ls())




ZtoP <- function(Zs){
  2*pnorm(-abs(Zs))
}

ZsigAdj <- function(Zs, cutoff){
  score = p.adjust(p = ZtoP(Zs), method = "BH")
  return(score <= cutoff)
}




specRef = "mm10"
specQue = "hg19"

load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/",specRef,".synthBin.RData", sep = ""))
noRepSynthBin.gr <- synthBin.gr
noRepSynthBin.df <- mcols(noRepSynthBin.gr)
noRepSynthBin.df <- noRepSynthBin.df[noRepSynthBin.df$missingGap + noRepSynthBin.df$seqGap < 20000,]

load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/",specRef,".synthBin.rep.RData", sep = ""))
repSynthBin.gr <- synthBin.gr
repSynthBin.df <- mcols(repSynthBin.gr)
repSynthBin.df <- repSynthBin.df[repSynthBin.df$missingGap + repSynthBin.df$seqGap < 20000,]







gapChoice <- c("refIns", "refDel", "queIns", "queDel")
mehtods = c("rep", "noRep")
nDistList <- NULL
sigRangeList <- NULL

#for(m in mehtods){

m = "rep"  

df <- data.frame(get(paste(m,"SynthBin.gr", sep = "")))
df$remainingBases <- (df$end - df$start + 1) - (df$missingGap + df$seqGap)

df <- df[df$remainingBases > 150001,]
df <- df[complete.cases(df),]

df[,c("refIns", "refDel", "queIns", "queDel")] <- (df[,c("refIns", "refDel", "queIns", "queDel")] * 200000 )/ df$remainingBases 

synthBinNorm.gr <- GRanges(df)

seqinfo(synthBinNorm.gr) <- seqinfo(synthBin.gr)



for(j in c(1:10, seq(12,20,2))){
  ol<-findOverlaps(synthBinNorm.gr, maxgap = j*width(synthBinNorm.gr)[1])
  ol <- ol[!(isRedundantHit(ol))]
  # remove NA hits
  ol <- ol[!is.na(df$refIns[queryHits(ol)])]
  ol <- ol[!is.na(df$refIns[subjectHits(ol)])]
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
  
  gapChoice <- c("refIns", "refDel", "queIns", "queDel")
  dfGscore <- data.frame(synthBinNorm.gr)
  
  for(i in 1:length(gapChoice)){
    score <- df[,gapChoice[i]]
    G <- localG(x = score,wMat)
    dfGscore[,gapChoice[i]] <- G
  }
  nDistList <- c(nDistList, list(dfGscore))
  
  
  sigRanges = NULL
  for(i in 1:length(gapChoice)){
    z <- ZsigAdj(dfGscore[,gapChoice[i]],cutoff = .05) & dfGscore[,gapChoice[i]] > 0
    sigRanges <- c(sigRanges , list((synthBinNorm.gr[z])))
  }
  names(sigRanges) <- gapChoice
  
  sigRangeList <- c(sigRangeList, list(sigRanges))
  
}

#assign(x = paste(m, "SigRanges", sep = ""), value = sigRanges)

#}


# we want to see how smooth we can make it before it gets out of hand

COR <- SMV <- NULL
for(i in 1:15){
  COR0 <- SMV0 <- NULL
  for(j in c("refIns", "refDel", "queIns", "queDel")){
    COR0 <- c(COR0,
             cor(df[,j], nDistList[[i]][[j]]))
    SMV0 <- c(SMV0,
             sd(diff(nDistList[[i]][[j]])))
  }
  names(SMV0) <- names(COR0) <- c("refIns", "refDel", "queIns", "queDel")
  COR <- c(COR, list(COR0))
  SMV <- c(SMV, list(SMV0))
}



pdf(file = paste("/home/buckley/Documents/RTN_domain_writing/manuscriptRound2/bioArxiv/sup/TexFigs/hotspotTest/",
                 specRef, "nSizeSmooth.pdf", sep = ""), width = 6, height = 10)
par(oma = c(0,1,0,0), mar = c(5,4,5,2))
layout(matrix(1:2, nrow = 2))
mainNames <- paste(specRef,c("gain", "loss"))
names(mainNames) <- c("refIns", "refDel")
for(j in c("refIns", "refDel")){
plot(c(sd(diff(scale(df[,j]))), sapply(SMV, "[[", j)), 
     c(1, sapply(COR, "[[", j)^2), 
     pch = 16, cex = .8, type = "b", xlab = "roughness", 
     ylab = parse(text = "R^2"), main = mainNames[j],
     ylim = c(0,1))
text(c(sd(diff(scale(df[,j]))), sapply(SMV, "[[", j)), 
     c(1, sapply(COR, "[[", j)^2), 
     labels = c(0:length(COR)), pos = 1, col = 2)
}

dev.off()

pdf(file = paste("/home/buckley/Documents/RTN_domain_writing/manuscriptRound2/bioArxiv/sup/TexFigs/hotspotTest/",
                 specRef, "nSizeWiggle.pdf", sep = ""), width = 6,height = 10)
layout(1:10)
par(oma = c(0,1,3,0), mar = c(1,4,1,2))
chrChoice <- "chr12"
xlim = c(80e6, 170e6)
for(j in c("refIns", "refDel")){
  plot((df$start + (df$width/2))[df$seqnames == chrChoice],
       scale(df[[j]]/df$remainingBases)[df$seqnames == chrChoice], 
       type = "l", main = mainNames[j],
       ylab = "Z", xlim = xlim)
  grid()
  abline(h = 0, lty = 2)
  
  for(i in c(1,3,5)){
    plot((df$start + (df$width/2))[df$seqnames == chrChoice],
         nDistList[[1]][[j]][df$seqnames == chrChoice], type = "l",
         ylab = "G*i", sub = chrChoice, xlim = xlim)
    grid()
    abline(h = 0, lty = 2)
    rect(xleft = start(sigRangeList[[i]][[j]])[as.character(seqnames(sigRangeList[[i]][[j]])) == chrChoice],
         xright = end(sigRangeList[[i]][[j]])[as.character(seqnames(sigRangeList[[i]][[j]])) == chrChoice],
         ytop = 10, ybottom = -10,
         col = scales::alpha(2,.2), border = NA)
  }
  plot.new()
  text(.5,.5,labels = chrChoice, cex = 1.2)
}
dev.off()



