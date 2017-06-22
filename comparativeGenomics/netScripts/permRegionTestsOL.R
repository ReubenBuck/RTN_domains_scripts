#### FInal idea is to look at significant regions and see if they share a significant overlap

### Both with high and low spots, 

# that reminds me, I don't think i split my data into high and low


# by randomizing their locations that way, it is probably likly that we will get a normal distribution
# the probability that we see an overlap as significant as that

# then we can alsp the local Z test. 


# is it possible to permute all the regions at once and look at all the potential overlaps?

# use regioner
rm(list = ls())



ZtoP <- function(Zs){
  2*pnorm(-abs(Zs))
}

ZsigAdj <- function(Zs, nb, cutoff){
  score = p.adjust(p = ZtoP(Zs), method = "BH")
  return(score <= cutoff)
}



library(igraph)
library(spdep)
library(regioneR)

load(file = "~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/mm10.synthBin.genome.variable.RData")


df<- data.frame(synthBin.gr)


colChoice <- c("dnasePeaks","ladCov","exon", "intron","ctcfMotif", 
               "L1Motif", "prdm9Motif", "ancient", "new_L1", "new_SINE", 
               "old_L1", "recHotSpot","cpgCov")

remainRefBases <- width(synthBin.gr)[1] - (df$queIns + df$refDel + df$seqGap)
df[remainRefBases > 0,colChoice] <- (df[remainRefBases > 0,colChoice]/remainRefBases[remainRefBases > 0]) * width(synthBin.gr)[1]

df <- df[df$seqGap + df$missingGap < 100001,]
df <- df[df$seqnames != "chrX",]
df <- df[complete.cases(df),]


synthBinNorm.gr <- GRanges(df)

seqinfo(synthBinNorm.gr) <- seqinfo(synthBin.gr)




ol<-findOverlaps(synthBinNorm.gr, maxgap = 3*width(synthBinNorm.gr)[1])
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

gapChoice <- c("refIns", "refDel", "queIns", "queDel",
               "gcContent", "dnasePeaks", "dnaseActivity",
               "exon", "intron", "ladCov", "ctcfMotif",
               "L1Motif", "prdm9Motif", "ancient", "new_L1",
               "new_SINE", "old_L1", "recombRate", "recHotSpot",
               "cpgCov")

dfGscore <- data.frame(synthBinNorm.gr)

for(i in 1:length(gapChoice)){
  score <- df[,gapChoice[i]]
  G <- localG(x = score,wMat)
  dfGscore[,gapChoice[i]] <- G
}


sigsRefDel <- ZsigAdj(dfGscore$refDel, nb = wMat$neighbours, cutoff = .05)
synthBinSigRefDel.gr <- synthBinNorm.gr[sigsRefDel & dfGscore$refDel >0]
sigRegion.gr <- reduce(synthBinSigRefDel.gr)

synthGenome.gr <- GRanges(seqinfo(synthBinNorm.gr))
mask.gr <- setdiff(synthGenome.gr, synthBinNorm.gr)

# pick random regions of same size as X
n = 10000
randRegion.grl <- replicate(n,expr = randomizeRegions(A = sigRegion.gr, 
                                                      genome = synthGenome.gr,
                                                      mask = mask.gr,
                                                      allow.overlaps = FALSE))
randRegionLen <- unlist(lapply(randRegion.grl, length))

n = sum(randRegionLen == length(sigRegion.gr))
randRegion.grl <- randRegion.grl[randRegionLen == length(sigRegion.gr)]

# should we apply a shift to get stuff in phase
randRegionNew.grl <- GRangesList(randRegion.grl)

start(randRegionNew.grl) <- IntegerList(round(start(GRangesList(randRegion.grl))/200000) * 200000 + 1)
end(randRegionNew.grl) <- IntegerList(round(end(GRangesList(randRegion.grl))/200000) * 200000)

# this means our regions are in phase and can overlap with any other sig region


ol <- findOverlaps(randRegionNew.grl, synthBinNorm.gr)



sigMat <- matrix(NA, nrow = nrow(dfGscore), ncol=length(gapChoice))
for(i in 1:length(gapChoice)){
  sigElement <- 
  sigMat[,i] <- ZsigAdj(dfGscore[,gapChoice[i]], cutoff = .05) & dfGscore[,gapChoice[i]] > 0
}



hitList <- split(x = data.frame(sigMat[subjectHits(ol),]), f = queryHits(ol))
hitTotals <- data.frame(t(as.data.frame(lapply(X = hitList, colSums))))
hitRates <- hitTotals/table(queryHits(ol))
colnames(hitRates) <- gapChoice

olSig <- overlapsAny(synthBinNorm.gr,synthBinSigRefDel.gr)
sigRates <- colSums(sigMat[olSig,])/length(synthBinSigRefDel.gr)
names(sigRates) <- gapChoice

pdf(file = "~/Desktop/mm10.repDel.sig.pdf", onefile = TRUE)
for(i in 1:length(gapChoice)){
varNo <- gapChoice[i]

sdDist <- 4 * sd(hitRates[,varNo]) 
lineDist <- abs(mean(hitRates[,varNo]) - sigRates[varNo])
if(sdDist > lineDist){
  xWidth <- c(mean(hitRates[,varNo]) - sdDist,
              mean(hitRates[,varNo]) + sdDist)
}else{
  xWidth <- c(mean(hitRates[,varNo]) - lineDist,
              mean(hitRates[,varNo]) + lineDist)
}


hist(hitRates[,varNo], breaks = 30, xlim = xWidth, main = varNo)
abline(v = sigRates[varNo], col = "purple", lwd = 3)
abline(v = mean(hitRates[,varNo]), col = 2, lwd = 2)
}

dev.off()

# count how many hits we got?

# build a matrix of all our significant bins








