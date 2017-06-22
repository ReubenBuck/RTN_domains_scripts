


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

load(file = "~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/hg19.synthBin.genome.variable.RData")


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

gapChoice <- c("refIns", "refDel", "queIns", "queDel")

dfGscore <- data.frame(synthBinNorm.gr)

for(i in 1:4){
  score <- df[,gapChoice[i]]
  G <- localG(x = score,wMat)
  dfGscore[,gapChoice[i]] <- G
}



# if we could work out our z.score cut off that would be good


sigs <- ZsigAdj(dfGscore$refDel, nb = wMat$neighbours, cutoff = .05)

synthBinSig.gr <- synthBinNorm.gr[sigs]
sigRegions.gr <- reduce(synthBinSig.gr)
synthGenome.gr <- GRanges(seqinfo(synthBinNorm.gr))
mask.gr <- setdiff(synthGenome.gr, synthBinNorm.gr)

# pick random regions of same size as X
n = 1000
randRegion.grl <- replicate(n,expr = randomizeRegions(A = sigRegions.gr, 
                                                      genome = synthGenome.gr,
                                                      mask = mask.gr,
                                                      allow.overlaps = FALSE))
# get regions lengths to remove small ones, part of pkg bug
randLengths <- unlist(lapply(randRegion.grl, length))

# overlap to get contigous regions of a similar size
qHitL <- findOverlaps(GRangesList(randRegion.grl), synthBinNorm.gr, minoverlap = 97000)

hist(table(queryHits(qHitL)), xlim = c(2000,3000), breaks = 1000)
abline(v = length(synthBinSig.gr))

qHitRanges <- data.frame(mcols(synthBinNorm.gr[subjectHits(qHitL)]))
qHitRanges <- split(qHitRanges, queryHits(qHitL))

qNonHitMeans <- NULL
for(i in 1:n){
  hits <- subjectHits(qHitL)[queryHits(qHitL) == i]
  qNonHitMeans <- rbind(qNonHitMeans, colMeans(data.frame(mcols(synthBinNorm.gr[-hits]))))
}




# get the means of each of our vairaibles
qHitMeans <- data.frame(t(as.data.frame(lapply(FUN = colMeans, qHitRanges, na.rm = T))))
rownames(qHitMeans ) <- NULL


qHitMeans  <- qHitMeans[randLengths == length(sigRegions.gr),]
qNonHitMeans  <- qNonHitMeans[randLengths == length(sigRegions.gr),]




# this histogram seems to be showing us that there may be some interesting things going on with the spacial aspects of the data



# mean of our sigRegions

sigMeans <- colMeans(data.frame(mcols(synthBinSig.gr)))

nonSigMeans <- colMeans(data.frame(mcols(synthBinNorm.gr[!sigs])))

genoMeans <- colMeans(data.frame(mcols(synthBinNorm.gr)))



varChoice <- "cpgCov"
for(i in names(sigMeans)){
  varChoice <- i
  xWidth <- max(abs(sigMeans - nonSigMeans)[varChoice], 6*sd((qHitMeans - qNonHitMeans)[,varChoice]) )
  hist((qHitMeans - qNonHitMeans)[,varChoice], 
       breaks = 40, 
       main = varChoice, 
       xlim = c(-xWidth, xWidth),
       xlab = "mean differences")
  abline(v = (sigMeans - nonSigMeans)[varChoice], col ="purple", lwd = 3)
}
# so there is some sort of inherant bias in my data


# the mean differences are not constant
varChoice <- "prdm9Motif"
hist(qHitMeans[,varChoice] - genoMeans[varChoice], breaks = 40, main = varChoice, add = TRUE)



# the other approach is to identify significant regions for everything 



sampDiffsMeans <-data.frame(t( t(rDf) - genoMeans))

sigDiffsMeans <- sigMeans - genoMeans



hist(qHitMeans$gcContent, breaks = 15)
m <- mean(qHitMeans$gcContent)
abline(v = m)
abline(v = genoMeans["gcContent"])
abline(v = sigMeans["gcContent"])


# we won't be able to see the contiribution to vairaition but will see if it is enriched.






hist(a[,9], breaks = 100)
abline(v = aMeans[9])
abline(v = sigMeans[9], col=2)

varNo = "gcContent"
dAll <- t(colMeans(data.frame(mcols(synthBinNorm.gr))) - t(rDf))[,varNo]
hist(dAll, breaks = 50)

d1 <- (sigMeans - colMeans(data.frame(mcols(synthBinNorm.gr))))[varNo]
abline(v = d1)

sum(dAll > d1)/100
sum(dAll < d1)/100

hist(dAll, breaks = 100)

hist()

sigs

varChoice <- "ladCov"
gcLevel <- replicate(1000, 
                     mean(df[,varChoice]) -  mean(sample(x = df[,varChoice], 
                                                       size = length(synthBinSig.gr),
                                                       replace = FALSE) 
                                                )
                     )

gcDiff <- mean(df[,varChoice]) - mean(mcols(synthBinSig.gr)[[varChoice]], na.rm = T)


hist(gcLevel,breaks = 100)
abline(v = gcDiff)





# most things seem significant

sigsPerm <- ZsigAdj(dfGscore$refIns - dfGscore$refDel, nb = wMat$neighbours, cutoff = .05)



y = df$ladCov


diff(by(y, sigs, mean))


dist0 <- replicate(1000, diff(by(y, 
                                 sample(sigsPerm, length(sigsPerm), replace = F), 
                                 mean)
)
)

# dist1 <- replicate(100, diff(by(y, 
#                                  sample(c(TRUE,FALSE), length(sigs), replace = T, prob = probS), 
#                                  mean)
#                               )
#                    )

hist(dist0, xlim=c(-12*sd(dist0), 12*sd(dist0)))
#hist(dist1, xlim=c(-6*sd(dist0), 6*sd(dist0)), add = T)

abline(v = diff(by(y, sigsPerm, mean)))

sum( abs(dist0) > abs(diff(by(y, sigs, mean))) )/1000



















# Check if chaning the sample size will give a skewed result

probS <- c(sum(sigs)/length(sigs) , 1 - (sum(sigs)/length(sigs)))

a <- sample(c(TRUE,FALSE), length(sigs), replace = T, prob = probS)










# we test wheater our sample is different to the genome,
# our sampling size is also likely to change too







col1 <- scales::alpha("black", .2)
col2 <- scales::alpha("red", .01)


gapChoice = "dnaseActivity"

breaks = range(mcols(synthBinNorm.gr)[[gapChoice]])
breaks = seq(from = floor(breaks[1]), to = ceiling(breaks[2]), length.out = 100)

hist(mcols(synthBinSig.gr)[[gapChoice]], col = col1, add = F, border = col1, freq = FALSE, breaks = breaks)

for(i in 1:length(randRegion.gr)){
  hist(mcols(synthBinNorm.gr[subjectHits(qHitL)[queryHits(qHitL) == i]])[[gapChoice]], 
       freq = FALSE,
       col = col2, border = col2, add= TRUE, breaks = breaks)
}

hist(mcols(synthBinSig.gr)[[gapChoice]], col = col1, border = col1, add = TRUE, freq = F, breaks = breaks)





length(sigSynthBin.gr)
length(qHits)
# 









samp.gr <- sample(synthBin.gr, size = 10)

randomizeRegions(A = samp.gr, genome = synthGenome.gr,allow.overlaps = FALSE)

region.gr <- (synthBin.gr[100:110])



# I can use the resample function to just pick regions randomly a bunch of times,
# I could have doen this anyway with the sample function and a dataframe

# alternativly, I could randomize my regions, that way I'll keep the genome structure intact
# since regions sapn a certain size, I will be able to maintain those sizes,
# I just hvae to get the setdiff from my genome to mask out the reigiosn that are out of bounds.


hist(table(queryHits(qHitL)), breaks = 500)

abline(v = length(synthBinSig.gr))


# why are we getting small ranges? 

hist(width(randRegion.gr[[31]]))

hist(width(sigRegions.gr))





reps <- replicate(colMeans(data.frame(mcols((resampleRegions(region.gr,universe = synthBin.gr))))),n = 1000)

reps <- as.data.frame(t(reps))

hist(reps$dnasePeaks)


