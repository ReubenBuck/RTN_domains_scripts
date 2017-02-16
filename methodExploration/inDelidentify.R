
#### test for indel locations 

# control for different chains. 
# If we get a gap somewhere it has to be part of the same chain.
# interupted chain is evidence of translocation or inversion.



### using dog chr22 as the reference

rm(list = ls())


library(GenomicRanges)

setwd("~/Desktop/RTN_domains/")


hg <- read.table(file = "data/chainAlignments/inDelTest/canFam2hg19chr22.txt", 
                 col.names = c("refChr", "refLen", "refStart", "refEnd", "refStrand", 
                               "refGap", "queChr", "queLen", "queStart", "queEnd", "queStrand", 
                               "queGap", "chainID"))

mm <- read.table(file = "data/chainAlignments/inDelTest/canFam2mm9chr22.txt",
                 col.names = c("refChr", "refLen", "refStart", "refEnd", "refStrand", 
                               "refGap", "queChr", "queLen", "queStart", "queEnd", "queStrand", 
                               "queGap", "chainID"))


hgRaw.gr <- GRanges(seqnames = Rle(hg$refChr), 
                    ranges = IRanges(start = hg$refStart, end = hg$refEnd))

mmRaw.gr <- GRanges(seqnames = Rle(mm$refChr), 
                    ranges = IRanges(start = mm$refStart, end = mm$refEnd))


hgWindow.gr <- GRanges(seqnames = Rle(hg$refChr[1]),
                       ranges = IRanges(start = hg$refStart[1], end = hg$refEnd[nrow(hg)])
                       )

mmWindow.gr <- GRanges(seqnames = Rle(mm$refChr[1]),
                       ranges = IRanges(start = mm$refStart[1], end = mm$refEnd[nrow(mm)])
)


intRange <- intersect(hgRaw.gr, mmRaw.gr)


# setting the start might be hard with multiple chromosomes
intRangeGaps <- gaps(intRange,start = end(intRange)[1])
start(intRangeGaps) = start(intRangeGaps) -1
end(intRangeGaps) = end(intRangeGaps) + 1




##### 
# Human
####
hgGapAll <- data.frame(refChr = hg$refChr[1:(nrow(hg)-1)], refStart = hg$refEnd[1:(nrow(hg)-1)],
                       refEnd = hg$refStart[2:(nrow(hg))], queChr = hg$queChr[1:(nrow(hg)-1)],
                       queStart = hg$queEnd[1:(nrow(hg)-1)], queEnd = hg$queStart[2:(nrow(hg))]
                       )
hgGapAll.gr <- GRanges(seqnames = Rle(hgGapAll$refChr), 
                       ranges = IRanges(start = hgGapAll$refStart, end = hgGapAll$refEnd),
                       queSeqnames = Rle(hgGapAll$queChr), 
                       queRanges = IRanges(start = hgGapAll$queStart, end = hgGapAll$queEnd)
)


hgQueGapElements <- (1:(nrow(hg)-1))[hg$refStart[2:nrow(hg)] == hg$refEnd[1:(nrow(hg)-1)]]

hgQueGap <- data.frame(refChr = hg$refChr[hgQueGapElements], refStart = hg$refEnd[hgQueGapElements],
                     refEnd = hg$refStart[hgQueGapElements + 1], queChr = hg$queChr[hgQueGapElements],
                     queStart = hg$queEnd[hgQueGapElements], queEnd = hg$queStart[hgQueGapElements + 1],
                     gap = "Q")

hist(log10(hgQueGap$queEnd - hgQueGap$queStart), breaks = 1000)
hist(hgQueGap$queEnd - hgQueGap$queStart, xlim = c(0,300), breaks = seq(0,200000,1))



hgRefGapElements <- (1:(nrow(hg)-1))[hg$queStart[2:nrow(hg)] == hg$queEnd[1:(nrow(hg)-1)]]

hgRefGap <- data.frame(refChr = hg$refChr[hgRefGapElements], refStart = hg$refEnd[hgRefGapElements],
                       refEnd = hg$refStart[hgRefGapElements + 1], queChr = hg$queChr[hgRefGapElements],
                       queStart = hg$queEnd[hgRefGapElements], queEnd = hg$queStart[hgRefGapElements + 1],
                       gap = "R")

hist(log10(hgRefGap$refEnd - hgRefGap$refStart), breaks = 1000)
hist(hgRefGap$refEnd - hgRefGap$refStart, xlim = c(0,700), breaks = seq(0,200000,1), ylim = c(0,200))


hgInDel <- rbind(hgRefGap,hgQueGap)


#####

#######
## mouse
####

mmGapAll <- data.frame(refChr = mm$refChr[1:(nrow(mm)-1)], refStart = mm$refEnd[1:(nrow(mm)-1)],
                       refEnd = mm$refStart[2:(nrow(mm))], queChr = mm$queChr[1:(nrow(mm)-1)],
                       queStart = mm$queEnd[1:(nrow(mm)-1)], queEnd = mm$queStart[2:(nrow(mm))]
)
mmGapAll.gr <- GRanges(seqnames = Rle(mmGapAll$refChr), 
                       ranges = IRanges(start = mmGapAll$refStart, end = mmGapAll$refEnd),
                       queSeqnames = Rle(mmGapAll$queChr), 
                       queRanges = IRanges(start = mmGapAll$queStart, end = mmGapAll$queEnd)
)




mmQueGapElements <- (1:(nrow(mm)-1))[mm$refStart[2:nrow(mm)] == mm$refEnd[1:(nrow(mm)-1)]]

mmQueGap <- data.frame(refChr = mm$refChr[mmQueGapElements], refStart = mm$refEnd[mmQueGapElements],
                       refEnd = mm$refStart[mmQueGapElements + 1], queChr = mm$queChr[mmQueGapElements],
                       queStart = mm$queEnd[mmQueGapElements], queEnd = mm$queStart[mmQueGapElements + 1],
                       gap = "Q")

hist(log10(mmQueGap$queEnd - mmQueGap$queStart), breaks = 1000)
hist(mmQueGap$queEnd - mmQueGap$queStart, xlim = c(0,700), breaks = seq(0,200000,1))

# how do we know these gaps are real
# there is definatly a concentration of small ones following a power curve. 


mmRefGapElements <- (1:(nrow(mm)-1))[mm$queStart[2:nrow(mm)] == mm$queEnd[1:(nrow(mm)-1)]]

mmRefGap <- data.frame(refChr = mm$refChr[mmRefGapElements], refStart = mm$refEnd[mmRefGapElements],
                       refEnd = mm$refStart[mmRefGapElements + 1], queChr = mm$queChr[mmRefGapElements],
                       queStart = mm$queEnd[mmRefGapElements], queEnd = mm$queStart[mmRefGapElements + 1],
                       gap = "R")

hist(log10(mmRefGap$refEnd - mmRefGap$refStart), breaks = 1000)
hist(mmRefGap$refEnd - mmRefGap$refStart, xlim = c(0,700), breaks = seq(0,200000,1), ylim = c(0,200))


mmInDel <- rbind(mmRefGap,mmQueGap)


### this is when we look to see if we can find stuff that is unique to a particular genome 
### in other words we can get the set diff and annotate using parsiomny the indels. 

hg.gr <-GRanges(seqnames = Rle(hgInDel$refChr), 
                ranges = IRanges(start = hgInDel$refStart, end = hgInDel$refEnd),
                queSeqnames = Rle(hgInDel$queChr), 
                queRanges = IRanges(start = hgInDel$queStart, end = hgInDel$queEnd),
                gap = hgInDel$gap)

mm.gr <-GRanges(seqnames = Rle(mmInDel$refChr),
                ranges = IRanges(start = mmInDel$refStart, end = mmInDel$refEnd),
                queSeqnames = Rle(mmInDel$queChr), 
                queRanges = IRanges(start = mmInDel$queStart, end = mmInDel$queEnd),
                gap = mmInDel$gap)


hg.gr <- subsetByOverlaps(hg.gr, 
                          GRanges(seqnames = seqnames(mm.gr)[1], 
                                  ranges = IRanges(start = start(mm.gr[1]), end = end(mm.gr[length(mm.gr)]))))

mm.gr <- subsetByOverlaps(mm.gr,
                          GRanges(seqnames = seqnames(hg.gr)[1], 
                                  ranges = IRanges(start = start(hg.gr[1]), end = end(hg.gr[length(hg.gr)]))))



hgRef.gr <- hg.gr[mcols(hg.gr)$gap == "R"]
olrefHG <- as.matrix(findOverlaps(hgRef.gr, intRangeGaps, type = "equal"))
hgRef.gr <- hgRef.gr[olrefHG[,1]]


mmRef.gr <- mm.gr[mcols(mm.gr)$gap == "R"]
olrefMM <- as.matrix(findOverlaps(mmRef.gr, intRangeGaps,  type = "equal"))
mmRef.gr <- mmRef.gr[olrefMM[,1]]


hgRefC <- subsetByOverlaps(hgRef.gr, mmRef.gr, type = "equal")
hgRefU <- hgRef.gr[!(hgRef.gr %in% hgRefC)]

mmRefC <- subsetByOverlaps(mmRef.gr, hgRef.gr, type = "equal")
mmRefU <- mmRef.gr[!(mmRef.gr %in% mmRefC)]

mmRefU <- subsetByOverlaps(mmRefU, reduce(gaps(hgGapAll.gr)), type = "within")
hgRefU <- subsetByOverlaps(hgRefU, reduce(gaps(mmGapAll.gr)), type = "within")



##### so got the ref regions that are comparable 

### now to get the query gaps 

# as long as they appear within the itersect it should be ok
# we cna remove the overlapping ones and then we have a clean set


hgQue.gr <- hg.gr[mcols(hg.gr)$gap == "Q"]
hgQue.gr <- subsetByOverlaps(hgQue.gr, intRange)


mmQue.gr <- mm.gr[mcols(mm.gr)$gap == "Q"]
mmQue.gr <- subsetByOverlaps(mmQue.gr, intRange)


hist(log10(width(mcols(mmQue.gr)$queRanges)), breaks = 100)
hist(log10(width(mcols(hgQue.gr)$queRanges)), breaks = 100, add = T, col = 2, density = 0)



hgQueC <- subsetByOverlaps(hgQue.gr, mmQue.gr)
hgQueU <- hgQue.gr[!(hgQue.gr %in% hgQueC)]

mmQueC <- subsetByOverlaps(mmQue.gr, hgQue.gr)
mmQueU <- mmQue.gr[!(mmQue.gr %in% mmQueC)]

hist(log10(width(mcols(mmQueC)$queRanges)), breaks = 100)
hist(log10(width(mcols(hgQueC)$queRanges)), breaks = 100)

qWd <- width(mcols(mmQueC)$queRanges) - width(mcols(hgQueC)$queRanges)


pdf(file = "plots/inDelIdentify/inDelSizeDist.pdf", onefile = TRUE)


hist(log10(width(mcols(mmQueC)$queRanges)), breaks = 50, main = "Overlapping gaps in que genomes", xlab = "gap size (log10 bp)")
hist(log10(width(mcols(hgQueC)$queRanges)), breaks = 50, col =2 , density = 0, add = TRUE)
legend("topright", legend = c("mm9", "hg19"), fill = c(1,2), title = "Que geneomes")


hist(log10(width(mmRefU)), breaks = 50, main = "Unique gaps in ref genome (Deletion)", xlab = "gap size (log10 bp)")
hist(log10(width(hgRefU)), breaks = 50, col =2 , density = 0, add = TRUE)
legend("topright", legend = c("mm9", "hg19"), fill = c(1,2), title = "Que geneomes")


hist(log10(width(mmRefC)), breaks = 50, main = "Overlapping gaps in ref genome", xlab = "gap size (log10 bp)")
hist(log10(width(hgRefC)), breaks = 50, col =2 , density = 0, add = TRUE)
legend("topright", legend = c("mm9", "hg19"), fill = c(1,2), title = "Que geneomes")

hist(log10(width(mcols(mmQueU)$queRanges)), breaks = 50, main = "Unique gaps in que genomes (insertions)", xlab = "gap size (log10 bp)")
hist(log10(width(mcols(hgQueU)$queRanges)), breaks = 50, col =2 , density = 0, add = TRUE)
legend("topright", legend = c("mm9", "hg19"), fill = c(1,2), title = "Que geneomes")


dev.off()

# so now I have the unique coordinates for each species. 
# What do they mean in regards to insertion and deletion. 

# these are primerily lineage specific insertions. 


sample(mmQueU[width(mcols(mmQueU)$queRanges) < 300 & width(mcols(mmQueU)$queRanges) > 250])

mmQueU[width(mcols(mmQueU)$queRanges) < 100][5]

hgQueU


# insertions aren't as clear cut as our deletions. 




mmQueC[qWd > 100]
hgQueC[qWd > 100]




# after finding these regions we need to go back and verify some things. 

test <- mmQueU[width(mcols(mmQueU)$queRanges) < 100][3]

subsetByOverlaps(hgGapAll.gr, test, maxgap = 300)


width(mcols(subsetByOverlaps(mmGapAll.gr, test, maxgap = 500))$queRanges)






MMwDiff <- data.frame(width(mmGapAll.gr), width(mcols(mmGapAll.gr)$queRanges))
MMwDiffmin <- apply(X = MMwDiff, MARGIN = 1, FUN = min)
MMwDiffmax <- apply(X = MMwDiff, MARGIN = 1, FUN = max)
# from our analysis alone it looks like the majority of events occur at our cutoff




HGwDiff <- data.frame(width(hgGapAll.gr), width(mcols(hgGapAll.gr)$queRanges))
HGwDiffmin <- apply(X = HGwDiff, MARGIN = 1, FUN = min)
HGwDiffmax <- apply(X = HGwDiff, MARGIN = 1, FUN = max)


pdf(file = "plots/inDelIdentify/gapSizes.pdf", onefile = T)
smoothScatter(log10(HGwDiffmin), log10(HGwDiffmax), xlim = c(0,6), ylim = c(0,6), nrpoints = 0, main = "hg19", xlab = "min gap (log10 bp)", ylab = "max gap (log10 bp)")
abline(a=0,b=1)
smoothScatter(log10(MMwDiffmin), log10(MMwDiffmax), xlim = c(0,6), ylim = c(0,6), nrpoints = 0, main = "mm9", xlab = "min gap (log10 bp)", ylab = "max gap (log10 bp)")
abline(a=0,b=1)


smoothScatter(log10(HGwDiffmin), log10(HGwDiffmax), xlim = c(1,6), ylim = c(0,6), nrpoints = 0, main = "hg19", xlab = "min gap (log10 bp)", ylab = "max gap (log10 bp)")
abline(a=0,b=1)
smoothScatter(log10(MMwDiffmin), log10(MMwDiffmax), xlim = c(1,6), ylim = c(0,6), nrpoints = 0, main = "mm9", xlab = "min gap (log10 bp)", ylab = "max gap (log10 bp)")
abline(a=0,b=1)

dev.off()


n = 10
hist(log10(HGwDiffmin[HGwDiffmax > n]), breaks = 100)
abline(v = log10(n))


hist(log10(HGwDiffmax[HGwDiffmin==1]), breaks = 100)
n = 1
hist(log10(HGwDiffmax[HGwDiffmin > n]), breaks = 100, add = T, col = 2, density = 0)
#abline(v = log10(n))


hist(log10(MMwDiffmax[MMwDiffmin==1]), breaks = 200, freq = F)
n = 1
m = 20
hist(log10(MMwDiffmax[MMwDiffmin > n & MMwDiffmin < m]), breaks = 200, add = T, col = 2, density = 0, freq = F)
abline(v = log10(m))


hist((MMwDiffmax[MMwDiffmin==1]), breaks = 10000, freq = T, xlim = c(0,400), ylim = c(0,200))
n = 1
m = 20
hist((MMwDiffmax[MMwDiffmin > n & MMwDiffmin < m]), breaks = 10000, add = T, col = 2, density = 0, freq = T)
abline(v = (m))



hist(log10(MMwDiffmax/MMwDiffmin), breaks = 100)
hist(log10(HGwDiffmax/HGwDiffmin), breaks = 100, col = 2, add = T, density = 0)

# dirty insertion/deletion

# as the size of a deletion increases the chance that it overlaps a complex region increases. 
# this is a problem becasue all complex regions are filtered out.


mmGapAll.gr[MMwDiffmax/MMwDiffmin > 3][width(mmGapAll.gr[MMwDiffmax/MMwDiffmin > 3])>1000]


# mouse turnover
# human bases not shared with mouse since divergence
(sum(width(subsetByOverlaps(hgRaw.gr, mmWindow.gr, type = "within"))) - sum(width(intRange))) / sum(width(intRange))

#human turnover
# mouse bases not shared with human since divergence
(sum(width(subsetByOverlaps(mmRaw.gr, hgWindow.gr, type = "within"))) - sum(width(intRange))) / sum(width(intRange))


# what cause bases to become unconserved? 
# becareful cause of the blocks we are using
# since we are only using a section on the genome

subset(mmRaw.gr, start = 12619626, end = 64321910)

windowInt <- intersect(hgWindow.gr, mmWindow.gr)

a <- subsetByOverlaps(hgRaw.gr, windowInt, type = "within")
b <- subsetByOverlaps(mmRaw.gr,windowInt, type = "within")

(sum(width(a))-sum(width(intRange))) / sum(width(a))
(sum(width(b))-sum(width(intRange))) / sum(width(b))

bins <- slidingWindows(x = windowInt,width = 50e3, step = 50e3)[[1]]

ola <- findOverlaps(a, bins)

olb <- findOverlaps(b, bins)

olInt <- findOverlaps(intRange, bins)

olHgRefU <- findOverlaps(hgRefU[width(hgRefU) > 20], bins)

olMmRefU <- findOverlaps(mmRefU[width(mmRefU) > 20], bins)

olHgQueU <- findOverlaps(hgQueU, bins)

olMmQueU <- findOverlaps(mmQueU, bins)

aggA <- aggregate(width(a[queryHits(ola)]), by = list(subjectHits(ola)), sum)
aggB <- aggregate(width(b[queryHits(olb)]), by = list(subjectHits(olb)), sum)

aggInt <- aggregate(width(intRange[queryHits(olInt)]), by = list(subjectHits(olInt)), sum)
aggHgRefU <- aggregate(width(hgRefU[width(hgRefU) > 20][queryHits(olHgRefU)]), by = list(subjectHits(olHgRefU)), sum)
aggMmRefU <- aggregate(width(mmRefU[width(mmRefU) > 20][queryHits(olMmRefU)]), by = list(subjectHits(olMmRefU)), sum)

aggHgQueU <- aggregate(width(mcols(hgQueU[queryHits(olHgQueU)])$queRanges), by = list(subjectHits(olHgQueU)), sum)
aggMmQueU <- aggregate(width(mcols(mmQueU[queryHits(olMmQueU)])$queRanges), by = list(subjectHits(olMmQueU)), sum)

df1 <- merge(merge(aggInt, aggHgRefU, by = 1), aggMmRefU , by = 1)

plot(df[,1], (df[,3]/df[,2]), type = "p", pch = 16, cex = .3, ylim = c(0,.1))
lines(smooth.spline(df[,1], (df[,3]/df[,2]),nknots = 800))

points(df[,1],(df[,4]/ df[,2]), col = 2, type = "p", pch = 16, cex = .3, ylim = c(0,1))
lines(smooth.spline(df[,1],(df[,4]/df[,2])))

df2 <- merge(merge(aggInt, aggHgQueU, by = 1), aggMmQueU , by = 1)

plot(df2[,1], (df2[,3]/df2[,2]), type = "p", pch = 16, cex = .3, ylim = c(0,1))
lines(smooth.spline(df2[,1], (df2[,3]/df2[,2])))

points(df2[,1],(df2[,4]/ df2[,2]), col = 2, type = "p", pch = 16, cex = .3, ylim = c(0,1))
lines(smooth.spline(df2[,1],(df2[,4]/df2[,2])))




plot((df2[,3]/df2[,2]) ,(df[,3]/df[,2]))



df <- merge(merge(merge(aggInt, aggA, by = 1, suffixes = c("int", "hgTO")), aggB, by = 1), aggHgRefU, by = 1)



library(zoo)

layout(c(1,2,3,4))
par(mar = c(1,1,1,1), oma = c(5,5,5,5))

plot(rollmedian(df$Group.1, k = 5) ,rollmedian((df[,2]/df[,3]), k = 5) , type = "l", ylim = c(0, 1))

plot(rollmedian((df1$Group.1), k=3), rollmedian((df1[,3]/df1[,2]), k = 3), type = "l", pch = 16, cex = .3, col = 2)

plot(rollmedian((df2$Group.1), k=11), rollmedian((df2[,4]/df2[,2]), k = 11), type = "l", ylim = c(0,.2), col = 3)

plot(rollmedian(df2$Group.1, k = 5), rollmedian(df2$x.x, k = 5), type = "l")

plot(rollmedian((df1$Group.1), k=5), rollmedian((df1[,3]/df1[,2]), k = 5), type = "l", pch = 16, cex = .3, col = 2)

# so the turn over rate is pretty similar to the alignment rate

plot(df1$Group.1, df1$x.y/df1$x.x, type = "l")
plot(df1$Group.1, df1$x/df1$x.x, type = "l")



sum(width(mmRefU))
length(mmRefU)

sum(width(hgRefU))
length(hgRefU)


width(windowInt)/
length(hgRefU)
sum(width(intRange))


# measuring turn over rates
# the prortion of 

# so if our insertion rates corre;ate to this, then we can map location of turnover between our species. 


# 50% of human bases aren't found in mouse

# 7 % of mouse bases arent found in human


