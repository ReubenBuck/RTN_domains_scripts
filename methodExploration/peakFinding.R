### read in both human and have a look at supported indels

setwd("~/Desktop/RTN_domains/")


library(GenomicRanges)

canFam3 <- read.table(file = "data/comparativeGenomics/inDel/canFam3/hg19_que.canFam3_ref.indel",header = TRUE)

susScr3 <- read.table(file = "data/comparativeGenomics/inDel/susScr3/hg19_que.susScr3_ref.indel", header = TRUE)



canFam3ins.gr <- GRanges(seqnames = Rle(canFam3$queRange.seqnames[canFam3$inDel == "ins"]), 
                      ranges = IRanges(start = canFam3$queRange.start[canFam3$inDel == "ins"],
                                       end = canFam3$queRange.end[canFam3$inDel == "ins"]))

susScr3ins.gr <- GRanges(seqnames = Rle(susScr3$queRange.seqnames[susScr3$inDel == "ins"]), 
                         ranges = IRanges(start = susScr3$queRange.start[susScr3$inDel == "ins"],
                                          end = susScr3$queRange.end[susScr3$inDel == "ins"]))

hist(log10(width(canFam3ins.gr)), breaks = 100)
hist(log10(width(susScr3ins.gr)), breaks = 100, add = TRUE, col = 2, density = 0)


olins <- findOverlaps(canFam3ins.gr, susScr3ins.gr, )


hist(log10(width(canFam3ins.gr[queryHits(olins)])), breaks = 100)
hist(log10(width(susScr3ins.gr[subjectHits(olins)])), breaks = 100, add = TRUE, col = 2, density = 0)




canFam3ins.gr[queryHits(ol)]





###### del 

canFam3del.gr <- GRanges(seqnames = Rle(canFam3$queRange.seqnames[canFam3$inDel == "del"]), 
                         ranges = IRanges(start = canFam3$queRange.start[canFam3$inDel == "del"] - (canFam3$width[canFam3$inDel == "del"]/2),
                                          width = canFam3$width[canFam3$inDel == "del"]))

susScr3del.gr <- GRanges(seqnames = Rle(susScr3$queRange.seqnames[susScr3$inDel == "del"]), 
                         ranges = IRanges(start = susScr3$queRange.start[susScr3$inDel == "del"] - (susScr3$width[susScr3$inDel == "del"]/2),
                                          width = susScr3$width[susScr3$inDel == "del"]))

hist(log10(width(canFam3del.gr)), breaks = 100)
hist(log10(width(susScr3del.gr)), breaks = 100, add = TRUE, col = 2, density = 0)


oldel <- findOverlaps(canFam3del.gr, susScr3del.gr, )


hist(log10(width(canFam3del.gr)), breaks = 50)
hist(log10(width(canFam3del.gr[subjectHits(oldel)])), breaks = 50, add = TRUE, col = 2, density = 0)

hist(log10(width(susScr3del.gr)), breaks = 50)
hist(log10(width(susScr3del.gr[subjectHits(oldel)])), breaks = 50, add = TRUE, col = 2, density = 0)

isDisjoint(susScr3ins.gr)
isDisjoint(canFam3ins.gr)

cov <- coverage(canFam3ins.gr)
sS1 <- slice(cov, lower = 1)
sS2 <- slice(cov, lower = 2)

# why are we geeting overlapping gaps?





covDel <- coverage(c(reduce(susScr3ins.gr), reduce(canFam3ins.gr)))

chr21 <- covDel$chr21
sl1 <- slice(covDel,lower = 1)
sl2 <- slice(covDel, lower = 2,upper=2)

sl3 <- slice(covDel, lower = 3)



g1 <- reduce(GRanges(sl1))

g2 <- GRanges(sl2)

ol <- findOverlaps(g1,g2)


a <- width(g1[queryHits(ol)]) / width(pintersect(g1[queryHits(ol)], g2[subjectHits(ol)]))


plot(log10(width(g1[queryHits(ol)])), log10(width(pintersect(g1[queryHits(ol)], g2[subjectHits(ol)]))))


sl1$chr1
sl2$chr1

?coverage

canFam3del.gr[queryHits(oldel)]
susScr3del.gr[subjectHits(oldel)]



# peak overlap as a percentage of individual peaks 




