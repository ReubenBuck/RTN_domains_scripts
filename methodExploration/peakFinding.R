### read in both human and have a look at supported indels

rm(list = ls())
setwd("~/Desktop/RTN_domains/")


library(GenomicRanges)

canFam3 <- read.table(file = "data/comparativeGenomics/inDel/canFam3/hg19_que.canFam3_ref.indel",header = TRUE)

susScr3 <- read.table(file = "data/comparativeGenomics/inDel/susScr3/hg19_que.susScr3_ref.indel", header = TRUE)



canFam3ins.gr <- reduce(GRanges(seqnames = Rle(canFam3$queRange.seqnames[canFam3$inDel == "ins"]), 
                      ranges = IRanges(start = canFam3$queRange.start[canFam3$inDel == "ins"],
                                       end = canFam3$queRange.end[canFam3$inDel == "ins"])))

susScr3ins.gr <- reduce(GRanges(seqnames = Rle(susScr3$queRange.seqnames[susScr3$inDel == "ins"]), 
                         ranges = IRanges(start = susScr3$queRange.start[susScr3$inDel == "ins"],
                                          end = susScr3$queRange.end[susScr3$inDel == "ins"])))


olins <- findOverlaps(canFam3ins.gr, susScr3ins.gr )




###### del 

canFam3del.gr <- reduce(GRanges(seqnames = Rle(canFam3$queRange.seqnames[canFam3$inDel == "del"]), 
                         ranges = IRanges(start = canFam3$queRange.start[canFam3$inDel == "del"] - (canFam3$width[canFam3$inDel == "del"]/2),
                                          width = canFam3$width[canFam3$inDel == "del"])))

susScr3del.gr <- reduce(GRanges(seqnames = Rle(susScr3$queRange.seqnames[susScr3$inDel == "del"]), 
                         ranges = IRanges(start = susScr3$queRange.start[susScr3$inDel == "del"] - (susScr3$width[susScr3$inDel == "del"]/2),
                                          width = susScr3$width[susScr3$inDel == "del"])))


oldel <- findOverlaps(canFam3del.gr, susScr3del.gr)



pdf(file = "plots/inDelIdentify/supportedInDelSize.pdf", height = 10, width = 12)
layout(matrix(1:4, byrow = FALSE, nrow = 2))
hist(log10(width(canFam3ins.gr)), breaks = seq(from =0, to = 10, by = .1), main = "hg19 insertions", 
     xlab = "width (bp)", xaxt = "n", yaxs = "i", xaxs = "i", xlim = c(0,5), ylim = c(0,2.8e5))
axis(1, at = 0:10,labels = as.expression(lapply(X = 0:10, FUN = function(y)bquote(10^.(y)))), line=0)
hist(log10(width(canFam3ins.gr[queryHits(olins)])), breaks = 50, add = TRUE, col = 2, density = 0)
legend("topright", c("canFam3 ref sites", "supported in susScr3"),  fill = c(1,2), bty = "n")

hist(log10(width(susScr3ins.gr)), breaks = seq(from =0, to = 10, by = .1), main = "hg19 insertions", 
     xlab = "width (bp)", xaxt = "n", yaxs = "i", xaxs = "i", xlim = c(0,5), ylim = c(0,2.8e5))
axis(1, at = 0:10,labels = as.expression(lapply(X = 0:10, FUN = function(y)bquote(10^.(y)))), line=0)
hist(log10(width(susScr3ins.gr[subjectHits(olins)])), breaks = 50, add = TRUE, col = 2, density = 0)
legend("topright", c("susScr3 ref sites", "supported in canFam3"),  fill = c(1,2), bty = "n")


hist(log10(width(canFam3del.gr)), breaks = seq(from =0, to = 10, by = .1), main = "hg19 deletions", 
     xlab = "width (bp)", xaxt = "n", yaxs = "i", xaxs = "i", xlim = c(0,3), ylim = c(0,1.6e5))
axis(1, at = 0:10,labels = as.expression(lapply(X = 0:10, FUN = function(y)bquote(10^.(y)))), line=0)
hist(log10(width(canFam3del.gr[queryHits(oldel)])), breaks = 50, add = TRUE, col = 2, density = 0)
legend("topright", c("canFam3 ref sites", "supported in susScr3"),  fill = c(1,2), bty = "n")

hist(log10(width(susScr3del.gr)), breaks = seq(from =0, to = 10, by = .1), main = "hg19 deletions", 
     xlab = "width (bp)", xaxt = "n", yaxs = "i", xaxs = "i", xlim = c(0,3), ylim = c(0,1.6e5))
axis(1, at = 0:10,labels = as.expression(lapply(X = 0:10, FUN = function(y)bquote(10^.(y)))), line=0)
hist(log10(width(susScr3del.gr[subjectHits(oldel)])), breaks = 50, add = TRUE, col = 2, density = 0)
legend("topright", c("susScr3 ref sites", "supported in canFam3"),  fill = c(1,2), bty = "n")


dev.off()






canFam3ins.gr <- GRanges(seqnames = Rle(canFam3$queRange.seqnames[canFam3$inDel == "ins"]), 
                                ranges = IRanges(start = canFam3$queRange.start[canFam3$inDel == "ins"],
                                                 end = canFam3$queRange.end[canFam3$inDel == "ins"]))

susScr3ins.gr <- GRanges(seqnames = Rle(susScr3$queRange.seqnames[susScr3$inDel == "ins"]), 
                                ranges = IRanges(start = susScr3$queRange.start[susScr3$inDel == "ins"],
                                                 end = susScr3$queRange.end[susScr3$inDel == "ins"]))




#### all the pieces are here to start removing our bad alignments


canFam3.gr <- GRanges(seqnames = Rle(canFam3$queRange.seqnames), 
                      ranges = IRanges(start = canFam3$queRange.start,
                                       end = canFam3$queRange.end),
                      refRange = GRanges(seqnames = canFam3$seqnames,
                                         ranges = IRanges(start = canFam3$start, end = canFam3$end)),
                      chainID = canFam3$chainID,
                      rowNum = 1:nrow(canFam3))


canFam3.gr <- canFam3.gr[order(mcols(canFam3.gr)$chainID)]
ol <- findOverlaps(canFam3.gr, type = "equal")
canFam3.gr <- canFam3.gr[-subjectHits(ol[(duplicated(queryHits(ol)))])]

# this will remove all the acounted for ranges






# taking too long, it might be best to use the overlaps function to split
# pull out the interesting ones to get the 
cU.gr <- canFam3.gr[order(width(canFam3.gr), decreasing = TRUE)]
ol <- findOverlaps(cU.gr)
cUOL.gr <- cU.gr[ subjectHits(ol[queryHits(ol) %in% unique(queryHits(ol[duplicated(queryHits(ol))]) )  ]) ]
mcols(cUOL.gr)$QhitId <- queryHits(ol[queryHits(ol) %in% unique(queryHits(ol[duplicated(queryHits(ol))]) )  ])


i = 1
while(length(unique(mcols(cUOL.gr)$QhitId)) > 0){
  cU.gr <- canFam3.gr[order(width(canFam3.gr), decreasing = TRUE)]
  ol <- findOverlaps(cU.gr)
  cUOL.gr <- cU.gr[ subjectHits(ol[queryHits(ol) %in% unique(queryHits(ol[duplicated(queryHits(ol))]) )  ]) ]
  mcols(cUOL.gr)$QhitId <- queryHits(ol[queryHits(ol) %in% unique(queryHits(ol[duplicated(queryHits(ol))]) )  ])
  
  
  cUOL.gr <- cUOL.gr[order(width(cUOL.gr), decreasing = TRUE)]
  sp <- split(x = mcols(cUOL.gr)$rowNum, f = as.factor(mcols(cUOL.gr)$QhitId))
  l <- unlist(lapply(X = sp, FUN = getElement, name = 1))
  canFam3.gr <- canFam3.gr[!(mcols(canFam3.gr)$rowNum %in% l)]
  i = i+1
  print(i)
}











canFam3.gr <- canFam3.gr[!(mcols(canFam3.gr)$rowNum %in% l)]

cU.gr <- canFam3.gr[order(width(canFam3.gr), decreasing = TRUE)]
ol <- findOverlaps(cU.gr)

cUOL.gr <- cU.gr[ subjectHits(ol[queryHits(ol) %in% unique(queryHits(ol[duplicated(queryHits(ol))]) )  ]) ]

mcols(cUOL.gr)$QhitId <- queryHits(ol[queryHits(ol) %in% unique(queryHits(ol[duplicated(queryHits(ol))]) )  ])

length(unique(mcols(cUOL.gr)$QhitId))

cUOL.gr <- cUOL.gr[order(width(cUOL.gr), decreasing = TRUE)]
sp <- split(x = mcols(cUOL.gr)$rowNum, f = as.factor(mcols(cUOL.gr)$QhitId))
l <- unlist(lapply(X = sp, FUN = getElement, name = 1))


# this chunk of code does the job 




sp <- split(x = cUOL.gr, f = as.factor(mcols(cUOL.gr)$QhitId))


sp <- split(x = mcols(cUOL.gr)$chainID, f = as.factor(mcols(cUOL.gr)$QhitId))

sort(sp[[4]])

a1 <- unlist(lapply(sp, getElement, name = 1))

a2 <- unlist(lapply(sp, getElement, name = 2))

all(a1 < a2)

# this could probably be easier if I quickly filter for gaps in query that overlap non gaps
# It seems a lot of there complex regions do not correspond to actual gaps in the reference



# order by chain Rank and identify unique the unique ranges, then we are just dealing with conflicting chains. 
# before that 

# so need to order by gap size
# group by overlap
# remove largest gap per group
# ungroup and repeat
# until all groups are size 1

# using this method we should end up with a clean set of non-overlapping gaps in our reference.

# it is almost like we want a complex number Rle object








library(dplyr)

dfl <- as.data.frame(sp)

a <- summarise(group_by(dfl, QhitId), max(chainID), row_number())$`row_number()`


lapply(X = sp, FUN = extractROWS)



elementNROWS(sp)


layout(1)


# number of bp that are outside the largest gap in conflict zones


# this will see if the ranges tend to overlap similarly
# width of them all overlapping

# actual size of zones - expected size of zones
pdf(file= "plots/inDelIdentify/genomeScaleError.pdf")

layout(mat = matrix(1:4, nrow = 2))
par(oma = c(5,5,5,5), mar = c(4,4,3,1))


hist(log10(sum(width(reduce(sp))) - (sum(width(sp))/ table(mcols(cUOL.gr)$QhitId)) + 1), 
     breaks = 50, main = "Gap region size\n(ovserved - expected)", xaxt = "n", xlab = "length (bp + 1)")
axis(1, at = seq(0,10,2),labels = as.expression(lapply(X = seq(0,10,2), FUN = function(y)bquote(10^.(y)))), line=0)

hist(log10(sum(width(reduce(sp))) - (sum(width(sp))/ table(mcols(cUOL.gr)$QhitId)) ), 
     breaks = 50, main = "", xaxt = "n", xlab = "length (bp)")
axis(1, at = seq(0,10,2),labels = as.expression(lapply(X = seq(0,10,2), FUN = function(y)bquote(10^.(y)))), line=0)


hist(log10( (sum(width(reduce(sp)))- max(width(sp))) + 1), breaks = 50, 
     main = "Non-nested gap regions", xaxt = "n", xlab = "length (bp + 1)")
axis(1, at = 0:10,labels = as.expression(lapply(X = 0:10, FUN = function(y)bquote(10^.(y)))), line=0)

hist(log10( (sum(width(reduce(sp)))- max(width(sp)))), breaks = 50, 
     main = "", xaxt = "n", xlab = "length (bp)")
axis(1, at = 0:10,labels = as.expression(lapply(X = 0:10, FUN = function(y)bquote(10^.(y)))), line=0)

title("Human overlaping insertion regions in canFam3", outer = TRUE)

dev.off()

# why are we geeting overlapping gaps?




unique(seqnames(cUOL.gr))



ref <- GRanges(seqnames = Rle(susScr3$seqnames), 
               ranges = IRanges(start = susScr3$start, susScr3$end))

# already have the chain IDs of all the gaps.
# just need to remove gaps with low chainIDs


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

cov <- coverage(c(reduce(canFam3ins.gr), reduce(susScr3ins.gr)), weight = c(width(reduce(canFam3ins.gr)), width(reduce(susScr3ins.gr))))


sl <- slice(cov, lower = 1)

hist(max(sl$chr10) / width(sl$chr10), breaks = 5000, ylim = c(0,4000), xlim = c(1,2))





canFam3del.gr <- (GRanges(seqnames = Rle(canFam3$queRange.seqnames[canFam3$inDel == "del"]), 
                                ranges = IRanges(start = canFam3$queRange.start[canFam3$inDel == "del"] - 5,
                                                 width = 10),
                          weight = canFam3$width[canFam3$inDel == "del"]))

susScr3del.gr <- (GRanges(seqnames = Rle(susScr3$queRange.seqnames[susScr3$inDel == "del"]), 
                                ranges = IRanges(start = susScr3$queRange.start[susScr3$inDel == "del"] - 5,
                                                 width = 10),
                          weight = susScr3$width[susScr3$inDel == "del"]))


canFam3del.gr <- canFam3del.gr[countOverlaps(canFam3del.gr, maxgap = 1) == 1]

susScr3del.gr <- susScr3del.gr[countOverlaps(susScr3del.gr, maxgap = 1) == 1]



cov <- coverage(c(canFam3del.gr, susScr3del.gr), weight = c(mcols(canFam3del.gr)$weight, mcols(susScr3del.gr)$weight))


sl <- slice(cov, lower = 1)

hist(max(sl$chr10) / width(sl$chr10), breaks = 1000)





