### read in both human and have a look at supported indels

# final stage of filtering would be to check if they actually corresponded to mouse/human gaps.

# for insertions they have to be missing from the mouse genome. 
# needs to overlap a human gap

# for deletions it has to be missing from the human genome, 
# needs to overlap a mouse sided gap in human.

# alternativly we could increase the overlap size to make sure our gaps are species specific.
# not all of our 


rm(list = ls())
setwd("~/Desktop/RTN_domains/")


library(GenomicRanges)

canFam3 <- read.table(file = "data/comparativeGenomics/inDel/canFam3/hg19_que.canFam3_ref.indel",header = TRUE)

susScr3 <- read.table(file = "data/comparativeGenomics/inDel/susScr3/hg19_que.susScr3_ref.indel", header = TRUE)



canFam3ins.gr <- GRanges(seqnames = Rle(canFam3$queRange.seqnames[canFam3$inDel == "ins"]), 
                      ranges = IRanges(start = canFam3$queRange.start[canFam3$inDel == "ins"],
                                       end = canFam3$queRange.end[canFam3$inDel == "ins"]))



susScr3ins.gr <- (GRanges(seqnames = Rle(susScr3$queRange.seqnames[susScr3$inDel == "ins"]), 
                         ranges = IRanges(start = susScr3$queRange.start[susScr3$inDel == "ins"],
                                          end = susScr3$queRange.end[susScr3$inDel == "ins"])))


olins <- findOverlaps(canFam3ins.gr, susScr3ins.gr )




###### del 

canFam3del.gr <- (GRanges(seqnames = Rle(canFam3$queRange.seqnames[canFam3$inDel == "del"]), 
                         ranges = IRanges(start = canFam3$queRange.start[canFam3$inDel == "del"] - (canFam3$width[canFam3$inDel == "del"]/2),
                                          width = canFam3$width[canFam3$inDel == "del"])))

susScr3del.gr <- (GRanges(seqnames = Rle(susScr3$queRange.seqnames[susScr3$inDel == "del"]), 
                         ranges = IRanges(start = susScr3$queRange.start[susScr3$inDel == "del"] - (susScr3$width[susScr3$inDel == "del"]/2),
                                          width = susScr3$width[susScr3$inDel == "del"])))


oldel <- findOverlaps(canFam3del.gr, susScr3del.gr)



pdf(file = "plots/inDelIdentify/supportedInDelSize.pdf", height = 10, width = 12)
layout(matrix(1:4, byrow = FALSE, nrow = 2))
hist(log10(width(canFam3ins.gr)), breaks = seq(from =0, to = 10, by = .05), main = "hg19 insertions", 
     xlab = "width (bp)", xaxt = "n", yaxs = "i", xaxs = "i", xlim = c(0,3), ylim = c(0,1.5e5))
axis(1, at = 0:10,labels = as.expression(lapply(X = 0:10, FUN = function(y)bquote(10^.(y)))), line=0)
hist(log10(width(canFam3ins.gr[queryHits(olins)])), breaks = seq(from =0, to = 10, by = .05), add = TRUE, col = 2, density = 0)
legend("topright", c("canFam3 ref sites", "supported in susScr3"),  fill = c(1,2), bty = "n")

hist(log10(width(susScr3ins.gr)), breaks = seq(from =0, to = 10, by = .05), main = "hg19 insertions", 
     xlab = "width (bp)", xaxt = "n", yaxs = "i", xaxs = "i", xlim = c(0,3), ylim = c(0,1.5e5))
axis(1, at = 0:10,labels = as.expression(lapply(X = 0:10, FUN = function(y)bquote(10^.(y)))), line=0)
hist(log10(width(susScr3ins.gr[subjectHits(olins)])), breaks = seq(from =0, to = 10, by = .05), add = TRUE, col = 2, density = 0)
legend("topright", c("susScr3 ref sites", "supported in canFam3"),  fill = c(1,2), bty = "n")


hist(log10(width(canFam3del.gr)),breaks = seq(from =0, to = 10, by = .05), main = "hg19 deletions", 
     xlab = "width (bp)", xaxt = "n", yaxs = "i", xaxs = "i", xlim = c(0,3), ylim = c(0,1.5e5))
axis(1, at = 0:10,labels = as.expression(lapply(X = 0:10, FUN = function(y)bquote(10^.(y)))), line=0)
hist(log10(width(canFam3del.gr[queryHits(oldel)])),breaks = seq(from =0, to = 10, by = .05), add = TRUE, col = 2, density = 0)
legend("topright", c("canFam3 ref sites", "supported in susScr3"),  fill = c(1,2), bty = "n")

hist(log10(width(susScr3del.gr)), breaks = seq(from =0, to = 10, by = .05), main = "hg19 deletions", 
     xlab = "width (bp)", xaxt = "n", yaxs = "i", xaxs = "i", xlim = c(0,3), ylim = c(0,1.5e5))
axis(1, at = 0:10,labels = as.expression(lapply(X = 0:10, FUN = function(y)bquote(10^.(y)))), line=0)
hist(log10(width(susScr3del.gr[subjectHits(oldel)])),breaks = seq(from =0, to = 10, by = .05), add = TRUE, col = 2, density = 0)
legend("topright", c("susScr3 ref sites", "supported in canFam3"),  fill = c(1,2), bty = "n")


dev.off()



### grouping and overlap

canFam3del.gr <- GRanges(seqnames = Rle(canFam3$queRange.seqnames[canFam3$inDel == "del"]), 
                          ranges = IRanges(start = canFam3$queRange.start[canFam3$inDel == "del"], width = 1),
                          refRange = GRanges(seqnames = Rle(canFam3$seqnames[canFam3$inDel == "del"]),
                                            ranges = IRanges(start = canFam3$start[canFam3$inDel == "del"],
                                                             end = canFam3$end[canFam3$inDel == "del"])),
                         refWidth = canFam3$width[canFam3$inDel == "del"],
                         genome = "canFam3")

susScr3del.gr <- GRanges(seqnames = Rle(susScr3$queRange.seqnames[susScr3$inDel == "del"]), 
                          ranges = IRanges(start = susScr3$queRange.start[susScr3$inDel == "del"], width = 1),
                         refRange = GRanges(seqnames = Rle(susScr3$seqnames[susScr3$inDel == "del"]),
                                            ranges = IRanges(start = susScr3$start[susScr3$inDel == "del"],
                                                             end = susScr3$end[susScr3$inDel == "del"])),
                         refWidth = susScr3$width[susScr3$inDel == "del"],
                         genome = "susScr3")

dels.gr <- c(canFam3del.gr,susScr3del.gr)
dels.gr <- sort.GenomicRanges(dels.gr)


oldel <- findOverlaps(dels.gr)
a <- oldel[mcols(dels.gr)$genome[queryHits(oldel)] != mcols(dels.gr)$genome[subjectHits(oldel)]]
layout(1)
smoothScatter(log10(mcols(dels.gr)$refWidth[queryHits(a)]), log10(mcols(dels.gr)$refWidth[subjectHits(a)]), 
              xlim = c(1,3), nrpoints = 0, ylim = c(1,3), xaxs = "i", yaxs = "i")




dels.gr[ queryHits(a[mcols(dels.gr)$refWidth[queryHits(a)] == mcols(dels.gr)$refWidth[subjectHits(a) ]])][30:42]


oldel <- findOverlaps(dels.gr, maxgap = 5)
a <- oldel[mcols(dels.gr)$genome[queryHits(oldel)] != mcols(dels.gr)$genome[subjectHits(oldel)]]
layout(1)
smoothScatter(log10(mcols(dels.gr)$refWidth[queryHits(a)]), log10(mcols(dels.gr)$refWidth[subjectHits(a)]), 
              xlim = c(1,3), nrpoints = 0, ylim = c(1,3), xaxs = "i", yaxs = "i")





# we can use the center coordinates 

canFam3ins.gr <- GRanges(seqnames = Rle(canFam3$queRange.seqnames[canFam3$inDel == "ins"]), 
                         ranges = IRanges(start = canFam3$queRange.start[canFam3$inDel == "ins"] + 
                                            (canFam3$queRange.width[canFam3$inDel == "ins"] / 2),
                                           width = 1),
                         refRange = GRanges(seqnames = Rle(canFam3$seqnames[canFam3$inDel == "ins"]),
                                            ranges = IRanges(start = canFam3$start[canFam3$inDel == "ins"],
                                                             end = canFam3$end[canFam3$inDel == "ins"])),
                         queWidth = canFam3$queRange.width[canFam3$inDel == "ins"],
                         genome = "canFam3")

susScr3ins.gr <- GRanges(seqnames = Rle(susScr3$queRange.seqnames[susScr3$inDel == "ins"]), 
                         ranges = IRanges(start = susScr3$queRange.start[susScr3$inDel == "ins"] + 
                                            (susScr3$queRange.width[susScr3$inDel == "ins"] / 2),
                                          width = 1),
                         refRange = GRanges(seqnames = Rle(susScr3$seqnames[susScr3$inDel == "ins"]),
                                            ranges = IRanges(start = susScr3$start[susScr3$inDel == "ins"],
                                                             end = susScr3$end[susScr3$inDel == "ins"])),
                         queWidth = susScr3$queRange.width[susScr3$inDel == "ins"],
                         genome = "susScr3")

ins.gr <- c(canFam3ins.gr,susScr3ins.gr)
ins.gr <- sort.GenomicRanges(ins.gr)


olins <- findOverlaps(ins.gr)
a <- olins[mcols(ins.gr)$genome[queryHits(olins)] != mcols(ins.gr)$genome[subjectHits(olins)]]
layout(1)
smoothScatter(log10(mcols(ins.gr)$queWidth[queryHits(a)]), log10(mcols(ins.gr)$queWidth[subjectHits(a)]), 
              xlim = c(1,5), nrpoints = 0, ylim = c(1,5), xaxs = "i", yaxs = "i")



ins.gr[ queryHits(a[mcols(ins.gr)$queWidth[queryHits(a)] == mcols(ins.gr)$queWidth[subjectHits(a) ]])][20:25]



smoothScatter((mcols(ins.gr)$queWidth[queryHits(a)]), (mcols(ins.gr)$queWidth[subjectHits(a)]),
              cex = .05, xaxs = "i", yaxs = "i", ylim = c(500,17000),xlim = c(500,17000), nrpoints = 0 )



olins <- findOverlaps(ins.gr, maxgap = 5)
a <- olins[mcols(ins.gr)$genome[queryHits(olins)] != mcols(ins.gr)$genome[subjectHits(olins)]]
layout(1)
smoothScatter(log10(mcols(ins.gr)$queWidth[queryHits(a)]), log10(mcols(ins.gr)$queWidth[subjectHits(a)]), 
              xlim = c(1,5), nrpoints = 0, ylim = c(1,5), xaxs = "i", yaxs = "i")



reduce(canFam3del.gr)

sort(canFam3del.gr[countOverlaps(canFam3del.gr, minoverlap = 0) > 1])






# Having some trouble identifying lineage specific indels
# we can find them when shared between a two species 
# for some reason though we are identifying insertions common to euchondrogleries.


# a potential problem with the insertoi side, gap is on query side, a gap may essentially overlap both species, bu

