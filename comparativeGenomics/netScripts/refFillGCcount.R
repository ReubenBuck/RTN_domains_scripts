

# now that we have the data, maybe we can compare it with other data,
# the question is how can we do this
# what data should we compare it to
# what should we be looking for


# There are a large numbre of possibilities at the moment that could all be very large levels of inquury
# What do I think is most important
# What is easiest to work with

# It would be helpful to place markers in our plots to know where we actually are





# lets get GC content

load("~/Desktop/RTN_domains/R_objects/mappedGaps/hg19.mm10.netData.RData")

library(BSgenome.Hsapiens.UCSC.hg19)

seqlengths(refFill.gr) <- seqlengths(refFill.gr) - 1

genoSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, GRangesList(refFill.gr[1:100], refFill.gr[101:200]))
genoSeq2 <- lapply(genoSeq,unlist)
lapply(FUN = letterFrequency, X = genoSeq2, "CG")


# how to shrink an expanded genome back to normal


# removing bases 
# what happens to the regions within 



specRef = "hg19"
specQue = "mm10"

load(paste("Desktop/RTN_domains/R_objects/mappedGaps/",specRef,".",specQue,".netData.RData",sep = ""))
load(paste("Desktop/RTN_domains/R_objects/mappedGaps/",specRef,".stretch.RData", sep = ""))
load(paste("Desktop/RTN_domains/R_objects/mappedGaps/",specRef,".expand.breaks.RData", sep = ""))



synthGenome <- GRanges(seqnames = seqlevels(stretchedRef.gr), 
                       ranges = IRanges(width = seqlengths(stretchedRef.gr), end = seqlengths(stretchedRef.gr)))

# select bin size
binSize <- 2e5

# bin synthetic genome
synthBin.gr <- unlist(slidingWindows(synthGenome, width = binSize, step = binSize))



a.gr <- genoExpandStretch(x.gr = newSynthRefShift[width(newSynthRefShift)>0], 
                          synthGenome = newSynthRefShift, 
                          expandedSeqlengths = seqlengths(stretchedRef.gr))

a.gr
g.gr <- gaps(a.gr)
g.gr <- g.gr[strand(g.gr) == "*"]


olA <- findOverlaps(synthBin.gr, a.gr, select = "last")
olG <- findOverlaps(synthBin.gr, g.gr, select = "last")

# do the shift and change the start

synthBinShrink.gr <- resize(head(synthBin.gr), width = a.gr[head(olA)]$shift, fix = "end")
synthBinShrink.gr <- shift(synthBinShrink.gr, shift = -a.gr[head(olA)]$shift)

# need more than a shift
# need more of a width change

hist((width(refFill.gr)), breaks = 100)
# its about 20000 bp long
# do we get the value for each range and then 
refFillSpecial.gr <- sort(refFill.gr)
refFillSpecial.gr$refRanges = granges(refFillSpecial.gr, use.mcols = FALSE)
refFillSpecial.gr$refRangeID <- 1:length(refFillSpecial.gr)

refFillSpecial.gr <- genoExpandBreak(refFillSpecial.gr, newSynthRefShift, seqlengths(stretchedRef.gr))
# sort them into bins 

# how do we do the bin sorting?

ol <- findOverlaps(refFillSpecial.gr, synthBin.gr)


df <- data.frame(refRangeID = refFillSpecial.gr$refRangeID[queryHits(ol)],
                 synthBinID = subjectHits(ol), 
                 olWidth = width(refFillSpecial.gr[queryHits(ol)]))
library(dplyr)
dfSum <- summarise(group_by(df, refRangeID, synthBinID),
                   olWidth = sum(olWidth))
head(dfSum)

# we can use the same approcah as before to fix the ties
# we can get all teh refFills and get their GC content
RefFillGC.gr <- resize(sort(refFill.gr)[dfSum$refRangeID], width = dfSum$olWidth, fix = "start")

RefFillGC.gr$synthBinID= dfSum$synthBinID

ol <- findOverlaps(RefFillGC.gr)
ol <- ol[!(isSelfHit(ol) | isRedundantHit(ol))]

gcShift <- end(RefFillGC.gr[queryHits(ol)]) - start(RefFillGC.gr[subjectHits(ol)]) + 1

gcShiftSum <- summarise(
  group_by(data_frame(subjectHit = subjectHits(ol), shift = gcShift), 
           subjectHit), shift = sum(shift))

RefFillGC.gr[gcShiftSum$subjectHit] <- shift(RefFillGC.gr[gcShiftSum$subjectHit], shift = gcShiftSum$shift)

# convert that to a list and extract the sequence

RefFillGC.grl <- GenomicRanges::split(RefFillGC.gr,f = RefFillGC.gr$synthBinID)



seqlengths(RefFillGC.grl) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
genoSeq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, RefFillGC.grl)
genoSeq2 <- lapply(genoSeq,unlist)
gcNumber <- lapply(FUN = letterFrequency, X = genoSeq2, "CG")
totalWidth <- lapply(genoSeq2, length)
gcContent <- unlist(gcNumber)/unlist(totalWidth)

plot(gcContent[1:2000], type = "l")


synthBin.gr$gcContent = NA
synthBin.gr$gcContent = as.numeric(synthBin.gr$gcContent)

mcols(synthBin.gr[as.integer(names(RefFillGC.grl))])$gcContent <- as.numeric(gcContent)


plot(start(synthBin.gr[seqnames(synthBin.gr)=="chr16"]), 
     synthBin.gr[seqnames(synthBin.gr)=="chr16"]$gcContent, 
     ylim = c(0.3,0.65) )


# can this be affected by the total number of fill bases



#the idea is return the relevant coordinates
# so we can extract the correct information


# changes to make to our earlier data.
# seqlengths are too long
# mayeb we don't need to add 1 to the end of everything

# it is probably worth learnign the correct convention for dealing with intervals 




dfU <- unique(df)



# which bins each piece belongs to.
# this will give us our GC content

# what sequnce do we consider?
# just fills?
# or more
# nonGap regions?
# nonGaps includes unmappable regions


# fill regions Only includes areas of DNA that align
# gap regions are those between fills in chains
# seq gaps are ignored
# nonRBH gaps
# queRanges Gaps


# if we do the fills first


x.gr <- refFillSpecial.gr

a.gr[head(ol)]
head(g.gr)

genoExpandStretch <- function(x.gr, synthGenome, expandedSeqlengths){
  x.gr <- sort(sortSeqlevels(x.gr))
  seqlengths(x.gr) <- seqlengths(synthGenome)
  olStart <- findOverlaps(x.gr, synthGenome, select = "first")
  olEnd <- findOverlaps(x.gr, synthGenome, select = "last")
  expanded.gr <- GRanges(seqnames = seqnames(x.gr), 
                         ranges = IRanges(start = start(x.gr) + synthGenome[olStart]$shift,
                                          end = end(x.gr) + synthGenome[olEnd]$shift))
  mcols(expanded.gr) <- mcols(x.gr)
  seqlengths(expanded.gr) <- expandedSeqlengths
  genome(expanded.gr) <- genome(x.gr)
  return(expanded.gr)
}

x.gr <- newSynthRefShift

# if it lands in the gap bring it back


# newSynthRefShift might be removing ranges with length 0
# because length has to be at least one


