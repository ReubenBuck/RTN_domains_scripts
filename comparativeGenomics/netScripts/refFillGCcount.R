

# now that we have the data, maybe we can compare it with other data,
# the question is how can we do this
# what data should we compare it to
# what should we be looking for


# There are a large numbre of possibilities at the moment that could all be very large levels of inquury
# What do I think is most important
# What is easiest to work with

# It would be helpful to place markers in our plots to know where we actually are





# lets get GC content




devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")


specRef = "hg19"
specQue = "mm10"


library(BSgenome.Hsapiens.UCSC.hg19)


load(paste("Desktop/RTN_domains/R_objects/mappedGaps/",specRef,".",specQue,".netData.RData",sep = ""))
load(paste("Desktop/RTN_domains/R_objects/mappedGaps/",specRef,".stretch.RData", sep = ""))
load(paste("Desktop/RTN_domains/R_objects/mappedGaps/",specRef,".expand.breaks.RData", sep = ""))



synthGenome <- GRanges(seqnames = seqlevels(stretchedRef.gr), 
                       ranges = IRanges(width = seqlengths(stretchedRef.gr), end = seqlengths(stretchedRef.gr)))

# select bin size
binSize <- 2e5

# bin synthetic genome
synthBin.gr <- unlist(slidingWindows(synthGenome, width = binSize, step = binSize))




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



# which bins each piece belongs to.
# this will give us our GC content


