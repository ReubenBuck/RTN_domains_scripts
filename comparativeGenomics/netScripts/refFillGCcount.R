

# we could try removing some of our seq length changes

rm(list = ls())

options(stringsAsFactors = FALSE)



library(dplyr)
library(GenomicRanges)
devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")



# lets get GC content
specRef = "hg19"
specQue = "mm10"


# binned Genome
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/",specRef,".synthBin.RData",sep = ""))
# fills and gaps
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/",specRef,".",specQue,".netData.RData",sep = ""))
#load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/",specRef,".stretch.RData", sep = ""))
# shift levels
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/shiftData/",specRef,".expand.breaks.RData", sep = ""))

# load in sequence data
if(specRef == "hg19"){
  library(BSgenome.Hsapiens.UCSC.hg19); wholeGenoSeq <- BSgenome.Hsapiens.UCSC.hg19
  
}
if(specRef == "mm10"){
  library(BSgenome.Mmusculus.UCSC.mm10); wholeGenoSeq <- BSgenome.Mmusculus.UCSC.mm10
}

# store original fill coordinates into metadata
refFillSpecial.gr <- sort(refFill.gr)
refFillSpecial.gr$refRanges = granges(refFillSpecial.gr, use.mcols = FALSE)
refFillSpecial.gr$refRangeID <- 1:length(refFillSpecial.gr)

# place fills into stretched genome
refFillSpecial.gr <- genoExpandBreak(refFillSpecial.gr, newSynthRefShift, seqlengths(synthBin.gr))
ol <- findOverlaps(refFillSpecial.gr, synthBin.gr)

# identify boundry fills and asign them correctly
df <- data.frame(refRangeID = refFillSpecial.gr$refRangeID[queryHits(ol)],
                 synthBinID = subjectHits(ol), 
                 olWidth = width(refFillSpecial.gr[queryHits(ol)]))
dfSum <- summarise(group_by(df, refRangeID, synthBinID),
                   olWidth = sum(olWidth))

# sort boundary overlapping regions
RefFillGC.gr <- resize(sort(refFill.gr)[dfSum$refRangeID], width = dfSum$olWidth, fix = "start")
RefFillGC.gr$synthBinID= dfSum$synthBinID
# sort ties
ol <- findOverlaps(RefFillGC.gr)
ol <- ol[!(isSelfHit(ol) | isRedundantHit(ol))]
gcShift <- end(RefFillGC.gr[queryHits(ol)]) - start(RefFillGC.gr[subjectHits(ol)]) + 1
gcShiftSum <- summarise(
  group_by(data_frame(subjectHit = subjectHits(ol), shift = gcShift), 
           subjectHit), shift = sum(shift))
RefFillGC.gr[gcShiftSum$subjectHit] <- shift(RefFillGC.gr[gcShiftSum$subjectHit], shift = gcShiftSum$shift)


# convert that to a list
RefFillGC.grl <- GenomicRanges::split(RefFillGC.gr,f = RefFillGC.gr$synthBinID)

# and extract the sequence and count GC content
#seqlengths(RefFillGC.grl) <- seqlengths(wholeGenoSeq)
genoSeq <- getSeq(wholeGenoSeq, RefFillGC.grl)
genoSeq2 <- lapply(genoSeq,unlist)
gcNumber <- lapply(FUN = letterFrequency, X = genoSeq2, "CG")
totalWidth <- lapply(genoSeq2, length)
gcContent <- unlist(gcNumber)/unlist(totalWidth)



synthBin.gr$gcContent = as.numeric(NA)
synthBin.gr$fill = 0

synthBin.gr[as.integer(names(RefFillGC.grl))]$gcContent <- as.numeric(gcContent)
synthBin.gr[as.integer(names(RefFillGC.grl))]$fill <- as.numeric(totalWidth)

save(synthBin.gr, file = paste("Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/",specRef,".synthBin.RData", sep = ""))



