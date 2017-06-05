

pkgs = names(sessionInfo()$otherPkgs)
pkgs = paste('package:', pkgs, sep = "")
lapply(pkg, detach, character.only = TRUE, unload = TRUE, force = TRUE)

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
if(specRef == "hg19")
  library(BSgenome.Hsapiens.UCSC.hg19); wholeGenoSeq <- BSgenome.Hsapiens.UCSC.hg19
if(specRef == "mm10")
  library(BSgenome.Mmusculus.UCSC.mm10); wholeGenoSeq <- BSgenome.Mmusculus.UCSC.mm10


# store original fill coordinates into metadata
refFillSpecial.gr <- sort(refFill.gr)
refFillSpecial.gr$refRanges = granges(refFillSpecial.gr, use.mcols = FALSE)
refFillSpecial.gr$refRangeID <- 1:length(refFillSpecial.gr)

# place fills into stretched genome
refFillSpecial.gr <- genoExpandBreak(refFillSpecial.gr, newSynthRefShift, seqlengths(stretchedRef.gr))


ol <- findOverlaps(refFillSpecial.gr, synthBin.gr)


df <- data.frame(refRangeID = refFillSpecial.gr$refRangeID[queryHits(ol)],
                 synthBinID = subjectHits(ol), 
                 olWidth = width(refFillSpecial.gr[queryHits(ol)]))


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


