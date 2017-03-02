#!/usr/bin/env Rscript


rm(list = ls())


options(stringsAsFactors = FALSE)


library("optparse")


# input files 
# 


option_list = list(
  make_option(c("-d", "--inputDir"), type="character", default=NA, 
              help="directory containing only ref directories", metavar="character"),
  make_option(c("-q", "--queryGaps1"), type="character", default=NA, 
              help="gaps in query 1 genome", metavar="character"),
  make_option(c("-Q", "--queryGaps2"), type="character", default=NA, 
              help="gaps in query 2 genome", metavar="character"),
  make_option(c("-n", "--queryName1"), type="character", default=NA, 
              help="name of query genome 1", metavar="character"),
  make_option(c("-N", "--queryName2"), type="character", default=NA, 
              help="name of query genome 2", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default="./", 
              help="out directory", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if(any(is.na(opt))){
  stop("argument missing")
}


opt$queryGaps1 <- "Desktop/RTN_domains/data/comparativeGenomics/queGenomes/gaps/hg19.mm10.que.indel"
opt$queryGaps2 <- "Desktop/RTN_domains/data/comparativeGenomics/queGenomes/gaps/mm10.hg19.que.indel"

opt$queryName1 <- "hg19"
opt$queryName1 <- "mm10"

opt$inputDir <- "Desktop/RTN_domains/data/comparativeGenomics/inDel/"




library(GenomicRanges)

queSpecie <- opt$queryName

genomeRefNames <- list.files(opt$inputDir)

allDel.gr <- GRanges()
allIns.gr <- GRanges()
for(ref in genomeRefNames){
  fileName <- paste(opt$inputDir,"/", ref, "/", queSpecie, "_que.", ref,"_ref.indel", sep = "")
  refInDel <- read.table(fileName, header = TRUE)
  # set up the genomic range and add the genome name
  
  refDel.gr <- GRanges(seqnames = Rle(refInDel$queRange.seqnames[refInDel$inDel == "del"]), 
                       ranges = IRanges(start = refInDel$queRange.start[refInDel$inDel == "del"], 
                                        end = refInDel$queRange.end[refInDel$inDel == "del"]),
                       refRange = GRanges(seqnames = Rle(refInDel$seqnames[refInDel$inDel == "del"]),
                                          ranges = IRanges(start = refInDel$start[refInDel$inDel == "del"],
                                                           end = refInDel$end[refInDel$inDel == "del"])),
                       gapWidth = refInDel$width[refInDel$inDel == "del"],
                       genome = ref,
                       indel = "del")
  refDel.gr <- resize(refDel.gr, width = 1, fix = "center")
  allDel.gr <- c(allDel.gr, refDel.gr)
  
  # might be something going on here
  refIns.gr <- GRanges(seqnames = Rle(refInDel$queRange.seqnames[refInDel$inDel == "ins"]), 
                       ranges = IRanges(start = refInDel$queRange.start[refInDel$inDel == "ins"],
                                        end = refInDel$queRange.end[refInDel$inDel == "ins"]),
                       refRange = GRanges(seqnames = Rle(refInDel$seqnames[refInDel$inDel == "ins"]),
                                          ranges = IRanges(start = refInDel$start[refInDel$inDel == "ins"],
                                                           end = refInDel$end[refInDel$inDel == "ins"])),
                       gapWidth = refInDel$queRange.width[refInDel$inDel == "ins"],
                       genome = ref, 
                       indel = "ins")
  refIns.gr <- resize(refIns.gr, width = 1, fix = "center")
  allIns.gr <- c(allIns.gr, refIns.gr)
  
  
}
allDel.gr <- sort.GenomicRanges(allDel.gr)
allIns.gr <- sort.GenomicRanges(allIns.gr)


##### Read in query genome and convert to GRange

qGaps <- read.table(opt$queryGaps, header = TRUE)
qGaps.gr <- GRanges(seqnames = qGaps$seqnames, 
                    ranges = IRanges(start = qGaps$start, end = qGaps$end),
                    queRange = GRanges(seqnames = qGaps$queRange.seqnames,
                                       ranges = IRanges(start = qGaps$queRange.start, end = qGaps$queRange.end)),
                    chainID = qGaps$chainID,
                    indel = qGaps$inDel,
                    queryGapID = qGaps$queryGapID
)

qIns.gr <- qGaps.gr[mcols(qGaps.gr)$indel == "ins"]
mcols(qIns.gr)$gapWidth = width(qIns.gr)
qIns.gr <- resize(qIns.gr, width = 1, fix = "center")

qDel.gr <- qGaps.gr[mcols(qGaps.gr)$indel == "del"]
mcols(qDel.gr)$gapWidth = width(mcols(qDel.gr)$queRange)
qDel.gr <- resize(qDel.gr, width = 1, fix = "center")






#### new

# how can we implement a percentage system?
# do a resize or get distances
distIns <- distanceToNearest(allIns.gr, qIns.gr)

distIns <- distIns[mcols(distIns)$distance/mcols(qIns.gr[subjectHits(distIns)])$gapWidth <= .1]
distIns <- as.matrix(distIns)
distIns <- distIns[complete.cases(distIns),]
olIns <- Hits(from = distIns[,2], to = distIns[,1], length(qIns.gr), length(allIns.gr))


# keep pairs where distances are within the width threshold
# so stuff that isn't too far off center, 
# stuff where widths are somewhat consitent


x <- mcols(qIns.gr)$gapWidth[queryHits(olIns)]
y <- mcols(allIns.gr)$gapWidth[subjectHits(olIns)]

#smoothScatter(log10(x), log10(y), xlim = c(1,6), ylim = c(1,6), xlab = "query wdith", ylab = "support width")
#lines(x = log10(seq(1,1e6,1000)), y = log10(seq(1,1e6,1000) - (.25*seq(1,1e6,1000))) , col = 2)
#lines(x = log10(seq(1,1e6,1000)), y = log10(seq(1,1e6,1000) + (.25*seq(1,1e6,1000))) , col = 2)

qInside <- unique(queryHits(olIns[abs(x - y) < .1 * x]))
qInsideOnly <- olIns[abs(x - y) < .1 * x]

qInsideOL <- olIns[(queryHits(olIns) %in% qInside)]
qOutsideOL <- olIns[!(queryHits(olIns) %in% qInside)]

qIns.gr[unique(queryHits(qInsideOL))]

# on qInside we get the support stats

# overlapping gaps 
# agreed gaps


sp <- split(x = mcols(allIns.gr)$genome[subjectHits(qInsideOnly)], f = as.factor(queryHits(qInsideOnly)))
supGenoName <- unlist(lapply(sp, paste0, collapse = ","))
supGenoNo <- unlist(lapply(sp, length))
supUniqGenoNo <- unlist(lapply(lapply(sp, unique), length))

# corresponding gap widths
sp <- split(x = mcols(allIns.gr)$gapWidth[subjectHits(qInsideOnly)], f = as.factor(queryHits(qInsideOnly)))
supGapWidth <- unlist(lapply(sp, paste0, collapse = ","))
supGapMean <- as.integer(unlist(lapply(sp, mean)))


sp <- split(x = paste("ins:",subjectHits(qInsideOnly), sep = ""), f = as.factor(queryHits(qInsideOnly)))
subjectIDs <- unlist(lapply(sp, paste0, collapse = ","))

df <- data.frame(supGenoNo, supGenoName, supGapWidth, supGapMean, subjectIDs)



supIns.gr <- allIns.gr[subjectHits(qInsideOnly[order(queryHits(qInsideOnly))])]
mcols(supIns.gr) <- data.frame(mcols(supIns.gr), subjectIDs = paste("ins:", subjectHits(qInsideOnly)[order(queryHits(qInsideOnly))], sep = ""))
qInsKeep.gr <- qIns.gr[sort(unique(queryHits(qInsideOnly)))]
mcols(qInsKeep.gr) <- data.frame(mcols(qInsKeep.gr), df)

# so finnaly have my high confidence insertions for human 




# hist(log10(mcols(qIns.gr)$gapWidth), breaks = 100, ylim = c(0,100000), main = paste("maxgap = ",maxGapWidth,sep = ""),
#      xlab = "width (log10 bp)")
# hist(log10(mcols(qIns.gr)$gapWidth[unique(queryHits(olIns))]), breaks = 100, add= TRUE, col = 4, density = 0)
# hist(log10(mcols(qIns.gr)$gapWidth[unique(queryHits(qInsideOL))]), breaks = 100, add= TRUE, col = 3, density = 0)
# hist(log10(mcols(qIns.gr)$gapWidth[unique(queryHits(qOutsideOL))]), breaks = 50, col = 2, add = TRUE, density = 0)
# legend("topright", legend = c("all human sided gaps", "overlap support gaps", "similar width", "different width"), fill = c(1,4,3,2))
# 
# sample(qIns.gr[unique(queryHits(qInsideOL))],10)

# where are the regions that do have overlaps but the sizes don't match in any species
# what do they look like

# all so tune it to suport levels too 
# what is the average level of support our 300 bp insertions get?
# 
# hist((mcols(qIns.gr)$gapWidth[unique(queryHits(qInsideOL))]), breaks = 10000, ylim = c(0,800), col = 1, density = 0, xlim = c(10,7000))
# 
# xOut <- mcols(qIns.gr)$gapWidth[queryHits(qOutsideOL)]
# yOut <- mcols(allIns.gr)$gapWidth[subjectHits(qOutsideOL)]
# smoothScatter(log10(xOut), log10(yOut), xlim = c(1,6), ylim = c(1,6), xlab = "query wdith", ylab = "support width")
# 
# xIn <- mcols(qIns.gr)$gapWidth[queryHits(qInsideOL)]
# yIn <- mcols(allIns.gr)$gapWidth[subjectHits(qInsideOL)]
# smoothScatter(log10(xIn), log10(yIn), xlim = c(1,6), ylim = c(1,6), xlab = "query wdith", ylab = "support width")
# 
# a <- identify(log10(xOut), log10(yOut))
# a
# 
# qIns.gr[queryHits(qOutsideOL[a])]
# allIns.gr[subjectHits(qOutsideOL[a])]











#########
######  Lets look at dels
########


distDel <- distanceToNearest(allDel.gr, qDel.gr)

distDel <- distDel[mcols(distDel)$distance/mcols(qDel.gr[subjectHits(distDel)])$gapWidth <= .1]
distDel <- as.matrix(distDel)
distDel <- distDel[complete.cases(distDel),]
olDel <- Hits(from = distDel[,2], to = distDel[,1], length(qDel.gr), length(allDel.gr))


x <- mcols(qDel.gr)$gapWidth[queryHits(olDel)]
y <- mcols(allDel.gr)$gapWidth[subjectHits(olDel)]

smoothScatter(log10(x), log10(y), xlim = c(1,6), ylim = c(1,6), xlab = "query wdith", ylab = "support width")
lines(x = log10(seq(1,1e6,1000)), y = log10(seq(1,1e6,1000) - (.25*seq(1,1e6,1000))) , col = 2)
lines(x = log10(seq(1,1e6,1000)), y = log10(seq(1,1e6,1000) + (.25*seq(1,1e6,1000))) , col = 2)

qInside <- unique(queryHits(olDel[abs(x - y) < .1 * x]))

qInsideOL <- olDel[(queryHits(olDel) %in% qInside)]
qOutsideOL <- olDel[!(queryHits(olDel) %in% qInside)]








# Once we get our sites that's it.

# what do our sites look like?
# do we need to get them supported in the same way as before. 

# it might be useful printing these plots 












##### old


# genomes
sp <- split(x = mcols(canDel.gr)$genome[subjectHits(ol)], f = as.factor(queryHits(ol)))
supGenName <- unlist(lapply(sp, paste0, collapse = ","))
supGenNo <- unlist(lapply(sp, length))
# corresponding gap widths
sp <- split(x = mcols(canDel.gr)$gapWidth[subjectHits(ol)], f = as.factor(queryHits(ol)))
supGapWidth <- unlist(lapply(sp, paste0, collapse = ","))
supGapMean <- as.integer(unlist(lapply(sp, mean)))

sp <- split(x = paste("del:",subjectHits(ol), sep = ""), f = as.factor(queryHits(ol)))
subjectIDs <- unlist(lapply(sp, paste0, collapse = ","))

df <- data.frame(supGenNo, supGenName, supGapWidth, supGapMean, subjectIDs, indel = "del")


mcols(ReduceCanDel.gr) <- df 



olIns <- findOverlaps(allIns.gr)

supIns <- olIns[mcols(allIns.gr)$genome[queryHits(olIns)] != mcols(allIns.gr)$genome[subjectHits(olIns)]]

olWidths <- data.frame(queryHitWidth = mcols(allIns.gr)$gapWidth[queryHits(supIns)],
                  subjectHitWidth = mcols(allIns.gr)$gapWidth[subjectHits(supIns)])

# diference between region size is < 20% of smallest regions
canIns.gr <- allIns.gr[queryHits(supIns)[ 
  abs(olWidths$queryHitWidth - olWidths$subjectHitWidth) < .2 * apply(X = olWidths,FUN = min,MARGIN = 1)]
  ]
canIns.gr <- canIns.gr[!(duplicated(mcols(canIns.gr)))]

mcols(canIns.gr)$subjectID <- paste("ins:", 1:length(canIns.gr), sep = "")
ReduceCanIns.gr <- reduce(resize(canIns.gr, width = 10, fix = "center"))

ol <- findOverlaps(ReduceCanIns.gr ,canIns.gr)

# supporting intervals

# genomes
sp <- split(x = mcols(canIns.gr)$genome[subjectHits(ol)], f = as.factor(queryHits(ol)))
supGenName <- unlist(lapply(sp, paste0, collapse = ","))
supGenNo <- unlist(lapply(sp, length))
# corresponding gap widths
sp <- split(x = mcols(canIns.gr)$gapWidth[subjectHits(ol)], f = as.factor(queryHits(ol)))
supGapWidth <- unlist(lapply(sp, paste0, collapse = ","))
supGapMean <- as.integer(unlist(lapply(sp, mean)))

sp <- split(x =  paste("ins:",subjectHits(ol), sep = ""), f = as.factor(queryHits(ol)))
subjectIDs <- unlist(lapply(sp, paste0, collapse = ","))

df <- data.frame(supGenNo, supGenName, supGapWidth, supGapMean, subjectIDs, indel = "ins")


mcols(ReduceCanIns.gr) <- df 


reducedCandidates <- c(ReduceCanIns.gr, ReduceCanDel.gr)
allCandidates <- c(canIns.gr, canDel.gr)
names(allCandidates) <- 1:length(allCandidates)


reducedCandidates.df <- as.data.frame(reducedCandidates)
allCandidates.df <- as.data.frame(allCandidates)


write.table(x = reducedCandidates.df, file = paste(opt$outDir,"/",opt$queryName,".supportedIndel.merged",sep = "" ), 
            sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE)

write.table(x = allCandidates.df, file = paste(opt$outDir,"/",opt$queryName,".supportedIndel.all",sep = "" ), 
            sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE)




