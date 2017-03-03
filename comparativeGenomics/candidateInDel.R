#!/usr/bin/env Rscript


rm(list = ls())


options(stringsAsFactors = FALSE)


library("optparse")


# input files 
# 


option_list = list(
  make_option(c("-d", "--inputDir"), type="character", default=NA, 
              help="directory containing only ref directories", metavar="character"),
  make_option(c("-q", "--queryGaps"), type="character", default=NA, 
              help="gaps in query genome", metavar="character"),
  make_option(c("-n", "--queryName"), type="character", default=NA, 
              help="name of query genome", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default="./", 
              help="out directory", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if(any(is.na(opt))){
  stop("argument missing")
}


opt$queryGaps <- "Desktop/RTN_domains/data/comparativeGenomics/queGenomes/gaps/mm10.hg19.que.indel"
opt$queryName <- "mm10"
opt$inputDir <- "Desktop/RTN_domains/data/comparativeGenomics/inDel/mappedIndels/"
opt$outDir <- "Desktop/RTN_domains/data/comparativeGenomics/supportedIndels/"



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

qInside <- unique(queryHits(olIns[abs(x - y) < .1 * x]))
qInsideOnly <- olIns[abs(x - y) < .1 * x]

qInsideOL <- olIns[(queryHits(olIns) %in% qInside)]
qOutsideOL <- olIns[!(queryHits(olIns) %in% qInside)]

qIns.gr[unique(queryHits(qInsideOL))]

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

qInside <- unique(queryHits(olDel[abs(x - y) < .1 * x]))
qInsideOnly <- olDel[abs(x - y) < .1 * x]

qInsideOL <- olDel[(queryHits(olDel) %in% qInside)]
qOutsideOL <- olDel[!(queryHits(olDel) %in% qInside)]


sp <- split(x = mcols(allDel.gr)$genome[subjectHits(qInsideOnly)], f = as.factor(queryHits(qInsideOnly)))
supGenoName <- unlist(lapply(sp, paste0, collapse = ","))
supGenoNo <- unlist(lapply(sp, length))
supUniqGenoNo <- unlist(lapply(lapply(sp, unique), length))

# corresponding gap widths
sp <- split(x = mcols(allDel.gr)$gapWidth[subjectHits(qInsideOnly)], f = as.factor(queryHits(qInsideOnly)))
supGapWidth <- unlist(lapply(sp, paste0, collapse = ","))
supGapMean <- as.integer(unlist(lapply(sp, mean)))

sp <- split(x = paste("del:",subjectHits(qInsideOnly), sep = ""), f = as.factor(queryHits(qInsideOnly)))
subjectIDs <- unlist(lapply(sp, paste0, collapse = ","))

df <- data.frame(supGenoNo, supGenoName, supGapWidth, supGapMean, subjectIDs)



supDel.gr <- allDel.gr[subjectHits(qInsideOnly[order(queryHits(qInsideOnly))])]
mcols(supDel.gr) <- data.frame(mcols(supDel.gr), subjectIDs = paste("del:", subjectHits(qInsideOnly)[order(queryHits(qInsideOnly))], sep = ""))
qDelKeep.gr <- qDel.gr[sort(unique(queryHits(qInsideOnly)))]
mcols(qDelKeep.gr) <- data.frame(mcols(qDelKeep.gr), df)



# combins insertions and dletions

candidates.gr <- c(qInsKeep.gr, qDelKeep.gr)
allCandidates <- c(supIns.gr, supDel.gr)
names(allCandidates) <- 1:length(allCandidates)


keepCandidates.df <- as.data.frame(candidates.gr)
allCandidates.df <- as.data.frame(allCandidates)




write.table(x = keepCandidates.df, file = paste(opt$outDir,"/",opt$queryName,".supportedIndel.que",sep = "" ), 
            sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE)


write.table(x = allCandidates.df, file = paste(opt$outDir,"/",opt$queryName,".supportedIndel.ref",sep = "" ), 
            sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE)









