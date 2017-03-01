#!/usr/bin/env Rscript


rm(list = ls())


options(stringsAsFactors = FALSE)


library("optparse")


# input files 
# 


option_list = list(
  make_option(c("-d", "--inputDir"), type="character", default=NA, 
              help="directory containing only ref directories", metavar="character"),
  make_option(c("-q", "--queryName"), type="character", default=NA, 
              help="name of query species", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default="./", 
              help="out directory", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if(any(is.na(opt))){
  stop("argument missing")
}


library(GenomicRanges)

queSpecie <- opt$queryName

genomeRefNames <- list.files(opt$inputDir)
#genomeRefNames <- genomeRefNames[-grep("inDel", genomeRefNames)]

allDel.gr <- GRanges()
allIns.gr <- GRanges()
for(ref in genomeRefNames){
  fileName <- paste(opt$inputDir,"/", ref, "/", queSpecie, "_que.", ref,"_ref.indel", sep = "")
  refInDel <- read.table(fileName, header = TRUE)
  # set up the genomic range and add the genome name
  
  refDel.gr <- GRanges(seqnames = Rle(refInDel$queRange.seqnames[refInDel$inDel == "del"]), 
                       ranges = IRanges(start = refInDel$queRange.start[refInDel$inDel == "del"], width = 1),
                       refRange = GRanges(seqnames = Rle(refInDel$seqnames[refInDel$inDel == "del"]),
                                          ranges = IRanges(start = refInDel$start[refInDel$inDel == "del"],
                                                           end = refInDel$end[refInDel$inDel == "del"])),
                       gapWidth = refInDel$width[refInDel$inDel == "del"],
                       genome = ref,
                       indel = "del")
  allDel.gr <- c(allDel.gr, refDel.gr)
  
  refIns.gr <- GRanges(seqnames = Rle(refInDel$queRange.seqnames[refInDel$inDel == "ins"]), 
                       ranges = IRanges(start = refInDel$queRange.start[refInDel$inDel == "ins"], width = 1),
                       refRange = GRanges(seqnames = Rle(refInDel$seqnames[refInDel$inDel == "ins"]),
                                          ranges = IRanges(start = refInDel$start[refInDel$inDel == "ins"],
                                                           end = refInDel$end[refInDel$inDel == "ins"])),
                       gapWidth = refInDel$queRange.width[refInDel$inDel == "ins"],
                       genome = ref, 
                       indel = "ins")
  allIns.gr <- c(allIns.gr, refIns.gr)
  
  
}
allDel.gr <- sort.GenomicRanges(allDel.gr)
allIns.gr <- sort.GenomicRanges(allIns.gr)


olDel <- findOverlaps(allDel.gr)

supDel <- olDel[mcols(allDel.gr)$genome[queryHits(olDel)] != mcols(allDel.gr)$genome[subjectHits(olDel)]]

olWidths <- data.frame(queryHitWidth = mcols(allDel.gr)$gapWidth[queryHits(supDel)],
                       subjectHitWidth = mcols(allDel.gr)$gapWidth[subjectHits(supDel)])

canDel.gr <- allDel.gr[queryHits(supDel)[ 
  abs(olWidths$queryHitWidth - olWidths$subjectHitWidth) < .2 * apply(X = olWidths,FUN = min,MARGIN = 1)]
  ]

# so we want only unique ranges, 
canDel.gr <- canDel.gr[!(duplicated(mcols(canDel.gr)))]

mcols(canDel.gr)$subjectID <- paste("del:", 1:length(canDel.gr), sep = "")
ReduceCanDel.gr <- reduce(resize(canDel.gr, width = 1, fix = "center"))

ol <- findOverlaps(ReduceCanDel.gr ,canDel.gr)

# supporting intervals

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


write.table(x = reducedCandidates.df, file = paste(opt$outDir,"/",opt$queryName,".supportedIndel.merged",sep = "" ), sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE)

write.table(x = allCandidates.df, file = paste(opt$outDir,"/",opt$queryName,".supportedIndel.all",sep = "" ), sep = "\t", quote = FALSE,row.names = FALSE, col.names = TRUE)




