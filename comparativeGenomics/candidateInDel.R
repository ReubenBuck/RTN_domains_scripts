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

# we can read them all in and get the interesting pairs



# should make it like 10% the size of the smaller element. 


olDel <- findOverlaps(allDel.gr)

supDel <- olDel[mcols(allDel.gr)$genome[queryHits(olDel)] != mcols(allDel.gr)$genome[subjectHits(olDel)]]
# layout(1)
# smoothScatter(log10(mcols(allDel.gr)$gapWidth[queryHits(supDel)]), log10(mcols(allDel.gr)$gapWidth[subjectHits(supDel)]), 
#               xlim = c(1,3.5), nrpoints = 0, ylim = c(1,3.5), xaxs = "i", yaxs = "i")
# lines(log10(seq(1,100000,by = 1000)),log10(seq(1,100000,by = 1000)+ (.2* seq(1,100000,by = 1000))), col = 2)
# lines(log10(seq(1,100000,by = 1000)),log10(seq(1,100000,by = 1000)- (.2* seq(1,100000,by = 1000))), col = 2)

olWidths <- data.frame(queryHitWidth = mcols(allDel.gr)$gapWidth[queryHits(supDel)],
                       subjectHitWidth = mcols(allDel.gr)$gapWidth[subjectHits(supDel)])

canDel.gr <- allDel.gr[queryHits(supDel)[ 
  abs(olWidths$queryHitWidth - olWidths$subjectHitWidth) < .2 * apply(X = olWidths,FUN = min,MARGIN = 1)]
  ]

canDel.gr <- GRanges(unique(as.data.frame(canDel.gr)))
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


# ReduceCanDel.gr[supGenNo == 4]
# 
# hist(((mcols(refDel.gr)$gapWidth)), breaks = 100)
# hist(((mcols(ReduceCanDel.gr[mcols(ReduceCanDel.gr)$supGenNo >1])$supGapMean)), breaks = 100, add = TRUE, col = 2, density = 0)
# 
# plot(density(log10((mcols(ReduceCanDel.gr[mcols(ReduceCanDel.gr)$supGenNo == 2])$supGapMean))))
# for(i in 3:9){
# lines(density(log10((mcols(ReduceCanDel.gr[mcols(ReduceCanDel.gr)$supGenNo == i])$supGapMean))))
# }
# ReduceCanDel.gr[overlapsAny(ReduceCanDel.gr,ReduceCanDel.gr[mcols(ReduceCanDel.gr)$subjectHits == "170770,170771,170772"], maxgap = 200)]




olIns <- findOverlaps(allIns.gr)

supIns <- olIns[mcols(allIns.gr)$genome[queryHits(olIns)] != mcols(allIns.gr)$genome[subjectHits(olIns)]]
# layout(1)
# # pair selection look like this
# smoothScatter(log10(mcols(allIns.gr)$gapWidth[queryHits(supIns)]), log10(mcols(allIns.gr)$gapWidth[subjectHits(supIns)]), 
#               xlim = c(1,5), nrpoints = 0, ylim = c(1,5), xaxs = "i", yaxs = "i")
# 
# lines(log10(seq(1,100000,by = 1000)),log10(seq(1,100000,by = 1000)+ (.2* seq(1,100000,by = 1000))), col = 2)
# lines(log10(seq(1,100000,by = 1000)),log10(seq(1,100000,by = 1000)- (.2* seq(1,100000,by = 1000))), col = 2)

olWidths <- data.frame(queryHitWidth = mcols(allIns.gr)$gapWidth[queryHits(supIns)],
                  subjectHitWidth = mcols(allIns.gr)$gapWidth[subjectHits(supIns)])

# diference between region size is < 10% of smallest regions
canIns.gr <- allIns.gr[queryHits(supIns)[ 
  abs(olWidths$queryHitWidth - olWidths$subjectHitWidth) < .2 * apply(X = olWidths,FUN = min,MARGIN = 1)]
  ]
canIns.gr <- GRanges(unique(as.data.frame(canIns.gr)))
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


write.table(x = reducedCandidates.df, file = paste(opt$outDir,"/",opt$queryName,".supportedIndel.merged",sep = "" ))

write.table(x = allCandidates.df, file = paste(opt$outDir,"/",opt$queryName,".supportedIndel.all",sep = "" ))


## combine obects
## Two outputs 
## Suprted ranges 
## merged ranges






# 
# hist(log10((mcols(refIns.gr)$gapWidth)), breaks = 100)
# hist(log10(df$supGapMean), breaks = 100, add = TRUE, col = 2, density = 0)
# hist(log10((mcols(refDel.gr)$gapWidth)), breaks = 100)
# hist(log10((mcols(ReduceCanDel.gr)$supGapMean)), breaks = 100, add = TRUE, col = 2, density = 0)
# 
# 
# 
# # get all ranges candidates that are sufficiently close to each other
# 
# hist(log10(mcols(distanceToNearest(reduce(allDel.gr)))$distance ), breaks = 100)
# 
# 
# refDel.gr[overlapsAny(mcols(refDel.gr)$refRange , mcols(canDel.gr[18])$refRange)]
# 
# refDel.gr[overlapsAny(refDel.gr, canDel.gr[18])]
# 
# 
# canDel.gr[overlapsAny(canDel.gr, refDel.gr[overlapsAny(refDel.gr, canDel.gr[18])] )]
# 
# 
# # now we are controlling for various effects.
# # its probably time to actually look to see if we can confirm any of our mouse and human gaps
# 
# 
# # also need to get the denomonator for the rate
# # that will be the part of the genome 
# 



# should we just merge both tables Ins and Dels, That way the subject hits will make sense




