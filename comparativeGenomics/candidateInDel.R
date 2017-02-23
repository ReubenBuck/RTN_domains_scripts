
# we can start by reading in species and then put them all together and search for 


rm(list = ls())

setwd("~/Desktop/RTN_domains/")

queSpecie <- "hg19"

genomeRefNames <- list.files("data/comparativeGenomics/inDel")
genomeRefNames <- genomeRefNames[-grep("inDel", genomeRefNames)]

allDel.gr <- GRanges()
allIns.gr <- GRanges()
for(ref in genomeRefNames[1:4]){
  fileName <- paste("data/comparativeGenomics/inDel/", ref, "/", queSpecie, "_que.", ref,"_ref.indel", sep = "")
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


olDel <- findOverlaps(allDel.gr, maxgap = 10)

supDel <- olDel[mcols(allDel.gr)$genome[queryHits(olDel)] != mcols(allDel.gr)$genome[subjectHits(olDel)]]
layout(1)
smoothScatter(log10(mcols(allDel.gr)$gapWidth[queryHits(supDel)]), log10(mcols(allDel.gr)$gapWidth[subjectHits(supDel)]), 
              xlim = c(1,3.5), nrpoints = 0, ylim = c(1,3.5), xaxs = "i", yaxs = "i")


canDel.gr <- allDel.gr[queryHits(supDel)[ 
  abs(mcols(allDel.gr)$gapWidth[queryHits(supDel)] - mcols(allDel.gr)$gapWidth[subjectHits(supDel)]) < 5]
  ]
canDel.gr <- GRanges(unique(as.data.frame(canDel.gr)))
ReduceCanDel.gr <- reduce(resize(canDel.gr, width = 10, fix = "center"))

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

sp <- split(x = subjectHits(ol), f = as.factor(queryHits(ol)))
subjectHits <- unlist(lapply(sp, paste0, collapse = ","))

df <- data.frame(supGenNo, supGenName, supGapWidth, supGapMean, subjectHits)


mcols(ReduceCanDel.gr) <- df 


ReduceCanDel.gr[supGenNo == 4]





olIns <- findOverlaps(allIns.gr)

supIns <- olIns[mcols(allIns.gr)$genome[queryHits(olIns)] != mcols(allIns.gr)$genome[subjectHits(olIns)]]
layout(1)
smoothScatter(log10(mcols(allIns.gr)$gapWidth[queryHits(supIns)]), log10(mcols(allIns.gr)$gapWidth[subjectHits(supIns)]), 
              xlim = c(1,5), nrpoints = 0, ylim = c(1,5), xaxs = "i", yaxs = "i")

canIns.gr <- allIns.gr[queryHits(supIns)[mcols(allIns.gr)$gapWidth[queryHits(supIns)] == mcols(allIns.gr)$gapWidth[subjectHits(supIns)]]]







g <- GRanges(seqnames = Rle("chr16"), 
             ranges = IRanges(start = 343659, end = 343659))
# so in the mcols, we can record our stats,
# number of supporters
# supporter names
# supporter gap Widths
# we can put ranges back to the center

# its possible to write a function to take care of putting all the names together. 
# and all the numbers. 



paste(c(1,2,3), c(4,5), collapse = ";", sep = ",")

g <- GRanges(seqnames = Rle(c("chr1", "chr1")),
        ranges = IRanges(start = c(1,10), end = c(14,17)), 
        deets = c("3;4", "1;5;6"))



GRanges(allIns.gr)
?GRangesList

as(unlist(allIns.gr), "GRanges")


# get all ranges candidates that are sufficiently close to each other

hist(log10(mcols(distanceToNearest(reduce(allDel.gr)))$distance ), breaks = 100)


refDel.gr[overlapsAny(mcols(refDel.gr)$refRange , mcols(canDel.gr[18])$refRange)]

refDel.gr[overlapsAny(refDel.gr, canDel.gr[18])]


canDel.gr[overlapsAny(canDel.gr, refDel.gr[overlapsAny(refDel.gr, canDel.gr[18])] )]

# the findOverlaps function might be exponentianly increasing overlaps of even an odd umber of things. 





