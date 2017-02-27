
# we can start by reading in species and then put them all together and search for 


rm(list = ls())

library(GenomicRanges)

setwd("~/Desktop/RTN_domains/")

queSpecie <- "hg19"

genomeRefNames <- list.files("data/comparativeGenomics/inDel")
genomeRefNames <- genomeRefNames[-grep("inDel", genomeRefNames)]

allDel.gr <- GRanges()
allIns.gr <- GRanges()
for(ref in genomeRefNames){
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


olDel <- findOverlaps(allDel.gr)

supDel <- olDel[mcols(allDel.gr)$genome[queryHits(olDel)] != mcols(allDel.gr)$genome[subjectHits(olDel)]]
layout(1)
smoothScatter(log10(mcols(allDel.gr)$gapWidth[queryHits(supDel)]), log10(mcols(allDel.gr)$gapWidth[subjectHits(supDel)]), 
              xlim = c(1,3.5), nrpoints = 0, ylim = c(1,3.5), xaxs = "i", yaxs = "i")
lines(log10(seq(1,100000,by = 1000)),log10(seq(1,100000,by = 1000)+ (.2* seq(1,100000,by = 1000))), col = 2)
lines(log10(seq(1,100000,by = 1000)),log10(seq(1,100000,by = 1000)- (.2* seq(1,100000,by = 1000))), col = 2)


olWidths <- data.frame(queryHitWidth = mcols(allDel.gr)$gapWidth[queryHits(supDel)],
                       subjectHitWidth = mcols(allDel.gr)$gapWidth[subjectHits(supDel)])

canDel.gr <- allDel.gr[queryHits(supDel)[ 
  abs(olWidths$queryHitWidth - olWidths$subjectHitWidth) < .2 * apply(X = olWidths,FUN = min,MARGIN = 1)]
  ]

canDel.gr <- GRanges(unique(as.data.frame(canDel.gr)))
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

sp <- split(x = subjectHits(ol), f = as.factor(queryHits(ol)))
subjectHits <- unlist(lapply(sp, paste0, collapse = ","))

df <- data.frame(supGenNo, supGenName, supGapWidth, supGapMean, subjectHits)


mcols(ReduceCanDel.gr) <- df 


ReduceCanDel.gr[supGenNo == 4]

hist(((mcols(refDel.gr)$gapWidth)), breaks = 100)
hist(((mcols(ReduceCanDel.gr[mcols(ReduceCanDel.gr)$supGenNo >1])$supGapMean)), breaks = 100, add = TRUE, col = 2, density = 0)

plot(density(log10((mcols(ReduceCanDel.gr[mcols(ReduceCanDel.gr)$supGenNo == 2])$supGapMean))))
for(i in 3:9){
lines(density(log10((mcols(ReduceCanDel.gr[mcols(ReduceCanDel.gr)$supGenNo == i])$supGapMean))))
}
ReduceCanDel.gr[overlapsAny(ReduceCanDel.gr,ReduceCanDel.gr[mcols(ReduceCanDel.gr)$subjectHits == "170770,170771,170772"], maxgap = 200)]

sample(ReduceCanDel.gr)

# we could also look at the kinds of support we're getting from differenr species

# frequency changes but probability density stays the same
# at least with deletion and amount of support

# looking at the effect of support on size 
# studying the effects of support 



olIns <- findOverlaps(allIns.gr)

supIns <- olIns[mcols(allIns.gr)$genome[queryHits(olIns)] != mcols(allIns.gr)$genome[subjectHits(olIns)]]
layout(1)
# pair selection look like this
smoothScatter(log10(mcols(allIns.gr)$gapWidth[queryHits(supIns)]), log10(mcols(allIns.gr)$gapWidth[subjectHits(supIns)]), 
              xlim = c(1,5), nrpoints = 0, ylim = c(1,5), xaxs = "i", yaxs = "i")

lines(log10(seq(1,100000,by = 1000)),log10(seq(1,100000,by = 1000)+ (.2* seq(1,100000,by = 1000))), col = 2)
lines(log10(seq(1,100000,by = 1000)),log10(seq(1,100000,by = 1000)- (.2* seq(1,100000,by = 1000))), col = 2)

olWidths <- data.frame(queryHitWidth = mcols(allIns.gr)$gapWidth[queryHits(supIns)],
                  subjectHitWidth = mcols(allIns.gr)$gapWidth[subjectHits(supIns)])

# diference between region size is < 10% of smallest regions
canIns.gr <- allIns.gr[queryHits(supIns)[ 
  abs(olWidths$queryHitWidth - olWidths$subjectHitWidth) < .2 * apply(X = olWidths,FUN = min,MARGIN = 1)]
  ]
canIns.gr <- GRanges(unique(as.data.frame(canIns.gr)))
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

sp <- split(x = subjectHits(ol), f = as.factor(queryHits(ol)))
subjectHits <- unlist(lapply(sp, paste0, collapse = ","))

df <- data.frame(supGenNo, supGenName, supGapWidth, supGapMean, subjectHits)


mcols(ReduceCanIns.gr) <- df 



hist(log10((mcols(refIns.gr)$gapWidth)), breaks = 100)
hist(log10(df$supGapMean), breaks = 100, add = TRUE, col = 2, density = 0)



hist(log10((mcols(refDel.gr)$gapWidth)), breaks = 100)
hist(log10((mcols(ReduceCanDel.gr)$supGapMean)), breaks = 100, add = TRUE, col = 2, density = 0)


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


# now we are controlling for various effects.
# its probably time to actually look to see if we can confirm any of our mouse and human gaps


# also need to get the denomonator for the rate
# that will be the part of the genome 









