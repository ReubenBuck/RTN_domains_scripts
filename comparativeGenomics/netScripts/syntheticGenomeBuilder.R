#!/usr/bin/env Rscript

# stretching chromosomes

# could write this as a scripts

# pkgs = names(sessionInfo()$otherPkgs)
# pkgs = paste('package:', pkgs, sep = "")
# lapply(pkg, detach, character.only = TRUE, unload = TRUE, force = TRUE)

rm(list = ls())

library(optparse)
option_list = list(
  make_option(c("-i", "--input"), default=NA, type='character',
              help="input Rdata file containting output from refGapFillStats"),
  make_option(c("-o", "--output"), default=NA, type='character',
              help="stretched genome output, svaed as RData"),
  make_option(c("-s", "--synthGenome"), default=NA, type='character',
              help="synthetic genome file, saved as RData")
)
opt = parse_args(OptionParser(option_list=option_list))

if(any(is.na(opt))){
  stop("missing options")
}



library(dplyr)
library(zoo)
library(GenomicRanges)


load(opt$input)


switchGenome <- function(queGenome.gr){
  GRanges(mcols(queGenome.gr)$queRanges, 
          queRanges = granges(queGenome.gr, use.mcols = FALSE), 
          sData = mcols(queGenome.gr)$sData, 
          chainID = mcols(queGenome.gr)$chainID)
}


## get the loss bases for the query genome
ol <- findOverlaps(queGap.gr,queAncDna.gr)
queGapAnc.gr <- pintersect(queGap.gr[queryHits(ol)], queAncDna.gr[subjectHits(ol)], drop.nohit.ranges=TRUE)

# get the gain bases for the query
ol <- findOverlaps(queGap.gr,
                   GenomicRanges::setdiff(queGap.gr,queAncDna.gr, ignore.strand = TRUE))
queGapNonAnc.gr <- pintersect(queGap.gr[queryHits(ol)], 
                              GenomicRanges::setdiff(queGap.gr,queAncDna.gr,ignore.strand = TRUE)[subjectHits(ol)], 
                              drop.nohit.ranges=TRUE)

# check seperation of ancestral and non ancestral
if(!all(!overlapsAny(queGapAnc.gr, queGapNonAnc.gr))){
  stop("ancestral and non andcetral did not seperate correctly for query")
}


## get the loss bases for the reference genome
ol <- findOverlaps(refGap.gr,refAncDna.gr)
refGapAnc.gr <- pintersect(refGap.gr[queryHits(ol)], refAncDna.gr[subjectHits(ol)], drop.nohit.ranges=TRUE)

# get the gain bases for the ref
ol <- findOverlaps(refGap.gr,
                   GenomicRanges::setdiff(refGap.gr,refAncDna.gr, ignore.strand = TRUE))
refGapNonAnc.gr <- pintersect(refGap.gr[queryHits(ol)], 
                              GenomicRanges::setdiff(refGap.gr,refAncDna.gr, ignore.strand = TRUE)[subjectHits(ol)], 
                              drop.nohit.ranges=TRUE)

# check seperation of ancestral and non ancestral
if(!all(!overlapsAny(refGapAnc.gr, refGapNonAnc.gr))){
  stop("ancestral and non andcetral did not seperate correctly for reference")
}



# set up genomes
queDel.gr <- refGapAnc.gr
queDel.gr$type = "queDel"

refIns.gr <- refGapNonAnc.gr
refIns.gr$type = "refIns"

queIns.gr <- switchGenome(queGapNonAnc.gr)
queIns.gr$type = "queIns"

refDel.gr <- switchGenome(queGapAnc.gr)
refDel.gr$type = "refDel"

mapped.gr <- c(refIns.gr,queDel.gr)
mapped.gr <- sort(sortSeqlevels(mapped.gr))

remapped.gr <- c(refDel.gr, queIns.gr)
remapped.gr <- sort(sortSeqlevels(remapped.gr))
# set widths for remapped bases
end(remapped.gr) = start(remapped.gr) + width(remapped.gr$queRanges)  -1

# need to make room for memmory
rm(refDel.gr, queIns.gr, queDel.gr, refIns.gr)
rm(refGapNonAnc.gr, refGapAnc.gr, queGapNonAnc.gr, queGapAnc.gr)
rm(queFill.gr, queGap.gr, refFill.gr, refGap.gr, refFillGaps.gr, queFillGaps.gr)

# calculate the total amount of shift required for each gap
agg <- aggregate(x = width(remapped.gr$queRanges), by = list(as.character(granges(remapped.gr, use.mcols = FALSE))), FUN = sum)


synthRefShift <- GRanges(agg$Group.1, shift = agg$x)
rm(agg)
seqlevels(synthRefShift) <- seqlevels(remapped.gr)
genome(synthRefShift) <- genome(remapped.gr)
seqlengths(synthRefShift) <- seqlengths(remapped.gr)
synthRefShift <- sort(sortSeqlevels(synthRefShift))
synthRefShift <- resize(synthRefShift, width = 1,fix = "start")

newDNA <- NULL
newSynthRefShift <- NULL
for(chr in seqlevels(synthRefShift)){
  newSynthRefShift0 <- synthRefShift[seqnames(synthRefShift) == chr]
  if(length(newSynthRefShift0) == 0){
    newSynthRefShift0 <- GRanges(seqnames = chr, ranges = IRanges(start = 1, width = 1))
  }else{
    startRange <- GRanges(seqnames = chr, 
                          ranges = IRanges(start = 1, end = start(newSynthRefShift0[1]) - 1), shift = 0)
    newSynthRefShift0 <- c(startRange, newSynthRefShift0)
  }
  newSynthRefShift0$shift[1] <- 0
  if(length(newSynthRefShift0) > 1){
    end(newSynthRefShift0) <- c(start(newSynthRefShift0)[2:length(newSynthRefShift0)] - 1, 
                                seqlengths(newSynthRefShift0)[chr])
  }else{
    end(newSynthRefShift0) <- seqlengths(synthRefShift)[chr]
  }
  newSynthRefShift0$shift <- cumsum(newSynthRefShift0$shift)
  newDNA <- c(newDNA,newSynthRefShift0$shift[length(newSynthRefShift0)] )
  newSynthRefShift <- c(newSynthRefShift, newSynthRefShift0)
}
names(newDNA) <- seqlevels(synthRefShift)
newSynthRefShift <- unlist(GRangesList(newSynthRefShift))

save(newSynthRefShift, file = opt$synthGenome)


# newSynthRefShift spans the entire reference genome containg intervals between gaps, 
# intervals with a specifc betweenGap regions are are shifted by that regions shift value
# This makes room for the gaped areas of the genome


# shift our mapped ref indels to make way for the query indels
ol <- findOverlaps(mapped.gr, newSynthRefShift)
newMapped.gr <- pintersect(mapped.gr[queryHits(ol)], newSynthRefShift[subjectHits(ol)],drop.nohit.ranges=TRUE)
seqlengths(newMapped.gr) <- seqlengths(newMapped.gr) + newDNA
newMapped.gr <- shift(newMapped.gr, shift = newSynthRefShift$shift[subjectHits(ol)])


# shift query regions into gap spaces made in reference
seqlengths(newSynthRefShift) <- seqlengths(newSynthRefShift) + 1
seqlengths(remapped.gr) <- seqlengths(remapped.gr) + 1
ol <- findOverlaps(resize(remapped.gr, width = 1, fix = "start"), shift(newSynthRefShift, 1))
newRemapped.gr <- remapped.gr
seqlengths(newRemapped.gr) <- seqlengths(newMapped.gr)
newRemapped.gr <- shift(newRemapped.gr[queryHits(ol)], newSynthRefShift$shift[subjectHits(ol)])

# make sure there is no overlap between query and ref regions

if(!all(!overlapsAny(newRemapped.gr, newMapped.gr))){
  stop("not sufficint room for remapped gaps in stretched genome")
}


# sort negative strand info
newRemapped.gr[newRemapped.gr$sData == "-"] <- sort(newRemapped.gr[newRemapped.gr$sData == "-"],
                                                    by = ~ seqnames(queRanges) + start(queRanges), decreasing = TRUE)


# identify ties
ol <- findOverlaps(newRemapped.gr)
ol <- ol[!(isSelfHit(ol) | isRedundantHit(ol))]

# shift ties so they become untied
untieShift <- abs(start(newRemapped.gr)[subjectHits(ol)] - end(newRemapped.gr)[queryHits(ol)]) + 1
aggShift <- dplyr::summarise(
  dplyr::group_by(
    dplyr::data_frame(range = subjectHits(ol), untieShift = untieShift),
    range),
  untieShift = sum(untieShift)
)
newRemapped.gr[aggShift$range] <- shift(newRemapped.gr[aggShift$range], aggShift$untieShift)

# make sure there is no overlaps
if(!all(!overlapsAny(newRemapped.gr, newMapped.gr))){
  stop("untied remmaped gaps overlap mapped gaps")
}

# save files 

stretchedRef.gr <- c(newRemapped.gr, newMapped.gr)
stretchedRef.gr <- sort(sortSeqlevels(stretchedRef.gr))

save(stretchedRef.gr, file = opt$output)



