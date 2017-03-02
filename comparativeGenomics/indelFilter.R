#!/usr/bin/env Rscript


rm(list = ls())


options(stringsAsFactors = FALSE)


library("optparse")


# input files 
# supported places and que gaps


option_list = list(
  make_option(c("-s", "--supportedIndels"), type="character", default=NA, 
              help="Merged sites in query with support from > 1 reference species", metavar="character"),
  make_option(c("-q", "--queryGaps"), type="character", default=NA, 
              help="gaps in query genome", metavar="character"),
  make_option(c("-n", "--queryName"), type="character", default=NA, 
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


opt$supportedIndels <- "Desktop/RTN_domains/data/comparativeGenomics/supportedIndels/hg19.supportedIndel.merged"

opt$queryGaps <- "Desktop/RTN_domains/data/comparativeGenomics/queGenomes/gaps/hg19.mm10.que.indel"

opt$queryName <- "hg19"


supIndel <- read.table(opt$supportedIndels, header = TRUE)

supIndel.gr <- GRanges(seqnames = Rle(supIndel$seqnames), 
                       ranges = IRanges(start = supIndel$start, end = supIndel$end)
                       )
mcols(supIndel.gr) <- supIndel[,c("supGenNo", "supGenName", "supGapWidth", "supGapMean", "subjectIDs", "indel")]


qGaps <- read.table(opt$queryGaps, header = TRUE)
qGaps.gr <- GRanges(seqnames = qGaps$seqnames, 
                    ranges = IRanges(start = qGaps$start, end = qGaps$end),
                    queRange = GRanges(seqnames = qGaps$queRange.seqnames,
                                       ranges = IRanges(start = qGaps$queRange.start, end = qGaps$queRange.end)),
                    chainID = qGaps$chainID,
                    indel = qGaps$inDel
                    )


supIns.gr <- supIndel.gr[mcols(supIndel.gr)$indel == "del"]
qGapsIns.gr <- qGaps.gr[mcols(qGaps.gr)$indel == "ins"]

#qGaps.gr <- resize(qGaps.gr, width = 1, fix = "center")


ol <- findOverlaps(supIns.gr, qGapsIns.gr, maxgap = 1)



length(supIns.gr[unique(queryHits(ol))])

length(qGaps.gr[unique(subjectHits(ol))])


a <- hist(log10(width(mcols(qGapsIns.gr[unique(subjectHits(ol))])$queRange)), breaks = 100, xlim = c(1,5))
b <- hist(log10(width(mcols(qGapsIns.gr)$queRange)), breaks = 100, xlim = c(1,5))

abline(v = 314, col = 2)
abline(v = 314 * 2, col = 2)
abline(v = 314 / 2, col = 2)

b <- hist((mcols(supIns.gr[unique(queryHits(ol))])$supGapMean), breaks = 100, ylim = c(0,800))

plot(a$mids, a$density, col = 1, type = "l")
lines(b$mids, b$density, col = 2, type = "l")
abline(v = log10(300), col = 2)

# a similar sort of profile however there is something funny going on with the ALus
# maybe they are interupted in some of our species. 

# is it better to use human mouse gaps as the candidates and attempt to support them with our references. 

# dels are a pretty solid match


# so what is going on with our sites we identified

# the proportion of suported places and the proportion of actual gaps. 

hist(log10(mcols(supIndel.gr[unique(queryHits(ol))])$supGapMean ))







# lets talk about our next goal 

# what would I need to be happy to say that ive found what I'm looking for!


# read in both mouse and human 

smoothScatter((width(mcols(qGapsIns.gr[subjectHits(ol)])$queRange)), 
     (mcols(supIns.gr[queryHits(ol)])$supGapMean), 
     xlim = c(80,1000), ylim = c(10,1000), cex = .3, xaxs = "i", yaxs = "i",
     xlab = "query species gaps", ylab = "supported indels")

smoothScatter(log10(width(mcols(qGapsIns.gr[subjectHits(ol)])$queRange)), 
     log10(mcols(supIns.gr[queryHits(ol)])$supGapMean), 
     cex = .3)







