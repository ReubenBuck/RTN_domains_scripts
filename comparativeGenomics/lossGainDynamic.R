
rm(list = ls())


library(optparse)


option_list = list(
  make_option(c("-q", "--query"), type="character", default=NA, 
              help="broken chains of UCSC liftover file containing query species 1", metavar="character"),
  make_option(c("-n", "--queryName"), type="character", default=NA, 
              help="query1 name for output file", metavar="character"),
  make_option(c("-b", "--baseRate"), type="character", default=NA, 
              help="union of outgroup species alignments in query speceis", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default="./", 
              help="out directory", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if(any(is.na(opt))){
  stop("argument missing")
}

#.all.chain.broken
opt$query = "~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/brokenChains/hg19/hg19.mm10.all.chain.broken"

opt$queryName = "hg19"

opt$baseRate = "~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/baseRateAllChain/hg19.base"



library(GenomicRanges)
library(RMySQL)

mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = opt$queryName)
chrInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")

chrInfo <- chrInfo[-(grep(pattern = "_", x = chrInfo$chrom)),]

queSpec <- read.table(file = opt$query,
                      col.names = c("refChr", "refLen", "refStart", "refEnd", "refStrand",
                                    "refGap", "queChr", "queLen", "queStart", "queEnd", "queStrand",
                                    "queGap", "chainID"),
                      colClasses = c("character", "numeric", "numeric", "numeric", "character", "numeric",
                                     "character", "numeric", "numeric", "numeric", "character", "numeric",
                                     "numeric")
)

queSpec <- queSpec[-(grep(pattern = "_", x = queSpec$refChr)),]


queSpec.gr <- GRanges(seqnames = Rle(queSpec$refChr),
                      ranges = IRanges(start = queSpec$refStart + 1, end = queSpec$refEnd + 1))
seqlevels(queSpec.gr) <- chrInfo$chrom
seqlengths(queSpec.gr) <- chrInfo$size + 1


baseRate <- read.table(opt$baseRate, header = TRUE)
baseRate <- baseRate[-(grep(pattern = "_", x = baseRate$seqnames)),]

baseRate.gr <- GRanges(seqnames = Rle(baseRate$seqnames),
                       ranges = IRanges(start = baseRate$start + 1, end = baseRate$end + 1))
seqlevels(baseRate.gr) <- chrInfo$chrom
seqlengths(baseRate.gr) <- chrInfo$size + 1


mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = opt$queryName)
gaps <- dbGetQuery(mychannel, "SELECT * FROM gap;")
gaps <- gaps[-(grep(pattern = "_", x = gaps$chrom)),]
gaps.gr <- GRanges(seqnames = Rle(gaps$chrom), 
                   ranges = IRanges(start = gaps$chromStart + 1, end = gaps$chromEnd + 1))
seqlevels(gaps.gr) <- chrInfo$chrom
seqlengths(gaps.gr) <- chrInfo$size + 1



u.gr <- union(baseRate.gr, queSpec.gr)

sum(width(u.gr))


gQue.gr <- gaps(reduce(queSpec.gr))
gQue.gr <- gQue.gr[strand(gQue.gr) == "*"]
sum(width(gQue.gr))

#gQue.gr <- gQue.gr[overlapsAny(gQue.gr, queSpec.gr, maxgap = 1)]
#sum(width(gQue.gr))

loss.gr <- intersect(gQue.gr, baseRate.gr)

sum(width(loss.gr))



gBR.gr <- gaps(reduce(baseRate.gr))
gBR.gr <- gBR.gr[strand(gBR.gr) == "*"]
sum(width(gBR.gr))

#gBR.gr <- gBR.gr[overlapsAny(gBR.gr, baseRate.gr, maxgap = 1)]
#sum(width(gBR.gr))

hist(log10(width(reduce(gQue.gr))), breaks = 100)


gain.gr <- intersect(gQue.gr, gBR.gr)
gain.gr <- setdiff(gain.gr, gaps.gr)


sum(as.numeric(width(gain.gr)))
sum(width(gain.gr))



gainHist <- hist((width(setdiff(gain.gr, gaps.gr))), breaks = seq(0,502000), xlim = c(0,1000), ylim = c(0,10000))


lossHist <- hist((width(loss.gr)), breaks = seq(0,502000), xlim = c(0,1000), ylim = c(0,10000))


plot(gainHist$mids, gainHist$counts, type = "l", xlim = c(0,1000), ylim = c(0,10000) )
lines(lossHist$mids, lossHist$counts, type = "l", xlim = c(0,1000), ylim = c(0,10000) , col = 2)


plot(gainHist$mids, gainHist$counts * gainHist$mids, type = "l", xlim = c(0,1000), ylim = c(0,1000000) )
lines(lossHist$mids, lossHist$counts * lossHist$mids, type = "l", xlim = c(0,1000), ylim = c(0,1000000) , col = 2)





# human deleted regions
# where they are located in mouse

# we can look at the gaps 

# need to get rid of unplaced chromosomes




# that data is already available

sum(width(gain.gr))
sum(width(loss.gr))

# gain overlap with repeats 







bins <- GRanges(seqnames = chrInfo$chrom, 
                ranges = IRanges(start = 1, end = chrInfo$size + 1))

bins.gr <- slidingWindows(bins, width = 10000, step = 10000)

bins.gr
bins.gr <- unlist(bins.gr)
bins.gr <- sort(sortSeqlevels(bins.gr))


ol <- findOverlaps(loss.gr, bins.gr)
pInt <- pintersect(bins.gr[subjectHits(ol)], loss.gr[queryHits(ol)])

agg <- aggregate(width(pInt), by = list(subjectHits(ol)), FUN = sum)
plot(agg$x[1:200], type = "l")

mcols(bins.gr)$delBases = 0
mcols(bins.gr)$delBases[agg$Group.1] = agg$x

mcols(bins.gr)
hist((mcols(bins.gr)$delBases), breaks = 100)



ol <- findOverlaps(gain.gr, bins.gr)
pInt <- pintersect(bins.gr[subjectHits(ol)], gain.gr[queryHits(ol)])

agg <- aggregate(width(pInt), by = list(subjectHits(ol)), FUN = sum)
plot(agg$x[1:200], type = "l")

mcols(bins.gr)$insBases = 0
mcols(bins.gr)$insBases[agg$Group.1] = agg$x

mcols(bins.gr)
hist( mcols(bins.gr)$delBases, breaks = 100)
hist((mcols(bins.gr)$insBases), breaks = 100, add = TRUE, col = 2, density = 0)


smoothScatter(mcols(bins.gr)$delBases, mcols(bins.gr)$insBases)



cols <- as.integer(seqnames(bins.gr))%%2

pdf(file = paste("~/Desktop/indelGenome",opt$queryName,".pdf", sep = ""), height = 6,width = 12)
layout(matrix(c(1:2),nrow = 2))
par(mar = c(3,6,0,0), oma = c(3,0,5,3))

plot(mcols(bins.gr)$insBases, type = "p", pch = 19, cex = .2, col = c(1,8)[cols + 1], ylim = c(0,1e5),
     ylab = "insertion (bp)", xaxt = "n")
#lines(smooth.spline(mcols(bins.gr)$insBases), col = 4)

axis(side = 1, labels = substring(unique(seqnames(bins.gr)), first = 4), 
     at = cumsum(table(seqnames(bins.gr))) - (table(seqnames(bins.gr))/2))

mtext(side = 3, text = opt$queryName, cex = 1.5, line = 2)

plot((mcols(bins.gr)$delBases), type = "p", pch = 19, cex = .2, col= c(1,8)[cols + 1], ylim = c(0,1e5),
     ylab = "deletion (bp)", xaxt = "n")
#lines(smooth.spline(mcols(bins.gr)$delBases), col = 2)

axis(side = 1, labels = substring(unique(seqnames(bins.gr)), first = 4), 
     at = cumsum(table(seqnames(bins.gr))) - (table(seqnames(bins.gr))/2))

mtext(side = 1, text = "Chromosome", cex = 1.2, line = 3)


dev.off()


layout(1)

hist(mcols(bins.gr)$insBases, breaks = 100)
hist(mcols(bins.gr)$delBases, breaks = 100, col = 2, density = 0, add = T)


# what are the important things to know 
# why do we want to calculate regional variation?


plot(mcols(bins.gr[seqnames(bins.gr) == "chr1"])$insBases, type = "p", pch = 19, cex = .2, col = c(1,8)[cols + 1], ylim = c(0,1e5),
     ylab = "insertion (bp)", xaxt = "n")
lines(smooth.spline(mcols(bins.gr[seqnames(bins.gr) == "chr1"])$insBases), col = 4)


plot(mcols(bins.gr[seqnames(bins.gr) == "chr1"])$delBases, type = "p", pch = 19, cex = .2, col = c(1,8)[cols + 1], ylim = c(0,1e5),
     ylab = "insertion (bp)", xaxt = "n")
lines(smooth.spline(mcols(bins.gr[seqnames(bins.gr) == "chr1"])$delBases), col = 2)


# proportion of insertions that overlap repeats
# maybe we could use a HMM to define the states
# HMM we can feed states and get a model and know where to draw the line
# We can use bin size and liklihood initially
# there is also the option of cross validation, however we may not have enough data

# At this poitn it is probably worth saving our results and doing humanisation.
# 


A <- as.data.frame(bins.gr)


head(A)


