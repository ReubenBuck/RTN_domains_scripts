
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


opt$query = "~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/brokenChains/hg19/hg19ToMm10.brokenChain"

opt$queryName = "hg19"

opt$baseRate = "~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/baseRates/hg19.base"



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


sum(as.numeric(width(gain.gr)))
sum(width(setdiff(gain.gr, gaps.gr)))



hist(log10(width(setdiff(gain.gr, gaps.gr))), breaks = 500)


hist(log10(width(loss.gr)), breaks = 200)



# need to get rid of unplaced chromosomes

# also remove Ns
mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = opt$queryName)
gaps <- dbGetQuery(mychannel, "SELECT * FROM gap;")
gaps <- gaps[-(grep(pattern = "_", x = gaps$chrom)),]
gaps.gr <- GRanges(seqnames = Rle(gaps$chrom), 
                   ranges = IRanges(start = gaps$chromStart + 1, end = gaps$chromEnd + 1))
seqlevels(gaps.gr) <- chrInfo$chrom
seqlengths(gaps.gr) <- chrInfo$size + 1


sum(width(setdiff(gain.gr, gaps.gr)))



# figured out genome size

# with the increase in size
# the insertion rate, there is likly to be a high amount of false positives and in the deletion side of things there is more likly to be false negatives


# however now we can sample a much larger range of species to get the ancestral DNA content





### run a test to see if as much mouse can be found in human as human can be found in mouse


hm <- read.table(file = "~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/brokenChains/hg19/hg19ToMm10.brokenChain",
              col.names = c("refChr", "refLen", "refStart", "refEnd", "refStrand",
                            "refGap", "queChr", "queLen", "queStart", "queEnd", "queStrand",
                            "queGap", "chainID"),
              colClasses = c("character", "numeric", "numeric", "numeric", "character", "numeric",
                             "character", "numeric", "numeric", "numeric", "character", "numeric",
                             "numeric")
)


mh <- read.table(file = "~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/brokenChains/mm10/mm10ToHg19.brokenChain",
                 col.names = c("refChr", "refLen", "refStart", "refEnd", "refStrand",
                               "refGap", "queChr", "queLen", "queStart", "queEnd", "queStrand",
                               "queGap", "chainID"),
                 colClasses = c("character", "numeric", "numeric", "numeric", "character", "numeric",
                                "character", "numeric", "numeric", "numeric", "character", "numeric",
                                "numeric")
)


hmH.gr <- reduce(GRanges(seqnames = Rle(hm$refChr), 
                 ranges = IRanges(start = hm$refStart, end = hm$refEnd)))
hmM.gr <- reduce(GRanges(seqnames = Rle(hm$queChr), 
                  ranges = IRanges(start = hm$queStart, end = hm$queEnd)))

mhM.gr <- reduce(GRanges(seqnames = Rle(mh$refChr), 
                         ranges = IRanges(start = mh$refStart, end = mh$refEnd)))
mhH.gr <- reduce(GRanges(seqnames = Rle(mh$queChr), 
                         ranges = IRanges(start = mh$queStart, end = mh$queEnd)))








# getting the direction right would improve both estimates
# the improvements in mouse would be better than human
# however, will they be enough to get confident deletion rates in human and insertion rates in mouse


sum(width(mhM.gr)) - sum(width(hmM.gr))
sum(width(hmH.gr)) - sum(width(mhH.gr))


sum(width(hmH.gr)) - sum(width(mhM.gr))
(sum(width(mhM.gr)) - sum(width(hmM.gr))) - (sum(width(hmH.gr)) - sum(width(mhH.gr)))


# we are failing to capture the ancestral mouse DNA 
# why





