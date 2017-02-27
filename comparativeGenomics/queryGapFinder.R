#!/usr/bin/env Rscript

# query gap finder, identify single ended gaps between two genomes

rm(list = ls())

options(stringsAsFactors = FALSE)


library("optparse")


# input files 
# 


option_list = list(
  make_option(c("-q", "--query"), type="character", default=NA, 
              help="broken chains of UCSC liftover file containing query species 1 and 2", metavar="character"),
  make_option(c("-qn1", "--queryName1"), type="character", default=NA, 
              help="query1 name for output file", metavar="character"),
  make_option(c("-qn2", "--queryName2"), type="character", default=NA, 
              help="query1 name for output file", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default="./", 
              help="out directory", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if(any(is.na(opt))){
  stop("argument missing")
}


#opt$query = "~/Desktop/RTN_domains/data/chainAlignments/hg19/mm10/hg19.mm10.brokenChain"
#opt$queryName1 = "hg19"
#opt$queryName2 = "mm10"



library("GenomicRanges")
library(RMySQL)



mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = opt$queryName1)
chrInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")

queFile = opt$query
queSpec <- read.table(file = queFile,
                      col.names = c("refChr", "refLen", "refStart", "refEnd", "refStrand",
                                    "refGap", "queChr", "queLen", "queStart", "queEnd", "queStrand",
                                    "queGap", "chainID"),
                      colClasses = c("character", "numeric", "numeric", "numeric", "character", "numeric",
                                     "character", "numeric", "numeric", "numeric", "character", "numeric",
                                     "numeric")
)

# lets do the switch here

queSpecRaw.gr <- GRanges(seqnames = Rle(queSpec$refChr), 
                         ranges = IRanges(start = queSpec$refStart, end = queSpec$refEnd),
                         queRange = GRanges(seqname = Rle(queSpec$queChr),
                                            ranges = IRanges(start = queSpec$queStart, end = queSpec$queEnd)
                         ),
                         chainID = queSpec$chainID,
                         refOrientation = queSpec$refStrand,
                         queOrientation = queSpec$queStrand,
                         refSeqlength = queSpec$refLen,
                         queSeqlength = queSpec$queLen
)

queSpecRaw.gr <- sortSeqlevels(queSpecRaw.gr)
queSpec <- queSpec[order(queSpecRaw.gr),]
queSpecRaw.gr <- sort(queSpecRaw.gr)

seqlevels(queSpecRaw.gr) <- chrInfo$chrom
seqlengths(queSpecRaw.gr) <- chrInfo$size 

# strand is probably important for when we have reversed things and we want to straighten out for compairson. 
# could be one of the last things we fix
# also important is the seqlength


# pulling out the query gaps of which some will be queSpec inertions
queSpecQueGapElements <- (1:(nrow(queSpec)-1))[queSpec$refStart[2:nrow(queSpec)] == queSpec$refEnd[1:(nrow(queSpec)-1)] &
                                                 queSpec$chainID[2:nrow(queSpec)] == queSpec$chainID[1:(nrow(queSpec)-1)]]

queSpecQueGap.gr <- queSpecRaw.gr[queSpecQueGapElements]
start(queSpecQueGap.gr) <- end(queSpecRaw.gr[queSpecQueGapElements])
end(queSpecQueGap.gr) <- start(queSpecRaw.gr[queSpecQueGapElements + 1])
start(mcols(queSpecQueGap.gr)$queRange) <- end(mcols(queSpecRaw.gr[queSpecQueGapElements])$queRange)
end(mcols(queSpecQueGap.gr)$queRange) <- start(mcols(queSpecRaw.gr[queSpecQueGapElements + 1])$queRange)


# pulling out the ref gaps of which some will be queSpec insertions
queSpecRefGapElements <- (1:(nrow(queSpec)-1))[queSpec$queStart[2:nrow(queSpec)] == queSpec$queEnd[1:(nrow(queSpec)-1)] & 
                                                 queSpec$chainID[2:nrow(queSpec)] == queSpec$chainID[1:(nrow(queSpec)-1)]]

queSpecRefGap.gr <- queSpecRaw.gr[queSpecRefGapElements]
start(queSpecRefGap.gr) <- end(queSpecRaw.gr[queSpecRefGapElements])
end(queSpecRefGap.gr) <- start(queSpecRaw.gr[queSpecRefGapElements + 1])
start(mcols(queSpecRefGap.gr)$queRange) <- end(mcols(queSpecRaw.gr[queSpecRefGapElements])$queRange)
end(mcols(queSpecRefGap.gr)$queRange) <- start(mcols(queSpecRaw.gr[queSpecRefGapElements + 1])$queRange)


gapsSepc <- list(Raw.gr = queSpecRaw.gr, queGap = queSpecQueGap.gr, refGap = queSpecRefGap.gr)


# gaps Spec contains gaps on either side






intRange <- reduce(gapsSepc$Raw.gr)
seqlevels(intRange) <- as.character(unique(seqnames(intRange)))
## apply a function to the ranges of each chromosome
startGap = unlist(lapply(end(split(intRange, seqnames(intRange))), min))
endGap = unlist(lapply(end(split(intRange, seqnames(intRange))), max))




# setting the start might be hard with multiple chromosomes
intRangeGaps <- gaps(intRange,
                     start = startGap[seqlevels(intRange)], 
                     end = endGap[seqlevels(intRange)])
start(intRangeGaps) = start(intRangeGaps) -1
end(intRangeGaps) = end(intRangeGaps) + 1

seqlevels(intRange) <- chrInfo$chrom
seqlengths(intRange) <- chrInfo$size 

print("intGap sorted")



# So finding gaps that aren't interupted by chains

olrefQue1 <- as.matrix(findOverlaps(gapsSepc$refGap, intRangeGaps, type = "equal"))
que1Ref.gr <- gapsSepc$refGap[olrefQue1[,1]]
mcols(que1Ref.gr)$inDel <- "del"



que1Que.gr <- subsetByOverlaps(gapsSepc$queGap, intRange)
mcols(que1Que.gr)$inDel <- "ins"



que1U.gr <- c(que1Que.gr, que1Ref.gr)


# adjust orientation

dfQue1 <- as.data.frame(mcols(que1U.gr)[c("queRange", "queOrientation", "queSeqlength")])
dfQue1[dfQue1$queOrientation == "-",c("queRange.start","queRange.end")] <- dfQue1$queSeqlength[dfQue1$queOrientation == "-"] - 
  dfQue1[dfQue1$queOrientation == "-",c("queRange.end","queRange.start")]

mcols(que1U.gr)$queRange <- GRanges(seqnames = Rle(dfQue1$queRange.seqnames),
                                    ranges = IRanges(start = dfQue1$queRange.start, 
                                                     end = dfQue1$queRange.end))




print("got indels")

# Remove gaps overlapping blocks in query
que1U.gr <- que1U.gr[ !overlapsAny(mcols(que1U.gr)$queRange, mcols(gapsSepc$Raw.gr)$queRange, minoverlap = 2) ]
mcols(que1U.gr)$rowNum <- 1:length(que1U.gr)

# remove both equal overlapping gaps and non-equal overlapping gaps until gaps are at depth 1

queU.gr <- get("que1U.gr")

queU.gr <- queU.gr[order(mcols(queU.gr)$chainID)]
ol <- findOverlaps(mcols(queU.gr)$queRange, type = "equal")
if(length(ol) != length(queU.gr)){
  queU.gr <- queU.gr[-subjectHits(ol[(duplicated(queryHits(ol)))])]
}

cU.gr <- queU.gr[order(width(mcols(queU.gr)$queRange), decreasing = TRUE)]
ol <- findOverlaps(mcols(cU.gr)$queRange)
cUOL.gr <- cU.gr[ subjectHits(ol[queryHits(ol) %in% unique(queryHits(ol[duplicated(queryHits(ol))]) )  ]) ]
mcols(cUOL.gr)$QhitId <- queryHits(ol[queryHits(ol) %in% unique(queryHits(ol[duplicated(queryHits(ol))]) )  ])

while(length(unique(mcols(cUOL.gr)$QhitId)) > 0){
  cU.gr <- queU.gr[order(width(mcols(queU.gr)$queRange), decreasing = TRUE)]
  ol <- findOverlaps(mcols(cU.gr)$queRange)
  cUOL.gr <- cU.gr[ subjectHits(ol[queryHits(ol) %in% unique(queryHits(ol[duplicated(queryHits(ol))]) )  ]) ]
  mcols(cUOL.gr)$QhitId <- queryHits(ol[queryHits(ol) %in% unique(queryHits(ol[duplicated(queryHits(ol))]) )  ])
  
  cUOL.gr <- cUOL.gr[order(width(cUOL.gr), decreasing = TRUE)]
  sp <- split(x = mcols(cUOL.gr)$rowNum, f = as.factor(mcols(cUOL.gr)$QhitId))
  l <- unlist(lapply(X = sp, FUN = getElement, name = 1))
  queU.gr <- queU.gr[!(mcols(queU.gr)$rowNum %in% l)]
}

queU.gr <- queU.gr[order(mcols(queU.gr)$rowNum)]
mcols(queU.gr) <- subset(mcols(queU.gr), select = -rowNum)
assign(x = "que1U.gr",value = queU.gr)



dfQue1 <- as.data.frame(que1U.gr)



print("df conversion")



#### Write out table

que1File <- paste(opt$outDir, "/",opt$queryName1,".", opt$queryName2,".", "que.indel", sep = "")
write.table(x = dfQue1 ,file = que1File, sep = "\t", quote = FALSE, row.names = FALSE)




