#!/usr/bin/env Rscript

rm(list = ls())


library(optparse)


option_list = list(
  make_option(c("-q", "--query"), type="character", default=NA, 
              help="broken UCSC chains between query and subject", metavar="character"),
  make_option(c("-n", "--queryName"), type="character", default=NA, 
              help="query name for output file", metavar="character"),
  make_option(c("-b", "--baseRate"), type="character", default=NA, 
              help="union of outgroup species alignments in query speceis", metavar="character"),
  make_option(c("-H", "--headChain"), type="character", default=NA, 
              help="header for chains", metavar="character"),
  make_option(c("-o", "--outGRangeDir"), type="character", default=NA, 
              help="output dir for GRanges opbject of mapped and unmapped gaps", metavar="character"),
  make_option(c("-O", "--outStatsDir"), type="character", default=NA, 
              help="output dir for summary statistics", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if(any(is.na(opt))){
  stop("argument missing")
}

#.all.chain.broken
opt$query = "~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/brokenChains/mm10/mm10.hg19.all.chain.broken"

opt$queryName = "mm10"

opt$baseRate = "~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/baseRateAllChain/mm10.base"

opt$headChain = "~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/unBrokenChains/mm10/mm10.hg19.all.chain.head"

library(GenomicRanges)
library(RMySQL)
library(dplyr)

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

# need to rotate negative strand 
queSpec[queSpec$queStrand == "-",c("queStart" ,"queEnd")] <- queSpec$queLen[queSpec$queStrand == "-"] - queSpec[queSpec$queStrand == "-",c("queEnd" ,"queStart")]


refSpec.gr <- GRanges(seqnames = Rle(queSpec$refChr),
                      ranges = IRanges(start = queSpec$refStart + 1, end = queSpec$refEnd + 1),
                      chainID = queSpec$chainID
                      )
seqlevels(refSpec.gr) <- chrInfo$chrom
seqlengths(refSpec.gr) <- chrInfo$size + 1

queSpec.gr <- GRanges(seqname = Rle(queSpec$queChr), 
                      ranges = IRanges(start = queSpec$queStart, end = queSpec$queEnd + 1)
                      )
mcols(refSpec.gr)$queRanges = queSpec.gr
refSpec.gr <- sort(sortSeqlevels(refSpec.gr))
aliDF <- data.frame(refSpec.gr)


## read in seqGaps
mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = opt$queryName)
seqGaps <- dbGetQuery(mychannel, "SELECT * FROM gap;")
seqGaps <- seqGaps[-(grep(pattern = "_", x = seqGaps$chrom)),]
seqGaps.gr <- GRanges(seqnames = Rle(seqGaps$chrom), 
                      ranges = IRanges(start = seqGaps$chromStart + 1, end = seqGaps$chromEnd + 1))
seqlevels(seqGaps.gr) <- chrInfo$chrom
seqlengths(seqGaps.gr) <- chrInfo$size + 1


## read in ancestral bases
baseRate <- read.table(opt$baseRate, header = TRUE)
baseRate <- baseRate[-(grep(pattern = "_", x = baseRate$seqnames)),]
baseRate.gr <- GRanges(seqnames = Rle(baseRate$seqnames),
                       ranges = IRanges(start = baseRate$start + 1, end = baseRate$end + 1))
seqlevels(baseRate.gr) <- chrInfo$chrom
seqlengths(baseRate.gr) <- chrInfo$size + 1

## read in header file of chains
headChain <- read.table(file = opt$headChain)
colnames(headChain) <- c("chain", "score", "refChr", "refLen", "refStrand", "refStart", "refEnd", 
                         "queChr", "queLen", "queStrand", "queStart", "queEnd", "chainID")
headChain.gr <- GRanges(seqnames = Rle(headChain$refChr[-grep("_",headChain$refChr)]), 
                        ranges = IRanges(start = headChain$refStart[-grep("_",headChain$refChr)]+1, 
                                         end = headChain$refEnd[-grep("_",headChain$refChr)]),
                        chainID = headChain$chainID[-grep("_",headChain$refChr)])
seqlevels(headChain.gr) <- chrInfo$chrom
seqlengths(headChain.gr) <- chrInfo$size + 1
headChainRed.gr <- reduce(headChain.gr)



# Reduce alignments and find gaps
rGaps <- gaps(reduce(refSpec.gr),start=1L, end=seqlengths(reduce(refSpec.gr)))
rGaps <- rGaps[strand(rGaps) == "*"]
rNoSeqGap <- GenomicRanges::setdiff(x = rGaps, y = seqGaps.gr)
rGaps = rNoSeqGap
rGaps <- rGaps[width(rGaps) > 9]

# annotate gap edges
sRgaps <- rGaps
width(sRgaps) <- 1
olS <- findOverlaps(sRgaps,refSpec.gr,maxgap = 1)

eRgaps <- rGaps
start(eRgaps) <- end(rGaps)
olE <- findOverlaps(eRgaps,refSpec.gr,maxgap = 1)

g.e.DF <- data.frame(queryHits = queryHits(olE), refSpec.gr[subjectHits(olE),])
r.e <- aggregate(g.e.DF$chainID, by = list(g.e.DF$queryHits), FUN = min)
r.e2 <- g.e.DF[paste(g.e.DF$queryHits,g.e.DF$chainID) %in% paste(r.e$Group.1,r.e$x),]

g.s.DF <- data.frame(queryHits = queryHits(olS), refSpec.gr[subjectHits(olS),])
r.s <- aggregate(g.s.DF$chainID, by = list(g.s.DF$queryHits), FUN = min)
r.s2 <- g.s.DF[paste(g.s.DF$queryHits,g.s.DF$chainID) %in% paste(r.s$Group.1,r.s$x),]


rGapsDF <- as.data.frame(rGaps)
rGapsDF$queryHits <- 1:nrow(rGapsDF)
mer <- merge(rGapsDF, r.s2[,c("chainID","queRanges.seqnames", "queRanges.start", "queRanges.end","queryHits")],by="queryHits", all.x = TRUE)
colnames(mer) <- c(colnames(mer)[1:6], paste("left.",colnames(mer)[7:ncol(mer)], sep = ""))
mer <- merge(mer, r.e2[,c("chainID","queRanges.seqnames", "queRanges.start", "queRanges.end","queryHits")],by="queryHits", all.x = TRUE)
colnames(mer) <- c(colnames(mer)[1:10], paste("right.", colnames(mer)[11:ncol(mer)], sep = ""))
allGaps <- GRanges(mer)
mcols(allGaps)$inChain <- overlapsAny(allGaps, headChainRed.gr, minoverlap = 2)
mcols(allGaps)$levelOneChain <- mcols(headChain.gr)$chainID[findOverlaps(allGaps, headChain.gr, minoverlap = 2, select = "first")]


# Annotating bases with ancestry and nonAncestry
brOL <- findOverlaps(allGaps, baseRate.gr)
brInt <- pintersect(allGaps[queryHits(brOL)], baseRate.gr[subjectHits(brOL)])
aggAncestry <- aggregate(width(brInt), by = list(queryHits(brOL)), FUN = sum)
mcols(allGaps)$ancestry[aggAncestry$Group.1] <- aggAncestry$x
mcols(allGaps)$ancestry[is.na(mcols(allGaps)$ancestry)] = 0
mcols(allGaps)$nonAncestry <- width(allGaps) - mcols(allGaps)$ancestry



df <- data.frame(qHit = mcols(allGaps)$queryHits, 
                 ancestry = mcols(allGaps)$ancestry, 
                 nonAncestry = mcols(allGaps)$nonAncestry,
                 chainEqual = mcols(allGaps)$right.chainID == mcols(allGaps)$left.chainID,
                 gapAdjacentOneSide = ( is.na(mcols(allGaps)$left.chainID) & !is.na(mcols(allGaps)$right.chainID) ) | 
                   ( is.na(mcols(allGaps)$right.chainID) & !is.na(mcols(allGaps)$left.chainID) ) ,
                 gapAdjacentTwoSide = is.na(mcols(allGaps)$left.chainID) & is.na(mcols(allGaps)$right.chainID), 
                 inChain = mcols(allGaps)$inChain, 
                 cisRelativeNear = ( abs(mcols(allGaps)$right.queRanges.start - mcols(allGaps)$left.queRanges.end) - width(allGaps) <= 2e4) & mcols(allGaps)$left.queRanges.seqnames == mcols(allGaps)$right.queRanges.seqnames, 
                 cisRelativeFar = ( abs(mcols(allGaps)$right.queRanges.start - mcols(allGaps)$left.queRanges.end) - width(allGaps) > 2e4) & mcols(allGaps)$left.queRanges.seqnames == mcols(allGaps)$right.queRanges.seqnames, 
                 trans = mcols(allGaps)$left.queRanges.seqnames != mcols(allGaps)$right.queRanges.seqnames,
                 deepNested = mcols(allGaps)$left.chainID != mcols(allGaps)$levelOneChain & mcols(allGaps)$right.chainID != mcols(allGaps)$levelOneChain & mcols(allGaps)$inChain)
df[is.na(df)] <- FALSE
for(i in 4:ncol(df)){
  df[,i] = factor(df[,i], levels = c(TRUE, FALSE))
}


dfGroup <- group_by(df, chainEqual, inChain, cisRelativeNear, cisRelativeFar, trans, gapAdjacentOneSide, gapAdjacentTwoSide)
dfSum <- summarise(dfGroup, 
          ancestralBases = sum(ancestry), nonAncestralBases = sum(nonAncestry), totalBases = sum(ancestry + nonAncestry), gapNo = n())
dfSum$ancestralPercent = round(dfSum$ancestralBases / sum(dfSum$ancestralBases) * 100, digits = 2)
dfSum$nonAncestralPercent = round(dfSum$nonAncestralBases / sum(dfSum$nonAncestralBases) * 100, digits = 2)
dfSum$totalPercent = round(dfSum$totalBases / sum(dfSum$totalBases) * 100, digits = 2)
dfSum$gapNoPercent = round(dfSum$gapNo/ sum(dfSum$gapNo) * 100, digits = 2)


for(i in 4:ncol(df)){
  df[,i] = as.logical.factor(df[,i])
}


mapGaps <- allGaps[df$chainEqual]
sumID <- group_indices(dfGroup)[df$chainEqual]
left <- data.frame(chr = mcols(mapGaps)$left.queRanges.seqnames,
                   start = mcols(mapGaps)$left.queRanges.end)
right <- data.frame(chr = mcols(mapGaps)$right.queRanges.seqnames,
                    start = mcols(mapGaps)$right.queRanges.end)
mcols(mapGaps) <- data.frame(sumID, ancestral = mcols(mapGaps)$ancestry, nonAncestral = mcols(mapGaps)$nonAncestry)
mcols(mapGaps)$leftQue = left
mcols(mapGaps)$rightQue = right
mappableGaps <- mapGaps





# left has higher rank
mapGaps <- allGaps[!df$chainEqual & df$inChain & !df$gapAdjacentTwoSide & !df$gapAdjacentOneSide]
sumID <- group_indices(dfGroup)[!df$chainEqual & df$inChain & !df$gapAdjacentTwoSide & !df$gapAdjacentOneSide]
sumID <- sumID[mcols(mapGaps)$left.chainID < mcols(mapGaps)$right.chainID]
mapGaps <- mapGaps[mcols(mapGaps)$left.chainID < mcols(mapGaps)$right.chainID]
leftCoord <- apply(MARGIN = 1,FUN = max,
  X = data.frame(mcols(mapGaps)$left.queRanges.start,
             mcols(mapGaps)$left.queRanges.end)
  )
left = data.frame(chr = mcols(mapGaps)$left.queRanges.seqnames, 
                  start = leftCoord)
right = left
mcols(mapGaps) <- data.frame(sumID, ancestral = mcols(mapGaps)$ancestry, nonAncestral = mcols(mapGaps)$nonAncestry)
mcols(mapGaps)$leftQue = left
mcols(mapGaps)$rightQue = right
mappableGaps <- c(mappableGaps,mapGaps)

#right has higher rank
mapGaps <- allGaps[!df$chainEqual & df$inChain & !df$gapAdjacentTwoSide & !df$gapAdjacentOneSide]
sumID <- group_indices(dfGroup)[!df$chainEqual & df$inChain & !df$gapAdjacentTwoSide & !df$gapAdjacentOneSide]
sumID <- sumID[mcols(mapGaps)$left.chainID > mcols(mapGaps)$right.chainID]
mapGaps <- mapGaps[mcols(mapGaps)$left.chainID > mcols(mapGaps)$right.chainID]
rightCoord <- apply(MARGIN = 1,FUN = min,
                   X = data.frame(mcols(mapGaps)$right.queRanges.start,
                                  mcols(mapGaps)$right.queRanges.end)
)
right = data.frame(chr = mcols(mapGaps)$right.queRanges.seqnames, 
                   start = rightCoord)
left = right
mcols(mapGaps) <- data.frame(sumID, ancestral = mcols(mapGaps)$ancestry, nonAncestral = mcols(mapGaps)$nonAncestry)
mcols(mapGaps)$leftQue = left
mcols(mapGaps)$rightQue = right
mappableGaps <- c(mappableGaps,mapGaps)


# left is not NA
mapGaps <- allGaps[df$gapAdjacentOneSide]
sumID <- group_indices(dfGroup)[df$gapAdjacentOneSide]
sumID <- sumID[!is.na(mcols(mapGaps)$left.chainID)]
mapGaps <- mapGaps[!is.na(mcols(mapGaps)$left.chainID)]
leftCoord <- apply(MARGIN = 1,FUN = max,
                   X = data.frame(mcols(mapGaps)$left.queRanges.start,
                                  mcols(mapGaps)$left.queRanges.end)
)
left = data.frame(chr = mcols(mapGaps)$left.queRanges.seqnames, 
                   start = leftCoord)
right = left
mcols(mapGaps) <- data.frame(sumID, ancestral = mcols(mapGaps)$ancestry, nonAncestral = mcols(mapGaps)$nonAncestry)
mcols(mapGaps)$leftQue = left
mcols(mapGaps)$rightQue = right
mappableGaps <- c(mappableGaps,mapGaps)


# right is not NA
mapGaps <- allGaps[df$gapAdjacentOneSide]
sumID <- group_indices(dfGroup)[df$gapAdjacentOneSide]
sumID <- sumID[!is.na(mcols(mapGaps)$right.chainID)]
mapGaps <- mapGaps[!is.na(mcols(mapGaps)$right.chainID)]
rightCoord <- apply(MARGIN = 1,FUN = min,
                    X = data.frame(mcols(mapGaps)$right.queRanges.start,
                                   mcols(mapGaps)$right.queRanges.end)
)
right = data.frame(chr = mcols(mapGaps)$right.queRanges.seqnames, 
                   start = rightCoord)
left = right
mcols(mapGaps) <- data.frame(sumID, ancestral = mcols(mapGaps)$ancestry, nonAncestral = mcols(mapGaps)$nonAncestry)
mcols(mapGaps)$leftQue = left
mcols(mapGaps)$rightQue = right
mappableGaps <- c(mappableGaps,mapGaps)


#out chain and cis near., this is a 50 50 one
mapGaps <- allGaps[!df$inChain & df$cisRelativeNear]
sumID <- group_indices(dfGroup)[!df$inChain & df$cisRelativeNear]

leftCoord <- apply(MARGIN = 1,FUN = max,
                   X = data.frame(mcols(mapGaps)$left.queRanges.start,
                                  mcols(mapGaps)$left.queRanges.end)
)

rightCoord <- apply(MARGIN = 1,FUN = min,
                   X = data.frame(mcols(mapGaps)$right.queRanges.start,
                                  mcols(mapGaps)$right.queRanges.end)
)
left = data.frame(chr = mcols(mapGaps)$left.queRanges.seqnames, 
                  start = leftCoord)
right = data.frame(chr = mcols(mapGaps)$right.queRanges.seqnames, 
                  start = rightCoord)
mcols(mapGaps) <- data.frame(sumID, ancestral = mcols(mapGaps)$ancestry, nonAncestral = mcols(mapGaps)$nonAncestry)
mcols(mapGaps)$leftQue = left
mcols(mapGaps)$rightQue = right
mappableGaps <- c(mappableGaps,mapGaps)




# select those with both side gaps
unMapGaps <- allGaps[df$gapAdjacentTwoSide]
sumID <- group_indices(dfGroup)[df$gapAdjacentTwoSide]
left = data.frame(chr =rep(NA,length(sumID)), 
                  start = rep(NA,length(sumID)))
right = data.frame(chr = rep(NA,length(sumID)), 
                   start = rep(NA,length(sumID)))
mcols(unMapGaps) <- data.frame(sumID, ancestral = mcols(unMapGaps)$ancestry, nonAncestral = mcols(unMapGaps)$nonAncestry)
mcols(unMapGaps)$leftQue = left
mcols(unMapGaps)$rightQue = right
unMappableGaps <- unMapGaps


unMapGaps <- allGaps[!df$inChain & (df$cisRelativeFar | df$trans)]
sumID <- group_indices(dfGroup)[!df$inChain & (df$cisRelativeFar | df$trans)]
left = data.frame(chr =rep(NA,length(sumID)), 
                  start = rep(NA,length(sumID)))
right = data.frame(chr = rep(NA,length(sumID)), 
                   start = rep(NA,length(sumID)))
mcols(unMapGaps) <- data.frame(sumID, ancestral = mcols(unMapGaps)$ancestry, nonAncestral = mcols(unMapGaps)$nonAncestry)
mcols(unMapGaps)$leftQue = left
mcols(unMapGaps)$rightQue = right
unMappableGaps <- c(unMappableGaps, unMapGaps)


allMapGap <- c(mappableGaps, unMappableGaps)

### Two things to save
## all Map Gaps
### and our summary data

save(allMapGap, file = paste(opt$outGRangeDir,"/",opt$queryName,"AllMapGap.RData",sep = ""))
write.csv(dfSum, file = paste(opt$outStatsDir,"/",opt$queryName,"sumData.csv",sep = ""))

# both species have been put through the pipeline and now they are ready to map over
# we can begin writing then




