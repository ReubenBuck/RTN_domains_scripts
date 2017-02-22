#!/usr/bin/env Rscript

# exacutable version of our program 

# read in the alignment blocks 

rm(list = ls())

options(stringsAsFactors = FALSE)

#!/usr/bin/env Rscript

library("optparse")


# input files 
# 


option_list = list(
  make_option(c("-q1", "--query1"), type="character", default=NA, 
              help="UCSC liftover file containing first query species", metavar="character"),
  make_option(c("-q2", "--query2"), type="character", default=NA, 
              help="UCSC liftover file containing seccond query species", metavar="character"),
  make_option(c("-qn1", "--queryName1"), type="character", default=NA, 
              help="query1 name for output file", metavar="character"),
  make_option(c("-qn2", "--queryName2"), type="character", default=NA, 
              help="query1 name for output file", metavar="character"),
  make_option(c("-rn", "--referenceName"), type="character", default=NA, 
              help="name of reference species will be used to name output", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default="./", 
              help="out directory", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if(any(is.na(opt))){
  stop("argument missing")
}



library("GenomicRanges")
#library("RMySQL")


### reading in data files

opt$query1 <- "~/Desktop/RTN_domains/data/chainAlignments/inDelTest/canFam2Hg19.out"
opt$query2 <- "~/Desktop/RTN_domains/data/chainAlignments/inDelTest/canFam2mm9.out"
opt$referenceName <- "canFam2"

#mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = opt$referenceName)


con <- gzcon(url(paste("http://hgdownload.soe.ucsc.edu/goldenPath/",opt$referenceName,"/database/chromInfo.txt.gz", sep="")))
txt <- readLines(con)
dat <- read.delim(textConnection(txt), header = FALSE)
colnames(dat) = c("chrom", "size", "fileName")
chrInfo <- dat

#chrInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")


queSpecies <- list(que1 = NA, que2 = NA )
for(i in 1:2){
  print(i)
  queFile = c(opt$query1, opt$query2)[i]
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
  queSpecies[[i]] = gapsSepc
  
}


# do we make them gRanges

# how do we remove overlapping blocks





intRange <- intersect(queSpecies$que1$Raw.gr, queSpecies$que2$Raw.gr)
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




olrefQue1 <- as.matrix(findOverlaps(queSpecies$que1$refGap, intRangeGaps, type = "equal"))
que1Ref.gr <- queSpecies$que1$refGap[olrefQue1[,1]]

olrefQue2 <- as.matrix(findOverlaps(queSpecies$que2$refGap, intRangeGaps, type = "equal"))
que2Ref.gr <- queSpecies$que2$refGap[olrefQue2[,1]]


que1RefC <- subsetByOverlaps(que1Ref.gr, que2Ref.gr, type = "equal")
que1RefU <- que1Ref.gr[!(que1Ref.gr %in% que1RefC)]

que2RefC <- subsetByOverlaps(que2Ref.gr, que1Ref.gr, type = "equal")
que2RefU <- que2Ref.gr[!(que2Ref.gr %in% que2RefC)]

# making sure that we have ref gaps that 
que1RefU <- subsetByOverlaps(que1RefU, reduce(queSpecies$que2$Raw.gr), type = "within")
mcols(que1RefU)$inDel <- "del"

que2RefU <- subsetByOverlaps(que2RefU, reduce(queSpecies$que1$Raw.gr), type = "within")
mcols(que2RefU)$inDel <- "del"

que1Que.gr <- subsetByOverlaps(queSpecies$que1$queGap, intRange)
que2Que.gr <- subsetByOverlaps(queSpecies$que2$queGap, intRange)


que1QueC <- subsetByOverlaps(que1Que.gr, que2Que.gr)
que1QueU <- que1Que.gr[!(que1Que.gr %in% que1QueC)]
mcols(que1QueU)$inDel <- "ins"

que2QueC <- subsetByOverlaps(que2Que.gr, que1Que.gr)
que2QueU <- que2Que.gr[!(que2Que.gr %in% que2QueC)]
mcols(que2QueU)$inDel <- "ins"


que1U.gr <- c(que1QueU, que1RefU)
que2U.gr <- c(que2QueU, que2RefU)


# adjust orientation

dfQue1 <- as.data.frame(mcols(que1U.gr)[c("queRange", "queOrientation", "queSeqlength")])
dfQue1[dfQue1$queOrientation == "-",c("queRange.start","queRange.end")] <- dfQue1$queSeqlength[dfQue1$queOrientation == "-"] - 
  dfQue1[dfQue1$queOrientation == "-",c("queRange.end","queRange.start")]

mcols(que1U.gr)$queRange <- GRanges(seqnames = Rle(dfQue1$queRange.seqnames),
                                    ranges = IRanges(start = dfQue1$queRange.start, 
                                                     end = dfQue1$queRange.end))


dfQue2 <- as.data.frame(mcols(que2U.gr)[c("queRange", "queOrientation", "queSeqlength")])
dfQue2[dfQue2$queOrientation == "-",c("queRange.start","queRange.end")] <- dfQue2$queSeqlength[dfQue2$queOrientation == "-"] - 
  dfQue2[dfQue2$queOrientation == "-",c("queRange.end","queRange.start")]

mcols(que2U.gr)$queRange <- GRanges(seqnames = Rle(dfQue2$queRange.seqnames),
                                    ranges = IRanges(start = dfQue2$queRange.start, 
                                                     end = dfQue2$queRange.end))


print("got indels")

# Remove gaps overlapping blocks in query
que1U.gr <- que1U.gr[ !overlapsAny(mcols(que1U.gr)$queRange, mcols(queSpecies$que1$Raw.gr)$queRange, minoverlap = 2) ]
mcols(que1U.gr)$rowNum <- 1:length(que1U.gr)
que2U.gr <- que2U.gr[ !overlapsAny(mcols(que2U.gr)$queRange, mcols(queSpecies$que2$Raw.gr)$queRange, minoverlap = 2) ]
mcols(que2U.gr)$rowNum <- 1:length(que2U.gr)

# remove both equal overlapping gaps and non-equal overlapping gaps until gaps are at depth 1
for(i in c("que1U.gr", "que2U.gr")){
  
  queU.gr <- get(i)
  
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
  assign(x = i,value = queU.gr)
  
}







dfQue1 <- as.data.frame(que1U.gr)

dfQue2 <- as.data.frame(que2U.gr)


print("df conversion")



#### Write out table

que1File <- paste(opt$outDir, "/",opt$queryName1,"_que.", opt$referenceName,"_ref.", "indel", sep = "")
write.table(x = dfQue1 ,file = que1File, sep = "\t", quote = FALSE, row.names = FALSE)

que2File <- paste(opt$outDir, "/",opt$queryName2,"_que.", opt$referenceName,"_ref.", "indel", sep = "")
write.table(x = dfQue2, file = que2File, sep = "\t", quote = FALSE, row.names = FALSE)


# 
# pdf(file = "plots/inDelIdentify/inDelSizeDist.pdf", onefile = TRUE)
# 
# 
# hist(log10(width(mcols(que2QueC)$queRange)), breaks = 50, main = "Overlapping gaps in que genomes", xlab = "gap size (log10 bp)")
# hist(log10(width(mcols(que1QueC)$queRange)), breaks = 50, col =2 , density = 0, add = TRUE)
# legend("topright", legend = c("mm9", "hg19"), fill = c(1,2), title = "Que geneomes")
# 
# 
# hist(log10(width(que2RefU)), breaks = 50, main = "Unique gaps in ref genome (Deletion)", xlab = "gap size (log10 bp)")
# hist(log10(width(que1RefU)), breaks = 50, col =2 , density = 0, add = TRUE)
# legend("topright", legend = c("mm9", "hg19"), fill = c(1,2), title = "Que geneomes")
# 
# 
# hist(log10(width(que2RefC)), breaks = 50, main = "Overlapping gaps in ref genome", xlab = "gap size (log10 bp)")
# hist(log10(width(que1RefC)), breaks = 50, col =2 , density = 0, add = TRUE)
# legend("topright", legend = c("mm9", "hg19"), fill = c(1,2), title = "Que geneomes")
# 
# hist(log10(width(mcols(que2QueU)$queRange)), breaks = 100, main = "Unique gaps in que genomes (insertions)", xlab = "gap size (log10 bp)")
# hist(log10(width(mcols(que1QueU)$queRange)), breaks = 50, col =2 , density = 0, add = TRUE)
# legend("topright", legend = c("mm9", "hg19"), fill = c(1,2), title = "Que geneomes")
# 
# 
# dev.off()
# 






