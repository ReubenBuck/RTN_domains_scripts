#!/usr/bin/env Rscript

## here we get our denomenator



rm(list = ls())


options(stringsAsFactors = FALSE)


library("optparse")


# input files 
# 


option_list = list(
  make_option(c("-d", "--inputDir"), type="character", default=NA, 
              help="directory containing only ref directories", metavar="character"),
  make_option(c("-n", "--queryName"), type="character", default=NA, 
              help="name of query genome", metavar="character"),
  make_option(c("-r", "--QueryAsRef"), type="character", default=NA, 
              help="specify 'TRUE' if query is in ref position", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default="./", 
              help="out directory", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if(any(is.na(opt))){
  stop("argument missing")
}


#opt$queryBlocks <- "~/Desktop/RTN_domains/data/chainAlignments/hg19/mm10/hg19.mm10.brokenChain"

#opt$queryName <- "hg19"
#opt$inputDir <- "Desktop/RTN_domains/data/chainAlignments/inDelTest/canFam2Hg19.out"
#opt$outDir <- "Desktop/RTN_domains/data/comparativeGenomics/supportedIndels/"



library(GenomicRanges)
library(RMySQL)
mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = opt$queryName)
chrInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")




genomeRefNames <- list.files(opt$inputDir)

proper=function(s) sub("(.)", ("\\U\\1"), tolower(s), pe=TRUE)
refSpecUnion.gr <- GRanges()
for(ref in genomeRefNames){
  fileName <- paste(opt$inputDir,"/",ref, sep = "")
  refSpec <- read.table(file = fileName,
                        col.names = c("refChr", "refLen", "refStart", "refEnd", "refStrand",
                                      "refGap", "queChr", "queLen", "queStart", "queEnd", "queStrand",
                                      "queGap", "chainID"),
                        colClasses = c("character", "numeric", "numeric", "numeric", "character", "numeric",
                                       "character", "numeric", "numeric", "numeric", "character", "numeric",
                                       "numeric")
  )
  
  # lets do the switch here
  refSpec[refSpec$queStrand == "-",c("queEnd", "queStart")] <- refSpec[refSpec$queStrand == "-", "queLen"] - refSpec[refSpec$queStrand == "-",c("queStart", "queEnd")]
  
  if(opt$QueryAsRef){
    refSpecRaw.gr <- GRanges(seqnames = Rle(refSpec$refChr), 
                             ranges = IRanges(start = refSpec$refStart, end = refSpec$refEnd)
    )
  }else if(!opt$QueryAsRef){
    refSpecRaw.gr <- GRanges(seqnames = Rle(refSpec$queChr), 
                           ranges = IRanges(start = refSpec$queStart, end = refSpec$queEnd)
    )
  }else{
    stop("specify if query is in ref position")
  }
  
  
  refSpecRaw.gr <- sortSeqlevels(refSpecRaw.gr)
  refSpec <- refSpec[order(refSpecRaw.gr),]
  refSpecRaw.gr <- sort(refSpecRaw.gr)
  
  seqlevels(refSpecRaw.gr) <- chrInfo$chrom
  seqlengths(refSpecRaw.gr) <- chrInfo$size 
  
  refSpecUnion.gr <- GenomicRanges::union(refSpecUnion.gr,reduce(refSpecRaw.gr))
}




denom.df <- as.data.frame(refSpecUnion.gr)

write.table(denom.df, file = paste(opt$outDir,"/",opt$queryName,".base",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

