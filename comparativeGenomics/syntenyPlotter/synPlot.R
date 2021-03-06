#!/usr/bin/env Rscript

# plots syntenic relationships between two species
# reads in seqGaps and chr lengths from ucsc
# there is also an option to include coordinates from another bedfile to plot above chromosomes
# bedfile must be three columns.


rm(list = ls())


library(optparse)


option_list = list(
  make_option(c("-r", "--ref"), type="character", default=NA, 
              help="reference genome", metavar="character"),
  make_option(c("-q", "--que"), type="character", default=NA, 
              help="query genome", metavar="character"),
  make_option(c("-s", "--synFile"), type="character", default=NA, 
              help="output from synBuildConvert.bash", metavar="character"),
  make_option(c("-b", "--bedFile1"), type="character", default="none", 
              help="bed file to plot above chromosomes in black, if no bed file use 'none'", metavar="character"),
  make_option(c("-B", "--bedFile2"), type="character", default="none", 
              help="bed file to plot above chromosomes in grey, if no bed file use 'none'", metavar="character"),
  make_option(c("-o", "--outPutPDF"), type="character", default="./out.pdf", 
              help="pdf output figure", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if(any(is.na(opt))){
  stop("argument missing")
}

# check output has correct extension, if not automaticly add one
ext <- rev(strsplit(x = opt$outPutPDF,split = "\\.")[[1]])[1]
if(ext != "pdf"){
  opt$outPutPDF <- paste(opt$outPutPDF,".pdf", sep = "")
}


library(RMySQL)
library(ggsci)

opt$ref = "hg19"
opt$que = "mm10"
opt$synFile = "~/Desktop/RTN_domains/data/synBuilder/hg19.mm10.highResSyn.txt"
opt$bedFile1="~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/missingData/hg19nonNet.bed"
opt$bedFile2="~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/missingData/hg19nonRBH.bed"
opt$outPutPDF = "~/Desktop/hg19syn.pdf"


ref = opt$ref
que = opt$que
synFile = opt$synFile


# get Ref chromosome data
# only placed chromosomes
# order according to name, will return warnings
mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = ref)
refChr <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")
refChr <- refChr[-(grep(pattern = "_", x = refChr$chrom)),]
refChr$chrom <- factor(refChr$chrom, 
                       levels =   c(
                         paste("chr",sort(as.integer(substring(refChr$chrom, 4))), sep = ""),
                         "chrX", "chrY", "chrM")
                       )

# repeat for query
mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = que)
queChr <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")
queChr <- queChr[-(grep(pattern = "_", x = queChr$chrom)),]
queChr$chrom <- factor(queChr$chrom, 
                       levels =   c(
                         paste("chr",sort(as.integer(substring(queChr$chrom, 4))), sep = ""),
                         "chrX", "chrY", "chrM")
)

# get syntenty data
synData <- read.table(synFile, header = FALSE, 
                      col.names = c("ref", "refChr", "refStart", "refEnd", "refStrand",
                                    "que", "queChr", "queStart", "queEnd", "queStrand"))
synData$refChr <- factor(synData$refChr, levels = levels(refChr$chrom))
synData$queChr <- factor(synData$queChr, levels = levels(queChr$chrom))

# get sequencing gaps in ref
# will retrun warnings
mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = ref)
seqGaps <- dbGetQuery(mychannel, "SELECT * FROM gap;")
seqGaps <- seqGaps[-(grep(pattern = "_", x = seqGaps$chrom)),]
seqGaps$chrom <- factor(seqGaps$chrom, levels = levels(refChr$chrom))


# get bed data
if(opt$bedFile1 != "none"){
bed1 = read.table(opt$bedFile1, header = FALSE)
  bed1 <- bed1[,1:3]
  colnames(bed1) <- c("chr", "start", "end")
  bed1$chr = factor(bed1$chr , levels = levels(refChr$chrom))
}

if(opt$bedFile2 != "none"){
  bed2 = read.table(opt$bedFile2, header = FALSE)
  bed2 <- bed2[,1:3]
  colnames(bed2) <- c("chr", "start", "end")
  bed2$chr = factor(bed2$chr , levels = levels(refChr$chrom))
}

# generate synteny plot
pdf(file = opt$outPutPDF, height = 9,width = 6)

#pdf("~/Desktop/gapsMouse.pdf", height = 9,width = 6)
plot.new()
plot.window(ylim = c(nrow(refChr)*3 - 3,2), xlim = c(-1,2.5e8), xaxs = "i")

grid()

# generate query colours
colFun <- pal_ucscgb()
queCols <- colFun(nrow(queChr))

# plto syntenty blocks
rect(ytop = (as.integer(synData$refChr) * 3) -1  ,
     ybottom  = (as.integer(synData$refChr)*3) -2, 
     xleft = synData$refStart, 
     xright = synData$refEnd,
     col = queCols[synData$queChr],
     border = NA
     )

# plot gaps
rect(ytop = (as.integer(seqGaps$chrom) * 3) -1  ,
     ybottom  = (as.integer(seqGaps$chrom)*3) -2, 
     xleft = seqGaps$chromStart, 
     xright = seqGaps$chromEnd,
     col = "black",
     border = NA
)

# plot bed data
if(opt$bedFile1 != "none"){
  bed1 <- bed1[bed1$end - bed1$start > 9999,]
  rect(ytop = (as.integer(bed1$chr) * 3) -2  ,
       ybottom  = (as.integer(bed1$chr)*3) -2.5, 
       xleft = bed1$start, 
       xright = bed1$end,
       col = "black",
       border = NA
  )
}

if(opt$bedFile2 != "none"){
  bed2 <- bed2[bed2$end - bed2$start > 9999,]
  rect(ytop = (as.integer(bed2$chr) * 3) -1  ,
       ybottom  = (as.integer(bed2$chr)*3) -.5, 
       xleft = bed2$start, 
       xright = bed2$end,
       col = "red",
       border = NA
  )
}

# plot chromosomes
rect(ytop = (as.integer(refChr$chrom) * 3) -1  ,
     ybottom  = (as.integer(refChr$chrom)*3) -2, 
     xleft = 0, 
     xright = refChr$size,
     lwd = .5)

# plot relative query position line
mer <- merge(x = synData, queChr[,c("chrom", "size")], by.x = "queChr", by.y = "chrom")
for(i in 1:nrow(mer)){
  x = c(mer$refStart[i], mer$refEnd[i])
  y = c(mer$queStart[i]/mer$size[i], mer$queEnd[i]/mer$size[i])
  y = (1 - y) + (as.integer(mer$refChr)[i] * 3) -2
  if(mer$queStrand[i] == "+"){
    x = rev(x)
  }
  lines(x, y)
}

mtext(side = 2, at = (as.integer(refChr$chrom) * 3) -1.5, text = refChr$chrom, las = 2, line = 1)
legend("bottomright", legend = levels(queChr$chrom), fill = queCols, title = que,bty = "n")
axis(side = 1, line = .1)
title(main = ref, xlab = "position (bp)")

dev.off()






