# collect statistics for our gap finding process and ancestral genome analysis

rm(list = ls())

options(stringsAsFactors = FALSE)

library(GenomicRanges)
library(RMySQL)


## read functions
rmSeqGapsFromNetOutput <- function(netOutput, seqGaps){
  # removes seq gaps from a net output file by taking the setdiff
  sDiff <- setdiff(netOutput, seqGaps)
  if(!all(overlapsAny(sDiff,netOutput))){
    stop("some netOutput ranges have been completely removed, check if net output file is correct")
  }
  ol <- findOverlaps(sDiff, netOutput)
  # this intersection makes sure no regions have been unnecesarily joined
  sDiff <- pintersect(sDiff[queryHits(ol)], netOutput[subjectHits(ol)])
  ol <- findOverlaps(sDiff, netOutput)
  if(!all(sDiff == sDiff[queryHits(ol)])){
    stop("not a single layer net output and orger has not been maintained, cannot correctly assign mcols")
  }
  mcols(sDiff) <- mcols(netOutput[subjectHits(ol)])
  return(sDiff)
}


## specify files

# genomes
genomes <- c(ref = "mm10",que = "hg19")

#gap file
gapFiles <- c(ref = "~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/processed/mm10.hg19.net.gaps",
              que = "~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/processed/hg19.mm10.net.gaps")


# fill file
fillFiles <- c(ref = "~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/processed/mm10.hg19.net.fills",
               que = "~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/processed/hg19.mm10.net.fills")

# ancestral DNA file
ancDNAfiles <- c(ref = "~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/ancestralGenome/mm10.merge.ancestral.bed",
                 que = "~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/ancestralGenome/hg19.merge.ancestral.bed")


# read in info
for(i in 1:length(genomes)){
  # get chr info
  mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = genomes[i])
  chrInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")
  #chrInfo <- chrInfo[-(grep(pattern = "_", x = chrInfo$chrom)),]
  
  # get seq gaps
  seqGaps <- dbGetQuery(mychannel, "SELECT * FROM gap;")
  #seqGaps <- seqGaps[-(grep(pattern = "_", x = seqGaps$chrom)),]
  seqGaps.gr <- GRanges(seqnames = Rle(seqGaps$chrom), 
                        ranges = IRanges(start = seqGaps$chromStart + 1, 
                                         end = seqGaps$chromEnd + 1)
  )
  seqlevels(seqGaps.gr) <- chrInfo$chrom
  seqlengths(seqGaps.gr) <- chrInfo$size + 1
  genome(seqGaps.gr) <- genomes[i]
  
  
  # get gap and fill data for ref and que
  
  for(j in 1:2){
    netOutputFile <- get(paste(c("gap", "fill")[j],"Files", sep = ""))
    
    netOutput <- read.table(file = netOutputFile[i], header = FALSE,
                            col.names = c("seqnames", "start", "end","que.seqnames","que.start", "que.end","strand"),
                            colClasses = c("character", "integer", "integer","character", "integer", "integer","character"))
    # convert to GRanges, move from 0 based to 1 based
    netOutput.gr <- GRanges(seqnames = netOutput$seqnames,
                            ranges = IRanges(start = netOutput$start + 1,
                                             end = netOutput$end + 1),
                            queRanges = GRanges(seqnames = netOutput$que.seqnames,
                                                ranges = IRanges(start = netOutput$que.start + 1,
                                                                 end = netOutput$que.end + 1)
                            )
    )
    seqlevels(netOutput.gr) <- chrInfo$chrom
    seqlengths(netOutput.gr) <- chrInfo$size + 1
    genome(netOutput.gr) <- genomes[i]
    
    assign(x = paste(names(genomes)[i],c("Gap","Fill")[j],".gr", sep = ""), value = netOutput.gr)
    
  }
  
  # get ancestral
  ancDna <- read.table(file = ancDNAfiles[i], header = FALSE,
                       col.names = c("seqnames", "start", "end"),
                       colClasses = c("character", "integer", "integer"))
  # convert to GRanges, move from 0 based to 1 based
  ancDna.gr <- GRanges(seqnames = ancDna$seqnames,
                       ranges = IRanges(start = ancDna$start  + 1,
                                        end = ancDna$end + 1)
  )
  seqlevels(ancDna.gr) <- chrInfo$chrom
  seqlengths(ancDna.gr) <- chrInfo$size + 1
  genome(ancDna.gr) <- genomes[i]
  
  
  assign(x = paste(names(genomes)[i],"AncDna.gr", sep = ""), value = ancDna.gr)
  assign(x = paste(names(genomes)[i],"ChrInfo", sep = ""), value = chrInfo)
  assign(x = paste(names(genomes)[i],"SeqGaps.gr", sep = ""), value = seqGaps.gr)
  
  rm(chrInfo, seqGaps.gr, seqGaps, netOutput, netOutput.gr, ancDna, ancDna.gr)
  
}

all_cons <- dbListConnections(MySQL())
for(con in all_cons) {
  dbDisconnect(con)
}

# have read in all the data of both query and species




refGap.gr <- rmSeqGapsFromNetOutput(netOutput = refGap.gr, seqGaps = refSeqGaps.gr)

refFill.gr <- rmSeqGapsFromNetOutput(netOutput = refFill.gr, seqGaps = refSeqGaps.gr)

refFillGaps.gr <- gaps(refFill.gr)
refFillGaps.gr <- refFillGaps.gr[strand(refFillGaps.gr) == "*"]
refFillGaps.gr <- rmSeqGapsFromNetOutput(netOutput = refFillGaps.gr, seqGaps = refSeqGaps.gr)

queGap.gr <- rmSeqGapsFromNetOutput(netOutput = queGap.gr, seqGaps = queSeqGaps.gr)

queFill.gr <- rmSeqGapsFromNetOutput(netOutput = queFill.gr, seqGaps = queSeqGaps.gr)

queFillGaps.gr <- gaps(queFill.gr)
queFillGaps.gr <- queFillGaps.gr[strand(queFillGaps.gr) == "*"]
queFillGaps.gr <- rmSeqGapsFromNetOutput(netOutput = queFillGaps.gr, seqGaps = queSeqGaps.gr)



# this will tell us the proportion of all gaps we can place
sum(width(refGap.gr))/sum(width(refFillGaps.gr))
sum(width(queGap.gr))/sum(width(queFillGaps.gr))


# lets see if we can explain the genome size discrepensies with our loss gain dynamics
# This should indicate if there are any bases unacounted for

refLoss <- sum(width(intersect(queFillGaps.gr, queAncDna.gr)))
queGain <- sum(width(queFillGaps.gr)) - refLoss

queLoss <- sum(width(intersect(refFillGaps.gr, refAncDna.gr)))
refGain <- sum(width(refFillGaps.gr)) - queLoss

refSize <- sum(as.numeric(seqlengths(refSeqGaps.gr))) - sum(width(refSeqGaps.gr))
queSize <- sum(as.numeric(seqlengths(queSeqGaps.gr))) - sum(width(queSeqGaps.gr))

((refGain - queGain) - (refLoss - queLoss)) / (refSize - queSize)


# what proportion of the gain and loss dynamics can be explined by our assignment of gained and lossed DNA


refLoss <- sum(width(intersect(queGap.gr, queAncDna.gr)))
queGain <- sum(width(queGap.gr)) - refLoss

queLoss <- sum(width(intersect(refGap.gr, refAncDna.gr)))
refGain <- sum(width(refGap.gr)) - queLoss

refSize <- sum(as.numeric(width(refGap.gr))) + sum(width(refFill.gr))
queSize <- sum(as.numeric(width(queGap.gr))) + sum(width(queFill.gr))

((queGain - refGain) - (queLoss - refLoss)) / (queSize - refSize)

# Two alternaitves, assembled chromosomes only, or consider the lineage specific TE contribution to identifying new DNA

# also think about the numbers that will lead to a smaller estimate and the numbers that will lead to a larger


## How much of our mouse human ancestral DNA can be seen as ancestral to other speceis
sum(width(intersect(refFill.gr,refAncDna.gr))) / sum(width(refFill.gr))
sum(width(intersect(queFill.gr,queAncDna.gr))) / sum(width(queFill.gr))

# these results suggest that more ancestral DNA from other species is missing in mouse than in human


# what if we look at just the mappable genome.

# what can we expect to explain based on our fills and thier level of ancestry. 




sum(width(intersect(refFill.gr, refAncDna.gr)))/(sum(width(refFill.gr)))

sum(width(intersect(refGap.gr, refAncDna.gr)))/(sum(width(refGap.gr)))


sum(as.numeric(seqlengths(refSeqGaps.gr))) - sum(width(refSeqGaps.gr))


sum(as.numeric(seqlengths(queSeqGaps.gr))) - sum(width(queSeqGaps.gr))

# this gets us total length of what we are interested in for the mouse
as.numeric(sum(width(refFillGaps.gr))) + as.numeric(sum(width(refFill.gr)))

# from here we migh be able to explain loss and gain and change in genome size.
# this can give us a margin of error to consider


### lets just do this equation 
# Mloss = Hun + Mun - Hloss - HMdiff

HMdiff <- (queSize - refSize)
Hun = sum(width(queGap.gr))
Mun = sum(width(refGap.gr))

HLoss = seq(0,Mun,length.out = 20)

MLoss = (as.numeric(Hun) + as.numeric(Mun) - (2*HLoss) - HMdiff)/2


plot(HLoss/1e6, MLoss/1e6, type = "b", xlab = "Human Loss", ylab = "Mosue Loss", 
     col = terrain.colors(30), xlim =c(0, 1800), ylim =c(0, 1800))
points(queLoss/1e6, refLoss/1e6)

plot((Hun - MLoss) /1e6, (Mun - HLoss)/1e6, type = "b", xlab = "Human Gain", 
     ylab = "Mosue Gain", col = terrain.colors(30), xlim =c(0, 1800), ylim =c(0, 1800))
points(queGain/1e6,refGain/1e6)

layout(c(1,2))
plot((Hun - MLoss) /1e6, HLoss/1e6, type = "b", xlab = "Gain", ylab = "Loss", col = terrain.colors(30), 
     pch = 16 , xlim =c(0, 1800), ylim =c(0, 1800), main = genomes[2])
points(queGain/1e6, queLoss/1e6, pch = 16)

plot((Mun - HLoss)/1e6, MLoss/1e6, type = "b", xlab = "Gain", ylab = "Loss", col = terrain.colors(30), 
     pch = 17, xlim =c(0, 1800), ylim =c(0, 1800), main = genomes[1])
points(refGain/1e6, refLoss/1e6, pch = 17)

# our estimates are consistent with a smaller ancestral genome between the two speceis, 
# however it is a bit hard to see how they wouldn't be. 

# based on the ancestry we discovered, is it possible to see what the ancestral size must have been 


# THe ancestral genome is about 2.3 - 2.4 Gb
sum(width(refFill.gr)) + as.numeric(queLoss + refLoss)

sum(width(queFill.gr)) + as.numeric(queLoss + refLoss)





sum(width(queGap.gr))



agg <- aggregate(width(refGap.gr), by = list(as.character(mcols(refGap.gr)$queRanges)), FUN = sum)

queAggGap.gr <- GRanges(as.character(agg$Group.1))



ol <- findOverlaps( queFill.gr, queAggGap.gr, type = "within")
int <- pintersect(queAggGap.gr[subjectHits(ol)], queFill.gr[queryHits(ol)])

y = width(queAggGap.gr)
y[subjectHits(ol)] <- y[subjectHits(ol)] - width(int)

# the idea here is we are comparing gap element sizes with the corresponding fill sizes removed 

# if we remove 
smoothScatter(log10(agg$x)[1:100000], log10(y)[1:100000], ylim = c(0,8), xlim = c(0,8) )
abline(a = 0, b = 1)




