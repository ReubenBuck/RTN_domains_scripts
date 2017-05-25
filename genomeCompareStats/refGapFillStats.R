# collect statistics for our gap finding process and ancestral genome analysis

rm(list = ls())

options(stringsAsFactors = FALSE)

library(GenomicRanges)
library(RMySQL)


## read functions
rmSeqGapsFromNetOutput <- function(netOutput, seqGaps){
  # removes seq gaps from a net output file by taking the setdiff
  sDiff <- GenomicRanges::setdiff(netOutput, seqGaps)
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
  
  seqGaps.gr <- sort(sortSeqlevels(seqGaps.gr))
  
  mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = genomes[genomes != genomes[i]])
  altChrInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")
  
  
  # get gap and fill data for ref and que
  
  for(j in 1:2){
    netOutputFile <- get(paste(c("gap", "fill")[j],"Files", sep = ""))
    
    netOutput <- read.table(file = netOutputFile[i], header = FALSE,
                            col.names = c("seqnames", "start", "end","que.seqnames","que.start", "que.end","strand","chainID"),
                            colClasses = c("character", "integer", "integer","character", "integer", "integer","character","integer"))
    # convert to GRanges, move from 0 based to 1 based
    netOutput.gr <- GRanges(seqnames = netOutput$seqnames,
                            ranges = IRanges(start = netOutput$start + 1,
                                             end = netOutput$end + 1),
                            queRanges = GRanges(seqnames = netOutput$que.seqnames,
                                                ranges = IRanges(start = netOutput$que.start + 1,
                                                                 end = netOutput$que.end + 1)
                            ),
                            strand = "*",
                            sData = netOutput$strand,
                            chainID = netOutput$chainID
    )
    seqlevels(netOutput.gr) <- chrInfo$chrom
    seqlengths(netOutput.gr) <- chrInfo$size + 1
    genome(netOutput.gr) <- genomes[i]
    
    seqlevels(netOutput.gr$queRanges) <- altChrInfo$chrom
    seqlengths(netOutput.gr$queRanges) <- altChrInfo$size + 1
    genome(netOutput.gr$queRanges) <- genomes[genomes != genomes[i]]
    
    netOutput.gr <- rmSeqGapsFromNetOutput(netOutput = netOutput.gr, seqGaps = seqGaps.gr)
    
    netOutput.gr <- sort(sortSeqlevels(netOutput.gr))
    netOutput.gr$queRanges <- sortSeqlevels(netOutput.gr$queRanges)
    
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
  
  ancDna.gr <- sort(sortSeqlevels(ancDna.gr))
  
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


refFillGaps.gr <- gaps(refFill.gr)
refFillGaps.gr <- refFillGaps.gr[strand(refFillGaps.gr) == "*"]
refFillGaps.gr <- rmSeqGapsFromNetOutput(netOutput = refFillGaps.gr, seqGaps = refSeqGaps.gr)


queFillGaps.gr <- gaps(queFill.gr)
queFillGaps.gr <- queFillGaps.gr[strand(queFillGaps.gr) == "*"]
queFillGaps.gr <- rmSeqGapsFromNetOutput(netOutput = queFillGaps.gr, seqGaps = queSeqGaps.gr)


# save RObjects

save("genomes","queChrInfo", "queAncDna.gr", "queFill.gr", "queFillGaps.gr", "queGap.gr", 
     "queSeqGaps.gr", "refChrInfo" ,"refAncDna.gr", "refFill.gr", "refFillGaps.gr", 
     "refGap.gr", "refSeqGaps.gr", 
     file = paste("~/Desktop/RTN_domains/R_objects/mappedGaps/", genomes[1],"." ,genomes[2], ".netData.RData", sep = ""))




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

queRemain  = sum(width(queFill.gr))
refReamin = sum(width(refFill.gr))

abs(queSize - refSize) 
abs(abs(queGain - refGain) - abs(refLoss - queLoss))


# so of course this works because we are acounting for all the bases
# so if we miss asigned in one species, will this lead to an incorrect level in the other speceis ?
as.numeric(queRemain) + as.numeric(queGain) + as.numeric(refLoss) - as.numeric(refReamin) - refGain - as.numeric(queLoss)

# so it adds up becasue all of these factors are different partitions of the same genome
# so is there some way to test if our partioning is actually effective?


# Two alternaitves, assembled chromosomes only, or consider the lineage specific TE contribution to identifying new DNA

# also think about the numbers that will lead to a smaller estimate and the numbers that will lead to a larger


## How much of our mouse human ancestral DNA can be seen as ancestral to other speceis
sum(width(intersect(refFill.gr,refAncDna.gr))) / sum(width(refFill.gr))
sum(width(intersect(queFill.gr,queAncDna.gr))) / sum(width(queFill.gr))


# what can we expect to explain based on our fills and thier level of ancestry. 
sum(width(intersect(refGap.gr, refAncDna.gr)))/(sum(width(refGap.gr)))
sum(width(intersect(queGap.gr, queAncDna.gr)))/(sum(width(queGap.gr)))


# genome size estimates
sum(as.numeric(seqlengths(refSeqGaps.gr))) - sum(width(refSeqGaps.gr))
sum(as.numeric(seqlengths(queSeqGaps.gr))) - sum(width(queSeqGaps.gr))


# remaining ancestral genome size
sum(width(refFill.gr))
sum(width(queFill.gr))
sum(width(refFill.gr)) - sum(width(queFill.gr))
# there is a 28 megabase difference in our remaining ancestral genomes
# maybe there has been more recent SDs in human than mouse,

# is there a way yo anlyse this
# we could take fills remapped and remove gaps and look for overlaps
# overlapping bases suggests there has been an SD event.

# these can be purged from the genome too





# from here we migh be able to explain loss and gain and change in genome size.
# this can give us a margin of error to consider


### lets just do this equation 
# Mloss = Hun + Mun - Hloss - HMdiff

# maybe we can factor in our estimates to get an upper and lower bound on these things
# based on the 28 megabase difference

# we can also ask what happens if we change our estimates by 5%
# we know a gain in one is equall to a loss in the other

HMdiff <- (queSize - refSize)
Hun = sum(width(queGap.gr))
Mun = sum(width(refGap.gr))

HLoss = seq(0,Mun,length.out = 20)

MLoss = (as.numeric(Hun) + as.numeric(Mun) - (2*HLoss) - HMdiff + as.numeric(sum(width(queFill.gr)) - sum(width(refFill.gr))))/2

# this is basicly predicting gain 

Hra <- 


mmAKloss = 1219
mmAKgain = 1007 + 66

hgAKloss = 650
hgAKgain = 815+84

mmAKinferedLoss <- (Hun/1e6) - hgAKgain
mmAKinferedGain <- (Mun/1e6) - hgAKloss

hgAKinferedLoss <- (Mun/1e6) - mmAKgain
hgAKinferedGain <- (Hun/1e6) - mmAKloss


((hgAKgain - mmAKgain) - (hgAKloss - mmAKloss)) * 1e6

# plot(HLoss/1e6, MLoss/1e6, type = "b", xlab = "Human Loss", ylab = "Mosue Loss", 
#      col = terrain.colors(30), xlim =c(0, 1800), ylim =c(0, 1800))
# points(queLoss/1e6, refLoss/1e6)
# plot((Hun - MLoss) /1e6, (Mun - HLoss)/1e6, type = "b", xlab = "Human Gain", 
#      ylab = "Mosue Gain", col = terrain.colors(30), xlim =c(0, 1800), ylim =c(0, 1800))
# points(queGain/1e6,refGain/1e6)

layout(c(1,2))
plot((Hun - MLoss) /1e6, HLoss/1e6, type = "b", xlab = "Gain (Mb)", ylab = "Loss (Mb)", col = terrain.colors(30), 
     pch = 16 , xlim =c(0, 1800), ylim =c(0, 1800), main = genomes[2])
points(queGain/1e6, queLoss/1e6, pch = 16)
points(hgAKgain, hgAKloss, col = 2, pch = 16)
#points(hgAKinferedGain, hgAKloss, col = 4, pch = 16)


plot((Mun - HLoss)/1e6, MLoss/1e6, type = "b", xlab = "Gain (Mb)", ylab = "Loss (Mb)", col = terrain.colors(30), 
     pch = 17, xlim =c(0, 1800), ylim =c(0, 1800), main = genomes[1])
points(refGain/1e6, refLoss/1e6, pch = 17)
points(mmAKgain, mmAKloss, col = 2, pch = 17)
#points(mmAKinferedGain, mmAKloss, col = 4, pch = 16)



# is there some way that our independant partioning can miss something?
# so far it might be that by using the left over genome as gain, we might always get the same numbers
# One thing is that we can put limits on our results, then we can know when those limits are exceeded.

# lets plot our partioning
HGain <- (Hun - MLoss)
MGain <- (Mun - HLoss)

layout(c(1,2))
plot(HGain/1e6, MLoss/1e6, type = "b", xlab = "Gain (Mb)", ylab = "Loss (Mb)", col = terrain.colors(30), 
     pch = 16 , xlim =c(0, 1800), ylim =c(0, 1800), main = genomes[2])
points(queGain/1e6, refLoss/1e6, pch = 16)
points(hgAKgain, mmAKloss, col = 2, pch = 16)
#points(hgAKinferedGain, hgAKloss, col = 4, pch = 16)


plot(MGain/1e6, HLoss/1e6, type = "b", xlab = "Gain (Mb)", ylab = "Loss (Mb)", col = terrain.colors(30), 
     pch = 17, xlim =c(0, 1800), ylim =c(0, 1800), main = genomes[1])
points(refGain/1e6, queLoss/1e6, pch = 17)
points(mmAKgain, hgAKloss, col = 2, pch = 17)

# because they care about total, they are off the line

# however, I'm not yet confinced our estimate is only an error becasue of the different begining amount



# the idea is that the total loss and gain along each lineage should be able to explain the discrepencies in genome size

# our estimates are consistent with a smaller ancestral genome between the two speceis, 
# however it is a bit hard to see how they wouldn't be. 

# based on the ancestry we discovered, is it possible to see what the ancestral size must have been 

# eutherians estimate a genome size of 2.8 GB
# THe ancestral genome is about 2.3 - 2.4 Gb
sum(width(refFill.gr)) + as.numeric(queLoss + refLoss)
sum(width(queFill.gr)) + as.numeric(queLoss + refLoss)



# 






# think about the erros of the estimate



# if we use the lineage specific TE estimates between human and mouse? 

mmAKgain = 1007 + 66
mmAKloss = 1219

hgAKgain = 815+84
hgAkloss = 650

points(mmAKgain, mmAKloss, col = 2, pch = 17)

points(mmAKgain, mmAKloss, col = 2, pch = 17)











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

### we need to assign the bases and put them together








