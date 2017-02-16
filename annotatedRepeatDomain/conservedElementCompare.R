rm(list = ls())

setwd("~/Desktop/RTN_domains/")
 

library(GenomicRanges)
library(dplyr)
library(RMySQL)
library(devtools)
devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")



repGroups = c("ancient", "new_SINE", "new_L1", "old_L1")
repCols = c("darkblue", "aquamarine3", "purple", "red")
snames <- c(s1name = "hg19", s2name = "mm9", s3name = "canFam3")


ints <- read.table(file = "data/repeatHotspot/intersect.txt", header= TRUE)
### Break it down to other names
## now that we can look at all three
## use those names to get the corresponding regions

# use the R api to get the UCSC info


# get all the genome data.
genomeInfo <- c(rep(list(NA),length(snames)))
names(genomeInfo) <- snames

for(s in snames){
  
  mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = s)
  
  seqInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")
  
  ensGene <- dbGetQuery(mychannel, "SELECT * FROM ensGene;")
  ensGene.gr <- GRanges(seqnames = Rle(ensGene$chrom), 
                        ranges = IRanges(start = ensGene$txStart, end = ensGene$txEnd, names = ensGene$name),
                        strand = Rle(ensGene$strand), ensGene[,7:ncol(ensGene)])
  
  seqlevels(ensGene.gr) = seqInfo[,1]
  isCircular(ensGene.gr) = seqlevels(ensGene.gr) == "chrM"
  seqlengths(ensGene.gr) = seqInfo[,2]
  genome(ensGene.gr) = s
  
  exon.gr <- GRanges(seqnames = Rle(do.call(rep, list(x = ensGene$chrom, times = ensGene$exonCount))),
                     ranges = IRanges(start = unlist(strsplitAsListOfIntegerVectors(x = ensGene$exonStarts,sep = ",")),
                                      end = unlist(strsplitAsListOfIntegerVectors(x = ensGene$exonEnds,sep = ",")),
                                      names = do.call(rep, list(x = ensGene$name, times = ensGene$exonCount))),
                     strand = Rle(do.call(rep, list(x = ensGene$strand, times = ensGene$exonCount)))
  )
  seqlevels(exon.gr) = seqInfo[,1]
  seqinfo(exon.gr) <- seqinfo(ensGene.gr)
  exon.gr <- reduce(exon.gr)
    
  file <- list.files(paste("data/genomeAno/consElement/", s, sep = ""))
  ce <- read.table(paste("data/genomeAno/consElement/", s, "/", file, sep = ""))
  if(s == "canFam3"){
    colnames(ce) <- c("chr", "start", "end", "name", "score")
  }else{
    colnames(ce) <- c("bin","chr", "start", "end", "name", "score")
  }
  
  ce.gr <- GRanges(seqnames = Rle(ce$chr),
                   ranges = IRanges(start = ce$start, end = ce$end))
  seqlevels(ce.gr) <- seqlevels(ensGene.gr)
  seqinfo(ce.gr) <- seqinfo(ensGene.gr)
  
  olexon <- as.matrix(findOverlaps(ce.gr, exon.gr )) 
  nonCodingCe.gr <- ce.gr[-unique(olexon[,1])]
  
  genomeInfo[[s]] = list(ensGene = ensGene.gr,exons = exon.gr, ce = ce.gr ,nonCodingCe = nonCodingCe.gr)
  
}


# we get the non-coding info and exon info from all species

# then we can cycle through and compare results

genomeData <- NULL

for(name1 in snames){
  
  # get the other genome ready for analysis
  otherNames <- snames[snames != name1]
  otherGenomes <- c(rep(list(NULL), length(otherNames)))
  names(otherGenomes) <- otherNames
  for(name2 in 1:length(otherNames)){
    dat <- read.table(file = paste("data/repeatHotspot/",name1,"/",name1,"_",otherNames[name2],"_conDif.txt", sep = ""), header = TRUE)
    datRef.gr <- GRanges(seqnames = Rle(dat$chr[dat$genome == "ref"]),
                         ranges = IRanges(start = dat$start[dat$genome == "ref"], end = dat$end[dat$genome == "ref"]))
    mcols(datRef.gr) <- dat[dat$genome == "ref", 4:ncol(dat)]
    seqlevels(datRef.gr) <- seqlevels(genomeInfo[[name1]]$ensGene)
    seqinfo(datRef.gr) <- seqinfo(genomeInfo[[name1]]$ensGene)
    
    refExon.gr <- genomeInfo[[name1]]$exons
    mcols(datRef.gr)$knownNonExon = mcols(datRef.gr)$known - overlapingBases(dat.gr = datRef.gr, overlap.gr = refExon.gr)
    
    refCe.gr <- genomeInfo[[name1]]$nonCodingCe
    mcols(datRef.gr)$nonCodingCe = overlapingBases(dat.gr = datRef.gr, overlap.gr = refCe.gr)

    mcols(datRef.gr)$refGenome <- name1
    mcols(datRef.gr)$queGenome <- otherNames[name2]
    
    datQue.gr <- GRanges(seqnames = Rle(dat$chr[dat$genome == "que"]),
                         ranges = IRanges(start = dat$start[dat$genome == "que"], end = dat$end[dat$genome == "que"]))
    mcols(datQue.gr) <- dat[dat$genome == "que", 4:ncol(dat)]
    seqlevels(datQue.gr) <- seqlevels(genomeInfo[[otherNames[name2]]]$ensGene)
    seqinfo(datQue.gr) <- seqinfo(genomeInfo[[otherNames[name2]]]$ensGene)
    
    
    queExon.gr <- genomeInfo[[otherNames[name2]]]$exons
    mcols(datQue.gr)$knownNonExon = mcols(datQue.gr)$known - overlapingBases(dat.gr = datQue.gr, overlap.gr = queExon.gr)
    
    queCe.gr <- genomeInfo[[otherNames[name2]]]$nonCodingCe
    mcols(datQue.gr)$nonCodingCe = overlapingBases(dat.gr = datQue.gr, overlap.gr = queCe.gr)
    
    mcols(datQue.gr)$refGenome <- name1
    mcols(datQue.gr)$queGenome <- otherNames[name2]
    
    
    genomeData <- rbind(genomeData, rbind(as.data.frame(datRef.gr), as.data.frame(datQue.gr)))
    
  }
  
}

# queGenome == otherNames[name2]

# there seems to be a problem with the number of known nonExon



for(name1 in snames){
   otherNames <- snames[snames != name1]
  #  for(name2 in 1:length(otherNames)){
  #    for(conState in c("con", "dif")){
  pdf(file = paste("plots/hotspotOverlap/ceRate/", name1, "_", otherNames[1],"_",otherNames[2],".pdf", sep = ""), height = 12, width = 12)
  
  layout(matrix(1:4, nrow = 2, byrow = T))
  par(xaxs="i", mar = c(5,3,1,1), oma = c(5,5,5,5))
  for(rep in repGroups){
    
    refInts <- filter(ints, genome == name1 & repGroup == rep)
    newDat <- filter(genomeData, refGenome == name1 & repGroup == rep & conState != "mid" & hotspotID %in% refInts$domains) %>%
      group_by(genome, conState, refGenome, queGenome, nonCodingCe, knownNonExon , hotspotID, hotspotGroup, known) %>%
      summarise(nonCodingCeSum = sum(nonCodingCe))
    
    newDat$genome <- factor(newDat$genome, levels = c("ref", "que"))
    
    newDat$conState <- droplevels(newDat$conState)
    boxplot( log10((nonCodingCeSum + 1)/knownNonExon * 100 ) ~ conState + genome+ queGenome ,
             data = newDat, las = 2, notch = TRUE, xlim = c(.5, 8.5), ylim = c(-2,2))
    legend("topleft", legend = c(rep), bty = "n")
    abline(v = 4.5, lwd = 3); abline(v = 2.5, lwd = 3, lty =2); abline(v = 6.5, lwd = 3, lty = 2)
    title(main = name1,outer = TRUE)
  }
  mtext(side = 2, text = "% Non-exonic conserved elements (log10)", outer = TRUE)
  dev.off()
  
}


for(name1 in snames){
  otherNames <- snames[snames != name1]
  #  for(name2 in 1:length(otherNames)){
  #    for(conState in c("con", "dif")){
  pdf(file = paste("plots/hotspotOverlap/insertionRateInt/", name1, "_", otherNames[1],"_",otherNames[2],".pdf", sep = ""), height = 12, width = 12)
  
  layout(matrix(1:4, nrow = 2, byrow = T))
  par(xaxs="i", mar = c(5,3,1,1), oma = c(5,5,5,5))
  for(rep in repGroups){
    
    refInts <- filter(ints, genome == name1 & repGroup == rep)
    newDat <- filter(genomeData, refGenome == name1 & repGroup == rep & conState != "mid" & hotspotID %in% refInts$domains) %>%
      group_by(genome, conState, refGenome, queGenome, insertionRate, knownNonExon , hotspotID, hotspotGroup, known) %>%
      summarise(nonCodingCeSum = sum(insertionRate))
    
    newDat$genome <- factor(newDat$genome, levels = c("ref", "que"))
    
    newDat$conState <- droplevels(newDat$conState)
    boxplot( (((nonCodingCeSum + 1)/known)*50e3 ) ~ conState + genome+ queGenome ,
             data = newDat, las = 2, notch = TRUE, xlim = c(.5, 8.5))
    legend("topleft", legend = c(rep), bty = "n")
    abline(v = 4.5, lwd = 3); abline(v = 2.5, lwd = 3, lty =2); abline(v = 6.5, lwd = 3, lty = 2)
    title(main = paste(name1, "insertionRate"),outer = TRUE)
  }
  mtext(side = 2, text = "Elements per 50kb", outer = TRUE)
  
  dev.off()
  
}




for(name1 in snames){
  otherNames <- snames[snames != name1]
  #  for(name2 in 1:length(otherNames)){
  #    for(conState in c("con", "dif")){
  pdf(file = paste("plots/hotspotOverlap/insertionRateInt/", name1, "_", otherNames[1],"_",otherNames[2],"noInt.pdf", sep = ""), height = 12, width = 12)
  
  layout(matrix(1:4, nrow = 2, byrow = T))
  par(xaxs="i", mar = c(5,3,1,1), oma = c(5,5,5,5))
  for(rep in repGroups){
    
    newDat <- filter(genomeData, refGenome == name1 & repGroup == rep & conState != "mid") %>%
      group_by(genome, conState, refGenome, queGenome, insertionRate, knownNonExon , hotspotID, hotspotGroup, known) %>%
      summarise(nonCodingCeSum = sum(insertionRate))
    
    newDat$genome <- factor(newDat$genome, levels = c("ref", "que"))
    
    newDat$conState <- droplevels(newDat$conState)
    boxplot( (((nonCodingCeSum + 1)/known)*50e3 ) ~ conState + genome+ queGenome ,
             data = newDat, las = 2, notch = TRUE, xlim = c(.5, 8.5))
    legend("topleft", legend = c(rep), bty = "n")
    abline(v = 4.5, lwd = 3); abline(v = 2.5, lwd = 3, lty =2); abline(v = 6.5, lwd = 3, lty = 2)
    title(main = paste(name1, "insertionRate"),outer = TRUE)
  }
  mtext(side = 2, text = "Elements per 50kb", outer = TRUE)
  dev.off()
  
}

# so far this is not getting the intersect. 





# so now I have dog coordinates


# not really seeing any differences between conserved and non conserved

# the real question would be what are genes doing.

# how to find the impact 

# maybe it is better to look at individual phylop scores


# There is a blind spot in our analysis. 
# it's what is happening in the other species.

# in the referecne everything is pretty similar


# there's a couple of interesting small indicators. 

# lets get these plots looking pretty

  
eLength <- genomeData[genomeData$knownNonExon < 0,]
  
e <- genomeInfo$mm9$exons[seqnames(genomeInfo$mm9$exons) == "chr17" & start(genomeInfo$mm9$exons) > 35050001 & end(genomeInfo$mm9$exons) < 35100000]

gd <- (genomeData[genomeData$knownNonExon < 0,])[2:6,]

50000 - sum(width(e)) * 5


# becasue multiple ref bins point to one query bin. 

genomeData[genomeData$hotspotID %in% gd$hotspotID & genomeData$genome == "ref" & genomeData$refGenome == "hg19" & genomeData$queGenome == "mm9",]




genomeData[genomeData$hotspotID == "new_SINE;4571" & genomeData$genome == "que" & genomeData$refGenome == "hg19" & genomeData$queGenome == "mm9",]


genomeData[genomeData$hotspotID == "new_SINE;1077"  & genomeData$refGenome == "hg19",]


# seems like the target missed.
# at some point during our overlap where we tried to assing bins , we missed. 





