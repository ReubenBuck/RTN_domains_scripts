
### hg19Script

rm(list = ls())

options(stringsAsFactors = FALSE)



library(dplyr)
library(reshape)
library(GenomicRanges)
devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")



# lets get GC content
specRef = "hg19"
specQue = "mm10"


# binned Genome
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/",specRef,".synthBin.RData",sep = ""))
# fills and gaps
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/",specRef,".",specQue,".netData.RData",sep = ""))
#load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/",specRef,".stretch.RData", sep = ""))
# shift levels
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/shiftData/",specRef,".expand.breaks.RData", sep = ""))



dnase <- read.table("~/Desktop/RTN_domains/data/UCSCtracks/hg19/wgEncodeAwgDnaseMasterSites.txt",
                    col.names = c("bin", "seqnames", "start", "end", "name", 
                                  "dispScore", "floatScore", "sourceCount", "sourceID"),
                    colClasses = c("integer", "character", "integer", "integer", 
                                   "integer", "integer", "double", "integer", "character"))
dnase$start = dnase$start + 1
dnase.gr <- GRanges(dnase)
seqlevels(dnase.gr) <- refChrInfo$chrom
seqlengths(dnase.gr) <- refChrInfo$size
genome(dnase.gr) <- genomes["ref"]

dnase.gr <- genoExpandStretch(dnase.gr, newSynthRefShift, seqlengths(synthBin.gr))


ol <- findOverlaps(synthBin.gr, dnase.gr)

dnaseSum <- summarise(group_by(data.frame(queryHits = queryHits(ol), 
                                          peaks = dnase.gr$sourceCount[subjectHits(ol)]), 
                               queryHits), 
                      peakNumbers = sum(peaks), peakSites = n())

synthBin.gr$dnasePeaks <- NA
synthBin.gr$dnasePeaks[dnaseSum$queryHits] <- dnaseSum$peakSites 

synthBin.gr$dnaseActivity <- NA
synthBin.gr$dnaseActivity[dnaseSum$queryHits] <- dnaseSum$peakNumbers 



#### Exons and Introns


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
exon.gr <- reduce(exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
intron.gr <- GenomicRanges::setdiff(unlist(intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)), exons(TxDb.Hsapiens.UCSC.hg19.knownGene))


exon.gr <- genoExpandBreak(exon.gr, newSynthRefShift, seqlengths(synthBin.gr))


ol <- findOverlaps(synthBin.gr, exon.gr)
pInt <- pintersect(synthBin.gr[queryHits(ol)], exon.gr[subjectHits(ol)])

exonSum <- summarise(group_by(data.frame(queryHits = queryHits(ol), 
                                          exon = width(pInt)), 
                               queryHits), 
                        exon = sum(exon))

synthBin.gr$exon <- NA
synthBin.gr$exon[exonSum$queryHits] <- exonSum$exon


intron.gr <- genoExpandBreak(intron.gr, newSynthRefShift, seqlengths(synthBin.gr))


ol <- findOverlaps(synthBin.gr, intron.gr)
pInt <- pintersect(synthBin.gr[queryHits(ol)], intron.gr[subjectHits(ol)])

intronSum <- summarise(group_by(data.frame(queryHits = queryHits(ol), 
                                            intron = width(pInt)), 
                                 queryHits), 
                        intron = sum(intron))

synthBin.gr$intron <- NA
synthBin.gr$intron[intronSum$queryHits] <- intronSum$intron



# laminaAssocaited domains 

lads <- read.table("~/Desktop/RTN_domains/data/UCSCtracks/hg19/laminB1Lads", 
                      col.names = c("bin", "seqnames", "start", "end"))
lads.gr <- GRanges(lads)
# from UCSC db need to convert to 1 base
start(lads.gr) <- start(lads.gr) + 1
seqlevels(lads.gr) <- refChrInfo$chrom
seqlengths(lads.gr) <- refChrInfo$size
genome(lads.gr) <- genomes["ref"]
lads.gr <- sort(sortSeqlevels(lads.gr))


lads.gr <- genoExpandBreak(lads.gr, newSynthRefShift, seqlengths(synthBin.gr))

ol <- findOverlaps(synthBin.gr, lads.gr)
pInt <- pintersect(synthBin.gr[queryHits(ol)], lads.gr[subjectHits(ol)])

ladSum <- data_frame( queryHit = queryHits(ol), width = width(pInt)) %>% 
  group_by(queryHit) %>% 
  summarise(width = sum(width))
  
synthBin.gr$ladCov <- 0
synthBin.gr$ladCov[ladSum$queryHit] <- ladSum$width

# replication Timing

# merge all rt data and just get the mean level per bin




# TE distribution 

load("~/Desktop/RTN_domains/R_objects/rmskTables/hg19.RData")

rep <- rep[complete.cases(rep),]

rep.gr <- GRanges(seqnames = rep$genoChr, 
                  ranges = IRanges(start = rep$genoStart + 1, end = rep$genoEnd),
                  repGroup = rep$repGroup,
                  repFamily = rep$repFamily)

rep.gr$repGroup[rep.gr$repFamily == "L1MA"] <- "old_L1"

seqlevels(rep.gr) <- refChrInfo$chrom
seqlengths(rep.gr) <- refChrInfo$size
genome(rep.gr) <- genomes["ref"]
rep.gr <- sort(sortSeqlevels(rep.gr))

repExpand.gr <- genoExpandBreak(x.gr = rep.gr, 
                                synthGenome = newSynthRefShift, 
                                expandedSeqlengths = seqlengths(synthBin.gr))

ol <- findOverlaps(synthBin.gr, repExpand.gr)

pInt <- pintersect(repExpand.gr[subjectHits(ol)], synthBin.gr[queryHits(ol)])

repSum <- data_frame( queryHit = queryHits(ol), width = width(pInt), repGroup = pInt$repGroup) %>% 
  group_by(queryHit, repGroup) %>% 
  summarise(width = sum(width))

repSum <- melt(repSum, c("queryHit", "repGroup"), "width")

repSum <- cast(repSum, queryHit ~ repGroup)

mcols(synthBin.gr)[,colnames(repSum)[2:ncol(repSum)]] <- as.integer(NA)
synthBin.gr[repSum$queryHit]$ancient <- repSum$ancient
synthBin.gr[repSum$queryHit]$new_L1 <- repSum$new_L1
synthBin.gr[repSum$queryHit]$new_SINE <- repSum$new_SINE
synthBin.gr[repSum$queryHit]$old_L1 <- repSum$old_L1

blueRes <- colorRampPalette(c("blue", "white", "red"))
heatmap(cor(as.matrix(mcols(synthBin.gr))[complete.cases(mcols(synthBin.gr)),]), scale = "none", col = blueRes(20), zlim = c(-1,1))

# Significant misclassification of mouse deletions in human, almost 10%
qDel <- GenomicRanges::intersect(refGap.gr, refAncDna.gr)
ol <- findOverlaps(rep.gr, refAncDna.gr)

pInt <- pintersect(rep.gr[queryHits(ol)], refAncDna.gr[subjectHits(ol)])

repSum <- data_frame( queryHit = queryHits(ol), width = width(pInt), repGroup = pInt$repGroup, repFamily = pInt$repFamily) %>%
  group_by(repGroup) %>%
  summarise(width = sum(width))

repSum2 <- data_frame( width = width(rep.gr), repGroup = rep.gr$repGroup, repFamily = rep.gr$repFamily) %>%
  group_by(repGroup) %>%
  summarise(width = sum(width))

hist(width(pInt[pInt$repFamily == "L1PA"]), breaks = 100)

sample(pInt[pInt$repFamily == "AluS" & width(pInt) > 200], 10)

repSum$width/sum(width(refAncDna.gr))

sum(width(GenomicRanges::intersect(refAncDna.gr, 
                                   GenomicRanges::union(refFill.gr, refGap.gr))))

# distance from end 

# 10 % false positive rate
# in the human there are 10% more insertions, and in the mouse there are 10% less deletions
# this may explain our results for genome shrinkage, espically oi GC rich areas
# damn
# at least now we know about it and have a stratergy to cope
# we can provide a quick patch up that takes the set difference of the ancestral regions as we bring them in
# this way they will be filtered from repeats


# get a new repeats table and fix this problem for good

s <- start(synthBin.gr)
sR <- seqlengths(synthBin.gr)[as.character(seqnames(synthBin.gr))] - s

synthBin.gr$distFromTelomere <- pmin(s,sR)





### recombination rate, may get these from else where

recombRate <- read.table("~/Desktop/RTN_domains/data/UCSCtracks/hg19/recombRate.txt",
                         col.names = c("seqnames",	"start",	"end",	"name",	
                                       "decodeAvg",	"decodeFemale",	"decodeMale",	
                                       "marshfieldAvg",	"marshfieldFemale",	"marshfieldMale",	
                                       "genethonAvg",	"genethonFemale",	"genethonMale"))


hist(recombRate$decodeAvg, breaks = 100)

layout(1:6)
chrChoice = "chr6"
plot(recombRate$start[recombRate$seqnames == chrChoice], 
     recombRate$decodeMale[recombRate$seqnames == chrChoice], type = "l")
plot(recombRate$start[recombRate$seqnames == chrChoice], 
      recombRate$marshfieldMale[recombRate$seqnames == chrChoice], type = "l")
plot(recombRate$start[recombRate$seqnames == chrChoice], 
      recombRate$genethonMale[recombRate$seqnames == chrChoice], type = "l")


plot(recombRate$start[recombRate$seqnames == chrChoice], 
     recombRate$decodeFemale[recombRate$seqnames == chrChoice], type = "l")
plot(recombRate$start[recombRate$seqnames == chrChoice], 
      recombRate$marshfieldFemale[recombRate$seqnames == chrChoice], type = "l")
plot(recombRate$start[recombRate$seqnames == chrChoice], 
      recombRate$genethonFemale[recombRate$seqnames == chrChoice], type = "l")

























# could set up synthetic data to see how it works


