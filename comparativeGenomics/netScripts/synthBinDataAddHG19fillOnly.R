
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
load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/",specRef,".synthBin.RData",sep = ""))
# fills and gaps
load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/",specRef,".",specQue,".netData.RData",sep = ""))
#load(paste("Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/",specRef,".stretch.RData", sep = ""))
# shift levels
load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/shiftData/",specRef,".expand.breaks.RData", sep = ""))

mapUnique <- read.table("~/Documents/dna_turnover/workStationDesktop/RTN_domains/data/mappability/uniqueMap/hg1936Unique.bed", 
                        col.names = c("seqnames", "start", "end", "score"))
mapUnique$start <- mapUnique$start +1
mapUnique.gr <- GRanges(mapUnique)
seqlevels(mapUnique.gr) <- refChrInfo$chrom
seqlengths(mapUnique.gr) <- refChrInfo$size
genome(mapUnique.gr) <- genomes["ref"]
mapUnique.gr <- genoExpandStretch(mapUnique.gr, newSynthRefShift, seqlengths(synthBin.gr))

ol <- findOverlaps(synthBin.gr, mapUnique.gr)
pInt <- pintersect(synthBin.gr[queryHits(ol)], mapUnique.gr[subjectHits(ol)])

mapUniqueSum <- summarise(group_by(data.frame(queryHits = queryHits(ol), 
                                         mapUnique = width(pInt)), 
                              queryHits), 
                     mapUnique = sum(mapUnique))

synthBin.gr$mapUnique <- 0
synthBin.gr$mapUnique[mapUniqueSum$queryHits] <- mapUniqueSum$mapUnique





dnase <- read.table("~/Documents/dna_turnover/workStationDesktop/RTN_domains/data/UCSCtracks/hg19/wgEncodeAwgDnaseMasterSites.txt",
                    col.names = c("bin", "seqnames", "start", "end", "name", 
                                  "dispScore", "floatScore", "sourceCount", "sourceID"),
                    colClasses = c("integer", "character", "integer", "integer", 
                                   "integer", "integer", "double", "integer", "character"))
dnase$start = dnase$start + 1
dnase.gr <- GRanges(dnase)
seqlevels(dnase.gr) <- refChrInfo$chrom
seqlengths(dnase.gr) <- refChrInfo$size
genome(dnase.gr) <- genomes["ref"]

#filter out non-fill elements
#dnase.gr <- dnase.gr[overlapsAny(dnase.gr,refFill.gr, minoverlap = 75)]

dnase.gr <- genoExpandStretch(dnase.gr, newSynthRefShift, seqlengths(synthBin.gr))


ol <- findOverlaps(synthBin.gr, dnase.gr)

dnaseSum <- summarise(group_by(data.frame(queryHits = queryHits(ol), 
                                          peaks = dnase.gr$sourceCount[subjectHits(ol)]), 
                               queryHits), 
                      peakNumbers = sum(peaks), peakSites = n())

synthBin.gr$dnasePeaks <- 0
synthBin.gr$dnasePeaks[dnaseSum$queryHits] <- dnaseSum$peakSites 

synthBin.gr$dnaseActivity <- 0
synthBin.gr$dnaseActivity[dnaseSum$queryHits] <- dnaseSum$peakNumbers 



#### Exons and Introns


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
exon.gr <- reduce(exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
#xon.gr <- intersect(exon.gr, refFill.gr)


intron.gr <- GenomicRanges::setdiff(unlist(intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)), exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
#intron.gr <- intersect(intron.gr, refFill.gr)

exon.gr <- genoExpandBreak(exon.gr, newSynthRefShift, seqlengths(synthBin.gr))


ol <- findOverlaps(synthBin.gr, exon.gr)
pInt <- pintersect(synthBin.gr[queryHits(ol)], exon.gr[subjectHits(ol)])

exonSum <- summarise(group_by(data.frame(queryHits = queryHits(ol), 
                                          exon = width(pInt)), 
                               queryHits), 
                        exon = sum(exon))

synthBin.gr$exon <- 0
synthBin.gr$exon[exonSum$queryHits] <- exonSum$exon


intron.gr <- genoExpandBreak(intron.gr, newSynthRefShift, seqlengths(synthBin.gr))


ol <- findOverlaps(synthBin.gr, intron.gr)
pInt <- pintersect(synthBin.gr[queryHits(ol)], intron.gr[subjectHits(ol)])

intronSum <- summarise(group_by(data.frame(queryHits = queryHits(ol), 
                                            intron = width(pInt)), 
                                 queryHits), 
                        intron = sum(intron))

synthBin.gr$intron <- 0
synthBin.gr$intron[intronSum$queryHits] <- intronSum$intron



# laminaAssocaited domains 

lads <- read.table("~/Documents/dna_turnover/workStationDesktop/RTN_domains/data/UCSCtracks/hg19/laminB1Lads", 
                      col.names = c("bin", "seqnames", "start", "end"))
lads.gr <- GRanges(lads)
# from UCSC db need to convert to 1 base
start(lads.gr) <- start(lads.gr) + 1
seqlevels(lads.gr) <- refChrInfo$chrom
seqlengths(lads.gr) <- refChrInfo$size
genome(lads.gr) <- genomes["ref"]
lads.gr <- sort(sortSeqlevels(lads.gr))
#lads.gr <- intersect(lads.gr, refFill.gr)

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






# motifs 
# 
# motifs <- read.table("~/Desktop/RTN_domains/data/UCSCtracks/hg19/hg19.motifs", header = TRUE, 
#                     colClasses = c("character", "integer", "integer", "character", "character"))
# 
# motifs.gr <- GRanges(motifs)
# 
# seqlevels(motifs.gr) <- refChrInfo$chrom
# seqlengths(motifs.gr) <- refChrInfo$size
# genome(motifs.gr) <- genomes["ref"]
# motifs.gr <- sort(sortSeqlevels(motifs.gr))
# 
# motifsExpand.gr <- genoExpandStretch(motifs.gr, newSynthRefShift, seqlengths(synthBin.gr))
# 
# ol <- findOverlaps(synthBin.gr, motifsExpand.gr)
# motifSum <- data_frame( queryHit = queryHits(ol), patternID = motifsExpand.gr[subjectHits(ol)]$patternID) %>% 
#   group_by(queryHit, patternID) %>% 
#   summarise(hits = n())
# 
# motifSum <- melt(motifSum, c("queryHit", "patternID"), "hits")
# motifSum <- cast(motifSum, queryHit ~ patternID)
# motifSum[is.na(motifSum)] <- 0
# 
# mcols(synthBin.gr)[,colnames(motifSum)[2:ncol(motifSum)]] <- 0
# 
# synthBin.gr[motifSum$queryHit]$ctcfMotif <- motifSum$ctcfMotif
# synthBin.gr[motifSum$queryHit]$L1Motif <- motifSum$L1Motif
# synthBin.gr[motifSum$queryHit]$prdm9Motif <- motifSum$prdm9Motif
# 



# TE distribution 

load("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/rmskTables/hg19.RData")

rep <- rep[complete.cases(rep),]

rep.gr <- GRanges(seqnames = rep$genoChr, 
                  ranges = IRanges(start = rep$genoStart + 1, end = rep$genoEnd),
                  repGroup = rep$repGroup,
                  repFamily = rep$repFamily)

# move L1MA to old group, because we are looking at a closer distance than our old paper
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

repSum <- melt(repSum, c("queryHit", "repGroup"), "width") %>% cast(queryHit ~ repGroup)
repSum[is.na(repSum)] <- 0


mcols(synthBin.gr)[,colnames(repSum)[2:ncol(repSum)]] <- 0
synthBin.gr[repSum$queryHit]$ancient <- repSum$ancient
synthBin.gr[repSum$queryHit]$new_L1 <- repSum$new_L1
synthBin.gr[repSum$queryHit]$new_SINE <- repSum$new_SINE
synthBin.gr[repSum$queryHit]$old_L1 <- repSum$old_L1


# # Significant misclassification of mouse deletions in human, almost 10%
# qDel <- GenomicRanges::intersect(refGap.gr, refAncDna.gr)
# ol <- findOverlaps(rep.gr, refAncDna.gr)
# 
# pInt <- pintersect(rep.gr[queryHits(ol)], refAncDna.gr[subjectHits(ol)])
# 
# repSum <- data_frame( queryHit = queryHits(ol), width = width(pInt), repGroup = pInt$repGroup, repFamily = pInt$repFamily) %>%
#   group_by(repGroup) %>%
#   summarise(width = sum(width))
# 
# repSum2 <- data_frame( width = width(rep.gr), repGroup = rep.gr$repGroup, repFamily = rep.gr$repFamily) %>%
#   group_by(repGroup) %>%
#   summarise(width = sum(width))
# 
# hist(width(pInt[pInt$repFamily == "AluJ"]), breaks = 1000)
# 
# sample(pInt[pInt$repFamily == "AluS" & width(pInt) > 200], 10)
# 
# repSum$width/sum(width(refAncDna.gr))
# 
# sum(width(GenomicRanges::intersect(refAncDna.gr, 
#                                    GenomicRanges::union(refFill.gr, refGap.gr))))
# 
# sum(width(GenomicRanges::intersect(refAncDna.gr, refFill.gr)))/ sum(width(refFill.gr))
# distance from end 

# 10 % false positive rate
# in the human there are 10% more insertions, and in the mouse there are 10% less deletions
# this may explain our results for genome shrinkage, espically oi GC rich areas
# damn
# at least now we know about it and have a stratergy to cope
# we can provide a quick patch up that takes the set difference of the ancestral regions as we bring them in
# this way they will be filtered from repeats


# get a new repeats table and fix this problem for good

# 
# library(RMySQL)
# mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = "hg19")
# seqGaps <- dbGetQuery(mychannel, "SELECT * FROM gap;")
# centromere.df <- seqGaps[seqGaps$type == "centromere",]
# telomere.df <- seqGaps[seqGaps$type == "telomere",]
# 
# centromere.gr <- GRanges(seqnames = centromere.df$chrom,
#                          ranges = IRanges(start = centromere.df$chromStart, end = centromere.df$chromEnd))
# seqlevels(centromere.gr) <- refChrInfo$chrom
# seqlengths(centromere.gr) <- refChrInfo$size
# genome(centromere.gr) <- genomes["ref"]
# centromere.gr <- sort(sortSeqlevels(centromere.gr))
# 
# centromere.gr <- genoExpandStretch(centromere.gr, newSynthRefShift, seqlengths(synthBin.gr))
# 
# 
# 
# telomere.gr <- GRanges(seqnames = telomere.df$chrom,
#                          ranges = IRanges(start = telomere.df$chromStart + 1, end = telomere.df$chromEnd))
# seqlevels(telomere.gr) <- refChrInfo$chrom
# seqlengths(telomere.gr) <- refChrInfo$size
# genome(telomere.gr) <- genomes["ref"]
# 
# seqlevels(telomere.gr)[table(seqnames(telomere.gr)) == 0]
# 
# extraTelomete.gr <- GRanges(seqnames <- rep(seqlevels(telomere.gr)[table(seqnames(telomere.gr)) == 0],2),
#                              ranges = IRanges(start = c(rep(1, sum(table(seqnames(telomere.gr)) == 0)),
#                                                         seqlengths(telomere.gr)[table(seqnames(telomere.gr)) == 0]),
#                                               width = 1))
# telomere.gr <- c(telomere.gr, extraTelomete.gr)                                    
# 
# 
# telomere.gr <- sort(sortSeqlevels(telomere.gr))
# 
# telomere.gr <- genoExpandStretch(telomere.gr, newSynthRefShift, seqlengths(synthBin.gr))
# 
# 
# distFromTelomere <- distanceToNearest(synthBin.gr, telomere.gr)
# distFromCentromere <- distanceToNearest(synthBin.gr, centromere.gr)
# 
# synthBin.gr$distFromTelomere <- 0
# synthBin.gr$distFromTelomere[queryHits(distFromTelomere)] <- mcols(distFromTelomere)$distance
# 
# synthBin.gr$distFromCentromere  <- 0
# synthBin.gr$distFromCentromere[queryHits(distFromCentromere)] <- mcols(distFromCentromere)$distance
# 


### recombination rate, may get these from else where

load("~/Documents/dna_turnover/workStationDesktop/RTN_domains/data/UCSCtracks/hg19/hg19.recombinationRate.RData")


recRate.gr <- GRanges(seqnames = RecombinationRate$Chromosome,
                      ranges = IRanges(start = RecombinationRate$Position.bp., width = 1),
                      Rate.cM.Mb. = RecombinationRate$Rate.cM.Mb.,
                      Map.cM. = RecombinationRate$Rate.cM.)

seqlevels(recRate.gr) <- refChrInfo$chrom
seqlengths(recRate.gr) <- refChrInfo$size
genome(recRate.gr) <- genomes["ref"]
recRate.gr <- sort(sortSeqlevels(recRate.gr))

recRateExpand.gr <- genoExpandBreak(x.gr = recRate.gr, 
                                synthGenome = newSynthRefShift, 
                                expandedSeqlengths = seqlengths(synthBin.gr))


ol <- findOverlaps(synthBin.gr, recRateExpand.gr)


recRateSum <- data_frame( queryHit = queryHits(ol), rate = recRateExpand.gr[subjectHits(ol)]$Rate.cM.Mb.) %>% 
  group_by(queryHit) %>% 
  summarise(rate = mean(rate))



synthBin.gr$recombRate <- 0

synthBin.gr[recRateSum$queryHit]$recombRate <- recRateSum$rate



### recmbination hotspot

recHot <- read.table("~/Documents/dna_turnover/workStationDesktop/RTN_domains/data/UCSCtracks/hg19/hg19.recombinationHotspots.bed",
                     col.names = c("seqnames", "start", "end"))

recHot.gr <- GRanges(recHot)

seqlevels(recHot.gr) <- refChrInfo$chrom
seqlengths(recHot.gr) <- refChrInfo$size
genome(recHot.gr) <- genomes["ref"]
recHot.gr <- sort(sortSeqlevels(recHot.gr))

recHotExpand.gr <- genoExpandBreak(x.gr = recHot.gr, 
                                synthGenome = newSynthRefShift, 
                                expandedSeqlengths = seqlengths(synthBin.gr))

ol <- findOverlaps(synthBin.gr, recHotExpand.gr)

pInt <- pintersect(recHotExpand.gr[subjectHits(ol)], synthBin.gr[queryHits(ol)])

recHotSum <- data_frame( queryHit = queryHits(ol), width = width(pInt)) %>% 
  group_by(queryHit) %>% 
  summarise(width = sum(width))

synthBin.gr$recHotSpot <- 0
synthBin.gr[recHotSum$queryHit]$recHotSpot <- recHotSum$width


## cpG island

cpg <- read.table("~/Documents/dna_turnover/workStationDesktop/RTN_domains/data/UCSCtracks/hg19/cpgIslandExtUnmasked.txt",
                  col.names = c("bin", "seqnames", "start" , "end","name", "cpgNum",
                                "length", "cpgNum", "gcNum", "perCpG" ,"perGc", "obsExp"))

cpg.gr <- GRanges(cpg)
start(cpg.gr) <- start(cpg.gr) + 1

seqlevels(cpg.gr) <- refChrInfo$chrom
seqlengths(cpg.gr) <- refChrInfo$size
genome(cpg.gr) <- genomes["ref"]
cpg.gr <- sort(sortSeqlevels(cpg.gr))

cpgExpand.gr <- genoExpandBreak(x.gr = cpg.gr, 
                                   synthGenome = newSynthRefShift, 
                                   expandedSeqlengths = seqlengths(synthBin.gr))

ol <- findOverlaps(synthBin.gr, cpgExpand.gr)

pInt <- pintersect(cpgExpand.gr[subjectHits(ol)], synthBin.gr[queryHits(ol)])

cpgSum <- data_frame( queryHit = queryHits(ol), width = width(pInt)) %>% 
  group_by(queryHit) %>% 
  summarise(width = sum(width))

synthBin.gr$cpgCov <- 0
synthBin.gr[cpgSum$queryHit]$cpgCov <- cpgSum$width




save(synthBin.gr, file = "~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/hg19.synthBin.genome.variable.fillOnly.RData")














