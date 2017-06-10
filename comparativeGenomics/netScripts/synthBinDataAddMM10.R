
# transposable elements, our four groups
# gene density
# coding density
# recombination rate
# replication timing

# DSP hotspot
# Genetic recombination is directed away from functional genomic elements in mice


# conserved elements





# gap pieces per gap bp

# the average size of the event, there is a distribution of events




laminB1 <- read.table("~/Desktop/RTN_domains/data/UCSCtracks/hg19/laminB1Lads", 
                      col.names = c("bin", "seqnames", "start", "end"))
# these are usually 100 thousand, means we can give a binary classification for lammina domain





rm(list = ls())

options(stringsAsFactors = FALSE)



library(dplyr)
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



if(specRef == "hg19"){
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  exon.gr <- reduce(exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
  intron.gr <- GenomicRanges::setdiff(unlist(intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)), exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
}

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

laminB1 <- read.table("~/Desktop/RTN_domains/data/UCSCtracks/hg19/laminB1Lads", 
                      col.names = c("bin", "seqnames", "start", "end"))






# distance from end 


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











cor(synthBin.gr$refIns, synthBin.gr$dnasePeaks, use = "complete.obs")

smoothScatter((synthBin.gr$refIns), log10(synthBin.gr$cdsExon) )

cor((synthBin.gr$queIns + synthBin.gr$queDel), (synthBin.gr$cdsExon) , use = "complete.obs", method = "spearman")


hist(log10(synthBin.gr$cdsExon), breaks = 100)

# so there may be some strange efects regarding rate of genome evolution and activity
synthBin.gr[synthBin.gr$cdsExon > 50000 & !is.na(synthBin.gr$cdsExon)]

# it might be worth holding on to the original coordinates for each bin 

plot(start(synthBin.gr[seqnames(synthBin.gr) == "chr2"]),synthBin.gr[seqnames(synthBin.gr) == "chr2"]$cdsExon )


