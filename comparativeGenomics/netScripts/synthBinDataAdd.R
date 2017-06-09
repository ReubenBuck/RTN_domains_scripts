# what other things can I measure that can assocaite with dependant variable

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



library(TxDb.Hsapiens.UCSC.hg19.knownGene)

ol <- findOverlaps(intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene), newSynthRefShift) 


sort(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))


# how to define house keeping genes


dnase <- read.table("~/Desktop/RTN_domains/data/UCSCtracks/hg19/wgEncodeAwgDnaseMasterSites.txt",
                    col.names = c("bin", "seqnames", "start", "end", "name", 
                                  "dispScore", "floatScore", "sourceCount", "sourceID"),
                    colClasses = c("integer", "character", "integer", "integer", 
                                   "integer", "integer", "double", "integer", "character"))
eFunc <- ecdf(dnase$sourceCount)
plot(ecdf(dnase$sourceCount))


df <- data.frame(x = 0:125, y = eFunc(seq(0,125)))

1/125
plot(df)
abline(v = which.max((df$y - seq(0,1,length.out = 126))))

# so we can use the vertex of the elbow curve to classify into nonActive/active
# then we can begin to unpack some of the regulatory relationships behind this

# total activity 

# active number 
# nonActive number 
# active:nonActive


hist(dnase$sourceCount, breaks = 125, ylim = c(0,100000))



lines(seq(0,1,length.out = 126))


plot(df$y - seq(0,1,length.out = 126))


?deriv

abline(v = 20)
# we could use the elbow of the curve to seperate into two groups




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



if(specRef = "hg19"){
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  exon.gr <- reduce(exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
  intron.gr <- setdiff(unlist(intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene)), exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
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



# distance from end 


s <- start(synthBin.gr)
sR <- seqlengths(synthBin.gr)[as.character(seqnames(synthBin.gr))] - s

synthBin.gr$distFromTelomere <- pmin(s,sR)

















### recombination rate

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



library(spdep)
library(igraph)



df <- data.frame(synthBin.gr)
df$chrType = "autosome"
df[df$seqnames == "chrX","chrType"] <- "sexChr"

df[df$seqGap + df$missingGap > 20000,] <- NA
df <- df[complete.cases(df),]

df$gcContent <- df$gcContent * 100

ol<-findOverlaps(synthBin.gr, maxgap = 30*width(synthBin.gr)[1])
ol <- ol[!(isRedundantHit(ol))]
# remove NA hits
ol <- ol[!is.na(df$refIns[queryHits(ol)])]
ol <- ol[!is.na(df$refIns[subjectHits(ol)])]
olMat <- data.frame(ol)
G <- graph.data.frame(d = olMat,directed=FALSE)
weight <- olMat$subjectHits - olMat$queryHits
weight <- -(weight - (max(weight) + 1)) / max(weight)
weight[isSelfHit(ol)] <- 0
E(G)$weight <- weight
A <- as_adjacency_matrix(G,type="both",names=FALSE,sparse=TRUE,edges = FALSE, attr = "weight")




dfChr <- df[df$seqnames == "chr1"| df$seqnames == "chr2",]

A <- A[df$seqnames == "chr1" | df$seqnames == "chr2" ,df$seqnames == "chr1" | df$seqnames == "chr2"]
wMat <- mat2listw(A,)


# con




spatMod <- errorsarlm(data = df, refDel ~ gcContent + dnasePeaks + dnaseActivity + exon + intron, listw = wMat)

summary(spatMod)


qqnorm(residuals(spatMod))

acf(residuals(spatMod),lag.max = 40)

plot(fitted(spatMod),residuals(spatMod))

# so now we have a modeling appraoch that might actually work


lm.morantest(spatMod, wMat)
moran.test(x = residuals(spatMod), wMat)
moran.test(x = residuals(modChr1), wMat)


# so now there is no autocorrelation in our residuals
# therefore we could succesfully remove the error

# next to get the R squered to determine our variance explained by the model.



modChr1 <- lm(data = dfChr, (queIns - queDel) ~ gcContent + dnasePeaks + dnaseActivity + exon + intron )

acf(residuals(modChr1),lag.max = 100)



h <- hist(scale(residuals(spatMod)), breaks = 1000)
plot( sort(rnorm(h$mids)), h$mids, xlim = c(-4,4), ylim = c(-4,4))


ESS <- sum((fitted(modChr1) - mean((df$queIns - df$queDel)))^2)
RSS <- sum((fitted(modChr1) - mean((df$queIns - df$queDel)))^2)

ts <- as.ts(modChr1)

acf(residuals(spatMod))
acf(residuals(modChr1))

# it is removing the auto correlation

# so we have a significant Morans I 
# and a significant autocorrelation
# when we apply the new spatial model we can factor that out.
lm.morantest(modQueIns, listw = wMat)



## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
g <- graph_from_data_frame(relations, directed=TRUE, vertices=actors)
print(g, e=TRUE, v=TRUE)

