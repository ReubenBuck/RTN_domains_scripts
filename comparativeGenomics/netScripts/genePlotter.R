

rm(list = ls())


# read in data and packages
library("Gviz")
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
load(file = "~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/mm10repNoRep.RData")
load(file = "~/Desktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/mm10.stretch.rep.RData")
load("~/Desktop/RTN_domains/R_objects/rmskTables/mm10/mm10New.RData")


mouseRange <- GRanges(seqnames = Rle("chr6"), 
                      ranges = IRanges(start = 80e6, end = 100e6))


# not all of our hotspots are all that hot.
# I don't think I correctly filtered them
# some of them are significant cold spots!!!!!
# shit
# thsi means I need to run a bunch of it again!!!!!

# place hotspots in origial cooridnates
sigG <- intersect(repSigRanges$refIns[seqnames(repSigRanges$refIns) == as.character(seqnames(mouseRange))], noRepSigRanges$refIns)
sigGgap <- gaps(sigG)
sigGgap <- sigGgap[strand(sigGgap) == "*" & seqnames(sigGgap) == as.character(seqnames(mouseRange))]
sigGgap$hotspot = 0
sigG <- granges(sigG,use.mcols = FALSE)
sigG$hotspot = 1
sigG <- sort(c(sigG,sigGgap))


downShift <- stretchedRef.gr[seqnames(stretchedRef.gr) == as.character(seqnames(mouseRange)) & (stretchedRef.gr$type == "queIns" | stretchedRef.gr$type == "refDel")]
ol <- findOverlaps(sigG, downShift)
pWid <- width(pintersect(sigG[queryHits(ol)], downShift[subjectHits(ol)]))
agg <- aggregate(x = pWid, by = list(queryHits(ol)), FUN = sum)
agg$x <- cumsum(agg$x) 

sigGShift <- sigG
start(sigGShift)[2:nrow(agg)] <- start(sigG)[2:nrow(agg)] - agg$x[1:(nrow(agg)-1)]
end(sigGShift) <- end(sigG)[agg$Group.1] - agg$x
hotspot <- sigGShift[overlapsAny(sigGShift,mouseRange)]
hotspot <- hotspot[hotspot$hotspot == 1]



# get TEs and make sure it is only the gap overlapping ones
newTE.gr <- newRep.gr[overlapsAny(newRep.gr, mouseRange)]
refIns <- stretchedRef.gr[stretchedRef.gr$type == "refIns"]
refIns <- refIns$refRanges[overlapsAny(refIns$refRanges, mouseRange)]
ol <- findOverlaps(newTE.gr, refIns)
newTE.gr <- pintersect(newTE.gr[queryHits(ol)], refIns[subjectHits(ol)])

sWindow <- unlist(slidingWindows(mouseRange,width = 100000,step = 100000))
ol <- findOverlaps(sWindow, refIns)
pWid <- width(pintersect(sWindow[queryHits(ol)], refIns[subjectHits(ol)]))
agg <- aggregate(pWid, list(queryHits(ol)), sum)
sWindow$newTransposon <- 0
sWindow$newTransposon[agg$Group.1] <- agg$x/1000
# range of interest


dT <- DataTrack(sWindow, type = "histogram")


# set up tracks
atrLINE <- AnnotationTrack(newTE.gr[grep("LINE", newTE.gr$repClass)], name="Recent LINEs",chromosome = as.character(seqnames(mouseRange)),genome = "mm10")
atrSINE <- AnnotationTrack(newTE.gr[grep("SINE", newTE.gr$repClass)], name="Recent SINEs",chromosome = as.character(seqnames(mouseRange)),genome = "mm10")
atrLTR <- AnnotationTrack(newTE.gr[grep("LTR", newTE.gr$repClass)], name="Recent LTRs",chromosome = as.character(seqnames(mouseRange)),genome = "mm10")

atrHS <- AnnotationTrack(hotspot, name="mm10 gain hotspot",chromosome = as.character(seqnames(mouseRange)),genome = "mm10")

atrNew <- AnnotationTrack(newTE.gr, name="Recent Transposons",chromosome = as.character(seqnames(mouseRange)),
                          genome = "mm10", fill=c(1,2,3), col.line=c(1,2,3))

gtr <- GenomeAxisTrack(range = mouseRange)
itr <- IdeogramTrack(genome="mm10", chromosome=as.character(seqnames(mouseRange)))


# getting gene objects
grtr <- GeneRegionTrack(txdb, name="Gene Model", showId=TRUE, chromosome = as.character(seqnames(mouseRange)), 
                        symbol = TRUE,stacking = "squish", start = start(mouseRange), end = end(mouseRange))

# chainging the names
gene(grtr)[gene(grtr) == "320189"] <- "231991"
x <- org.Mm.egSYMBOL
mappedGene <- mappedkeys(org.Mm.egSYMBOL)
symbol(grtr) <- unlist(as.list(x[gene(grtr)]))

plotTracks(list(itr, gtr,grtr , atrHS,dT), from = 80e6, to = 100e6, collapseTranscripts=FALSE, shape="arrow")







# chr7:24560000-29720000

library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

humanRange <- GRanges(seqnames = Rle("chr7"), 
                      ranges = IRanges(start = 24560000, end = 29720000))



load("~/Desktop/RTN_domains/R_objects/rmskTables/hg19/hg19New.RData")
newTE.gr <- newRep.gr[overlapsAny(newRep.gr, humanRange)]



# set up tracks
atrLINE <- AnnotationTrack(newTE.gr[grep("LINE", newTE.gr$repClass)], name="Recent LINEs",chromosome = "chr7",genome = "mm10")
atrSINE <- AnnotationTrack(newTE.gr[grep("SINE", newTE.gr$repClass)], name="Recent SINEs",chromosome = "chr7",genome = "mm10")
atrLTR <- AnnotationTrack(newTE.gr[grep("LTR", newTE.gr$repClass)], name="Recent LTRs",chromosome = "chr7",genome = "mm10")


gtr <- GenomeAxisTrack(range = humanRange)
itr <- IdeogramTrack(genome="hg19", chromosome="chr7")


# getting gene objects
grtr <- GeneRegionTrack(txdb, name="Gene Model", showId=TRUE, chromosome = "chr7", symbol = TRUE,stacking = "squish", start = 24560000, end = 29720000)

# chainging the names
x <- org.Hs.egSYMBOL
mappedGene <- mappedkeys(org.Hs.egSYMBOL)
symbol(grtr) <- unlist(as.list(x[gene(grtr)]))

plotTracks(list(itr, gtr,grtr ,atrLINE, atrSINE, atrLTR), from = 24560000, to = 29720000, collapseTranscripts=TRUE, shape="arrow")










# place hotspots in origial cooridnates
load(file = "~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/mm10repNoRep.RData")
load(file = "~/Desktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/mm10.stretch.rep.RData")

sigG <- intersect(repSigRanges$refIns[seqnames(repSigRanges$refIns) == "chr6"], noRepSigRanges$refIns)
sigGgap <- gaps(sigG)
sigGgap <- sigGgap[strand(sigGgap) == "*" & seqnames(sigGgap) == "chr6"]
sigGgap$hotspot = 0
sigG <- granges(sigG,use.mcols = FALSE)
sigG$hotspot = 1
sigG <- sort(c(sigG,sigGgap))


downShift <- stretchedRef.gr[seqnames(stretchedRef.gr) == "chr6" & (stretchedRef.gr$type == "queIns" | stretchedRef.gr$type == "refDel")]
ol <- findOverlaps(sigG, downShift)
pWid <- width(pintersect(sigG[queryHits(ol)], downShift[subjectHits(ol)]))
agg <- aggregate(x = pWid, by = list(queryHits(ol)), FUN = sum)
agg$x <- cumsum(agg$x) 

sigGShift <- sigG
start(sigGShift)[2:nrow(agg)] <- start(sigG)[2:nrow(agg)] - agg$x[1:(nrow(agg)-1)]
end(sigGShift) <- end(sigG)[agg$Group.1] - agg$x
hotspot <- sigGShift[overlapsAny(sigGShift,mouseRange)]
hotspot <- hotspot[hotspot$hotspot == 1]



sigGShift





