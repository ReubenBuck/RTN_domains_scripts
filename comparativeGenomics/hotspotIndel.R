## Produce all our plots across one script

## so lets go back a couple of steps 

## what are we actually measuring with our 4 stats

## human bases that have been lost since divergence with mouse

## human bases that have been gained since divergence with mouse

## mouse bases that have been gaines since divergence with human 

## mouse bases that have been lost since divergence with human

## We are looking to see how the independant gain and loss of DNA along each lineage associates with changes in TE content

## in shared regions, we are expecting to see similar levels of gain 

## in new TE lineage specific, lineage specific amounts of gain

## how these two things relate 
## need to correctly assign values too 


library(GenomicRanges)
library(dplyr)
### 
# regions identified



# differnces in TE content

hgGap <- read.table("Desktop/RTN_domains/data/comparativeGenomics/queGenomes/gaps/hg19.mm10.que.indel", header = TRUE)
hgGap.gr <- GRanges(hgGap)
mcols(hgGap.gr)$gapWidth <- mcols(hgGap.gr)$queRange.width
mcols(hgGap.gr[mcols(hgGap.gr)$inDel == "ins"])$gapWidth <- width(hgGap.gr[mcols(hgGap.gr)$inDel == "ins"])


hgHot <- read.table("Desktop/RTN_domains/data/repeatHotspot/hg19/hg19_mm9_conDif.txt", header=TRUE)
hgHot.gr <- GRanges(hgHot)
hgHot.gr <- hgHot.gr[mcols(hgHot.gr)$genome == "ref"]

#head(hg)


# how do we get the lineage specific 


hgIndel <- read.table("~/Desktop/RTN_domains/data/comparativeGenomics/supportedIndels/hg19.supportedIndel.que", header = TRUE)
hgIndel.gr <- GRanges(hgIndel)


hgBr <- read.table("~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/baseRates/hg19.base", header = TRUE)
hgBr.gr <- GRanges(hgBr)


mmIndel <- read.table("~/Desktop/RTN_domains/data/comparativeGenomics/supportedIndels/mm10.supportedIndel.que", header = TRUE)
#mmBr <- read.table("~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/baseRates/mm10.base")
mmIndel.gr <- GRanges(seqnames = Rle(mmIndel$queRange.seqnames),
                      ranges = IRanges(start = mmIndel$queRange.start, end = mmIndel$queRange.end))
                      
mcols(mmIndel.gr) <- mmIndel[,11:ncol(mmIndel)]




# 1 for 1 overlap conflicting
ol <- findOverlaps(hgIndel.gr, mmIndel.gr, minoverlap = 1)
smoothScatter(mcols(hgIndel.gr[queryHits(ol)])$supGenoNo, mcols(mmIndel.gr[subjectHits(ol)])$supGenoNo)

# is there a bad species 
a1 <- table(unlist(strsplit(as.character(mcols(hgIndel.gr)$supGenoName), split = ",")))
a2 <- table(unlist(strsplit(as.character(mcols(hgIndel.gr[queryHits(ol)])$supGenoName), split = ",")))

plot(a1/sum(a1))
points(a2/sum(a2), col = 2, type = "p")

b1 <- table(unlist(strsplit(as.character(mcols(mmIndel.gr)$supGenoName), split = ",")))
b2 <- table(unlist(strsplit(as.character(mcols(mmIndel.gr[subjectHits(ol)])$supGenoName), split = ",")))

plot(b1/sum(b1))
points(b2/sum(b2), col = 2, type = "p")


echTel <- grep("echTel2",as.character(mcols(hgIndel.gr[queryHits(ol)])$supGenoName))
hist(mcols(hgIndel.gr[queryHits(ol)])$supGenoNo, breaks = 20, freq = FALSE, ylim = c(0,2))
hist(mcols(hgIndel.gr[queryHits(ol)][echTel])$supGenoNo, add = T, col = 2, density = 0, breaks = 20, freq = FALSE)

hist(mcols(hgIndel.gr[queryHits(ol)][echTel])$supGenoNo, col = 2, density = 0, breaks = 20)



olSamp <- ol[mcols(mmIndel.gr[subjectHits(ol)])$supGenoNo == mcols(hgIndel.gr[queryHits(ol)])$supGenoNo]
a3 <- table(unlist(strsplit(as.character(mcols(hgIndel.gr[queryHits(olSamp)])$supGenoName), split = ",")))
b3 <- table(unlist(strsplit(as.character(mcols(mmIndel.gr[subjectHits(olSamp)])$supGenoName), split = ",")))

mmIndel.gr[subjectHits(olSamp)]

plot(a1/sum(a1))
points(a2/sum(a2), col = 2, type = "p")
points(a3/sum(a3), col = 3, type = "p")


plot(b1/sum(b1))
points(b2/sum(b2), col = 2, type = "p")
points(b3/sum(b3), col = 3, type = "p")

smoothScatter(mcols(hgIndel.gr[queryHits(ol)])$supGenoNo,mcols(mmIndel.gr[subjectHits(ol)])$supGenoNo )


# resolve by counting votes

length(ol)






olIndel <- (findOverlaps(hgHot.gr, hgIndel.gr))
dfIndel <- data.frame(queryHits = queryHits(olIndel), mcols(hgIndel.gr[subjectHits(olIndel)])[c("gapWidth","indel")] )
dfIndel <- summarise(group_by(dfIndel,queryHits, indel = indel), gapWidth = sum(gapWidth))

dfIndel$genome = "hg19"

olMindel <- findOverlaps(hgHot.gr, mmIndel.gr)
dfMindel <- data.frame(queryHits = queryHits(olMindel), mcols(mmIndel.gr[subjectHits(olMindel)])[c("gapWidth","indel")])
dfMindel <- summarise(group_by(dfMindel, queryHits, indel = indel), gapWidth = sum(gapWidth))

dfMindel$genome = "mm10"

dfIndel <- rbind(dfIndel, dfMindel)


olBr <- findOverlaps(hgHot.gr, hgBr.gr)
dfBr <- data.frame(queryHits = queryHits(olBr), baseRate = width(hgBr.gr[subjectHits(olBr)]))
dfBr <- summarise(group_by(dfBr, queryHits), baseRate = sum(baseRate))
mer <- merge(dfIndel, dfBr)

olGap <- findOverlaps(hgHot.gr, hgGap.gr)
dfGap <- data.frame(queryHits = queryHits(olGap), mcols(hgGap.gr[subjectHits(olGap)])[c("inDel","gapWidth")])
dfGap <- summarise(group_by(dfGap, queryHits, GapIndel = inDel), HMgapWidth = sum(gapWidth))

mer <- merge(mer, dfGap)


dfAll <- data.frame(mcols(hgHot.gr[mer$queryHits]), mer)


s <- summarise(group_by(dfAll, repGroup, hotspotID, genome.1, conState, indel, GapIndel), 
               gapWidth = sum(gapWidth), baseRate = sum(baseRate), known = sum(known), HMgapWidth = sum(HMgapWidth))

s <- s[s$conState!="mid",]


s <- filter(s,(genome.1 == "hg19" & indel == "del" & GapIndel == "del") |
              (genome.1 == "hg19" & indel == "ins" & GapIndel == "ins") |
              (genome.1 == "mm10" & indel == "ins" & GapIndel == "ins") |
              (genome.1 == "mm10" & indel == "del" & GapIndel == "del"))


layout(1)
par(mar = c(10,5,5,5))
boxplot((s$gapWidth/s$baseRate) ~ s$conState + s$repGroup + s$indel + s$genome.1, las = 2, notch = TRUE, log = "y")

boxplot((s$gapWidth/s$known) ~ s$repGroup + s$indel + s$genome.1, las = 2, notch = TRUE, log = "y")


abline(v = 12, lty = 2)
abline(v = 24, lty = 2)
abline(v = 36, lty = 2)





boxplot(log10(s$HMgapWidth/s$known) ~ s$repGroup + s$GapIndel, las = 2, notch = TRUE)
stripchart((s$HMgapWidth/s$known) ~ s$conState + s$repGroup + s$GapIndel,
           method= "jitter",jitter = .5, pch = 16, cex = .1, vert = TRUE, 
           ylim = c(0,.4), las = 2)


boxplot(log10(s$gapWidth/s$HMgapWidth) ~ s$conState + s$repGroup + s$GapIndel + s$genome.1, las = 2, notch = TRUE)


# relative contribution of insertions and deleations is similar


# we are trying to look at associated factors of what the gaps could be due to
# proportion of human sided gaps and proportion of mouse sided gaps 

# so this will help us identify particular classes of TE assocaited evolution 
# so really there is little variation between regions 


# it might mean TEs ahve little impact



# if we consider it as a proportion of gapped sequence. 








head(mmIndel)






# so now i have the full details of hat happned in human and mouse in human regions. 
# that might mean no more awkward matching up


layout(matrix(c(1,2), nrow = 1))
par(mar = c(5,4,4,2))

hist(log10(width(hgGap.gr[mcols(hgGap.gr)$inDel == "ins"])), breaks = 100, main = "Human side gaps",
     xlab = "width (log10 bp)", ylim = c(0,6e5), xlim = c(1,4.5))
hist(log10(mcols(hgIndel.gr[mcols(hgIndel.gr)$indel == "ins"])$gapWidth), breaks = 100, add = TRUE, col = 2, density = 0)
hist(log10(mcols(mmIndel.gr[mcols(mmIndel.gr)$indel == "del"])$gapWidth), breaks = 50, add = TRUE, col = 3, density = 0)
legend("topright",legend = c("all", "human insertion", "mouse deletion"), fill = c(1,2,3))



hist(log10(mcols(hgGap.gr[mcols(hgGap.gr)$inDel == "del"])$gapWidth), breaks = 100, main = "Mouse side gaps", 
     xlab = "width (log10 bp)", ylim = c(0,6e5),xlim = c(1,4.5))
hist(log10(mcols(mmIndel.gr[mcols(mmIndel.gr)$indel == "ins"])$gapWidth), breaks = 50, add = TRUE, col = 2, density = 0)
hist(log10(mcols(hgIndel.gr[mcols(hgIndel.gr)$indel == "del"])$gapWidth), breaks = 50, add = TRUE, col = 3, density = 0)
legend("topright",legend = c("all", "mouse insertion", "human deletion"), fill = c(1,2,3))


