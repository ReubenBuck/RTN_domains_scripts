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


rm(list = ls())


inter <- read.table("~/Desktop/RTN_domains/data/repeatHotspot/intersect.txt", header = TRUE)

# differnces in TE content

# hgGap <- read.table("~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/gaps/hg19.mm10.que.indel", header = TRUE)
# hgGap.gr <- GRanges(hgGap)
# mcols(hgGap.gr)$gapWidth <- mcols(hgGap.gr)$queRange.width
# mcols(hgGap.gr[mcols(hgGap.gr)$inDel == "ins"])$gapWidth <- width(hgGap.gr[mcols(hgGap.gr)$inDel == "ins"])

hgHot <- read.table("~/Desktop/RTN_domains/data/repeatHotspot/hg19/hg19_mm10_conDif.txt", header=TRUE)
#hgHot <- hgHot[hgHot$hotspotID %in% inter$domains[inter$genome == "hg19"],]
hgHot.gr <- GRanges(hgHot)
hgHot.gr <- hgHot.gr[mcols(hgHot.gr)$genome == "ref"]

# how do we get the lineage specific 
hgIndel <- read.table("~/Desktop/RTN_domains/data/comparativeGenomics/supportedIndels/hg19.supportedIndel.que", header = TRUE)
hgIndel.gr <- GRanges(hgIndel)
mcols(hgIndel.gr)$queRange <- GRanges(seqnames = Rle(hgIndel$queRange.seqnames),
                                                   ranges = IRanges(start = hgIndel$queRange.start, end = hgIndel$queRange.end))

hgBr <- read.table("~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/baseRates/hg19.base", header = TRUE)
hgBr.gr <- GRanges(hgBr)



#mmGap <- read.table("~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/gaps/mm10.hg19.que.indel", header = TRUE)
#mmGap.gr <- GRanges(mmGap)
#mcols(mmGap.gr)$gapWidth <- mcols(mmGap.gr)$queRange.width
#mcols(mmGap.gr[mcols(mmGap.gr)$inDel == "ins"])$gapWidth <- width(mmGap.gr[mcols(mmGap.gr)$inDel == "ins"])

mmHot <- read.table("~/Desktop/RTN_domains/data/repeatHotspot/mm10/mm10_hg19_conDif.txt", header=TRUE)
#mmHot <- mmHot[mmHot$hotspotID %in% inter$domains[inter$genome == "mm10"],]
mmHot.gr <- GRanges(mmHot)
mmHot.gr <- mmHot.gr[mcols(mmHot.gr)$genome == "ref"]

mmIndel <- read.table("~/Desktop/RTN_domains/data/comparativeGenomics/supportedIndels/mm10.supportedIndel.que", header = TRUE)
mmIndel.gr <- GRanges(mmIndel)
mcols(mmIndel.gr)$queRange = GRanges(seqnames = Rle(mmIndel$queRange.seqnames),
                                     ranges = IRanges(start = mmIndel$queRange.start, end = mmIndel$queRange.end))

mmBr <- read.table("~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/baseRates/mm10.base", head = TRUE)
mmBr.gr <- GRanges(mmBr)






# 
# # 1 for 1 overlap conflicting
# ol <- findOverlaps(hgIndel.gr, mmIndel.gr, minoverlap = 1)
# smoothScatter(mcols(hgIndel.gr[queryHits(ol)])$supGenoNo, mcols(mmIndel.gr[subjectHits(ol)])$supGenoNo)
# 
# # is there a bad species 
# a1 <- table(unlist(strsplit(as.character(mcols(hgIndel.gr)$supGenoName), split = ",")))
# a2 <- table(unlist(strsplit(as.character(mcols(hgIndel.gr[queryHits(ol)])$supGenoName), split = ",")))
# 
# plot(a1/sum(a1))
# points(a2/sum(a2), col = 2, type = "p")
# 
# b1 <- table(unlist(strsplit(as.character(mcols(mmIndel.gr)$supGenoName), split = ",")))
# b2 <- table(unlist(strsplit(as.character(mcols(mmIndel.gr[subjectHits(ol)])$supGenoName), split = ",")))
# 
# plot(b1/sum(b1))
# points(b2/sum(b2), col = 2, type = "p")
# 
# 
# echTel <- grep("echTel2",as.character(mcols(hgIndel.gr[queryHits(ol)])$supGenoName))
# hist(mcols(hgIndel.gr[queryHits(ol)])$supGenoNo, breaks = 20, freq = FALSE, ylim = c(0,2))
# hist(mcols(hgIndel.gr[queryHits(ol)][echTel])$supGenoNo, add = T, col = 2, density = 0, breaks = 20, freq = FALSE)
# 
# hist(mcols(hgIndel.gr[queryHits(ol)][echTel])$supGenoNo, col = 2, density = 0, breaks = 20)
# 
# 
# 
# olSamp <- ol[mcols(mmIndel.gr[subjectHits(ol)])$supGenoNo == mcols(hgIndel.gr[queryHits(ol)])$supGenoNo]
# a3 <- table(unlist(strsplit(as.character(mcols(hgIndel.gr[queryHits(olSamp)])$supGenoName), split = ",")))
# b3 <- table(unlist(strsplit(as.character(mcols(mmIndel.gr[subjectHits(olSamp)])$supGenoName), split = ",")))
# 
# mmIndel.gr[subjectHits(olSamp)]
# 
# plot(a1/sum(a1))
# points(a2/sum(a2), col = 2, type = "p")
# points(a3/sum(a3), col = 3, type = "p")
# 
# 
# plot(b1/sum(b1))
# points(b2/sum(b2), col = 2, type = "p")
# points(b3/sum(b3), col = 3, type = "p")
# 
# smoothScatter(mcols(hgIndel.gr[queryHits(ol)])$supGenoNo,mcols(mmIndel.gr[subjectHits(ol)])$supGenoNo )
# 
# 
# # resolve by counting votes
# 
# length(ol)




#hgHot <- hgHot[hgHot$hotspotID %in% inter$domains[inter$genome == "hg19"],]


olIndel <- (findOverlaps(hgHot.gr, hgIndel.gr))
dfIndel <- data.frame(queryHits = queryHits(olIndel), mcols(hgIndel.gr[subjectHits(olIndel)])[c("gapWidth","indel")] )
dfIndel <- summarise(group_by(dfIndel,queryHits, indel = indel), gapWidth = sum(gapWidth))

dfIndel$genome = "hg19"

olMindel <- findOverlaps(hgHot.gr, mcols(mmIndel.gr)$queRange)
dfMindel <- data.frame(queryHits = queryHits(olMindel), mcols(mmIndel.gr[subjectHits(olMindel)])[c("gapWidth","indel")])
dfMindel <- summarise(group_by(dfMindel, queryHits, indel = indel), gapWidth = sum(gapWidth))

dfMindel$genome = "mm10"

dfIndel <- rbind(dfIndel, dfMindel)


olBr <- findOverlaps(hgHot.gr, hgBr.gr)
dfBr <- data.frame(queryHits = queryHits(olBr), baseRate = width(hgBr.gr[subjectHits(olBr)]))
dfBr <- summarise(group_by(dfBr, queryHits), baseRate = sum(baseRate))
mer <- merge(dfIndel, dfBr)

#olGap <- findOverlaps(hgHot.gr, hgGap.gr)
#dfGap <- data.frame(queryHits = queryHits(olGap), mcols(hgGap.gr[subjectHits(olGap)])[c("inDel","gapWidth")])
#dfGap <- summarise(group_by(dfGap, queryHits, GapIndel = inDel), HMgapWidth = sum(gapWidth))

#mer <- merge(mer, dfGap)


dfAll <- data.frame(mcols(hgHot.gr[mer$queryHits]), mer)


s <- summarise(group_by(dfAll, repGroup, hotspotID, genome.1, conState, indel), 
               gapWidth = sum(gapWidth), baseRate = sum(baseRate), known = sum(known))

s <- s[s$conState!="mid",]

s$conState <- factor(s$conState, levels = c("con","dif"))

# s <- filter(s,(genome.1 == "hg19" & indel == "del" & GapIndel == "del") |
#               (genome.1 == "hg19" & indel == "ins" & GapIndel == "ins") |
#               (genome.1 == "mm10" & indel == "ins" & GapIndel == "ins") |
#               (genome.1 == "mm10" & indel == "del" & GapIndel == "del"))


pdf(file = "~/Desktop/RTN_domains/plots/inDelIdentify/HG19indelRates.pdf", width = 12, height = 8, onefile = TRUE)
par(mar = c(10,5,4,4))
stripchart((s$gapWidth/s$baseRate) ~ + s$repGroup + s$indel + s$genome.1,
           method= "jitter",jitter = .3, pch = 16, cex = .3, vert = TRUE, 
           las = 2, main = "hg19 all repeat enriched regions", xaxs = "i",
           col = c("darkblue", "purple", "aquamarine3", "red"), log = "y")
boxplot((s$gapWidth/s$baseRate) ~   s$repGroup  + s$indel + s$genome.1,
        las = 2, main = "", xaxs = "i", outline = FALSE, add = TRUE, col = NA, 
        border = c("darkblue", "purple", "aquamarine3", "red"), log = "y")
for(i in 0:16){abline(v = (i)+ .5, lty = 2, lwd = 1)};for(i in 0:4){abline(v = (i * 4) + .5, lty = 1, lwd = 2)};abline(v = 8.5, lty = 1, lwd = 3, col= 2)



stripchart(log10(s$gapWidth/s$baseRate) ~ s$conState + s$repGroup  + s$indel + s$genome.1,
           method= "jitter",jitter = .3, pch = 16, cex = .3, vert = TRUE, 
           las = 2, main = "hg19 high/low in query", xaxs = "i",
           col = c("darkblue", "darkblue","purple", "purple","aquamarine3", "aquamarine3","red", "red"))
boxplot(log10(s$gapWidth/s$baseRate) ~  s$conState + s$repGroup  + s$indel + s$genome.1,
        las = 2, main = "", xaxs = "i", outline = FALSE, add = TRUE, col = NA, 
        border =  c("darkblue", "darkblue","purple", "purple","aquamarine3", "aquamarine3","red", "red"))
for(i in 0:16){abline(v = (i * 2) + .5, lty = 2)};for(i in 0:4){abline(v = (i * 8) + .5, lty = 1, lwd = 2)};abline(v = 16.5, lty = 1, lwd = 3, col= 2)




r1 <- s[s$genome.1 == "hg19" & s$hotspotID %in% inter$domains[inter$genome == "hg19"],]
r2 <- s[s$genome.1 == "mm10" & s$hotspotID %in% inter$domains[inter$genome == "mm10"],]

r <- rbind(r1,r2)

stripchart(log10(r$gapWidth/r$baseRate) ~ r$conState + r$repGroup  + r$indel + r$genome.1,
           method= "jitter",jitter = .3, pch = 16, cex = .3, vert = TRUE, 
           las = 2, main = "hg19 shared/lineage specific", xaxs = "i",
           col=  c("darkblue", "darkblue","purple", "purple","aquamarine3", "aquamarine3","red", "red"))
boxplot(log10(r$gapWidth/r$baseRate) ~  r$conState + r$repGroup  + r$indel + r$genome.1,
        las = 2, main = "", xaxs = "i", outline = FALSE, add = TRUE, col = NA, 
        border =  c("darkblue", "darkblue","purple", "purple","aquamarine3", "aquamarine3","red", "red"))
for(i in 0:16){abline(v = (i * 2) + .5, lty = 2)};for(i in 0:4){abline(v = (i * 8) + .5, lty = 1, lwd = 2)};abline(v = 16.5, lty = 1, lwd = 3, col= 2)


dev.off()



head(inter)




# relative contribution of insertions and deleations is similar


# we are trying to look at associated factors of what the gaps could be due to
# proportion of human sided gaps and proportion of mouse sided gaps 

# so this will help us identify particular classes of TE assocaited evolution 
# so really there is little variation between regions 


# it might mean TEs ahve little impact



# if we consider it as a proportion of gapped sequence. 








olIndel <- (findOverlaps(mmHot.gr, mmIndel.gr))
dfIndel <- data.frame(queryHits = queryHits(olIndel), mcols(mmIndel.gr[subjectHits(olIndel)])[c("gapWidth","indel")] )
dfIndel <- summarise(group_by(dfIndel,queryHits, indel = indel), gapWidth = sum(gapWidth))

dfIndel$genome = "mm10"

olMindel <- findOverlaps(mmHot.gr, mcols(hgIndel.gr)$queRange)
dfMindel <- data.frame(queryHits = queryHits(olMindel), mcols(hgIndel.gr[subjectHits(olMindel)])[c("gapWidth","indel")])
dfMindel <- summarise(group_by(dfMindel, queryHits, indel = indel), gapWidth = sum(gapWidth))

dfMindel$genome = "hg19"

dfIndel <- rbind(dfIndel, dfMindel)


olBr <- findOverlaps(mmHot.gr, mmBr.gr)
dfBr <- data.frame(queryHits = queryHits(olBr), baseRate = width(mmBr.gr[subjectHits(olBr)]))
dfBr <- summarise(group_by(dfBr, queryHits), baseRate = sum(baseRate))
mer <- merge(dfIndel, dfBr)

#olGap <- findOverlaps(mmHot.gr, mmGap.gr)
#dfGap <- data.frame(queryHits = queryHits(olGap), mcols(mmGap.gr[subjectHits(olGap)])[c("inDel","gapWidth")])
#dfGap <- summarise(group_by(dfGap, queryHits, GapIndel = inDel), HMgapWidth = sum(gapWidth))

#mer <- merge(mer, dfGap)


dfAll <- data.frame(mcols(mmHot.gr[mer$queryHits]), mer)


s <- summarise(group_by(dfAll, repGroup, hotspotID, genome.1, conState, indel), 
               gapWidth = sum(gapWidth), baseRate = sum(baseRate), known = sum(known))

s <- s[s$conState!="mid",]

s$conState <- factor(s$conState, levels = c("con","dif"))

# s <- filter(s,(genome.1 == "hg19" & indel == "del" & GapIndel == "del") |
#               (genome.1 == "hg19" & indel == "ins" & GapIndel == "ins") |
#               (genome.1 == "mm10" & indel == "ins" & GapIndel == "ins") |
#               (genome.1 == "mm10" & indel == "del" & GapIndel == "del"))



pdf(file = "~/Desktop/RTN_domains/plots/inDelIdentify/MM10indelRate.pdf", width = 12, height = 8, onefile = TRUE)
par(mar = c(10,4,4,4))

stripchart((s$gapWidth/s$baseRate) ~ s$repGroup + s$indel + s$genome.1,
           method= "jitter",jitter = .3, pch = 16, cex = .3, vert = TRUE, 
           las = 2, main = "mm10 all repeat enriched regions", xaxs = "i",
           col = c("darkblue", "purple", "aquamarine3", "red"), ylim = c(0,.1))
boxplot((s$gapWidth/s$baseRate) ~   s$repGroup  + s$indel + s$genome.1,
        las = 2, main = "", xaxs = "i", outline = FALSE, add = TRUE, col = NA,
        border = c("darkblue", "purple", "aquamarine3", "red"))
for(i in 0:16){abline(v = (i)+ .5, lty = 2, lwd = 1)};for(i in 0:4){abline(v = (i * 4) + .5, lty = 1, lwd = 2)};abline(v = 8.5, lty = 1, lwd = 3, col= 2)


stripchart(log10(s$gapWidth/s$baseRate) ~ s$conState + s$repGroup  + s$indel + s$genome.1,
           method= "jitter",jitter = .3, pch = 16, cex = .3, vert = TRUE, 
           las = 2, main = "mm10 high/low in query", xaxs = "i",
           col = c("darkblue","darkblue", "purple", "purple", "aquamarine3", "aquamarine3", "red", "red"))
boxplot(log10(s$gapWidth/s$baseRate) ~ s$conState + s$repGroup  + s$indel + s$genome.1,
        las = 2, main = "", xaxs = "i", outline = FALSE, add = TRUE, col = NA, 
        border = c("darkblue","darkblue", "purple", "purple", "aquamarine3", "aquamarine3", "red", "red"))
for(i in 0:16){abline(v = (i * 2) + .5, lty = 2)};for(i in 0:4){abline(v = (i * 8) + .5, lty = 1, lwd = 2)};abline(v = 16.5, lty = 1, lwd = 3, col= 2)



r1 <- s[s$genome.1 == "hg19" & s$hotspotID %in% inter$domains[inter$genome == "hg19"],]
r2 <- s[s$genome.1 == "mm10" & s$hotspotID %in% inter$domains[inter$genome == "mm10"],]

r <- rbind(r1,r2)



stripchart(log10(r$gapWidth/r$baseRate) ~ r$conState + r$repGroup  + r$indel + r$genome.1,
           method= "jitter",jitter = .3, pch = 16, cex = .3, vert = TRUE, 
           las = 2, main = "mm10 shared/lineage specific", xaxs = "i", add = FALSE,
           col = c("darkblue","darkblue", "purple", "purple", "aquamarine3", "aquamarine3", "red", "red"))

boxplot(log10(r$gapWidth/r$baseRate) ~ r$conState + r$repGroup  + r$indel + r$genome.1,
        las = 2, main = "", xaxs = "i", outline = FALSE, add = TRUE, col = NA, border = c("darkblue","darkblue", "purple", "purple", "aquamarine3", "aquamarine3", "red", "red"))

for(i in 0:16){abline(v = (i * 2) + .5, lty = 2)};for(i in 0:4){abline(v = (i * 8) + .5, lty = 1, lwd = 2)};abline(v = 16.5, lty = 1, lwd = 3, col= 2)

dev.off()



# so now i have the full details of hat happned in human and mouse in human regions. 
# that might mean no more awkward matching up


### overlap of human delations in mouse returned nothing 
### I wonder why










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


