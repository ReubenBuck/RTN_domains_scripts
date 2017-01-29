rm(list = ls())

setwd("~/Desktop/RTN_domains/")


s1name = "hg19"
s2name = "canFam3"


data0 <- read.table(file = paste("data/repeatHotspot/", s1name, "/", s1name,"_",s2name,"_conDif.txt", sep = ""), header = T)


load(paste("R_objects/rmskMapTables/binSizes/", s1name, "/repData_", s1name, "_50000.RData", sep = ""))
s1Dat <- repDataList

load(paste("R_objects/rmskMapTables/binSizes/", s2name, "/repData_", s2name, "_50000_size.RData", sep = ""))
s2Dat <- repDataList



# so what are we going to test between these regions
# levels of conservation 
# density of conserved elements

# test conservation levels of orthologous genes.




head(data)

ce <- read.table("data/genomeAno/hg19/phastConsElements46way.txt")

ce.gr <- GRanges(seqnames = Rle(ce[,2]), 
                 ranges = IRanges(start = ce[,3], end = ce[,4]))



data1 = data0[data0$genome == "ref",]

data1.gr <- GRanges(seqnames = Rle(as.character(data1$chr)),
                   ranges = IRanges(start = data1$start, end = data1$end))



ol <- as.matrix(findOverlaps(data1.gr, ce.gr))



data3 <- data.frame(data1[ol[,1],], ce = width(ce.gr[ol[,2]]))


#data3 <- data2[data2$repGroup == "new_SINE",]

df <- summarise(group_by(data3, hotspotGroup, conState, hotspotID, repGroup), sum(ce))

layout(1)
par(mar = c(10,5,5,5))
boxplot((df$`sum(ce)`) ~ df$conState+ df$repGroup, las = 2, outline = F)
stripchart((df$`sum(ce)`) ~ df$conState + df$repGroup, method = "jitter", pch = 16, cex = .3, vert = T, add = T, jitter = .35)


# look at conserved sequence in noncoding regions





library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


olexon <- as.matrix(findOverlaps(ce.gr, exons(txdb))) 

nonCodingCe.gr <- ce.gr[-unique(olexon[,1])]



ol <- as.matrix(findOverlaps(data1.gr, nonCodingCe.gr))

data3 <- data.frame(data1[ol[,1],], ce = width(ce.gr[ol[,2]]))

df <- summarise(group_by(data3, hotspotGroup, conState, hotspotID, repGroup), sum(ce))
df <- df[df$conState != "mid",]
df$conState <- droplevels(df$conState)

df1 <- summarise(group_by(df, hotspotGroup, conState, repGroup), mean = mean(`sum(ce)`), n())


layout(1)
par(mar = c(10,5,5,5))
boxplot((df$`sum(ce)`) ~ df$conState + df$repGroup , las = 2, outline = F, notch  = T, ylim = c(1.5, (max(df1$mean))))
stripchart((df1$mean) ~ df1$conState + df1$repGroup, 
           pch = 16, cex = sqrt(df1$`n()`/5), vert = T, add = T, jitter = .35)









