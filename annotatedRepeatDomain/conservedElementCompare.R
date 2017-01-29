rm(list = ls())

setwd("~/Desktop/RTN_domains/")


s1name = "hg19"
s2name = "mm9"
s3name = "canFam3"

datas2 <- read.table(file = paste("data/repeatHotspot/", s1name, "/", s1name,"_",s2name,"_conDif.txt", sep = ""), header = T)
datas3 <- read.table(file = paste("data/repeatHotspot/", s1name, "/", s1name,"_",s3name,"_conDif.txt", sep = ""), header = T)


cefile <- list.files(paste("data/genomeAno/consElement/", s1name, sep = ""))
ce <- read.table(paste("data/genomeAno/consElement/", s1name, "/", cefile, sep = ""))
ce.gr <- GRanges(seqnames = Rle(ce[,2]), 
                 ranges = IRanges(start = ce[,3], end = ce[,4]))

if(s1name == "mm9"){
  library("TxDb.Mmusculus.UCSC.mm9.knownGene")
  txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
}else if(s1name == "hg19"){
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
}else{
  print("no txdb object")
}


olexon <- as.matrix(findOverlaps(ce.gr, exons(txdb))) 
nonCodingCe.gr <- ce.gr[-unique(olexon[,1])]

conInts <- unique(intersect(as.character(datas2$hotspotID[datas2$genome == "ref" & datas2$conState == "con" ]), 
                            as.character(datas3$hotspotID[datas3$genome == "ref" & datas3$conState == "con" ])
)
)
difInts <- unique(intersect(as.character(datas2$hotspotID[datas2$genome == "ref" & datas2$conState == "dif"]), 
                            as.character(datas3$hotspotID[datas3$genome == "ref" & datas3$conState == "dif"])
)
)

ints <- c(conInts, difInts)

data0 <- datas2


data1 = data0[data0$genome == "ref" & (data0$conState == "con" | data0$conState == "dif"),]

#data1 = data1[(as.character(data1$hotspotID) %in% ints),]

data1.gr <- GRanges(seqnames = Rle(as.character(data1$chr)),
                   ranges = IRanges(start = data1$start, end = data1$end))


oldomain <- as.matrix(findOverlaps(data1.gr, exons(txdb)))
a <- aggregate(width(GenomicRanges::pintersect(data1.gr[oldomain[,1]], exons(txdb)[oldomain[,2]])), by = list(oldomain[,1]), FUN = sum)

data1$known[a$Group.1] = data1$known[a$Group.1] - a$x


ol <- as.matrix(findOverlaps(data1.gr, nonCodingCe.gr))

data3 <- data.frame(data1[ol[,1],], ce = width(ce.gr[ol[,2]]))

df <- summarise(group_by(data3, hotspotGroup, conState, hotspotID, repGroup, known), sum(ce))
df <- df[df$conState != "mid",]
df$conState <- droplevels(df$conState)
df$ceRate <- df$`sum(ce)`/df$known

df1 <- summarise(group_by(df, hotspotGroup, conState, repGroup), mean = mean(ceRate), n())


layout(1)
par(mar = c(10,5,5,5))
boxplot((df$ceRate) ~ df$conState + df$repGroup , las = 2, outline = T, notch  = T, log = "y")
stripchart((df1$mean) ~ df1$conState + df1$repGroup, method = "jitter",
           pch = 16, cex = sqrt(df1$`n()`/5), 
           vert = T, add = T, jitter = .35)


# not really seeing any differences between conserved and non conserved

# the real question would be what are genes doing.

# how to find the impact 

# maybe it is better to look at individual phylop scores


# There is a blind spot in our analysis. 
# it's what is happening in the other species.

# in the referecne everything is pretty similar



plot(density((df$ceRate[df$repGroup == "ancient" & df$conState == "con"])))
lines(density((df$ceRate[df$repGroup == "ancient" & df$conState == "dif"])),col = 2)


a <- summarise(group_by(data1, repGroup, conState),n())
a[,3] <- a[,3]/20
