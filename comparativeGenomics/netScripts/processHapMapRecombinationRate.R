

# human hotspots

fileList <- list.files("~/Desktop/RTN_domains/data/UCSCtracks/hg19/genetic_map_HapMapII_GRCh37/")

fileList <- fileList[grep("genetic_map", fileList)]

AllRates <- NULL
for(i in fileList){
  r <- read.table(paste("~/Desktop/RTN_domains/data/UCSCtracks/hg19/genetic_map_HapMapII_GRCh37/",i, sep = ""), header = TRUE)
  AllRates <- rbind(AllRates, r)
  }

head(AllRates)

table(AllRates[grep("_", AllRates$Chromosome),"Chromosome"])

AllRates[grep("_", AllRates$Chromosome),]$Chromosome <- "chrX"

RecombinationRate <- AllRates

save(RecombinationRate, file = "~/Desktop/RTN_domains/data/UCSCtracks/hg19/hg19.recombinationRate.RData")




# mouse data 

rate <- read.csv("~/Desktop/RTN_domains/data/UCSCtracks/mm9/recombinationHotspot/mm9_recombRate.csv")

rate$Kb = rate$Kb*1000

rate <- data.frame(seqnames = rate$Chr, start = rate$Kb, end = rate$Kb+1, rate = rate$X4Ner.kb)
write.table(rate, file = "~/Desktop/RTN_domains/data/UCSCtracks/mm9/recombinationHotspot/mm9_recombRate.bed", 
            sep = "\t",quote = FALSE,row.names = FALSE, col.names = FALSE)


# hotspot
rate <- read.csv("~/Desktop/RTN_domains/data/UCSCtracks/mm9/recombinationHotspot/mm9_recombHotspot.csv")

rate$Chr <- paste("chr",rate$Chr, sep = "")

write.table(rate, file = "~/Desktop/RTN_domains/data/UCSCtracks/mm9/recombinationHotspot/mm9_recombHotSpot.bed", 
            sep = "\t",quote = FALSE,row.names = FALSE, col.names = FALSE)


##### peaks


fileList <- list.files("~/Desktop/RTN_domains/data/mm9DNAse1Peaks/")

AllRates <- NULL
for(i in fileList){
  
  r <- read.table(paste("~/Desktop/RTN_domains/data/mm9DNAse1Peaks/",i, sep = ""), header = FALSE)
  
  AllRates <- rbind(AllRates, r)
  
}

AllRates.gr <- GRanges(seqnames = AllRates[,1],
                       ranges = IRanges(start = AllRates[,2], end = AllRates[,3]))


redRates.gr <- reduce(AllRates.gr)

redRates.gr$activity <- countOverlaps(redRates.gr, AllRates.gr)

dfRate <- data.frame(redRates.gr)


write.table(dfRate, file = "~/Desktop/RTN_domains/data/UCSCtracks/mm9/DHS_peaks/mm9_dnase1cluster.bed", 
            sep = "\t",col.names = F,row.names = F,quote = F)
