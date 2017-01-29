library(GenomicRanges)
library(devtools)
library(dplyr)
setwd("~/Desktop/RTN_domains/")

options(stringsAsFactors = TRUE)

rm(list = ls())

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")



s1name <- "hg19"
s2name <- "canFam3"

repGroups <- c("ancient", "new_SINE", "new_L1", "old_L1")
repCols <- c("darkblue", "aquamarine3", "purple", "red")


load(file = paste("R_objects/rmskMapTables/",s1name,"/repData_",s1name,"_50000.RData", sep = ""))
s1DataList <- repDataList


load(file = paste("R_objects/rmskMapTables/",s2name,"/repData_",s2name,"_50000.RData", sep = ""))
s2DataList <- repDataList


datS1.gr <- GRanges(seqnames = Rle(s1DataList$bin$chr),
                    ranges = IRanges(start = s1DataList$bin$start, end = s1DataList$bin$end),
                    binID = s1DataList$repSummary$binID,
                    known = s1DataList$bin$Known
)
for(i in 1:length(repGroups)){
  elementMetadata(datS1.gr)[,repGroups[i]] <- ((s1DataList$repSummary[,repGroups[i]]/s1DataList$bin$Known)*50000)
}

datS2.gr <- GRanges(seqnames = Rle(s2DataList$bin$chr),
                    ranges = IRanges(start = s2DataList$bin$start, end = s2DataList$bin$end),
                    binID = s2DataList$repSummary$binID,
                    known = s2DataList$bin$Known
)
for(i in 1:length(repGroups)){
  elementMetadata(datS2.gr)[,repGroups[i]] <- ((s2DataList$repSummary[,repGroups[i]]/s2DataList$bin$Known)*50000)
}



s1 <- read.table(paste("data/repeatHotspot/",s1name,"/",s1name,"Hotspots.bed",sep = "" ), 
                 col.names = c("chr", "start", "end", "domainID"),
                 colClasses = c("character", "integer", "integer", "character"))

s1_s2 <- read.table(paste("data/repeatHotspot/",s1name,"/","hg19_lift_canFam3.bed",sep = "" ), 
                    col.names = c("chr", "start", "end", "domainID"),
                    colClasses = c("character", "integer", "integer", "character"))


s2 <- read.table(paste("data/repeatHotspot/",s2name,"/",s2name,"Hotspots.bed",sep = "" ), 
                 col.names = c("chr", "start", "end", "domainID"),
                 colClasses = c("character", "integer", "integer", "character"))

s2_s1 <- read.table(paste("data/repeatHotspot/",s2name,"/","canFam3_lift_hg19.bed",sep = "" ), 
                    col.names = c("chr", "start", "end", "domainID"),
                    colClasses = c("character", "integer", "integer", "character"))


dataSets <- c("s1", "s1_s2", "s2", "s2_s1")

for(d in dataSets){
  data <- get(d)
  dataAdomain <- t(
  as.data.frame(
    strsplit(
      data$domainID[grep("ancient", data$domainID)],
      "_"
      )
    )
  )

dataNdomain <- t(
  as.data.frame(
    strsplit(
      data$domainID[-grep("ancient", data$domainID)],
      "_"
    )
  )
)

data$repGroup <- c(dataAdomain[,1], paste(dataNdomain[,1], dataNdomain[,2], sep = "_"))
data$hotspotID <- c(dataAdomain[,2], dataNdomain[,3])
data$hotspotID <- paste(data$repGroup, data$hotspotID, sep = ";")
data$hotspotGroup <- c(dataAdomain[,3], dataNdomain[,4])
data$hotspotGroup <- paste(data$repGroup, data$hotspotGroup, sep = ";")


datagr <- GRanges(seqnames = Rle(data$chr),
                     ranges = IRanges(start = data$start, end = data$end),
                     domainID = data$domainID,
                     repGroup = data$repGroup,
                     hotspotID = data$hotspotID,
                     hotspotGroup = data$hotspotGroup)
assign(x = d, data)
assign(x = paste(d,".gr", sep = ""),value = datagr)
}



# be careful not to use s2_s1 as a proxy for s2
# we can get our s2 meausres from s2 directly

pdf(file = paste("plots/hotspotOverlap/", s2name, "_",s1name,"_overlap.pdf", sep = ""), height = 12, width = 6, onefile = TRUE)
layout(c(1,2))
barplot(olHotspotSummary(s2.gr, s1_s2.gr, repGroups), beside = TRUE, col = repCols, ylim = c(0,1),
        xlab = "reference repeat groups", ylab = "overlaping hotspots (proportion of reference)",
        main = paste(s2name, "overlapping", s1name, "repeat hotspots"))
legend("topright", legend = repGroups, fill = repCols, title = "query repeat groups", bty = "n",cex = .75)

barplot(olHotspotSummary(s2.gr, s2.gr, repGroups), beside = TRUE, col = repCols, ylim = c(0,1),
        xlab = "reference repeat groups", ylab = "overlaping hotspots (proportion of reference)",
        main = paste(s2name, "self overlap"))

barplot(olHotspotSummary(s1.gr, s2_s1.gr, repGroups), beside = TRUE, col = repCols, ylim = c(0,1),
        xlab = "reference repeat groups", ylab = "overlaping hotspots (proportion of reference)",
        main = paste(s1name, "overlapping", s2name, "repeat hotspots"))
legend("topright", legend = repGroups, fill = repCols, title = "query repeat groups", bty = "n", cex = .75)

barplot(olHotspotSummary(s1.gr, s1.gr, repGroups), beside = TRUE, col = repCols, ylim = c(0,1),
        xlab = "reference repeat groups", ylab = "overlaping hotspots (proportion of reference)",
        main = paste(s1name, "self overlap"))
dev.off()



# maybe get a function to pull regions. 

# so groups should be maintained in the ref
# pull out groups that have a high proportion of self overlap


# can now pull out corrsponding groups now its time to turn it into a function
# also it helps to supply scaled data

s1Corresponding <- extractCorrespondingHotspots(ref.gr = s1.gr,que_ref.gr = s2_s1.gr, ref_que.gr = s1_s2.gr,que.gr = s2.gr,repGroups = repGroups)

s2Corresponding <- extractCorrespondingHotspots(ref.gr = s2.gr,que_ref.gr = s1_s2.gr, ref_que.gr = s2_s1.gr,que.gr = s1.gr,repGroups = repGroups)




s1_s2_insertionRate <- extractInsertionRates(refGenome.gr = datS1.gr, queGenome.gr = datS2.gr, 
                                             refCorresponding = s1Corresponding, repGroups = repGroups, 
                                             minoverlap = 10e3)
s2_s1_insertionRate <- extractInsertionRates(refGenome.gr = datS2.gr, queGenome.gr = datS1.gr, 
                                             refCorresponding = s2Corresponding, repGroups = repGroups, 
                                             minoverlap = 10e3)


pdf(file = paste("plots/hotspotOverlap/", s1name,"_",s2name,"_insertionRates.pdf", sep = ""),
    width = 12,height = 12, onefile = TRUE)

layout(matrix(c(1,2,3,4), nrow = 2, byrow = T))
par(mar = c(3,3,0,0), oma = c(5,5,5,5))
for(i in 1:length(repGroups)){
  dat <- filter(s1_s2_insertionRate, repGroup == repGroups[i] & (conState == "con" | conState == "dif"))
  dat$conState <- droplevels(dat$conState)
  dat1 <- summarise(group_by(dat, hotspotID, repGroup, conState, genome), mean(insertionRate))
  boxplot(`mean(insertionRate)` ~ genome + conState, data = dat1, las = 2, outline = FALSE, 
          notch = FALSE, xaxt = "n", col = c("grey50","grey90"), ylim = c(0, max(dat1$`mean(insertionRate)`)))
  axis(side = 1,at = c(1.5,3.5), labels = c("conserved hotspots", "non-conserved hotspots"))
  legend("topright", legend = repGroups[i], bty = "n")
  
  abline(h = mean(elementMetadata(datS1.gr)[[repGroups[i]]], na.rm = T), lwd = 2)
  abline(h = mean(elementMetadata(datS2.gr)[[repGroups[i]]], na.rm = T), lwd = 2, col = "grey80")
  
  boxplot(`mean(insertionRate)` ~ genome + conState, data = dat1, las = 2, outline = FALSE, 
          notch = FALSE, xaxt = "n", col = c("grey50","grey90"), add = T)
  
  dat2 <- summarise(group_by(dat, genome, conState, hotspotGroup, repGroup), mean(insertionRate), n())
  stripchart(`mean(insertionRate)` ~ genome + conState, data = dat2, add = TRUE, vert = TRUE,
             method = "jitter", jitter = .35, pch = 16, cex = sqrt(dat2$`n()`/10), col = repCols[i])
  
}
mtext("Insertion rate (Z score)", side = 2, line = 2.5, outer = T)
title(main = paste(s1name,"hotspots mapped to" ,s2name),outer = TRUE)

legend("topleft",legend = c("100kb", "500kb", "1Mb"), pch = c(16,16,16), pt.cex = sqrt(c(2/10,10/10,20/10)),
       bty = "n", title = "hotspot size")
dev.off()


pdf(file = paste("plots/hotspotOverlap/", s2name,"_",s1name,"_insertionRates.pdf", sep = ""), 
    width = 12,height = 12, onefile = TRUE)


layout(matrix(c(1,2,3,4), nrow = 2, byrow = T))
par(mar = c(3,3,0,0), oma = c(5,5,5,5))
for(i in 1:length(repGroups)){
  dat <- filter(s2_s1_insertionRate, repGroup == repGroups[i] & (conState == "con" | conState == "dif"))
  dat$conState <- droplevels(dat$conState)
  dat1 <- summarise(group_by(dat, hotspotID, repGroup, conState, genome), mean(insertionRate))
  
  boxplot(`mean(insertionRate)` ~ genome + conState, data = dat1, las = 2, outline = FALSE, 
          notch = FALSE, xaxt = "n", col = c("grey50","grey90"), ylim = c(0, max(dat1$`mean(insertionRate)`)))
  axis(side = 1,at = c(1.5,3.5), labels = c("conserved hotspots", "non-conserved hotspots"))
  legend("topright", legend = repGroups[i], bty = "n")
  
  abline(h = mean(elementMetadata(datS2.gr)[[repGroups[i]]], na.rm = T), lwd = 2)
  abline(h = mean(elementMetadata(datS1.gr)[[repGroups[i]]], na.rm = T), lwd = 2, col = "grey80")
  
  boxplot(`mean(insertionRate)` ~ genome + conState, data = dat1, las = 2, outline = FALSE, 
          notch = FALSE, xaxt = "n", col = c("grey50","grey90"), add= T)
  
  dat2 <- summarise(group_by(dat, genome, conState, hotspotGroup, repGroup), mean(insertionRate), n())
  stripchart(`mean(insertionRate)` ~ genome + conState, data = dat2, add = TRUE, vert = TRUE,
             method = "jitter", jitter = .35, pch = 16, cex = sqrt(dat2$`n()`/10), col = repCols[i])
  
  
}
mtext("Insertion rate (Z score)", side = 2, line = 2.5, outer = T)
title(main = paste(s2name,"hotspots mapped to" ,s1name),outer = TRUE)

legend("topleft",legend = c("100kb", "500kb", "1Mb"), pch = c(16,16,16), pt.cex = sqrt(c(2/10,10/10,20/10)), 
       bty = "n", title = "hotspot size")

dev.off()



write.table(x = s1_s2_insertionRate, file = paste("data/repeatHotspot/",s1name,"/", s1name, "_", s2name, "_conDif.txt", sep = ""),
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

write.table(x = s2_s1_insertionRate, file = paste("data/repeatHotspot/",s2name ,"/", s2name, "_", s1name, "_conDif.txt", sep = ""),
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)



