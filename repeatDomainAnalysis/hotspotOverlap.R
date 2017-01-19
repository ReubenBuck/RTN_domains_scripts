library(GenomicRanges)
library(devtools)
library(dplyr)
setwd("~/Desktop/RTN_domains/")

options(stringsAsFactors = TRUE)

rm(list = ls())

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")


# might be worth also downloadign the seq info for each species.
# this will stop the false positive warning messages.


# seperate out IDs

# so we want to look at exclusive groups
# Those with minimal overlap
# also domains that are large. 

s1name <- "canFam3"
s2name <- "mm9"

repGroups <- c("ancient", "new_SINE", "new_L1", "old_L1")
repCols <- c("darkblue", "aquamarine3", "purple", "red")


load(file = paste("R_objects/rmskMapTables/",s1name,"/repData_",s1name,"_50000.RData", sep = ""))
s1DataList <- repDataList
s1bin.gr <- GRanges(seqnames = Rle(s1DataList$bin$chr),
                    ranges = IRanges(start = s1DataList$bin$start, end = s1DataList$bin$end),
                    binID = s1DataList$bin$binID)

load(file = paste("R_objects/rmskMapTables/",s2name,"/repData_",s2name,"_50000.RData", sep = ""))
s2DataList <- repDataList
s2bin.gr <- GRanges(seqnames = Rle(s2DataList$bin$chr),
                    ranges = IRanges(start = s2DataList$bin$start, end = s2DataList$bin$end),
                    binID = s2DataList$bin$binID)


s1 <- read.table("data/repeatHotspot/ canFam3 Hotspots.bed", 
                 col.names = c("chr", "start", "end", "domainID"),
                 colClasses = c("character", "integer", "integer", "character"))

s1_s2 <- read.table("data/repeatHotspot/canFam3_lift_mm10_lift_mm9.bed", 
                    col.names = c("chr", "start", "end", "domainID"),
                    colClasses = c("character", "integer", "integer", "character"))


s2 <- read.table("data/repeatHotspot/ mm9 Hotspots.bed", 
                 col.names = c("chr", "start", "end", "domainID"),
                 colClasses = c("character", "integer", "integer", "character"))

s2_s1 <- read.table("data/repeatHotspot/mm9_lift_mm10_lift_canFam3.bed", 
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

# now the group information is more accesible


s1.gr <- s1.gr[elementMetadata(s1.gr)$domainID  %in% elementMetadata(s1_s2.gr)$domainID]
s2.gr <- s2.gr[elementMetadata(s2.gr)$domainID  %in% elementMetadata(s2_s1.gr)$domainID]


sum(width(s1.gr))
sum(width(s1_s2.gr))

#hist(width(s1_s2.gr), breaks = 1000, xlim = c(0,2e5))
#hist(width(s2_s1.gr), breaks = 1000, xlim = c(0,2e5))

par(mar = c(5,5,5,5), oma = c(0,0,0,0))
plot(density(width(s1_s2.gr)), xlim = c(0,2e5), main = "width of lifted hotspots", xlab = "width (bp)", lwd = 3)
lines(density(width(s2_s1.gr)), col = 2, lwd = 3)
abline(v = 5e4, col = 4, lwd = 3)
legend("topright", legend = c(paste(s1name, "to", s2name), paste(s2name, "to", s1name)), col= c(1,2), lwd = 3)

# have got data in the form I want now It is time to begin analyzing the degree of conservation between hotspots

# when we've got our groups 
# get the actuall repeat insertion rate 
# make sure it is significantly different


# so we could merge

# we could see if the groups were maintained after the self alignment.


# be careful not to use s2_s1 as a proxy for s2
# we can get our s2 meausres from s2 directly

pdf(file = paste("plots/hotspotOverlap/", s2name, "_",s1name,"_overlap.pdf"), height = 12, width = 6, onefile = TRUE)
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

s1Corresponding <- extractCorrespondingHotspots(ref.gr = s1.gr,que_ref.gr = s2_s1.gr, ref_que.gr = s1_s2.gr,repGroups = repGroups)

s2Corresponding <- extractCorrespondingHotspots(ref.gr = s2.gr,que_ref.gr = s1_s2.gr, ref_que.gr = s2_s1.gr,repGroups = repGroups)



hist(log10(width(s1Corresponding$que$dif$new_SINE)), breaks = 100); abline(v = 5e4)
hist(log10(width(s1Corresponding$que$con$new_SINE)), breaks = 100); abline(v = 5e4)

boxplot(width(s2Corresponding$que$con$new_SINE),
        width(s2Corresponding$que$dif$new_SINE), 
        outline = T, names = c("con", "dif"),
        ylab = "width (bp)")


# how to extract repeat content.

datS1 <- s1DataList
datS2 <- s2DataList

datS1.gr <- GRanges(seqnames = Rle(s1DataList$bin$chr),
                    ranges = IRanges(start = s1DataList$bin$start, end = s1DataList$bin$end),
                    binID = s1DataList$repSummary$binID
                    )
for(i in 1:length(repGroups)){
  elementMetadata(datS1.gr)[,repGroups[i]] <- scale((s1DataList$repSummary[,repGroups[i]]/s1DataList$bin$Known)*10000)
}

datS2.gr <- GRanges(seqnames = Rle(s2DataList$bin$chr),
                    ranges = IRanges(start = s2DataList$bin$start, end = s2DataList$bin$end),
                    binID = s2DataList$repSummary$binID
)
for(i in 1:length(repGroups)){
  elementMetadata(datS2.gr)[,repGroups[i]] <- scale((s2DataList$repSummary[,repGroups[i]]/s2DataList$bin$Known)*10000)
}




s1_s2_insertionRate <- extractInsertionRates(ref.gr = datS1.gr, que.gr = datS2.gr, refCorresponding = s1Corresponding, repGroups = repGroups)
s2_s1_insertionRate <- extractInsertionRates(ref.gr = datS2.gr, que.gr = datS1.gr, refCorresponding = s2Corresponding, repGroups = repGroups)


pdf(file = paste("plots/hotspotOverlap/", s1name,"_",s2name,"_insertionRates.pdf", sep = ""), width = 5,height = 10, onefile = TRUE)
layout(matrix(c(1,2,3,4), nrow = 4))
par(mar = c(3,0,0,0), oma = c(5,5,5,5))
for(i in 1:length(repGroups)){
  dat <- filter(s1_s2_insertionRate, repGroup == repGroups[i] & (conState == "con" | conState == "dif"))
  dat$conState <- droplevels(dat$conState)
  boxplot(Z_score ~ genome + conState, data = dat, las = 2, outline = FALSE, 
          notch = FALSE, xaxt = "n", col = c("white","grey80"), ylim = c(-2,7))
  axis(side = 1,at = c(1.5,3.5), labels = c("conserved hotspots", "differential hotspots"))
  legend("topright", legend = repGroups[i], bty = "n")
  stripchart(Z_score ~ genome + conState, data = dat, add = TRUE, vert = TRUE,
             method = "jitter", jitter = .35, pch = 16, cex = .3, col = repCols[i])
  
}
mtext("Insertion rate (Z score)", side = 2, line = 2.5, outer = T)
title(main = paste(s1name,"hotspots mapped to" ,s2name),outer = TRUE)


layout(matrix(c(1,2,3,4), nrow = 4))
par(mar = c(3,0,0,0), oma = c(5,5,5,5))
for(i in 1:length(repGroups)){
  dat <- filter(s2_s1_insertionRate, repGroup == repGroups[i] & (conState == "con" | conState == "dif"))
  dat$conState <- droplevels(dat$conState)
  boxplot(Z_score ~ genome + conState, data = dat, las = 2, outline = FALSE, 
          notch = FALSE, xaxt = "n", col = c("white","grey80"), ylim = c(-2,7))
  axis(side = 1,at = c(1.5,3.5), labels = c("conserved", "differential"))
  legend("topright", legend = repGroups[i], bty = "n")
  
  stripchart(Z_score ~ genome + conState, data = dat, add = TRUE, vert = TRUE,
             method = "jitter", jitter = .35, pch = 16, cex = .3, col = repCols[i])
  
}
mtext("Insertion rate (Z score)", side = 2, line = 2.5, outer = T)
title(main = paste(s2name,"hotspots mapped to" ,s1name),outer = TRUE)
dev.off()



# plan now is to make the results cleaner
# when we pull regions make sure we are getting the right ones and not the ones beside them.


# at least here there is something now to compare


# may be good to get percentile rank scores

# the problem is the whole apples and oranges thing

# get the reference query stuff sorted. 






