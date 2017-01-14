library(GenomicRanges)
library(devtools)
setwd("~/Desktop/RTN_domains/")

rm(list = ls())

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")




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

barplot(olHotspotSummary(s2.gr, s1_s2.gr, repGroups), beside = TRUE, col = repCols, ylim = c(0,1),
        xlab = "reference repeat groups", ylab = "overlaping hotspots (proportion of reference)",
        main = paste(s2name, "overlapping", s1name, "repeat hotspots"))
legend("topright", legend = repGroups, fill = repCols, title = "query repeat groups", bty = "n")

barplot(olHotspotSummary(s2.gr, s2.gr, repGroups), beside = TRUE, col = repCols, ylim = c(0,1),
        xlab = "reference repeat groups", ylab = "overlaping hotspots (proportion of reference)",
        main = paste(s2name, "self overlap"))
legend("topright", legend = repGroups, fill = repCols, title = "query repeat groups", bty = "n")

barplot(olHotspotSummary(s1.gr, s2_s1.gr, repGroups), beside = TRUE, col = repCols, ylim = c(0,1),
        xlab = "reference repeat groups", ylab = "overlaping hotspots (proportion of reference)",
        main = paste(s1name, "overlapping", s2name, "repeat hotspots"))
legend("topright", legend = repGroups, fill = repCols, title = "query repeat groups", bty = "n")

barplot(olHotspotSummary(s1.gr, s1.gr, repGroups), beside = TRUE, col = repCols, ylim = c(0,1),
        xlab = "reference repeat groups", ylab = "overlaping hotspots (proportion of reference)",
        main = paste(s1name, "self overlap"))
legend("topright", legend = repGroups, fill = repCols, title = "query repeat groups", bty = "n")



# maybe get a function to pull regions. 

# so groups should be maintained in the ref
# pull out groups that have a high proportion of self overlap



ref.gr <- s1.gr
que.gr <- s2_s1.gr

ol <- as.matrix(findOverlaps(ref.gr, que.gr))

dfRepGroup <- data.frame(elementMetadata(ref.gr)$repGroup[ol[,1]], 
                         elementMetadata(que.gr)$repGroup[ol[,2]])

#for(i in 1:length(repGroups))
  
i = 2
 
conHotspot.gr <- ref.gr[unique(ol[dfRepGroup[,1] == repGroups[i] & dfRepGroup[,1] == repGroups[i], 1])]

tabRef <- table(elementMetadata(ref.gr[elementMetadata(ref.gr)$repGroup == repGroups[i]])$hotspotGroup)
tabCon <- table(elementMetadata(conHotspot.gr)$hotspotGroup)

tabConScores <- tabCon/tabRef[names(tabRef) %in% names(tabCon)][names(tabCon)]

conGroups <- names(tabConScores[tabConScores >= .7])
midGroups <- names(tabConScores[tabConScores < .7 & tabConScores >= .3])
difGroups <- c(names(tabConScores[tabConScores < .3]), names(tabRef[!(names(tabRef) %in% names(tabCon))]))


conDomainQue <- s1_s2.gr[elementMetadata(s1_s2.gr)$hotspotGroup %in% conGroups]
midDomainQue <- s1_s2.gr[elementMetadata(s1_s2.gr)$hotspotGroup %in% midGroups]
difDomainQue <- s1_s2.gr[elementMetadata(s1_s2.gr)$hotspotGroup %in% difGroups]


conDomainRef <- s1.gr[elementMetadata(s1.gr)$hotspotGroup %in% conGroups]
midDomainRef <- s1.gr[elementMetadata(s1.gr)$hotspotGroup %in% midGroups]
difDomainRef <- s1.gr[elementMetadata(s1.gr)$hotspotGroup %in% difGroups]


s2RepsCon <- data.frame(domains = "s2RepsCon", 
                        repInsertion = s2DataList$repSummary[s2DataList$repSummary$binID %in% elementMetadata(subsetByOverlaps(s2bin.gr, conDomainQue))$binID, repGroups[i]])
s1RepsCon <- data.frame(domains = "s1RepsCon",
                        repInsertion = s1DataList$repSummary[s1DataList$repSummary$binID %in% elementMetadata(subsetByOverlaps(s1bin.gr, conDomainRef))$binID, repGroups[i]])

s2RepsMid <- data.frame(domains = "s2RepsMid",
                        repInsertion = s2DataList$repSummary[s2DataList$repSummary$binID %in% elementMetadata(subsetByOverlaps(s2bin.gr, midDomainQue))$binID, repGroups[i]])
s1RepsMid <- data.frame(domains = "s1RepsMid",
                        repInsertion = s1DataList$repSummary[s1DataList$repSummary$binID %in% elementMetadata(subsetByOverlaps(s1bin.gr, midDomainRef))$binID, repGroups[i]])

s2RepsDif <- data.frame(domains = "s2RepsDif",
                        repInsertion = s2DataList$repSummary[s2DataList$repSummary$binID %in% elementMetadata(subsetByOverlaps(s2bin.gr, difDomainQue))$binID, repGroups[i]])
s1RepsDif <- data.frame(domains = "s1RepsDif",
                        repInsertion = s1DataList$repSummary[s1DataList$repSummary$binID %in% elementMetadata(subsetByOverlaps(s1bin.gr, difDomainRef))$binID, repGroups[i]])

dfSummary <- rbind(s1RepsDif, s2RepsDif, s1RepsMid, s2RepsMid, s1RepsCon, s2RepsCon)

boxplot(dfSummary$repInsertion ~ dfSummary$domains, outline=FALSE)
stripchart(dfSummary$repInsertion ~ dfSummary$domains, vertical = TRUE, method = "jitter", pch = 16, cex = .3, jitter = .35, add = TRUE)


# plan now is to make the results cleaner
# when we pull regions make sure we are getting the right ones and not the ones beside them.


# at least here there is something now to compare


# may be good to get percentile rank scores

# the problem is the whole apples and oranges thing

# get the reference query stuff sorted. 





ol[elementMetadata(ref.gr)$hotspotGroup[ol[,1]] %in% midGroups, 2]


que.gr[ol[,2]]



# from these results it seems the distribtuion of newL1s is almost random


# the conserved one has to be the same family of element overlapping. 

# the sizes of the groups 

df <- data.frame(s1 = elementMetadata(s1.gr)$repGroup[ol[,1]], s2_s1 = elementMetadata(s2_s1.gr)$repGroup[ol[,2]])
table(df)

# 
table(elementMetadata(s1.gr)$repGroup)
table(elementMetadata(s2_s1.gr)$repGroup)

# 
table(elementMetadata(s2_s1.gr[-ol[,2]])$repGroup)
table(elementMetadata(s1.gr[-ol[,1]])$repGroup)

# are we counting the number of regions or the number of overlaps
# might be better to do the number of regions

#find conserved regions

canFam3_mm9canFameOl <- as.matrix(findOverlaps(canFam3gr, mm9_canFam3gr))








