
rm(list=ls())

setwd("~/Desktop/RTN_domains/")
library(GenomicRanges)


# seperate out IDs

# so we want to look at exclusive groups
# Those with minimal overlap
# also domains that are large. 

s1name <- "canFam3"
s2name <- "mm9"

repGroups <- c("ancient", "new_SINE", "new_L1", "old_L1")

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

olSummary <- function(ref.gr, que.gr, repGroups){
  ol <- as.matrix(findOverlaps(ref.gr, que.gr))
  dfRepGroup <- data.frame(elementMetadata(ref.gr)$repGroup[ol[,1]], 
                           elementMetadata(que.gr)$repGroup[ol[,2]])
  dfDomainID <- data.frame(elementMetadata(ref.gr)$domainID[ol[,1]], 
                           elementMetadata(que.gr)$domainID[ol[,2]])
  
  # go through each combination and count how many unique regions
  olSumS1 <- matrix(nrow = 4, ncol = 4, dimnames = list(repGroups,repGroups))
  olTotals <- matrix(nrow = 4, ncol = 1, dimnames = list(repGroups,"totalOL"))

  for(i in 1:length(repGroups)){
    for(j in 1:length(repGroups)){
      olSumS1[i,j] <- length(unique(dfDomainID[dfRepGroup[,1] == repGroups[i] & dfRepGroup[,2] == repGroups[j],1])) / 
        length(ref.gr[elementMetadata(ref.gr)$repGroup == repGroups[i]])
    }
    olTotals[i] <- length(unique(dfDomainID[dfRepGroup[,1] == repGroups[i], 1] )) / length(ref.gr[elementMetadata(ref.gr)$repGroup == repGroups[i]])
  }
  return(cbind(olSumS1, olTotals))
}


olSummary(s2.gr, s1_s2.gr, repGroups)





ol <- as.matrix(findOverlaps(s1.gr, s2_s1.gr))
dfRepGroup <- data.frame(elementMetadata(s1.gr)$repGroup[ol[,1]], elementMetadata(s2_s1.gr)$repGroup[ol[,2]])
dfDomainID <- data.frame(elementMetadata(s1.gr)$domainID[ol[,1]], elementMetadata(s2_s1.gr)$domainID[ol[,2]])

# go through each combination and count how many unique regions
olSumS1 <- matrix(nrow = 4, ncol = 4, dimnames = list(repGroups,repGroups))
#olSumS2_S1 <- matrix(nrow = 4, ncol = 4, dimnames = list(repGroups,repGroups))

for(i in 1:length(repGroups)){
  for(j in 1:length(repGroups)){
    olSumS1[i,j] <- length(unique(dfDomainID[dfRepGroup[,1] == repGroups[i] & dfRepGroup[,2] == repGroups[j],1])) / length(s1.gr[elementMetadata(s1.gr)$repGroup == repGroups[i]])
#    olSumS2_S1[i,j] <- length(unique(dfDomainID[dfRepGroup[,1] == repGroups[i] & dfRepGroup[,2] == repGroups[j],2])) / length(s2_s1.gr[elementMetadata(s2_s1.gr)$repGroup == repGroups[i]])
  }
}


barplot(olSumS1, beside = TRUE, ylim = c(0,1))




olSumS2_S1
barplot(olSumS2_S1, beside = TRUE, ylim = c(0,1))


# the proportion of things that are overlapped
olTotals <- data.frame(s1 = rep(NA,4), s2_s1 = rep(NA,4), row.names = repGroups)
for(i in 1:length(repGroups)){
  olTotals$s1[i] <- length(unique(dfDomainID[dfRepGroup[,1] == repGroups[i], 1] )) / length(s1.gr[elementMetadata(s1.gr)$repGroup == repGroups[i]])
  
  olTotals$s2_s1[i] <- length(unique(dfDomainID[dfRepGroup[,2] == repGroups[i], 2] )) / length(s2_s1.gr[elementMetadata(s2_s1.gr)$repGroup == repGroups[i]])
  
}


barplot(t(olTotals), beside = TRUE, ylim = c(0,1) )

### lets look at self overlap

ol <- as.matrix(findOverlaps(s1.gr, s1.gr))
dfRepGroup <- data.frame(elementMetadata(s1.gr)$repGroup[ol[,1]], elementMetadata(s1.gr)$repGroup[ol[,2]])
dfDomainID <- data.frame(elementMetadata(s1.gr)$domainID[ol[,1]], elementMetadata(s1.gr)$domainID[ol[,2]])


for(i in 1:length(repGroups)){
  for(j in 1:length(repGroups)){
    olSumS1[i,j] <- length(unique(dfDomainID[dfRepGroup[,1] == repGroups[i] & dfRepGroup[,2] == repGroups[j],1])) / length(s1.gr[elementMetadata(s1.gr)$repGroup == repGroups[i]])
 #   olSumS2_S1[i,j] <- length(unique(dfDomainID[dfRepGroup[,1] == repGroups[i] & dfRepGroup[,2] == repGroups[j],2])) / length(s2_s1.gr[elementMetadata(s2_s1.gr)$repGroup == repGroups[i]])
  }
}
olSumS1
barplot(olSumS1, beside = TRUE)


ol <- as.matrix(findOverlaps(s2_s1.gr, s2_s1.gr))
dfRepGroup <- data.frame(elementMetadata(s2_s1.gr)$repGroup[ol[,1]], elementMetadata(s2_s1.gr)$repGroup[ol[,2]])
dfDomainID <- data.frame(elementMetadata(s2_s1.gr)$domainID[ol[,1]], elementMetadata(s2_s1.gr)$domainID[ol[,2]])


for(i in 1:length(repGroups)){
  for(j in 1:length(repGroups)){
      olSumS2_S1[i,j] <- length(unique(dfDomainID[dfRepGroup[,1] == repGroups[i] & dfRepGroup[,2] == repGroups[j],2])) / length(s2_s1.gr[elementMetadata(s2_s1.gr)$repGroup == repGroups[i]])
  }
}
olSumS2_S1










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








