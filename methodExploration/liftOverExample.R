

#### good to know that this tool is available in the R environment.

### It can be helpful in building the pipeline

#if (!exists("cur")) load("cur.rda")

## ----lkOne,eval=FALSE----------------------------------------------------
## library(gwascat)

## ----lkcur---------------------------------------------------------------


## ----getch---------------------------------------------------------------

setwd("Desktop/RTN_domains/")

library(rtracklayer)
library(dplyr)

ch = import.chain("data/chainAlignments/liftOver/hg19/mm9/hg19ToMm9.over.chain")
ch
str(ch[[1]])

## ----dolift--------------------------------------------------------------
#seqlevelsStyle(cur) = "UCSC"  # necessary
test.df <- read.table("data/repeatHotspot/hg19/hg19Hotspots.bed", header = FALSE, 
                      col.names = c("chr", "start", "end", "hotspotID"))
test.gr <- GRanges(seqnames = Rle(test.df$chr),
                   ranges = IRanges(start = test.df$start, end = test.df$end),
                   hotspotID = test.df$hotspotID
)
cur19 = liftOver(test.gr, ch)
#unlist(cur19)
red <- reduce(x = resize((cur19), width = 1500, fix = "center"))
s <- cur19



tList <- cur19


# could just use a for loop to assign groups 

names(tList) <- unique(as.data.frame(tList)$hotspotID)

red <- reduce(x = resize((tList), width = 50000, fix = "center"))
tList <- reduce(tList)


tDF <- as.data.frame((tList))

for(i in 1:length(tList)){
  tDF[tDF$group == i, "rangeGroup"] <- subjectHits(findOverlaps(tList[[i]], red[[i]]))
}


# maybe do something on tDF to remove any ranges that are two far from the one in front of them
# as an iterative technique
# then look at how the alignmnet has worked

# 


# by converting to data frame we can use dplyr methods to get information


a <- tDF

# a$distFromLast = c(NaN, a$start[2:nrow(a)] - a$end[1:(nrow(a) - 1)])
# a$distFromLast[c(FALSE, a$group_name[1:(nrow(a)-1)] != a$group_name[2:(nrow(a))] )] = NaN
# a$distFromLast[c(FALSE, a$seqnames[1:(nrow(a)-1)] != a$seqnames[2:(nrow(a))] )] = NaN
# 
# a$distToNext = c(a$start[2:nrow(a)] - a$end[1:(nrow(a) - 1)],NaN)
# a$distToNext[ a$group_name[1:(nrow(a)-1)] != a$group_name[2:(nrow(a))]] = NaN
# a$distToNext[ a$seqnames[1:(nrow(a)-1)] != a$seqnames[2:(nrow(a))]] = NaN
# 
# a$rangeGroup = subjectHits(ol)




a[c(51408,51409,51410,51411, 51412),]

#filter(a, (distToNext > 30000 & distFromLast > 30000) | (is.nan(distToNext) & distFromLast > 30000) | (is.nan(distFromLast) & distToNext > 30000))




# so this gets rid of any lonly regions
# what if we want to get rid of more than one region 
# if there are two regions near each other but far from the others. 

# we could group things below a certain level then get the group with the biggest weight and split it again
# that way we'll remove stuff that is a certain distance away


# a trimming method


# build group IDs for stuff that is close together

tDF <- a


a <- summarise(group_by(tDF, group_name, seqnames, rangeGroup, group), width = sum(width), start = min(start), end = max(end)) %>% 
  arrange(group_name, desc(width)) 

a <- filter(group_by(a, group_name) ,width == max(width)) %>%
  mutate(queWidth = end - start + 1)



hist((a$queWidth), breaks = 50, xlim = c(0,100000))


nrow(a)
nrow(a[a$width/50000 > .1,])

a2 <- a[a$width/a$queWidth > .1 & a$width/50000 > .1,]






hist((a2$width), breaks = 100)
plot(a2$queWidth, pch = 16, cex = .5)
id <- identify(a2$queWidth)

a2[id,]

cov <- coverage(tList[["ancient_1946_501"]])
sl <- IRanges::slice(cov, lower = 1)
gr <- GRanges(sl)
gr <- gr[seqnames(gr) == "chr12"]


# the width has to be big enough 

# the query width has to be small enough

# are we going to worry about rouge small sections that could be messing up the bigger picture






lapply(lapply(lapply(tList, seqnames),table), names)

# list elements with 











cur19[1:2]


aligningWidth <- unlist(lapply(lapply(reduce(cur19), width), sum))

testStart = unlist(lapply(lapply(cur19, start), min))
testEnd = unlist(lapply(lapply(cur19, end), max))

absoluteRanges(cur19[[2]])

# inf is returned if there are missing vlues

# what do we do if part of it is on another chromosome








coveredWidth <- testEnd - testStart

# the ranges so far map over



# what does it look like when they don't



# shit, it is possible to make an absolute ranges object 
# this merges all the chromosomes into one long chromosome
# this would be helpful for plotting 
# may even be helpful for random sampling

gr <- GRanges(Rle(c("chr2", "chr1", "chr3", "chr1"), 4:1),
              IRanges(1:10, width=5),
              seqinfo=Seqinfo(c("chr1", "chr2", "chr3"), c(100, 50, 20)))

ar <- absoluteRanges(gr)
ar

gr2 <- relativeRanges(ar, seqlengths(gr))
gr2

## Sanity check:
stopifnot(all(gr == gr2))


# by switching between the two types of ranges we can pick a random position in an unbiased manner. 
# not really mammalian genomes are two large for our max integer size


# there seems to be no other good ways of randomizing stuff
# another way is to make randomized bins according to the size dsitrbutions of the bins we are sampling and then pick some.

# so we can guarentee no overlap
# however the bin sizes wont be exactly the same 








