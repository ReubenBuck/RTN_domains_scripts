rm(list = ls())




DNaseCovEdge <- function(regionsRefTE, dnase.gr){
  regionsRefTE5.gr <- GRanges(seqnames = regionsRefTE$chr,
                                ranges = IRanges(start = regionsRefTE$start, end = regionsRefTE$start + 500))
  
  regionsRefTE3.gr <- GRanges(seqnames = regionsRefTE$chr,
                                ranges = IRanges(start = regionsRefTE$end - 500, end = regionsRefTE$end))
  
  dnase5.gr <- GenomicRanges::intersect(dnase.gr, regionsRefTE5.gr)
  dnase3.gr <- GenomicRanges::intersect(dnase.gr, regionsRefTE3.gr)
  
  ol5 <- as.matrix(findOverlaps(regionsRefTE5.gr,dnase5.gr))
  ol3 <- as.matrix(findOverlaps(regionsRefTE3.gr,dnase3.gr))
  
  range5 <- IRanges(start = start(dnase5.gr[ol5[,2]]) - start(regionsRefTE5.gr[ol5[,1]]), 
                    end =  end(dnase5.gr[ol5[,2]]) - start(regionsRefTE5.gr[ol5[,1]]) )
  
  range3 <- IRanges(start = start(dnase3.gr[ol3[,2]]) - start(regionsRefTE3.gr[ol3[,1]]), 
                    end =  end(dnase3.gr[ol3[,2]]) - start(regionsRefTE3.gr[ol3[,1]]) )
  out <- list(prime5 = coverage(range5)[1:500]/nrow(regionsRefTE), prime3 = coverage(range3)[1:500]/nrow(regionsRefTE))
  
}


# I wonder what mouse looks like
# if we see the same level of overlap



DNaseCovAll <- function(regionsRefTE, dnase.gr){
  regionsRefTE.gr <- GRanges(seqnames = regionsRefTE$chr,
                               ranges = IRanges(start = regionsRefTE$start, end = regionsRefTE$end))
  
 # dnase.gr <- GenomicRanges::intersect(dnase.gr,regionsRefTE.gr)
  ol <- as.matrix(findOverlaps(regionsRefTE.gr,dnase.gr))
  
  range <- IRanges(start = start(dnase.gr[ol[,2]]) - start(regionsRefTE.gr[ol[,1]]), 
                   end =  end(dnase.gr[ol[,2]]) - start(regionsRefTE.gr[ol[,1]]) )
  
  out <- coverage(range)
  
}


setwd("~/Desktop/RTN_domains/")




gapInfoRef <- read.table(file = "data/chainAlignments/test_insertionHumanVsMouse.txt")

# gaps in mouse are equal to zero
queNoGap <- (1:nrow(gapInfoRef))[gapInfoRef[2:(nrow(gapInfoRef)),9] -  gapInfoRef[1:(nrow(gapInfoRef)-1),10] == 0]


plot(log10(gapInfoRef[queNoGap+1,3] -  gapInfoRef[queNoGap,4]))
abline(h = log10(300))
abline(h = log10(6000))


mHIST <- hist(log10(gapInfoRef[queNoGap+1,3] -  gapInfoRef[queNoGap,4]),plot=F, breaks = 100)
layout(matrix(c(1,2), nrow = 1))
plot(sort(log10(gapInfoRef[queNoGap+1,3] -  gapInfoRef[queNoGap,4])))
barplot(mHIST$counts, horiz = T)

queNoGap2 <- queNoGap[((gapInfoRef[queNoGap+1,3] -  gapInfoRef[queNoGap,4])) < 330 & ((gapInfoRef[queNoGap+1,3] -  gapInfoRef[queNoGap,4])) > 310 ]


gapInfoQue<- read.table(file = "data/chainAlignments/test_insertionMouse.txt")
refNoGap <- (1:nrow(gapInfoQue))[gapInfoQue[2:(nrow(gapInfoQue)),9] -  gapInfoQue[1:(nrow(gapInfoQue)-1),10] == 0]

plot(log10(gapInfoQue[refNoGap+1,3] -  gapInfoQue[refNoGap,4]))
abline(h = log10(300))
abline(h = log10(6000))


hHIST <- hist(log10(gapInfoQue[refNoGap+1,3] -  gapInfoQue[refNoGap,4]),plot=F, breaks = 200)
layout(matrix(c(1,2), nrow = 1))
plot(sort(log10(gapInfoQue[refNoGap+1,3] -  gapInfoQue[refNoGap,4])))
barplot(hHIST$counts, horiz = T)


refNoGap2 <- refNoGap[((gapInfoQue[refNoGap+1,3] -  gapInfoQue[refNoGap,4])) < 330 & ((gapInfoQue[refNoGap+1,3] -  gapInfoQue[refNoGap,4])) > 310 ]







# pull out flanking regions 
regionsRefTE <- data.frame(chr = gapInfoRef[queNoGap2,1], 
                             start = gapInfoRef[queNoGap2,4] - 0, 
                             end = gapInfoRef[queNoGap2 + 1,3] + 0, 
                             ID = paste("refTE","Ref", 1:length(queNoGap2), sep = "_"))
regionsRefTE.gr <- GRanges(seqnames = Rle(regionsRefTE$chr),
                             ranges = IRanges(start = regionsRefTE$start,
                                              end = regionsRefTE$end),
                             ID = regionsRefTE$ID)


regionsRefNoTE <- data.frame(chr = gapInfoQue[refNoGap2,7], 
                               start = gapInfoQue[refNoGap2,10] - 50, 
                               end = gapInfoQue[refNoGap2 + 1,9] + 50, 
                               ID = paste("refNoTE","Ref", 1:length(refNoGap2), sep = "_"))
regionsRefNoTE <- regionsRefNoTE[complete.cases(regionsRefNoTE),]
regionsRefNoTE.gr <- GRanges(seqnames = Rle(regionsRefNoTE$chr),
                              ranges = IRanges(start = regionsRefNoTE$start,
                                                end = regionsRefNoTE$end),
                              ID = regionsRefNoTE$ID)


refRegions.gr <- c(regionsRefTE.gr, regionsRefNoTE.gr)


regionsQueTE <- data.frame(chr = gapInfoQue[refNoGap2,1], 
                           start = gapInfoQue[refNoGap2,4] - 0, 
                           end = gapInfoQue[refNoGap2 + 1,3] + 0, 
                           ID = paste("refNoTE","Que",1:length(refNoGap2), sep = "_"))
regionsQueTE <- regionsQueTE[complete.cases(regionsQueTE),]
regionsQueTE.gr <- GRanges(seqnames = Rle(regionsQueTE$chr),
                           ranges = IRanges(start = regionsQueTE$start,
                                            end = regionsQueTE$end),
                           ID = regionsQueTE$ID)


regionsQueNoTE <- data.frame(chr = gapInfoRef[queNoGap2,7], 
                             start = gapInfoRef[queNoGap2,10] - 50, 
                             end = gapInfoRef[queNoGap2 + 1,9] + 50, 
                             ID = paste("refTE","Que",1:length(queNoGap2), sep = "_"))
regionsQueNoTE.gr <- GRanges(seqnames = Rle(regionsQueNoTE$chr),
                             ranges = IRanges(start = regionsQueNoTE$start,
                                              end = regionsQueNoTE$end),
                             ID = regionsQueNoTE$ID)

queRegions.gr <- c(regionsQueTE.gr, regionsQueNoTE.gr)

# from eyeballing the data isn't very clean
# some gaps don't correspond to new insertions
# this can be fixed though.



dnaseClust <- read.table(file = "~/Downloads/permissive_enhancers.bed",header = F, skip = 1)
dnaseClust.gr <- GRanges(seqnames = dnaseClust[,1], ranges = IRanges(start = dnaseClust[,2],
                                                                  end = dnaseClust[,3]),
                         number = dnaseClust[,5])

#dnaseClust.gr <- subsetByOverlaps(dnaseClust.gr, refRegions)

#dnaseClust.gr <- dnaseClust.gr[seqnames(dnaseClust.gr) == "chr1"]

# dnaseSites <- NULL
# for(i in 1:length(dnaseClust.gr)){
#   repNum <- elementMetadata(dnaseClust.gr)$number[i]
# #  roW <- data.frame(chr = rep(dnaseClust[i,2], repNum), 
# #                    start = rep(dnaseClust[i,3], repNum),
# #                    end = rep(dnaseClust[i,4], repNum))
#   stack <- as.data.frame(rep(dnaseClust.gr[i], repNum))
#   dnaseSites <- rbind(dnaseSites, stack)
# }
# 
# dnaseSites.gr <- GRanges(seqnames = Rle(dnaseSites$seqnames), 
#                          ranges = IRanges(end = dnaseSites$end, start = dnaseSites$start))
# 
# 




# noTE <- DNaseCovAll(regionsRefTE = regionsRefNoTE[complete.cases(regionsRefNoTE),], dnase.gr = dnaseClust.gr)
# 
# 
# TE <- DNaseCovAll(regionsRefTE = regionsRefTE[complete.cases(regionsRefTE),], dnase.gr = dnaseClust.gr )
# 
# layout(c(1,2))
# plot(noTE,type= "l")
# plot(TE,type = "l")


# Not an open site but it's not closed either

# not really sure how to interpret the TE associated dip though.

# it happen right where the breakpoint is

# the other idea is to turn around the alignment fil to get cleaner gaps
# didn't really work


# I t might be worth identifying conserved non-gaps between human and an out group
# get an enhancer set of conserved non gaps
# then test to see if there is some sort of enrichment of gaps in the other species. 
# the question becomes much more easier, are you more likly to see a gap in an enhancer than by pure chance.

length(dnaseClust.gr)
range(dnaseClust$V2)
range(gapInfoRef[,3])
range(gapInfoQue[,9])

gapInfoRef_ref.gr <- GRanges(seqnames = Rle(gapInfoRef[,1]),
                             ranges = IRanges( start = gapInfoRef[,3], end = gapInfoRef[,4])
                             )


layout(1)
plot(c(0,50), c(0,.2), type = "n")
for(i in 1:50){
clustNo <- i
points(i,sum(overlapsAny(regionsRefNoTE.gr, dnaseClust.gr[elementMetadata(dnaseClust.gr)$number > clustNo]))/length(regionsRefNoTE.gr), col = 1)

points(i,sum(overlapsAny(regionsRefTE.gr, dnaseClust.gr[elementMetadata(dnaseClust.gr)$number > clustNo]))/length(regionsRefTE.gr), col = 2)

points(i,sum(width(GenomicRanges::intersect(gapInfoRef_ref.gr, dnaseClust.gr[elementMetadata(dnaseClust.gr)$number > clustNo])))/sum(width(gapInfoRef_ref.gr)), col = 3)
}

# so lets set up a fair permutation test
# this will tell us whats going on.
# so we might also want to think about tidying up our dataset


# increase the strength of the cluster see if we get a change



# this isn't working





# so approximatly 18% of alginable sequence doubles as a DNase1 cluster

df <- data.frame(clusts = c(18,35568746), nonClusts = c(152, 162491720), row.names = c("site", "nonSite"))
fisher.test(df)

# from our initail analysis, it appears a TE site in chimp is unlikly to be an open DNA cluster in human

# so maybe this percentage would be higher if our selection was more refined.

# also we might only be looking at something that avoids multi mapping

sum(width(gapInfoRef_ref.gr))








