
library(GenomicRanges)


specRef = "hg19"
specQue = "mm10"

load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/",specRef,".",specQue,".netData.RData",sep = ""))



resGenomeStat <- matrix(data = NA, ncol = 9, nrow = 2)
rownames(resGenomeStat) <- genomes
colnames(resGenomeStat) <- c("genome_size", "fills","ancestral" ,
                             "fills_percent_ancestral", "all_gaps", "placed_gaps", 
                             "placed_gaps_percentage","gain", "loss")
resGenomeStat <- as.data.frame(resGenomeStat)

# genome size
## mm10
resGenomeStat$genome_size[1] <- round( (sum(as.numeric(seqlengths(refSeqGaps.gr))) - sum(width(refSeqGaps.gr)))/1e6 )
## hg19
resGenomeStat$genome_size[2] <-round((sum(as.numeric(seqlengths(queSeqGaps.gr))) - sum(width(queSeqGaps.gr)))/1e6)

# all non gaps
## mm10
resGenomeStat$fills[1] <- round(sum(width(refFill.gr))/1e6)
## hg19
resGenomeStat$fills[2] <- round(sum(width(queFill.gr))/1e6)

#ancestral
resGenomeStat$ancestral[1] <- round(sum(width(refAncDna.gr))/1e6)
resGenomeStat$ancestral[2] <- round(sum(width(queAncDna.gr))/1e6)


# ancestral %
# mm10
resGenomeStat$fills_percent_ancestral[1] <- round(sum(width(GenomicRanges::intersect(refFill.gr, refAncDna.gr))) / sum(width(refFill.gr)) * 100, digits = 1)
# hg19
resGenomeStat$fills_percent_ancestral[2] <- round(sum(width(GenomicRanges::intersect(queFill.gr, queAncDna.gr))) / sum(width(queFill.gr)) * 100, digits = 1)


# all gaps
#mm10
resGenomeStat$all_gaps[1] <- round(sum(width(refFillGaps.gr))/1e6)
# hg19
resGenomeStat$all_gaps[2] <- round(sum(width(queFillGaps.gr))/1e6)


# placeable gaps 
#mm10
resGenomeStat$placed_gaps[1] <- round(sum(width(refGap.gr))/1e6)
#hg19
resGenomeStat$placed_gaps[2] <- round(sum(width(queGap.gr))/1e6)



# placeable gaps percentage
#mm10
resGenomeStat$placed_gaps_percentage[1] <- round(sum(width(refGap.gr))/sum(width(refFillGaps.gr))*100,digits = 1)
#hg19
resGenomeStat$placed_gaps_percentage[2] <- round(sum(width(queGap.gr))/sum(width(queFillGaps.gr))*100, digits = 1)


# gain
## mm10
resGenomeStat$gain[1] <- round(sum(width(setdiff(refGap.gr, refAncDna.gr)   ))/1e6)
# hg19
resGenomeStat$gain[2] <- round(sum(width(setdiff(queGap.gr, queAncDna.gr)   ))/1e6)

# loss
# mm10 
resGenomeStat$loss[1] <- round(sum(width(intersect(queGap.gr, queAncDna.gr)   ))/1e6)
# hg19
resGenomeStat$loss[2] <- round(sum(width(intersect(refGap.gr, refAncDna.gr)   ))/1e6)


# we're ignoring segmental duplication
# we're ignoring seq gaps
# we're ignoring 

# ref seqGaps is not included in genome totals
sum(width(refSeqGaps.gr))/1e6

(sum(width(refNonRBHFill.gr))/1e6)

(sum(width(refNonRBHGap.gr))/1e6)

((sum(width(refFillGaps.gr))/1e6) - (sum(width(refGap.gr))/1e6) - (sum(width(refNonRBHGap.gr))/1e6))




resGenomeStat$genome_size[1] - (sum(width(refFill.gr))/1e6 + sum(width(refGap.gr))/1e6)
# 416 Mb of missing sequence, We need to document its source




filterTable = matrix(data = NA, nrow = 6, ncol = 2)
colnames(filterTable) <- c("hg19", "mm10")
rownames(filterTable) <- c("sequencedBases",
                           "nonNetGaps",
                           "nonRBHGaps",
                           "nonRBHFills",
                           "remainingFills",
                           "remainingGaps"
                           )
filterTable <- data.frame(t(filterTable))

# total genome

filterTable$sequencedBases[1] = resGenomeStat$genome_size[1] 
filterTable$sequencedBases[2] = resGenomeStat$genome_size[2] 


# minus nonFill gaps
refNonNetGap = setdiff(refFillGaps.gr, union(refGap.gr, refNonRBHGap.gr))
queNonNetGap = setdiff(queFillGaps.gr, union(queGap.gr, queNonRBHGap.gr))

filterTable$nonNetGaps[1] = sum(width(refNonNetGap))/1e6
filterTable$nonNetGaps[2] = sum(width(queNonNetGap))/1e6

# minus nonRBHGaps
filterTable$nonRBHGaps[1] = (sum(width(refNonRBHGap.gr))/1e6)
filterTable$nonRBHGaps[2] = (sum(width(queNonRBHGap.gr))/1e6)

# minus nonRBHFills 
filterTable$nonRBHFills[1] = (sum(width(refNonRBHFill.gr))/1e6)
filterTable$nonRBHFills[2] = (sum(width(queNonRBHFill.gr))/1e6)

# fills
filterTable$remainingFills[1] = sum(width(refFill.gr))/1e6
filterTable$remainingFills[2] = sum(width(queFill.gr))/1e6

# gaps
filterTable$remainingGaps[1] = sum(width(refGap.gr))/1e6
filterTable$remainingGaps[2] = sum(width(queGap.gr))/1e6


t(round(filterTable, digits = 1))

# we have 2.5 Gb of gapped sequnce we can assign to four catagories


refNonRBH <- data.frame(union(refNonRBHFill.gr, refNonRBHGap.gr))[,1:3]
queNonRBH <- data.frame(union(queNonRBHFill.gr, queNonRBHGap.gr))[,1:3]



write.table(refNonRBH,file = "~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/missingData/hg19nonRBH.bed", 
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

write.table(queNonRBH,file = "~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/missingData/mm10nonRBH.bed", 
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)


write.table(data.frame(refNonNetGap)[,1:3],
            file = "~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/missingData/hg19nonNet.bed", 
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

write.table(data.frame(queNonNetGap)[,1:3],
            file = "~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/missingData/mm10nonNet.bed", 
            quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)







union(union(queNonNetGap, queNonRBHFill.gr), queNonRBHGap.gr)







sum(width(GenomicRanges::intersect(refNonNetGap, refAncDna.gr)))/1e6
sum(width(GenomicRanges::intersect(queNonNetGap, queAncDna.gr)))/1e6



sum(width(intersect(refFill.gr, refAncDna.gr)))/1e6

sum(width(intersect(queFill.gr, queAncDna.gr)))/1e6




load(paste("~/Desktop/RTN_domains/R_objects/rmskTables/",genomes[1],"/",genomes[1],".RData",sep = ""))

refRep <- rep

refRepAll.gr <- GRanges(seqnames = refRep$genoChr, 
                        ranges = IRanges(start = refRep$genoStart + 1, end = refRep$genoEnd),
                        repName = refRep$repName,
                        repClass = refRep$repClass,
                        perDiv = refRep$perDiv)

refRep <- rep[complete.cases(rep),]

refRep.gr <- GRanges(seqnames = refRep$genoChr, 
                     ranges = IRanges(start = refRep$genoStart + 1, end = refRep$genoEnd),
                     repGroup = refRep$repGroup,
                     repFamily = refRep$repFamily)

# move L1MA to old group, because we are looking at a closer distance than our old paper
refRep.gr$repGroup[refRep.gr$repFamily == "L1MA"] <- "old_L1"


ol <- findOverlaps(refRep.gr, refFill.gr)
pInt <- GenomicRanges::pintersect(refRep.gr[queryHits(ol)], refFill.gr[subjectHits(ol)])
repSumFillRef <- data_frame( queryHit = queryHits(ol), width = width(pInt), repGroup = pInt$repGroup, repFamily = pInt$repFamily) %>%
  group_by(repGroup) %>%
  summarise(width = sum(width))


ol <- findOverlaps(refRep.gr, GenomicRanges::intersect(refAncDna.gr, refGap.gr))
pInt <- GenomicRanges::pintersect(refRep.gr[queryHits(ol)], GenomicRanges::intersect(refAncDna.gr, refGap.gr)[subjectHits(ol)])
repSumDelRef <- data_frame( queryHit = queryHits(ol), width = width(pInt), repGroup = pInt$repGroup, repFamily = pInt$repFamily) %>%
  group_by(repGroup) %>%
  summarise(width = sum(width))

ol <- findOverlaps(refRep.gr, GenomicRanges::setdiff(refGap.gr,refAncDna.gr))
pInt <- GenomicRanges::pintersect(refRep.gr[queryHits(ol)], GenomicRanges::setdiff(refGap.gr,refAncDna.gr)[subjectHits(ol)])
repSumInsRef <- data_frame( queryHit = queryHits(ol), width = width(pInt), repGroup = pInt$repGroup, repFamily = pInt$repFamily) %>%
  group_by(repGroup) %>%
  summarise(width = sum(width))







# query genome

load(paste("~/Desktop/RTN_domains/R_objects/rmskTables/",genomes[2],"/",genomes[2],".RData",sep = ""))

queRep <- rep

queRepAll.gr <- GRanges(seqnames = queRep$genoChr, 
                        ranges = IRanges(start = queRep$genoStart + 1, end = queRep$genoEnd),
                        repName = queRep$repName,
                        repClass = queRep$repClass,
                        perDiv = queRep$perDiv)




queRep <- rep[complete.cases(rep),]

queRep.gr <- GRanges(seqnames = queRep$genoChr, 
                  ranges = IRanges(start = queRep$genoStart + 1, end = queRep$genoEnd),
                  repGroup = queRep$repGroup,
                  repFamily = queRep$repFamily)

# move L1MA to old group, because we are looking at a closer distance than our old paper
queRep.gr$repGroup[queRep.gr$repFamily == "L1MA"] <- "old_L1"


ol <- findOverlaps(queRep.gr, queFill.gr)
pInt <- GenomicRanges::pintersect(queRep.gr[queryHits(ol)], queFill.gr[subjectHits(ol)])
repSumFillQue <- data_frame( queryHit = queryHits(ol), width = width(pInt), repGroup = pInt$repGroup, repFamily = pInt$repFamily) %>%
  group_by(repGroup) %>%
  summarise(width = sum(width))


ol <- findOverlaps(queRep.gr, GenomicRanges::intersect(queAncDna.gr, queGap.gr))
pInt <- GenomicRanges::pintersect(queRep.gr[queryHits(ol)], GenomicRanges::intersect(queAncDna.gr, queGap.gr)[subjectHits(ol)])
repSumDelQue <- data_frame( queryHit = queryHits(ol), width = width(pInt), repGroup = pInt$repGroup, repFamily = pInt$repFamily) %>%
  group_by(repGroup) %>%
  summarise(width = sum(width))

ol <- findOverlaps(queRep.gr, GenomicRanges::setdiff(queGap.gr,queAncDna.gr))
pInt <- GenomicRanges::pintersect(queRep.gr[queryHits(ol)], GenomicRanges::setdiff(queGap.gr,queAncDna.gr)[subjectHits(ol)])
repSumInsQue <- data_frame( queryHit = queryHits(ol), width = width(pInt), repGroup = pInt$repGroup, repFamily = pInt$repFamily) %>%
  group_by(repGroup) %>%
  summarise(width = sum(width))



repSumInsQue
repSumDelQue
repSumFillQue


repSumInsRef
repSumDelRef
repSumFillRef
# all the basic data we need 


repeatOLtable <- data.frame(matrix(NA, ncol = 4, nrow = 6))
colnames(repeatOLtable) <- c("regionType", "repeatType", "hg19", "mm10")


repeatOLtable$regionType <- c("fill", "fill", "Del", "Del", "Ins", "Ins")
repeatOLtable$repeatType <- rep(c("new_SINE", "ancient"),3)


repeatOLtable$hg19 <- round(c(repSumFillRef$width[3], repSumFillRef$width[1], 
                        repSumDelRef$width[3], repSumDelRef$width[1], 
                        repSumInsRef$width[3], repSumInsRef$width[1])/1e6,digits = 2)

repeatOLtable$mm10 <- round(c(repSumFillQue$width[3], repSumFillQue$width[1], 
                        repSumDelQue$width[3], repSumDelQue$width[1], 
                        repSumInsQue$width[3], repSumInsQue$width[1])/1e6,digits = 2)


# finding a significant amount of ancient repeats classified as insertions
# it is as though we can't identify where they are coming from, since a few of them have no ancestral bases



# so for our gapped areas, there is stuff that maybe should be a deletion in human


repeatOLtable

t(resGenomeStat)

t(round(filterTable, digits = 1))



sum(width(refRep.gr[refRep.gr$repGroup == "new_SINE"]))/1e6


table(queRepAll.gr$repClass)


queRepAll.gr <- queRepAll.gr[!(queRepAll.gr$repClass == "Simple_repeat" | 
                                 queRepAll.gr$repClass == "Satellite" |  
                                 queRepAll.gr$repClass == "Low_complexity" )]

queRepGap.gr <- queRepAll.gr[overlapsAny(queRepAll.gr, GenomicRanges::setdiff(queGap.gr, queAncDna.gr), type = "within")]
queRepFill.gr <- queRepAll.gr[overlapsAny(queRepAll.gr, queAncDna.gr, type = "within")]


dfQueGap <- data.frame(queRepGap.gr)
dfQueGap <- dfQueGap[ complete.cases(dfQueGap), ]
dfQueGap$type = "gap"


dfQueFill <- data.frame(queRepFill.gr)
dfQueFill <- dfQueFill[ complete.cases(dfQueFill), ]
dfQueFill$type = "fill"


dfQueAll <- rbind(dfQueGap, dfQueFill)


n = table(dfQueAll$repName)
dfQueAll$totalN = as.integer(n[dfQueAll$repName])

n = table(dfQueAll$type)
dfQueAll$typeN = as.integer(n[dfQueAll$type])

dfRepQue <- dfQueAll %>% 
  group_by(repName, type, totalN, typeN) %>% 
  summarise(divMean = mean(perDiv), divSD = sd(perDiv), number = n(), width = sum(width), 
            bot25 = quantile(perDiv, probs = .25),top25 = quantile(perDiv, probs = .75))

dfRepQue <- dfRepQue[order(dfRepQue$divMean),]


refRepAll.gr <- refRepAll.gr[!(refRepAll.gr$repClass == "Simple_repeat" | 
                                 refRepAll.gr$repClass == "Satellite" |  
                                 refRepAll.gr$repClass == "Low_complexity" )]

refRepGap.gr <- refRepAll.gr[overlapsAny(refRepAll.gr, GenomicRanges::setdiff(refGap.gr, refAncDna.gr), type = "within")]
refRepFill.gr <- refRepAll.gr[overlapsAny(refRepAll.gr, refAncDna.gr, type = "within")]


dfRefGap <- data.frame(refRepGap.gr)
dfRefGap <- dfRefGap[ complete.cases(dfRefGap), ]
dfRefGap$type = "gap"


dfRefFill <- data.frame(refRepFill.gr)
dfRefFill <- dfRefFill[ complete.cases(dfRefFill), ]
dfRefFill$type = "fill"


dfRefAll <- rbind(dfRefGap, dfRefFill)


n = table(dfRefAll$repName)
dfRefAll$totalN = as.integer(n[dfRefAll$repName])

n = table(dfRefAll$type)
dfRefAll$typeN = as.integer(n[dfRefAll$type])

dfRepRef <- dfRefAll %>% 
  group_by(repName, type, totalN, typeN) %>% 
  summarise(divMean = mean(perDiv), divSD = sd(perDiv), number = n(), width = sum(width), 
            bot25 = quantile(perDiv, probs = .25),top25 = quantile(perDiv, probs = .75))

dfRepRef <- dfRepRef[order(dfRepRef$divMean),]


layout(1:2)
par(mar=c(5,1,1,1), oma = c(2,2,2,2))
for(i in c("Ref","Que")){
  
  repSum <- get(paste("dfRep",i, sep = ""))
  
  plot(repSum$divMean[repSum$type == "gap"],
       (repSum$number/repSum$totalN)[repSum$type == "gap"],
       ylab = " gap elements per family",
       xlab = "percentage divergence from consensus", type = "n", xlim = c(0,35))
  
  #points(repSum$divMean[repSum$type == "gap"],
  #       (repSum$number/repSum$totalN)[repSum$type == "gap"], pch = 16, cex = .1)
  
  rect(xleft = repSum$bot25[repSum$type == "gap"],
       xright = repSum$top25[repSum$type == "gap"],
       ytop = (repSum$number/repSum$totalN)[repSum$type == "gap"] + repSum$width[repSum$type == "gap"]/5e8,
       ybottom = (repSum$number/repSum$totalN)[repSum$type == "gap"] - repSum$width[repSum$type == "gap"]/5e8,
       col = scales::alpha("black", .2), border = NA)
  
  #points(repSum$divMean[repSum$type == "gap"],
  #       (repSum$number/repSum$totalN)[repSum$type == "gap"], pch = 16, cex = .1, col = 2)

  rect(xleft = repSum$bot25[repSum$type == "fill"],
       xright = repSum$top25[repSum$type == "fill"],
       ytop = (1 - (repSum$number/repSum$totalN)[repSum$type == "fill"]) + repSum$width[repSum$type == "fill"]/5e8,
       ybottom = (1 - (repSum$number/repSum$totalN)[repSum$type == "fill"]) - repSum$width[repSum$type == "fill"]/5e8,
       col = scales::alpha("red", .2), border = NA)
  
  legend("bottomleft", pch = 15, legend = c("gap", "fill"),
         col = c(scales::alpha("black", .2),
                 scales::alpha("red", .2)),
         bty = "n", pt.cex = 2)
  
}

# a mixture of insertion and deletion,
# families that overlap species divergence

# families in regions deleted in mouse

data.frame((repSum[repSum$type == "fill",][order(repSum$width[repSum$type == "fill"], decreasing = TRUE),]))

# shared families between the species
# difficult 

# now we likely have our new repeat families
# devise a plan


# asign as either insertion or deletion

# we purge areas of the most recent repeats









head(dfRepQue)

int <- setdiff(dfRepQue$repName, dfRepRef$repName)

plot(sort(dfRepQue[dfRepQue$repName %in% int & dfRepQue$type == "fill",]$number / dfRepQue[dfRepQue$repName %in% int & dfRepQue$type == "fill",]$totalN))
plot(sort(dfRepRef[dfRepRef$repName %in% int & dfRepRef$type == "fill",]$number / dfRepRef[dfRepRef$repName %in% int & dfRepRef$type == "fill",]$totalN))

# even our common families are not exclusivly distributed across fills and gaps

# we just need all new DNA, the idea is that everything else is old


data.frame(dfRepRef[1 - (dfRepRef$number/dfRepRef$totalN) > .9,][order( dfRepRef$width[1 - (dfRepRef$number/dfRepRef$totalN) > .9], decreasing = TRUE ),])


# fills can go in gaps but gaps can't go in fills,
# this is how families should be split into new/old DNA
# unless that family spans the age gap


# again what is the chances that this DNA in a gap in new in ref or del in query
# are there particular sequnces that are favoured


layout(1:2)
plot(dfRepQue$divMean[dfRepQue$type == "gap"],
     (dfRepQue$number/dfRepQue$totalN)[dfRepQue$type == "gap"],
     ylab = " gap elements per family",
     xlab = "percentage divergence from consensus", type = "p", xlim = c(0,35))


ids <- identify(dfRepQue$divMean[dfRepQue$type == "gap"],
                (dfRepQue$number/dfRepQue$totalN)[dfRepQue$type == "gap"])
idNames <- repSum[repSum$type == "gap",][ids,]$repName
dfRepRef[dfRepRef$repName %in% idNames,]


col = rep(1, sum(dfRepRef$type == "gap"))
col[dfRepRef$repName[dfRepRef$type == "gap"] %in% idNames] = 2
plot(dfRepRef$divMean[dfRepRef$type == "gap"],
     (dfRepRef$number/dfRepRef$totalN)[dfRepRef$type == "gap"],
     ylab = " gap elements per family",
     xlab = "percentage divergence from consensus", type = "p", xlim = c(0,35),
     col = col, cex = col/2, pch = 16)



repSum[repSum$type == "gap",][ids,]


# the venn diagram approach
# how many families in one catagory overlap with families from another catagory 

# catagories, human new, mouse new, human old, mouse old

# the only group that should have any overlap is human old and mouse old

# human new and mouse new may have some overlap

# should that be a weighted intersect ? 
# number of bp from overlapping families 

# we should also melt and recast the dataframes so we can analyse the data properly

library(reshape2)

head(repSum)

motifSum <- melt(repSum, c("repName", "divMean"), "type")[1:100,]
motifSum <- cast(motifSum, divMean ~ repName)
motifSum[is.na(motifSum)] <- 0









