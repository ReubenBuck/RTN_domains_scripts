#### alignments and repeats
library(reshape)
library(dplyr)
library(GenomicRanges)





ldaLine <- function(fit){
  gmean <- fit$prior%*%fit$means
  const <- drop(gmean%*%fit$scaling)
  
  slope <- -fit$scaling[1]/fit$scaling[2]
  intercept <- const/fit$scaling[2]
  
  res <- c(intercept, slope)
  names(res) <- c("intercept", "slope")
  return(res)
}


specRef = "hg19"
specQue = "mm10"

load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/",specRef,".",specQue,".netData.RData",sep = ""))


load(paste("~/Desktop/RTN_domains/R_objects/rmskTables/",genomes[1],"/",genomes[1],".RData",sep = ""))

refRep <- rep

refRepAll.gr <- GRanges(seqnames = refRep$genoChr, 
                        ranges = IRanges(start = refRep$genoStart + 1, end = refRep$genoEnd),
                        repName = refRep$repName,
                        repClass = refRep$repClass,
                        perDiv = refRep$perDiv)


# query genome

load(paste("~/Desktop/RTN_domains/R_objects/rmskTables/",genomes[2],"/",genomes[2],".RData",sep = ""))

queRep <- rep

queRepAll.gr <- GRanges(seqnames = queRep$genoChr, 
                        ranges = IRanges(start = queRep$genoStart + 1, end = queRep$genoEnd),
                        repName = queRep$repName,
                        repClass = queRep$repClass,
                        perDiv = queRep$perDiv)





queRepAll.gr <- queRepAll.gr[!(queRepAll.gr$repClass == "Simple_repeat" | 
                                 queRepAll.gr$repClass == "Satellite" |  
                                 queRepAll.gr$repClass == "Low_complexity" )]

queRepGap.gr <- queRepAll.gr[overlapsAny(queRepAll.gr, GenomicRanges::setdiff(queGap.gr, queAncDna.gr), type = "within")]
queRepFill.gr <- queRepAll.gr[overlapsAny(queRepAll.gr, GenomicRanges::union(queAncDna.gr,queFill.gr), type = "within")]


dfQueGap <- data.frame(queRepGap.gr)
dfQueGap <- dfQueGap[ complete.cases(dfQueGap), ]
dfQueGap$type = "gap"
dfQueFill <- data.frame(queRepFill.gr)
dfQueFill <- dfQueFill[ complete.cases(dfQueFill), ]
dfQueFill$type = "fill"


dfQueAll <- rbind(dfQueGap, dfQueFill)

dfRepQue <- dfQueAll %>% 
  group_by(repName, type) %>% 
  summarise(divMean = mean(perDiv), number = n(), width = sum(width), 
            bot25 = quantile(perDiv, probs = .25),top25 = quantile(perDiv, probs = .75),
            mid50 = quantile(perDiv, probs = .5))

dfRepQue <- melt(setDT(dfRepQue), id=1:2)
dfRepQue <- cast(dfRepQue, repName  ~ type + variable)
dfRepQue[is.na(dfRepQue)] <- 0
dfRepQue$totalNo <- dfRepQue$fill_number + dfRepQue$gap_number
dfRepQue$totalWidth <- dfRepQue$fill_width + dfRepQue$gap_width



refRepAll.gr <- refRepAll.gr[!(refRepAll.gr$repClass == "Simple_repeat" | 
                                 refRepAll.gr$repClass == "Satellite" |  
                                 refRepAll.gr$repClass == "Low_complexity" )]

refRepGap.gr <- refRepAll.gr[overlapsAny(refRepAll.gr, GenomicRanges::setdiff(refGap.gr, refAncDna.gr), type = "within")]
refRepFill.gr <- refRepAll.gr[overlapsAny(refRepAll.gr, GenomicRanges::union(refAncDna.gr,refFill.gr), type = "within")]


dfRefGap <- data.frame(refRepGap.gr)
dfRefGap <- dfRefGap[ complete.cases(dfRefGap), ]
dfRefGap$type = "gap"
dfRefFill <- data.frame(refRepFill.gr)
dfRefFill <- dfRefFill[ complete.cases(dfRefFill), ]
dfRefFill$type = "fill"

dfRefAll <- rbind(dfRefGap, dfRefFill)


dfRepRef <- dfRefAll %>% 
  group_by(repName, type) %>% 
  summarise(divMean = mean(perDiv), number = n(), width = sum(width), 
            bot25 = quantile(perDiv, probs = .25),top25 = quantile(perDiv, probs = .75), 
            mid50 = quantile(perDiv, probs = .5))

dfRepRef <- melt(setDT(dfRepRef), id=1:2)
dfRepRef <- cast(dfRepRef, repName  ~ type + variable)
dfRepRef[is.na(dfRepRef)] <- 0
dfRepRef$totalNo <- dfRepRef$fill_number + dfRepRef$gap_number
dfRepRef$totalWidth <- dfRepRef$fill_width + dfRepRef$gap_width




dfQueAllLda <- dfQueAll
familyGapOL <- dfRepQue$gap_width / dfRepQue$totalWidth
names(familyGapOL) <- dfRepQue$repName
dfQueAllLda$gapOL <-  familyGapOL[dfQueAllLda$repName]

dfRefAllLda <- dfRefAll
familyGapOL <- dfRepRef$gap_width / dfRepRef$totalWidth
names(familyGapOL) <- dfRepRef$repName
dfRefAllLda$gapOL <-  familyGapOL[dfRefAllLda$repName]


fit <- lda(data = rbind(dfQueAllLda, dfRefAllLda), type ~ perDiv + gapOL)

fitRef <- lda(data = dfRefAllLda, type ~ perDiv + gapOL)
fitQue <- lda(data = dfQueAllLda, type ~ perDiv + gapOL)










# we want to see if there is a difference between the divergence levels for each repeat
layout(matrix(1:2, nrow=1))
par(mar=c(5,1,1,1), oma = c(2,5,2,2))
for(i in c("Ref","Que")){
repSum <- get(paste("dfRep",i, sep = ""))

plot(repSum$gap_mid50, repSum$gap_width/repSum$totalWidth, pch = 16, cex = .3,
     ylab = "gap percentage overlap",
     xlab = "percentage divergence from consensus", type = "n", xlim = c(0,40),
     main = genomes[tolower(i)])


rect(xleft = repSum$gap_bot25,
     xright = repSum$gap_top25,
     ytop = (repSum$gap_width/repSum$totalWidth) + (repSum$gap_width/5e8),
     ybottom = (repSum$gap_width/repSum$totalWidth) - (repSum$gap_width/5e8),
     col = scales::alpha("black", .2), border = NA)


rect(xleft = repSum$fill_bot25,
     xright = repSum$fill_top25,
     ytop = ( (repSum$gap_width/repSum$totalWidth)) + (repSum$fill_width/5e8),
     ybottom = ( (repSum$gap_width/repSum$totalWidth)) - (repSum$fill_width/5e8),
     col = scales::alpha("red", .2), border = NA)

res <- ldaLine(fit)
abline( res[c("intercept", "slope")], lty = 2, col = 3)
res <- ldaLine(get(paste("fit",i, sep = "")))
abline( res[c("intercept", "slope")], lty = 2)

}


# compare our speceis specifc classifier to general classifier by looking at overlapping repeat names
# what ever gets the best result we can go with

# maybe repeat name should have been a variable 
# especially in the shared dataset 

# get the repeat family mean



# based on family name and percent divergence, we can place repeats into new/old groups

# count the number of new sequences that overlap

# joint and indivdual
dfRepRefBoth <- dfRefAll %>% 
  group_by(repName) %>% 
  summarise(perDiv = mean(perDiv), width = sum(width))

familyGapOL <- dfRepRef$gap_width / dfRepRef$totalWidth
names(familyGapOL) <- dfRepRef$repName
dfRepRefBoth$gapOL <-  familyGapOL[dfRepRefBoth$repName]

dfRepQueBoth <- dfQueAll %>% 
  group_by(repName) %>% 
  summarise(perDiv = mean(perDiv), width = sum(width))

familyGapOL <- dfRepQue$gap_width / dfRepQue$totalWidth
names(familyGapOL) <- dfRepQue$repName
dfRepQueBoth$gapOL <-  familyGapOL[dfRepQueBoth$repName]


# now we can classify families based on lda

# the total sequnce assigned to each group is pretty constant in the query
# maybe not so much in the ref

refClass <- data.frame(class = predict(fitRef, dfRepRefBoth)$class, width = dfRepRefBoth$width)
rownames(refClass) <- dfRepRefBoth$repName

queClass <- data.frame(class = predict(fitQue, dfRepQueBoth)$class, width = dfRepQueBoth$width)
rownames(queClass) <- dfRepQueBoth$repName



queSharedName <- queClass[intersect(rownames(queClass), rownames(refClass)), ]
refSharedName <- refClass[intersect(rownames(queClass), rownames(refClass)), ]

sum(queSharedName$class == refSharedName$class)/nrow(queSharedName)

sum(queSharedName$width[queSharedName$class == refSharedName$class])
sum(queSharedName$width[queSharedName$class != refSharedName$class])



queClass[intersect(rownames(queClass), rownames(refClass)), ] %>%
  group_by(class) %>%
  summarise( width = sum(width))


refClass[intersect(rownames(queClass), rownames(refClass)), ] %>%
  group_by(class) %>%
  summarise( width = sum(width))


queClass[setdiff(rownames(queClass), rownames(refClass)), ] %>%
  group_by(class) %>%
  summarise( width = sum(width))


refClass[setdiff(rownames(refClass), rownames(queClass)), ] %>%
  group_by(class) %>%
  summarise( width = sum(width))


# a significan tnumber of familes active during divergence
# families that are not shared but are 

# shared famileis that overlap gaps are new


refShare <- refClass[intersect(rownames(queClass), rownames(refClass)), ]
refNonShare <- refClass[setdiff(rownames(refClass), rownames(queClass)), ]

queShare <- queClass[intersect(rownames(queClass), rownames(refClass)), ]
queNonShare <- queClass[setdiff(rownames(queClass), rownames(refClass)), ]


plot(density(dfRefAll$perDiv[
  dfRefAll$repName %in% rownames(refShare[refShare$class == "gap",])]))
lines(density(dfRefAll$perDiv[
  dfRefAll$repName %in% rownames(refShare[refShare$class == "fill",])]), col= 2)
lines(density(dfRefAll$perDiv[
  dfRefAll$repName %in% rownames(refNonShare[refNonShare$class == "gap",])]), col= 3)
#lines(density(dfRefAll$perDiv[
#  dfRefAll$repName %in% rownames(refNonShare[refNonShare$class == "fill",])]), col= 4)


plot(density(dfQueAll$perDiv[
  dfQueAll$repName %in% rownames(queShare[queShare$class == "gap",])]))
lines(density(dfQueAll$perDiv[
  dfQueAll$repName %in% rownames(queShare[queShare$class == "fill",])]), col= 2)
lines(density(dfQueAll$perDiv[
  dfQueAll$repName %in% rownames(queNonShare[queNonShare$class == "gap",])]), col= 3)
#lines(density(dfQueAll$perDiv[
#  dfQueAll$repName %in% rownames(queNonShare[queNonShare$class == "fill",])]), col= 4)




# there are four classes of family
# new
# shared name, new distribtuion
# old
# new name, old distribution

# we need to define wheater a particular repeat belongs to a specifc family or not
hist((dfRefAll$perDiv[dfRefAll$repName %in% rownames(refShare[refShare$class == "fill",])]), col= 2, breaks = 50)




# now we have divergence from consensus
# we have a few inbetween families 

# now to split families into new and old and make comparison



# how to draw a straight line through the data



refNewShare.gr <- refRepAll.gr[refRepAll.gr$repName %in% rownames(refShare[refShare$class == "gap",])]
refNew.gr <- refRepAll.gr[refRepAll.gr$repName %in% rownames(refNonShare[refNonShare$class == "gap",])]



sum(width(intersect(refNewShare.gr, 
                    GenomicRanges::setdiff(refGap.gr, refAncDna.gr))))

sum(width(intersect(refNewShare.gr, 
                    GenomicRanges::intersect(refGap.gr, refAncDna.gr))))


sum(width(intersect(refNew.gr, 
                    GenomicRanges::setdiff(refGap.gr, refAncDna.gr))))

sum(width(intersect(refNew.gr, 
                    GenomicRanges::intersect(refGap.gr, refAncDna.gr))))


# now it is something like a 5% chance our data is misasigned 

newRep.gr <- refRepAll.gr[refRepAll.gr$repName %in% rownames(refClass[refClass$class == "gap",])]
repIns <- GenomicRanges::intersect(refGap.gr, newRep.gr)
repDel <- GenomicRanges::setdiff(refGap.gr, newRep.gr)

ancIns <- GenomicRanges::setdiff(refGap.gr, refAncDna.gr)
ancDel <- GenomicRanges::intersect(refGap.gr, refAncDna.gr)


gapSum <- sum(width(refGap.gr))

conTable <- matrix(data = NA, nrow = 2, ncol = 2, 
                   dimnames = list(c("repIns", "repDel"),
                                   c("ancIns", "ancDel")))

conTable["repIns", "ancIns"] <- sum(width(GenomicRanges::intersect(repIns, ancIns)))
conTable["repDel", "ancIns"] <- sum(width(GenomicRanges::intersect(repDel, ancIns)))
conTable["repIns", "ancDel"] <- sum(width(GenomicRanges::intersect(repIns, ancDel)))
conTable["repDel", "ancDel"] <- sum(width(GenomicRanges::intersect(repDel, ancDel)))

# looks like methods agree 80% of the time
# many of the "insertions" we found based on ancestral method, do not overlap repeats
# it is probably likly that they are small, ie 10bp
# in addition many of the deleted regions did not overlap repeats. 

# Of the shared families, are the same ones getting the same classification




library(MASS)
data(iris)
# fit model
fit <- lda(Species~., data=iris)

# summarize the fit
summary(fit)
# make predictions
predictions <- predict(fit, iris[,1:4])$class
# summarize accuracy
table(predictions, iris$Species)



# to use lda we have two groups
# but they are weighted 
dfQueAllLda <- dfQueAll
familyGapOL <- dfRepQue$gap_width / dfRepQue$totalWidth
names(familyGapOL) <- dfRepQue$repName
dfQueAllLda$gapOL <-  familyGapOL[dfQueAllLda$repName]

dfRefAllLda <- dfRefAll
familyGapOL <- dfRepRef$gap_width / dfRepRef$totalWidth
names(familyGapOL) <- dfRepRef$repName
dfRefAllLda$gapOL <-  familyGapOL[dfRefAllLda$repName]


fit <- lda(data = rbind(dfQueAllLda, dfRefAllLda), type ~ perDiv + gapOL)

# we do it based on families, 
# fit them according to their percent overlap
# or add the percent family gap overlap to the data




ldaLine <- function(fit){
  gmean <- fit$prior%*%fit$means
  const <- drop(gmean%*%fit$scaling)
  
  slope <- -fit$scaling[1]/fit$scaling[2]
  intercept <- const/fit$scaling[2]
  
  res <- c(intercept, slope)
  names(res) <- c("intercept", "slope")
  return(res)
}


res <- ldaLine(fit)






