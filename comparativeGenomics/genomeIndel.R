#!/usr/bin/env Rscript

## here we get our denomenator

rm(list = ls())


options(stringsAsFactors = FALSE)


library("optparse")


# input files 
# 


option_list = list(
  make_option(c("-b", "--baseRate"), type="character", default=NA, 
              help="base alignmnet blocks", metavar="character"),
  make_option(c("-q", "--queryIndel"), type="character", default=NA, 
              help="confident indels in query species", metavar="character"),
  make_option(c("-n", "--queryName"), type="character", default=NA, 
              help="name of query genome", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default="./", 
              help="out directory", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if(any(is.na(opt))){
  stop("argument missing")
}



opt$baseRate = "~/Desktop/RTN_domains/data/comparativeGenomics/queGenomes/baseRates/hg19.base"

opt$queryName = "hg19"

opt$queryIndel = "~/Desktop/RTN_domains/data/comparativeGenomics/supportedIndels/hg19.supportedIndel.que"


library(GenomicRanges)
library(RMySQL)

br <- read.table(opt$baseRate, header = TRUE)
br.gr <- GRanges(br)




queIndel <- read.table(opt$queryIndel, header = TRUE)
queIndel.gr <- GRanges(queIndel)



# about 300 human lineage specifc indel per Mb of alignable sequence



sum(mcols(queIndel.gr[mcols(queIndel.gr)$indel == "del"])$gapWidth)/ sum(width(br.gr)) * 1000000
sum(mcols(queIndel.gr[mcols(queIndel.gr)$indel == "ins"])$gapWidth)/ sum(width(br.gr)) * 1000000

## for every human Mb that shares ancestry with mouse and some other mammal, we find approxiamtly 2 kb human lineage specific deletion since divergence form mouse

## for every human Mb that shares ancestry with mouse and some other mammal, we find approximatly 60 kb of human lineage specific insertion since divergence from mouse.


length(queIndel.gr[mcols(queIndel.gr)$indel == "ins"])
length(queIndel.gr[mcols(queIndel.gr)$indel == "del"])

# 200 k human specific insertion events

# 100 k human specific deletion events


# is there any spatial correlation between the types of events. 

# let's see how these events are distributed 



mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = opt$queryName)
chrInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")

chrInfo.gr <- GRanges(seqnames = chrInfo$chrom, ranges = IRanges(start = 1,width = chrInfo$size))

slWindow <- unlist(slidingWindows(chrInfo.gr, width = 1000000, step = 1000000))

olbr <- findOverlaps(slWindow, br.gr)
intBr <- pintersect(slWindow[queryHits(olbr)],
                    br.gr[subjectHits(olbr)])
aggBr <- aggregate(x = width(intBr), by = list(queryHits(olbr)), FUN = sum)



olins <- findOverlaps(slWindow, queIndel.gr[mcols(queIndel.gr)$indel == "ins"])

aggIns <- aggregate(x = mcols(queIndel.gr[mcols(queIndel.gr)$indel == "ins"])$gapWidth[subjectHits(olins)], 
                    by = list(queryHits(olins)), FUN = sum)


plot(aggIns$x[1:500]/aggBr$x[1:500], type = "l")
plot(aggBr$x[1:500], type = "l")


merIns <- merge(aggIns, aggBr, by = 1)
smoothScatter(merIns$x.y, merIns$x.x)
smoothScatter(merIns$x.y, merIns$x.x/mer$x.y)
plot(merIns$Group.1,(merIns$x.x/merIns$x.y), type = "p", col = cols, xlim = c(500,1000))
spIns <- smooth.spline(x = merIns$Group.1, y = (merIns$x.x/merIns$x.y))
colsIns <- as.integer(as.factor(seqnames(slWindow)[merIns$Group.1]))
lines(spIns$x, spIns$y, type = "p", col = cols)
points(merIns$Group.1, (merIns$x.x/merIns$x.y), pch = 16, cex = .3)




oldel <- findOverlaps(slWindow, queIndel.gr[mcols(queIndel.gr)$indel == "del"])

aggDel <- aggregate(x = mcols(queIndel.gr[mcols(queIndel.gr)$indel == "del"])$gapWidth[subjectHits(oldel)], 
                    by = list(queryHits(oldel)), FUN = sum)


plot(aggDel$x[1:500]/aggBr$x[1:500], type = "l")
plot(aggBr$x[1:500], type = "l")


merDel <- merge(aggDel, aggBr, by = 1)
smoothScatter(merDel$x.y, merDel$x.x)
smoothScatter(merDel$x.y, merDel$x.x/merDel$x.y)
plot(merDel$Group.1,(merDel$x.x/merDel$x.y), type = "p", col = cols, xlim = c(1,500))
spDel <- smooth.spline(x = merDel$Group.1, y = (merDel$x.x/merDel$x.y))
colsDel <- as.integer(as.factor(seqnames(slWindow)[merDel$Group.1]))
lines(spDel$x, spDel$y, type = "p", col = cols)
points(merDel$Group.1, (merDel$x.x/merDel$x.y), pch = 16, cex = .3)


spDel <- smooth.spline(x = merDel$Group.1, y = scale(merDel$x.x/merDel$x.y),all.knots = TRUE, cv = TRUE)
spIns <- smooth.spline(x = merIns$Group.1, y = scale(merIns$x.x/merIns$x.y), all.knots = TRUE)


layout(c(1,2))
xlim = c(1500,2000)
plot(spDel$x, spDel$y, type = "p", col = colsDel, xlim = xlim, ylim = c(-3,3),pch = 16)
points(merDel$Group.1, scale(merDel$x.x/merDel$x.y), pch = 16, cex = .3)
plot(spIns$x, spIns$y, col = colsIns, pch = 16, xlim = xlim, ylim = c(-3,3))
points(merIns$Group.1, scale(merIns$x.x/merIns$x.y), pch = 16, cex = .3)


plot(spDel$x, spDel$y, type = "p", col = colsDel, xlim = xlim, ylim = c(-3,3),pch = 16)
lines(c(0,spDel$x,0),spDel$fit$coef)


lines(spIns$x, spIns$y, col = colsIns, pch = 16, xlim = xlim, ylim = c(-3,3))

library(zoo)
layout(1)
merAll <- merge(merIns, merDel, by = 1)
smoothScatter( rollmean(merAll$x.x.x/merAll$x.y.x, k = 10), rollmean(merAll$x.x.y/merAll$x.y.y, k = 10))


hist(merDel$x.y, breaks = 100)
hist(merDel$x.x, breaks = 100)
hist(merDel$x.x/merDel$x.y, breaks = 100)
smoothScatter(merDel$x.y,merDel$x.x/merDel$x.y)
modDel <- lm(merDel$x.x/merDel$x.y ~ merDel$x.y) 


hist(merIns$x.y, breaks = 100)
hist(merIns$x.x, breaks = 100)
hist(merIns$x.x/merIns$x.y, breaks = 100)
smoothScatter(merIns$x.y,merIns$x.x/merIns$x.y)
modIns <- lm((merIns$x.x/merIns$x.y)~ merIns$x.y) 

layout(1)
plot(modDel)
plot(modIns)

hist(merIns$x.y, breaks = 100)







l <- loess.smooth(x = merDel$Group.1,  y = scale(merDel$x.x/merDel$x.y), degree = 10,evaluation = 100, span = .01)
plot(l$x, l$y)

l <- loess(scale(merDel$x.x/merDel$x.y) ~ merDel$Group.1, span = .05)
plot(l$fitted)
predict(l)

spDel <- smooth.spline(x = merDel$Group.1, y = scale(merDel$x.x/merDel$x.y),all.knots = TRUE, cv = TRUE)
spIns <- smooth.spline(x = merIns$Group.1, y = scale(merIns$x.x/merIns$x.y), all.knots = TRUE)



layout(c(1,2))
xlim = c(1500,2000)
plot(spDel$x, spDel$y, type = "p", col = colsDel, xlim = xlim, ylim = c(-3,3),pch = 16)
points(merDel$Group.1, scale(merDel$x.x/merDel$x.y), pch = 16, cex = .3)
plot(spIns$x, spIns$y, col = colsIns, pch = 16, xlim = xlim, ylim = c(-3,3))
points(merIns$Group.1, scale(merIns$x.x/merIns$x.y), pch = 16, cex = .3)






fit <- smooth.spline(merDel$Group.1, scale(merDel$x.x/merDel$x.y), all.knots = TRUE)      # smooth.spline fit
res <- (fit$yin - fit$y)/(1-fit$lev)      # jackknife residuals
sigma <- sqrt(var(res))                     # estimate sd

upper <- fit$y + 2.0*sigma*sqrt(fit$lev)   # upper 95% conf. band
lower <- fit$y - 2.0*sigma*sqrt(fit$lev)   # lower 95% conf. band
matplot(fit$x, cbind(upper, fit$y, lower), type="plp", pch=".", xlim = c(1500,2000),col = c(1,2,1), ylim = c(-3,3))

fit <- smooth.spline(merIns$Group.1, scale(merIns$x.x/merIns$x.y), all.knots = TRUE)      # smooth.spline fit
res <- (fit$yin - fit$y)/(1-fit$lev)      # jackknife residuals
sigma <- sqrt(var(res))                     # estimate sd

upper <- fit$y + 2.0*sigma*sqrt(fit$lev)   # upper 95% conf. band
lower <- fit$y - 2.0*sigma*sqrt(fit$lev)   # lower 95% conf. band
matplot(fit$x, cbind(upper, fit$y, lower), type="plp", pch=".", xlim = c(1500,2000), add = TRUE, col = c(3,4,3))



# what is the correaltion of smoothed data 
# remove some of the randomness to see if there is a real difference




# so can we start comparing species 


# so now that I have some regions lets do some biology






