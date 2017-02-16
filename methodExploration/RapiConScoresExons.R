### now we have genes with conserved intron/exon number 
### all the orthologous gene information is there too

### it is possible now to look at these regions in reference to their conservation scores
### are there any changes in conservation levels between the regions and are these changes constant. 

### maybe we could provide a count of the data instead. 

### this may in some ways be more accurate than a rate. 
rm(list = ls())

setwd("~/Desktop/RTN_domains/")

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(RMySQL)
library(devtools)
devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")


repGroups = c("ancient", "new_SINE", "new_L1", "old_L1")
repCols = c("darkblue", "aquamarine3", "purple", "red")
snames <- c(s1name = "hg19", s2name = "mm9", s3name = "canFam3")


ints <- read.table(file = "data/repeatHotspot/intersect.txt", header= TRUE)


load(file = "R_objects/ensGene/genomeInfoAll.RData")



mySession <- browserSession()
genome(mySession) <- "hg19"
track.names <- trackNames(ucscTableQuery(mySession))
#head(track.names)

track <- track.names[grep("Cons 46", names(track.names))]

table.names = tableNames(ucscTableQuery(mySession))

table.names[grep("46way",table.names )]

table = "phastCons46wayPlacental"

range <- genomeInfoAll$ensGene$hg19$exons


cons <- ucscTableQuery(mySession, track=track,
                       range=range[1:250], table=table)

a <- track(cons)

head(a)


hist(score(a))

ol <- as.matrix(findOverlaps(a, genomeInfoAll$ensGene$hg19$ensGene))

plot(score(a[ol[ol[,2] == 1,1]]), type = "l", ylim = c(0,1), col = 8)


for(i in 1:100){
  par(new = T)
plot(score(a[ol[ol[,2] == i,1]]), type = "l", col = 8, ylim = c(0,1))
}

mat <- matrix(NA, nrow = 29, ncol = 100)
ourMeans <- NULL
ourMedians <- NULL
for(i in 1:29){
x <- aggregate( x = score(a[ol[ol[,2] == i,1]]), by = list(cut(1:length(ol[ol[,2] == i,1]), breaks = 100)), mean)$x
mat[i,1:length(x)] <- x

ourMeans <- c(ourMeans, mean(score(a[ol[ol[,2] == i,1]])))
ourMedians <- c(ourMedians, median(score(a[ol[ol[,2] == i,1]])))
}

image(t(mat))

hist(ourMeans, freq = F, breaks = 20)
hist(ourMedians, col = 2, density = 0, add = F, freq = F, breaks = 20)

# when it comes to exons it seems that bases are conserved or not conserved.
plot(score(a[ol[ol[,2] == 8,1]]))

# funny we seem to be picking up different scores in the UTRs.

# really need to write something in GO that can get us the data we need. 


plot(score(a[ol[ol[,2] == 8,1]]), xlim=c(2500,2600), type = "l")
abline(v=(seq(2,3000,3)))




