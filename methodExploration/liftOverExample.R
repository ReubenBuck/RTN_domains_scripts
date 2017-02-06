

#### good to know that this tool is available in the R environment.

### It can be helpful in building the pipeline

#if (!exists("cur")) load("cur.rda")

## ----lkOne,eval=FALSE----------------------------------------------------
## library(gwascat)

## ----lkcur---------------------------------------------------------------


## ----getch---------------------------------------------------------------
library(rtracklayer)
ch = import.chain("data/chainAlignments/liftOver/hg19/mm9/hg19ToMm9.over.chain")
ch
str(ch[[1]])

## ----dolift--------------------------------------------------------------
#seqlevelsStyle(cur) = "UCSC"  # necessary
test.gr <- exon.gr[10:1000]
cur19 = liftOver(test.gr, ch)
#unlist(cur19)

#lapply(cur19, length)



class(cur19)


# we cna map regions across within the R environment
# however they are subject to multi mapping
# we do findOut when ther is no overlap

plot(log10(width(test.gr)), log10(sapply(lapply(cur19, width), sum) ))

# functions to get specific regions.
s <- sapply(lapply(X = cur19, FUN = start), min)
e <- sapply(lapply(X = cur19, FUN = end), max)

points(log10(width(test.gr)), log10(e-s), col = 2)
abline(a = 0, b=1)

# there are methods we can implement to make sure we get a good mapping ratio.
# size stays similar
# percentage of mapped bases across the interval. 
# we can give our own states and actually detect things that are lost. 



start(cur19)

## ----ul------------------------------------------------------------------
cur19 = unlist(cur19)
genome(cur19) = "mm9"
cur19 = new("gwaswloc", cur19)
cur19

## ----lkloss--------------------------------------------------------------
length(cur)-length(cur19)
setdiff(cur$SNPs, cur19$SNPs)

