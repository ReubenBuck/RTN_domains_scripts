# so with this script will be used to bin genomes and create repeat data objects with binned information 

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")

library(optparse)
library(dplyr)
library('spdep')
library(GenomicRanges)

option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path to genome", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="out.txt", 
              help="genome name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

genome <- opt$genome
path <- opt$path
print(opt$genome)
print(opt$path)
print(opt$out)

path = "~/Desktop/RTN_domains/R_objects/"
genome = "hg19"

load(file = paste(path,genome,".RData", sep =""))
load(file = "Desktop/RTN_domain/R_objects/hg19.RData")


head(rep)
tail(rep)


#Next we need to generate bin data 
# use the function to generate reasonable sized bins

repGR <- GRanges(seqnames = Rle(rep$genoChr),
                 ranges = IRanges(start = rep$genoStart, end = rep$genoEnd),
                 repUID = rep$repUID)


# genreate a table of repUIDs and binIDs for bins at different sizes
# so appropriate bin sizes 
sizes <- as.integer(c(50000, 100000, 200000,500000, 1000000, 1500000,2000000))

binList <- binned.genome.reader(genome = "hg19", bin.size = sizes, keep.rate = .5)

for(i in 1:length(binList)){binList[[i]]$binID <- paste(binList[[i]]$chr,":",binList[[i]]$start, "-",binList[[i]]$end, sep = "" )}

# cool now we have our bins
head(binList$hg19_bin.size_2000000)


# lets bin the data first 


bin <- binList$hg19_bin.size_200000

head(bin)

binGR <- GRanges(seqnames = Rle(bin$chr), 
                 ranges = IRanges(start = bin$start, end = bin$end),
                 binID = bin$binID)


OL <- as.matrix(findOverlaps(repGR, binGR, select = "first"))

OL[duplicated(OL[,1]),]

binningRep <- cbind(rep$repUID, bin$binID[OL])

#rownames(binningRep) = rep$repUID

# there will be two kinds of edge cases
# how to we want to assign those
# majority assignment seems like the most logical
# however first match seems easyest


# how can we get the binning rep into the bin data we need.
# Once we have the weight list we just need a single vector of X
# basically all the insertions of that TE


a <- rep %>% 
  group_by(binned = factor(x = as.character(binningRep[,2]),levels = unique(binningRep[,2])), repGroup) %>%
  summarise(insertionN = n_distinct(repID)) %>%
  dcast(formula = binned ~ repGroup, value.var = "insertionN") %>%
  collect

b <- a

b[is.na(b)] <- 0

region <- "chr2:"
xlim = c(100,300)

plot(smooth(b$ancient[grep(region, b$binned)]), type = "l", xlim = xlim)
par(new=TRUE)
plot(smooth(b$new_SINE[grep(region, b$binned)]), col = 2, type = "l", xlim = xlim)
par(new=TRUE)
plot(smooth(b$new_L1[grep(region, b$binned)]), col = 4, type = "l", xlim = xlim)
par(new=TRUE)
plot(smooth(b$old_L1[grep(region, b$binned)]), col = 3, type = "l", xlim = xlim)




a[[1]]
# we can use the start position for each chromosome to build a distance matrix
# every numbr over a certain value we can remove
# from there we chan build our neighborhood lists

bin <- binList$hg19_bin.size_2000000

bin <- bin[bin$chr == "chr1",]
dim(bin)
newStart <- (bin$start + 50000 - 1)/50000
newStartDist <- as.matrix(dist(newStart))
newStartDist[1:10,1:10]
newStartDist[newStartDist > 6] = 0
newStartDist[1:10,1:10]

a <- newStartDist
a[newStartDist > 0] <- 1/newStartDist[newStartDist > 0]
a[newStartDist > 0] <- 2/(2^newStartDist[newStartDist > 0])


newStartDist[newStartDist > 0] = abs(newStartDist[newStartDist > 0] - 7)
newStartDist[1:10,1:10]


weightMaker <- function(x,k,model){
  # x is a matrix with the distnaces neightbors
  # select the number of neighbors willing to be considered
  # apply one of several pre-selected weight models
  
  # 
  
  if(model == "nearestNeightbor"){
    
  }
  
  
}



# at least we have a way to generate the distance matrix 
# however we are not sure what kind of distance matrix would be best
# we are also not sure what outcomes we could expect form a better distance matrix
# the idea is that we capture hot spots at the best resolution. 

# therefore the best matrix is the one that can explian the most at the smaller bin sizes
# get the best shifts at the boundaries

# what is the plan do we want to gi



m2 <- Matrix(0, nrow = 1000, ncol = 1000, sparse = TRUE)
m2[1:10,1:10] <- a[1:10,1:10]
m2[1:10,1:10]

lw <- mat2listw(m2[1:10,1:10])



