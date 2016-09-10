# so with this script will be used to bin genomes and create repeat data objects with binned information 

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")

library(optparse)
library(dplyr)

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

path = "~/Desktop/RTN_domains/R_objects/rmskTables/"
genome = "hg19"

load(file = paste(path,genome, sep =""))



head(rep)
tail(rep)


#Next we need to generate bin data 
# use the function to generate reasonable sized bins


# genreate a table of repUIDs and binIDs for bins at different sizes
# so appropriate bin sizes 
sizes <- c(50000, 100000, 200000,500000, 1000000, 150000,2000000)

binList <- binned.genome.reader(genome = "hg19", bin.size = sizes, keep.rate = .5)

# cool now we have our bins

tail(binList$hg19_bin.size_50000)

# we can use the start position for each chromosome to build a distance matrix
# every numbr over a certain value we can remove
# from there we chan build our neighborhood lists

bin <- binList$hg19_bin.size_50000

bin <- bin[bin$chr == "chr1",]
dim(bin)
newStart <- (bin$start + 50000 - 1)/50000
newStartDist <- as.matrix(dist(newStart))
newStartDist[newStartDist > 10] = 0
newStartDist[newStartDist > 0] = abs(newStartDist[newStartDist > 0] - 11)


# at least we have a way to generate the distance matrix 
# however we are not sure what kind of distance matrix would be best
# we are also not sure what outcomes we could expect form a better distance matrix
# the idea is that we capture hot spots at the best resolution. 

# therefore the best matrix is the one that can explian the most at the smaller bin sizes
# get the best shifts at the boundaries


