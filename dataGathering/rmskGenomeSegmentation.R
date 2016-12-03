# so with this script will be used to bin genomes and create repeat data objects with binned information 


rm(list = ls())

library(optparse)
option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path to genome", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NA, 
              help="genome name", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NA, 
              help="distance from bin or size of bin, dist or size", metavar="character"),
  make_option(c("-o", "--outPathRobject"), type="character", 
              help="output path for processed repeats, default = current dir", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(!(opt$type == "size" | opt$type == "dist")){
  stop("must specify type of segmentation as size or dist")
}

genome <- opt$genome
path <- opt$path
outPath <- opt$outPathRobject
print(opt$genome)
print(opt$path)
print(opt$outPathRobject)

if(is.na(genome)){stop("need to specify genome")}
if(rev(strsplit(opt$path,split = "")[[1]])[1] != "/"){
  opt$inPath = paste(opt$inPath,"/",sep = "")
}
if(rev(strsplit(opt$outPathRobject,split = "")[[1]])[1] != "/"){
  opt$outPathRobject = paste(opt$outPathRobject,"/",sep = "")
}

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")


library(dplyr)
#library('spdep')
library(Matrix)
library(GenomicRanges)
library(reshape2)

# Comment these back after testing
# path = "~/Desktop/RTN_domains/R_objects/rmskTables/"
# outPath="~/Desktop/RTN_domains/R_objects/repMapData/canFam3"
# genome = "canFam3"

 
 
load(file = paste(path,genome,".RData", sep =""))


# so to make a rep obfect we need to have m


#Next we need to generate bin data 
# use the function to generate reasonable sized bins

repGR <- GRanges(seqnames = Rle(rep$genoChr),
                 ranges = IRanges(start = rep$genoStart, end = rep$genoEnd),
                 repUID = rep$repUID)


# genreate a table of repUIDs and binIDs for bins at different sizes
# so appropriate bin sizes 
if(opt$type == "size"){
  sizes <- as.integer(c(20000,50000,100000,250000, 500000, 750000, 1000000, 1500000,2000000))
  binDistances <- 1
  
}

if(opt$type == "dist"){
  sizes <- as.integer(50000)
  binDistances <- c(1,2,5,7,10,15)
}

binList <- binned.genome.reader(genome = genome, bin.size = sizes, keep.rate = .5)
for(i in 1:length(binList)){binList[[i]]$binID <- paste(binList[[i]]$chr,":",binList[[i]]$start, "-",binList[[i]]$end, sep = "" )}

# cool now we have our bins



#sizes = sizes[1:2]
for(i in 1:length(sizes)){
  
  sizeSelect <- sizes[i]
  binObject <- paste(genome,"_bin.size_", as.integer(sizeSelect), sep = "")
  
  bin <- binList[[binObject]]
  
  
  #segmenterSummeriser <- function(bin, rep, groupingFactor)
  
  # some quick tests
  # make sure objects are of the correct types, a good way to make them as S4 classes 
  # check they have the correct rownames
  
  binGR <- GRanges(seqnames = Rle(bin$chr), 
                   ranges = IRanges(start = bin$start, end = bin$end),
                   binID = bin$binID)
  OL <- as.matrix(findOverlaps(repGR, binGR, select = "first"))
  
  
  repMap <- data.frame(repUID = rep$repUID, binID = bin$binID[OL])

  repSummary <- rep %>% 
    group_by(binID = factor(x = as.character(bin$binID[OL]),levels = unique(bin$binID)), repGroup) %>%
    summarise(insertionN = n_distinct(repID)) %>%
    dcast(formula = binID ~ repGroup, value.var = "insertionN") %>%
    collect
  
  repSummary <- repSummary[1:(nrow(repSummary)-1),1:(ncol(repSummary)-1)]
  
  if(nrow(repSummary) != nrow(bin)){
    binAdd <- data.frame(
      matrix(data = NA,
             nrow = sum(!(bin$binID %in% repSummary$binID)), 
             ncol = ncol(repSummary),
             dimnames = list(1:sum(!(bin$binID %in% repSummary$binID)), colnames(repSummary)))
    )
    binAdd$binID <- factor(bin[!(bin$binID %in% repSummary$binID),"binID"], levels = bin$binID)
    repSummary <- rbind(repSummary, binAdd)
    repSummary <- repSummary[order(repSummary$binID),]
    rownames(repSummary) <- 1:nrow(bin)
  }
  repSummary[is.na(repSummary)] <- 0
  
  neighborList <- NULL
  for(matDistance in binDistances){
    rownames(bin) <- 1:nrow(bin)
    mat <- Matrix(data = 0, nrow = nrow(bin), ncol= nrow(bin))
    for(i in 1:length(unique(bin$chr))){
      binChr <- bin[bin$chr == unique(bin$chr)[i],]
      neighborDist <- as.matrix(dist((binChr$start + sizeSelect - 1)/sizeSelect))
      neighborDist[neighborDist > matDistance] = 0
      neighborDist[neighborDist > 0] = 1
      mat[as.integer(rownames(binChr)), as.integer(rownames(binChr))] <- neighborDist
    }
    diag(mat) <- 1
    
    neighborList <- c(neighborList, list(mat))
  }
  names(neighborList) <- paste("binDistance", binDistances, sep = "")
  #  lw <- mat2listw(mat)
 
  
  
  repDataList <- list(repSummary = repSummary,bin = bin, repMap = repMap ,neighborMat = neighborList)

saveName = paste("repData_",genome,"_",sizeSelect,"_",opt$type,".RData", sep = "")
save(repDataList, file = paste(outPath, saveName,sep = "/"))

}

