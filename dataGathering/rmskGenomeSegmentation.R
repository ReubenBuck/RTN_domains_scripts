# so with this script will be used to bin genomes and create repeat data objects with binned information 


rm(list = ls())

library(optparse)
option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path to genome", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NA, 
              help="genome name", metavar="character"),
  make_option(c("-o", "--outPathRobject"), type="character", default=getwd(), 
              help="output path for processed repeats, default = current dir", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

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

# 
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
sizes <- as.integer(c(50000,100000,250000, 500000, 750000, 1000000, 1500000,2000000))

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
  
  
  rownames(bin) <- 1:nrow(bin)
  mat <- Matrix(data = 0, nrow = nrow(bin), ncol= nrow(bin))
  for(i in 1:length(unique(bin$chr))){
    binChr <- bin[bin$chr == unique(bin$chr)[i],]
    neighborDist <- as.matrix(dist((binChr$start + sizeSelect - 1)/sizeSelect))
    neighborDist[neighborDist > 1] = 0
    mat[as.integer(rownames(binChr)), as.integer(rownames(binChr))] <- neighborDist
  }
  diag(mat) <- 1
  
#  lw <- mat2listw(mat)
 
  
  
  repDataList <- list(repSummary = repSummary,bin = bin, repMap = repMap ,neighborMat = mat)

saveName = paste("repData_",genome,"_",sizeSelect,".RData", sep = "")
save(repDataList, file = paste(outPath, saveName,sep = "/"))

}

# objects to save is the mapping, the matrix, the summaries


# 
# 
# 
# lw <- mat2listw(mat)
# 
# 
# moran.plot(x = (repSummary$new_L1/bin$Known) * 1e6, listw = lw,labels = FALSE, pch = 16, cex = .3)
# moran(x = repSummary$new_L1/(bin$Known)* 1e6, listw = lw,n = nrow(b)-1, S0 = sum(mat))
# moran.mc(x = repSummary$new_SINE/(bin$Known)* 1e6, listw = lw, nsim = 100)
# moran.mc(x = repSummary$ancient/(bin$Known)* 1e6, listw = lw, nsim = 100)
# 
# 
# lmAncient <- localmoran(x = repSummary$ancient, listw = lw)
# lmNew_SINE <- localmoran(x = repSummary$new_SINE, listw = lw)
# 
# lgAncient <- localG(x = repSummary$ancient/(bin$Known)* 1e6, listw = lw)
# lgNew_SINE <- localG(x = repSummary$new_SINE/(bin$Known)* 1e6, listw = lw)
# lgNew_L1 <- localG(x = repSummary$new_L1/(bin$Known)* 1e6, listw = lw)
# lgOld_L1 <- localG(x = repSummary$old_L1/(bin$Known)* 1e6, listw = lw)
# 
# 
# layout(matrix(c(1,2,3,4),nrow = 4))
# par(mar=c(2,5,2,5))
# chr.choice = "chr1"
# xlim = c(10e6,220e6)
# plot(bin$start[bin$chr == chr.choice] + (sizeSelect/2),lgAncient[bin$chr == chr.choice], col = "darkblue",type = "l", xlim = xlim);abline(h=3);grid()
# plot(bin$start[bin$chr == chr.choice] + (sizeSelect/2),lgNew_SINE[bin$chr == chr.choice], col = "darkgreen", type = "l", xlim = xlim);abline(h=3);grid()
# plot(bin$start[bin$chr == chr.choice] + (sizeSelect/2),lgNew_L1[bin$chr == chr.choice],col = "purple",type = "l", xlim = xlim);abline(h=3);grid()
# plot(bin$start[bin$chr == chr.choice] + (sizeSelect/2),lgOld_L1[bin$chr == chr.choice],col = "red",type = "l", xlim = xlim);abline(h=3);grid()
# 
# layout(matrix(c(1,2),nrow = 2))
# chr.choice = "chr20"
# xlim = c(30e6,40e6)
# plot(bin$start[bin$chr == chr.choice] + (sizeSelect/2),lgAncient[bin$chr == chr.choice], col = "darkblue",type = "l", xlim = xlim);abline(h=3);grid()
# lines(bin$start[bin$chr == chr.choice] + (sizeSelect/2), lmAncient[bin$chr == chr.choice,4], col = "orange",type = "l", xlim = xlim)
# plot(bin$start[bin$chr == chr.choice] + (sizeSelect/2),lgNew_SINE[bin$chr == chr.choice], col = "darkgreen",type = "l", xlim = xlim);abline(h=3);grid()
# lines(bin$start[bin$chr == chr.choice] + (sizeSelect/2), lmNew_SINE[bin$chr == chr.choice,4], col = "orange",type = "l", xlim = xlim)
# 
# layout(matrix(c(1,2),nrow = 2))
# chr.choice = "chr3"
# xlim = c(10e6,20e6)
# plot(bin$start[bin$chr == chr.choice] + (sizeSelect/2),lgAncient[bin$chr == chr.choice], col = "darkblue",type = "l", xlim = xlim);abline(h=3);grid()
# lines(bin$start[bin$chr == chr.choice] + (sizeSelect/2),lgNew_SINE[bin$chr == chr.choice], col = "darkgreen",type = "l", xlim = xlim)
# plot(bin$start[bin$chr == chr.choice] + (sizeSelect/2), lmAncient[bin$chr == chr.choice,4], col = "darkblue",type = "l", xlim = xlim);abline(h=3);grid()
# lines(bin$start[bin$chr == chr.choice] + (sizeSelect/2), lmNew_SINE[bin$chr == chr.choice,4], col = "darkgreen",type = "l", xlim = xlim)
# 
# 
# # the idea here is that clusters 
# bluered <- colorRampPalette(colors = c("darkblue","blue","green","orange","red"))
# 
# h2 <- hist2d((repSummary[,c("new_SINE", "ancient")]/bin$Known * 10e6), nbins=50, 
#              col=bluered(30), xlab = "Alu per Mb", ylab = "MIR/L2 per Mb", xlim = c(0,20000), ylim = c(0,20000),
#              main = paste("bin size =",sizeSelect), zlim = c(0,log10(nrow(repSummary)/50)),
#              FUN=function(x) log10(length(x)))
# abline(v = sort(repSummary$new_SINE/bin$Known * 10e6)[.9 * nrow(bin)])
# abline(h = sort(repSummary$ancient/bin$Known * 10e6)[.9 * nrow(bin)])
# legend("topright",legend = "90th percentile", lty = 1, bty = "n")
# 
# 
# # so the idea is for each element to use the spatial data that best describes the spatial relationships
# # this is determined by where these paterns most emerge. 
# 
# # maybe remove the X chromasome form the analysis because it is a speceial case
# 
# 
# newStartDist[newStartDist > 0] = abs(newStartDist[newStartDist > 0] - 7)
# newStartDist[1:10,1:10]
# 
# 
# weightMaker <- function(x,k,model){
#   # x is a matrix with the distnaces neightbors
#   # select the number of neighbors willing to be considered
#   # apply one of several pre-selected weight models
#   
#   # 
#   
#   if(model == "nearestNeightbor"){
#     
#   }
#   
#   
# }
# 
# 

# at least we have a way to generate the distance matrix 
# however we are not sure what kind of distance matrix would be best
# we are also not sure what outcomes we could expect form a better distance matrix
# the idea is that we capture hot spots at the best resolution. 

# therefore the best matrix is the one that can explian the most at the smaller bin sizes
# get the best shifts at the boundaries

# what is the plan do we want to gi

# we can get pretty good hotspots at reasonable resolutions.
# we go down to a point where clusters no longer overlap.

# there's probably a way to measure the optimal at different scales for different elements


