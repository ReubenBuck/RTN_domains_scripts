binned.genome.reader <- function(genome, bin.size, keep.rate){
#  library(rtracklayer)
  
  # getting chrom info 

  my_db <- src_mysql(genome,
                     host = "genome-mysql.cse.ucsc.edu", 
                     user = "genomep", 
                     password = "password")
  
  chrom_info <- data.frame(tbl(my_db, "chromInfo"))
  chrom_info$size <- as.integer(chrom_info$size)
  
  # gaps may not work depending on the database
  # gaps are sometimes split across chromosome 
  
  if(any(src_tbls(my_db) == "gap")){
    gaps <- data.frame(collect(tbl(my_db, "gap")))
  }else{
    gaps <- NULL
    for(i in 1:nrow(chrom_info)){
      gaps <- rbind(gaps, data.frame(tbl(my_db, paste(chrom_info$chrom[i], "gap", sep = "_"))))
    }
  }
  
  gaps$chromStart = as.integer(gaps$chromStart)
  gaps$chromEnd = as.integer(gaps$chromEnd)
  gaps$size = as.integer(gaps$size)
  

  gaps.gr <- GRanges(seqnames = Rle(gaps$chrom),
                     ranges = IRanges(start = gaps$chromStart, end = gaps$chromEnd)
  )		 
  
  Final <- NULL
  for(z in 1:length(bin.size)){
    chrom_info1 <- chrom_info[ chrom_info$size > bin.size[z],]
    
    bins <- NULL
    for( i in 1:dim(chrom_info1)[1]){
      start <- seq(from = 1, to = chrom_info1[i,2], by = bin.size[z])
      end <- c(seq(from = bin.size[z], to = chrom_info1[i,2], by = bin.size[z]), chrom_info1[i,2])
      end <- end[1:length(start)]
      pre.bin <- data.frame(chr = chrom_info1[i,1], start = start, end = end)
      bins <- rbind(bins, pre.bin)
    }
    
    bins$Known <- bins$end - bins$start + 1
    bins.gr <- GRanges(seqnames = Rle(bins$chr), 
                       ranges = IRanges(start = bins$start, end = bins$end - 1)
    )
    
    gap.int.gr <- intersect(gaps.gr, bins.gr)	   
    gap.bin.ol <- as.matrix(
      findOverlaps(bins.gr, gap.int.gr)
    )	   
    gap.size <- data.frame(bin = gap.bin.ol[,1], width = width(gap.int.gr[gap.bin.ol[,2]]))
    gap.size.bin <- aggregate(gap.size[,2],list(bin = gap.size$bin), sum)
    bins$Known[gap.size.bin$bin] <- bins$Known[gap.size.bin$bin] - gap.size.bin$x + 1
    # so now we have all the known bp for bins for a genome of any bin size
    bins <- bins[bins$Known >= keep.rate*bin.size[z],1:4]
    rownames(bins) <- 1:nrow(bins)
    Final <- c(Final, list(bins))
  }
  names(Final) <- paste(genome,"bin.size", bin.size, sep = "_")
  return(Final)
}


# convert squere matrix into another squere matrix with a triangle in the middle
triangulate <- function(mat){
  mat[1:nrow(mat),] <- mat[nrow(mat):1,]
  triangle <- matrix(NA,nrow = nrow(mat),ncol = nrow(mat)*2)
  for(i in nrow(mat):2){
    dat <- (diag(mat[i:1,1:i]))
    triangle[(nrow(mat):1)[i],seq((nrow(mat):1)[i],by = 2,length.out = i)] <- (dat)
    triangle[(nrow(mat):1)[i],seq((nrow(mat):1)[i] + 1,by = 2,length.out = i)] <- (dat)
  }
  triangle[(nrow(mat):1)[1],seq((nrow(mat):1)[1] ,by = 2,length.out = 1)] <- rep(mat[1,1])
  triangle[(nrow(mat):1)[1],seq((nrow(mat):1)[1] + 1,by = 2,length.out = 1)] <- rep(mat[1,1])
  return(t(triangle))
}



# calculate insulation score for hi-c data
insulationScore <- function(mat, squareSize){
  mat[is.infinite(mat)] <- NA
  r <- c(1:squareSize)
  c <- c((squareSize + 1):(squareSize + squareSize))
  insulate = NULL
  for(i in 1:(nrow(mat)-(squareSize + squareSize - 1))){
    insulate <- c(insulate, mean(mat[r,c], na.rm = T))
    r <- r + 1
    c <- c + 1
  }
  return(insulate)
}

# report percentage of particular repeat hotspots in reference that overlap particular repeat hotspots in query
olHotspotSummary <- function(ref.gr, que.gr, repGroups){
  ol <- as.matrix(findOverlaps(ref.gr, que.gr))
  dfRepGroup <- data.frame(elementMetadata(ref.gr)$repGroup[ol[,1]], 
                           elementMetadata(que.gr)$repGroup[ol[,2]])
  dfDomainID <- data.frame(elementMetadata(ref.gr)$domainID[ol[,1]], 
                           elementMetadata(que.gr)$domainID[ol[,2]])
  
  # go through each combination and count how many unique regions
  olSumS1 <- matrix(nrow = 4, ncol = 4, dimnames = list(repGroups,repGroups))
  olTotals <- matrix(nrow = 4, ncol = 1, dimnames = list(repGroups,"totalOL"))
  
  for(i in 1:length(repGroups)){
    for(j in 1:length(repGroups)){
      olSumS1[i,j] <- length(unique(dfDomainID[dfRepGroup[,1] == repGroups[i] & dfRepGroup[,2] == repGroups[j],1])) / 
        length(ref.gr[elementMetadata(ref.gr)$repGroup == repGroups[i]])
    }
    olTotals[i] <- length(unique(dfDomainID[dfRepGroup[,1] == repGroups[i], 1] )) / length(ref.gr[elementMetadata(ref.gr)$repGroup == repGroups[i]])
  }
  return(cbind(olSumS1, olTotals))
}



# so it can be 


