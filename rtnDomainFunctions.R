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


# Function designed to extract corresponding coordinates between two species according to a reference
# the reference hotspots
# the query hotspots in the reference genome
# the reference hotspots mapped to the query genome

extractCorrespondingHotspots <- function(ref.gr, que_ref.gr, ref_que.gr, que.gr, repGroups){ 
  
  repList <- NULL
  for(i in 1:length(repGroups)){
    repList <- c(repList, list(NULL))
  }
  names(repList) = repGroups
  
  queList <- list(con = repList, mid = repList, dif= repList)
  refList <- list(con = repList, mid = repList, dif= repList)
  
  hotspotList <- list(ref = refList, que = queList)
  
  
  for(i in 1:length(repGroups)){
    sampRef.gr <- ref.gr[elementMetadata(ref.gr)$repGroup == repGroups[i]]
    sampQue_Ref.gr <- que_ref.gr[elementMetadata(que_ref.gr)$repGroup == repGroups[i]]
    
    sampQue.gr <- que.gr[elementMetadata(que.gr)$repGroup == repGroups[i]]
    sampRef_Que.gr <- ref_que.gr[elementMetadata(ref_que.gr)$repGroup == repGroups[i]]
    
    
    olRef <- as.matrix(findOverlaps(sampRef.gr,sampQue_Ref.gr ))
  
    olRef_Que <- as.matrix(findOverlaps(sampRef_Que.gr,sampQue.gr))
   # dfRepGroup <- data.frame(elementMetadata(sampRef.gr)$repGroup[ol[,1]], 
    #                       elementMetadata(sampQue_Ref.gr)$repGroup[ol[,2]])
    
    conID <- unique(c(elementMetadata(sampRef_Que.gr[olRef_Que[,1]])$hotspotID, 
      elementMetadata(sampRef.gr[olRef[,1]])$hotspotID))
    
    conHotspot.gr <- sampRef.gr[elementMetadata(sampRef.gr)$hotspotID %in% conID]
    
   # conHotspot.gr <- ref.gr[unique(ol[dfRepGroup[,1] == repGroups[i] & dfRepGroup[,1] == repGroups[i], 1])]
    
    tabRef <- table(elementMetadata(ref.gr[elementMetadata(ref.gr)$repGroup == repGroups[i]])$hotspotGroup)
    tabCon <- table(elementMetadata(unique(conHotspot.gr))$hotspotGroup)
    
    tabConScores <- tabCon/tabRef[names(tabRef) %in% names(tabCon)][names(tabCon)]
    
    # here we get the names of groups
    conGroups <- names(tabConScores[tabConScores >= .7])
    midGroups <- names(tabConScores[tabConScores < .7 & tabConScores >= .3])
    difGroups <- c(names(tabConScores[tabConScores < .3]), names(tabRef[!(names(tabRef) %in% names(tabCon))]))
    
    # here we extract the coordinates in the que species
    conDomainQue <- ref_que.gr[elementMetadata(ref_que.gr)$hotspotGroup %in% conGroups]
    midDomainQue <- ref_que.gr[elementMetadata(ref_que.gr)$hotspotGroup %in% midGroups]
    difDomainQue <- ref_que.gr[elementMetadata(ref_que.gr)$hotspotGroup %in% difGroups]
    
    hotspotList$que$con[[repGroups[i]]] <- conDomainQue
    hotspotList$que$mid[[repGroups[i]]] <- midDomainQue
    hotspotList$que$dif[[repGroups[i]]] <- difDomainQue
    
    
    # here we extract the coordinates in the ref species
    conDomainRef <- ref.gr[elementMetadata(ref.gr)$hotspotGroup %in% conGroups]
    midDomainRef <- ref.gr[elementMetadata(ref.gr)$hotspotGroup %in% midGroups]
    difDomainRef <- ref.gr[elementMetadata(ref.gr)$hotspotGroup %in% difGroups]
    
    hotspotList$ref$con[[repGroups[i]]] <- conDomainRef
    hotspotList$ref$mid[[repGroups[i]]] <- midDomainRef
    hotspotList$ref$dif[[repGroups[i]]] <- difDomainRef
    
    
    # filtering step to make sure that groups are the correct size in both ref and que
    # difference between actuall and assumed width must be either 150kb or 20% within actuall
    for(conState in c("con", "mid", "dif")){
      KeepList <- NULL
      for(genome in c("que","ref")){
        df <- as.data.frame(hotspotList[[genome]][[conState]][[repGroups[i]]])
        groupSum <- summarise(group_by(df, hotspotGroup), n(), min(start), max(end))
        groupSum$sizeDif = abs((groupSum$`n()` * 5e4) - (groupSum$'max(end)' - groupSum$'min(start)'))
        keepers <- filter(groupSum, sizeDif < 1.5e5 )$hotspotGroup 
        KeepList <- c(KeepList, keepers)
      }
      KeepList <- KeepList[duplicated(KeepList)]
      keepRef <- elementMetadata(hotspotList[["ref"]][[conState]][[repGroups[i]]])$hotspotGroup%in%KeepList
      keepQue <- elementMetadata(hotspotList[["que"]][[conState]][[repGroups[i]]])$hotspotGroup%in%KeepList
      
      
      hotspotList[["ref"]][[conState]][[repGroups[i]]] = hotspotList[["ref"]][[conState]][[repGroups[i]]][keepRef]
      hotspotList[["que"]][[conState]][[repGroups[i]]] = hotspotList[["que"]][[conState]][[repGroups[i]]][keepQue]
    }
    

  }
  
  return(hotspotList)
  
}

# may need to do a final check to see if any of the non conserved hot spots are actually just hot spots that couldn't map.




# for extracting insertion rates, 
# the next important step is to make sure we only get stuff that can be directly compared from ref to que

extractInsertionRates <- function(refGenome.gr, queGenome.gr, refCorresponding, repGroups, minoverlap = 1){
  # input reference and query binned granges with TE content in metadata
  insertionRate <- NULL
  for(rep in repGroups){
    for(conState in c("con", "mid", "dif")){
      conStateInsertionRate <- NULL
      for(genome in c("ref", "que")){
        genome.gr <- get(paste(genome, "Genome.gr", sep = ""))
        refCor <- refCorresponding[[genome]][[conState]][[rep]]
        ol <- as.matrix(findOverlaps(genome.gr,refCor, minoverlap = minoverlap))
        pOl <- width(pintersect(genome.gr[ol[,1]], refCor[ol[,2]]))/width(refCor[ol[,2]])
        ol <- ol[pOl > .5,]
        mD <- data.frame(chr = seqnames(genome.gr[ol[,1]]),
                         start = start(genome.gr[ol[,1]]), end = end(genome.gr[ol[,1]]),
                         elementMetadata(genome.gr[ol[,1]])[c(rep, "binID", "known")],
                         elementMetadata(refCor[ol[,2]])[c("hotspotID", "hotspotGroup")] )
        colnames(mD)[4] <- c("insertionRate")
        df <- data.frame(mD,repGroup = rep, genome = genome, conState = conState)
        conStateInsertionRate <- rbind(conStateInsertionRate, df)
      }
      keep <- conStateInsertionRate$hotspotID[duplicated(conStateInsertionRate$hotspotID)]
      insertionRate <- rbind(insertionRate, conStateInsertionRate[(conStateInsertionRate$hotspotID %in% keep),])
    }
  }
  return(insertionRate)
}



# could write a second funciton for group.
# get the groups first, do the overlap 
# then aggregate the overlapped scores

#refCor <- refCorresponding[[genome]][[conState]][[rep]]
#aggregate(x = start(refCor), by = elementMetadata(refCor)
#df <- data.frame(insertionRate = as.matrix(elementMetadata(genome.gr[ol[,1]])[rep]), elementMetadata(refCor[ol[,2]])[c("hotspotID", "hotspotGroup")] )

# the number of overlapping bases between a dat.gr and overlap.gr
# the results are in reference to dat.gr
overlapingBases <- function(dat.gr, overlap.gr){
  m <- as.matrix(findOverlaps(dat.gr, overlap.gr)) 
  p <- pintersect(dat.gr[m[,1]], overlap.gr[m[,2]])
  m2 <- as.matrix(findOverlaps(dat.gr, p)) 
  agg <- aggregate(width(p[m2[,2]]), by = list(m2[,1]), FUN= sum)
  res = rep(0, length(dat.gr))
  res[agg$Group.1] = agg$x
  return(res)
}

