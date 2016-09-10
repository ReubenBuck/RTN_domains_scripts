binned.genome.reader <- function(genome, bin.size, keep.rate){
  
  # getting chrom info 
  web <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/", genome, "/database/chromInfo.txt.gz", sep = "")
  con <- gzcon(url(web))
  txt <- readLines(con)
  chrom_info <- read.delim(textConnection(txt), header = FALSE)
  
  
  ses <- browserSession("UCSC")
  genome(ses) <- genome
  
  data.types <- c("gaps")
  track.name <- c("gap")
  table.name <- c("gap")
  
  
  for(i in 1:length(data.types)){
    dat <- getTable(
      ucscTableQuery(
        ses, 
        track = track.name[i],
        table = table.name[i]
      )
    )
    assign(data.types[i], dat)    		
  }
  
  
  
  gaps.gr <- GRanges(seqnames = Rle(gaps$chrom),
                     ranges = IRanges(start = gaps$chromStart, end = gaps$chromEnd)
  )		 
  
  Final <- NULL
  for(z in 1:length(bin.size)){
    chrom_info1 <- chrom_info[chrom_info[,2] > bin.size[z],]
    
    bins <- NULL
    for( i in 1:dim(chrom_info1)[1]){
      start <- seq(from = 1, to = chrom_info1[i,2], by = bin.size[z])
      end <- c(seq(from = bin.size[z], to = chrom_info1[i,2], by = bin.size[z]), chrom_info1[i,2])
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


