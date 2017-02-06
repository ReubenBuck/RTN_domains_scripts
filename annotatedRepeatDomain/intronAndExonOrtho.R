###
### Using the mysql Api and our species we can download the repeat 
### information and we can look at the genes
###
###

rm(list = ls())

setwd("~/Desktop/RTN_domains/")


library(GenomicRanges)
library(dplyr)
library(RMySQL)


repGroups = c("ancient", "new_SINE", "new_L1", "old_L1")
repCols = c("darkblue", "aquamarine3", "purple", "red")
snames <- c(s1name = "hg19", s2name = "mm9", s3name = "canFam3")


ints <- read.table(file = "data/repeatHotspot/intersect.txt", header= TRUE)

orthoDf <- read.table(file = "data/orthoAnalysis/orhtoList/ortholist.txt", header = TRUE)

# get all the genome data.
genomeInfo <- c(rep(list(NA),length(snames)))
names(genomeInfo) <- snames

for(s in snames){
  
  mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = s)
  
  seqInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")
  
  ensGene <- dbGetQuery(mychannel, "SELECT * FROM ensGene;")
  ensGene <- ensGene[ensGene$name %in% orthoDf[,s],]
  ensGene.gr <- GRanges(seqnames = Rle(ensGene$chrom), 
                        ranges = IRanges(start = ensGene$txStart, end = ensGene$txEnd, names = ensGene$name),
                        strand = Rle(ensGene$strand), ensGene[,7:ncol(ensGene)])
  
  seqlevels(ensGene.gr) = seqInfo[,1]
  isCircular(ensGene.gr) = seqlevels(ensGene.gr) == "chrM"
  seqlengths(ensGene.gr) = seqInfo[,2]
  genome(ensGene.gr) = s
  
  exon.gr <- GRanges(seqnames = Rle(do.call(rep, list(x = ensGene$chrom, times = ensGene$exonCount))),
                     ranges = IRanges(start = unlist(strsplitAsListOfIntegerVectors(x = ensGene$exonStarts,sep = ",")),
                                      end = unlist(strsplitAsListOfIntegerVectors(x = ensGene$exonEnds,sep = ",")) ),
                     names = do.call(rep, list(x = ensGene$name, times = ensGene$exonCount)),
                     strand = Rle(do.call(rep, list(x = ensGene$strand, times = ensGene$exonCount)))
  )
  seqlevels(exon.gr) = seqInfo[,1]
  seqinfo(exon.gr) <- seqinfo(ensGene.gr)
  
  #exonList.gr <- makeGRangesListFromDataFrame(df = as.data.frame(exon.gr),split.field = "names")
  
  ce <- cumsum(ensGene$exonCount)
  cs <- ce + 1 - ensGene$exonCount
  intron.gr <- GRanges(seqnames = Rle(do.call(rep, list(x = ensGene$chrom, times = ensGene$exonCount))[-cs]),
                     ranges = IRanges(start = unlist(strsplitAsListOfIntegerVectors(x = ensGene$exonEnds,sep = ","))[-ce] ,
                                      end = unlist(strsplitAsListOfIntegerVectors(x = ensGene$exonStarts,sep = ","))[-cs] ),
                     names = do.call(rep, list(x = ensGene$name, times = ensGene$exonCount))[-cs],
                     strand = Rle(do.call(rep, list(x = ensGene$strand, times = ensGene$exonCount))[-cs])
  )
  seqlevels(intron.gr) = seqInfo[,1]
  seqinfo(intron.gr) <- seqinfo(ensGene.gr)
  
  # one way to look at is, we can keep all the intro and exon info for each ortho
  
 # exon.gr <- reduce(exon.gr)
  
  
  
  
  genomeInfo[[s]] = list(ensGene = ensGene.gr,exons = exon.gr, introns = intron.gr)
  
}

genomeInfoAll = list(orthoInfo = orthoDf, ensGene = genomeInfo)

save(genomeInfoAll, file = "R_objects/ensGene/genomeInfoAll.RData")

### next time I can load the R object and there's all the ortho data to do some analysis.
### the next easiest thing to compare would be exonic phastcon scores. 


### plot these results to bring to the meeting tomorrow

scols <- c(1,2,4)
names(scols) <- snames

pdf(file = "plots/hotspotStats/intronSize/intronIntersect.pdf")
for(s in snames){
  name1 = s
  
  
  otherNames = snames[snames != name1]
  #name2 = 1
  
  dat <- read.table(file = paste("data/repeatHotspot/",name1, "/",name1, "_", otherNames[1], "_conDif.txt", sep = ""), header = T)
  
  
  datRef.gr <- GRanges(seqnames = Rle(dat$chr[dat$genome == "ref"]),
                       ranges = IRanges(start = dat$start[dat$genome == "ref"], end = dat$end[dat$genome == "ref"]))
  mcols(datRef.gr) <- dat[dat$genome == "ref", 4:ncol(dat)]
  seqlevels(datRef.gr) <- seqlevels(genomeInfo[[name1]]$ensGene)
  seqinfo(datRef.gr) <- seqinfo(genomeInfo[[name1]]$ensGene)
  
  
  olRef <- as.matrix(findOverlaps(query = genomeInfo[[name1]]$ensGene, subject = datRef.gr))
  dfRef <- data.frame(ensGene = names(genomeInfo[[name1]]$ensGene[olRef[,1]]), mcols(datRef.gr[olRef[,2]]) )
  
  
  
  m <- merge(x = dfRef, y = orthoDf, by.x = 1, by.y = (1:4)[colnames(orthoDf) == name1] )
  
  
  for(r in repGroups){
    repChoice = r
    m2 <- m[ m$hotspotID %in% ints$domains[ints$genome == name1 & ints$repGroup == repChoice],]
    
    
    layout(matrix(1:2, nrow = 2))
    par(mar = c(2,2,2,2), oma = c(5,5,5,5))
    
    if(length(genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% 
                                          unique(m2[m2$repGroup == repChoice & m2$conState == "con", "ensGene"])]) > 2){
      
      plot(density(log10(width( genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(m2[m2$repGroup == repChoice & m2$conState == "con", "ensGene"])])), 
      ), xlim = c(1,6), ylim= c(0,.8),main = "", lwd = 2, lty = 1, col = scols[name1])
      lines(density(log10(width( genomeInfo[[otherNames[1]]]$introns[mcols(genomeInfo[[otherNames[1]]]$introns)$name %in% unique(m2[m2$repGroup == repChoice & m2$conState == "con", otherNames[1] ] ) ])), 
      ), col = scols[otherNames[1]], lwd = 2, lty = 2)
      lines(density(log10(width( genomeInfo[[otherNames[2]]]$introns[mcols(genomeInfo[[otherNames[2]]]$introns)$name %in% unique(m2[m2$repGroup == repChoice & m2$conState == "con", otherNames[2]])])), 
      ), col = scols[otherNames[2]], lwd = 2, lty = 2)
      
      legend("topright", legend = c(name1, otherNames), col = scols[c(name1, otherNames)], lwd = 2, bty = "n", lty = c(1,2,2))
      title(main = "shared")
      
      
      a <- genomeInfo[[otherNames[2]]]$introns[mcols(genomeInfo[[otherNames[2]]]$introns)$name %in% unique( m2[m2$repGroup == repChoice & m2$conState == "con", otherNames[2]]   ) ]
      legend("topleft", legend = paste("n =", length(a)), bty = "n")
    }else{
      plot.new()
    }
    
    if(length(genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% 
                                          unique(m2[m2$repGroup == repChoice & m2$conState == "dif", "ensGene"])]) > 2){
      
      plot(density(log10(width( genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(m2[m2$repGroup == repChoice & m2$conState == "dif", "ensGene"]) ])), 
      ), xlim = c(1,6), main = "", ylim= c(0,.8), lty = 1, lwd = 2, col = scols[name1])
      lines(density(log10(width( genomeInfo[[otherNames[1]]]$introns[mcols(genomeInfo[[otherNames[1]]]$introns)$name %in% unique(m2[m2$repGroup == repChoice & m2$conState == "dif", otherNames[1] ])])), 
      ), col = scols[otherNames[1]], lty = 2, lwd = 2)
      lines(density(log10(width( genomeInfo[[otherNames[2]]]$introns[mcols(genomeInfo[[otherNames[2]]]$introns)$name %in% unique(m2[m2$repGroup == repChoice & m2$conState == "dif", otherNames[2]]) ])), 
      ), col = scols[otherNames[2]], lty = 2, lwd = 2)
      legend("topright", legend = c(name1, otherNames), col = scols[c(name1, otherNames)], lwd = 2, bty = "n", lty = c(1,2,2))
      
      
      a <- genomeInfo[[otherNames[2]]]$introns[mcols(genomeInfo[[otherNames[2]]]$introns)$name %in% unique( m2[m2$repGroup == repChoice & m2$conState == "dif", otherNames[2]]   ) ]
      legend("topleft", legend = paste("n =", length(a)), bty = "n")
      title(main = "lineage specific")
    }else{
      plot.new()
    }
    title(main = paste(s,repChoice), outer = T,cex.main =2)
    mtext(text = "intron size (log10 bp)",side = 1,outer = TRUE, line = 1.25, cex = 1.5)
    mtext(text = "kernal density",side = 2,outer = TRUE, line = 1.25, cex = 1.5)
    
    
  }
  
}
dev.off()


# intron size is pretty conserved. 

# how do we compare the non-intersect


#### Non intersect







pdf(file = "plots/hotspotStats/intronSize/intronNonIntersect.pdf")

for(s in snames){
  name1 = s
  
  
  otherNames = snames[snames != name1]
  
  dat1 <- read.table(file = paste("data/repeatHotspot/",name1, "/",name1, "_", otherNames[1], "_conDif.txt", sep = ""), header = T)
  dat2 <- read.table(file = paste("data/repeatHotspot/",name1, "/",name1, "_", otherNames[2], "_conDif.txt", sep = ""), header = T)
  
  datRef1.gr <- GRanges(seqnames = Rle(dat1$chr[dat1$genome == "ref"]),
                        ranges = IRanges(start = dat1$start[dat1$genome == "ref"], end = dat1$end[dat1$genome == "ref"]))
  mcols(datRef1.gr) <- dat1[dat1$genome == "ref", 4:ncol(dat1)]
  seqlevels(datRef1.gr) <- seqlevels(genomeInfo[[name1]]$ensGene)
  seqinfo(datRef1.gr) <- seqinfo(genomeInfo[[name1]]$ensGene)
  
  olRef1 <- as.matrix(findOverlaps(query = genomeInfo[[name1]]$ensGene, subject = datRef1.gr))
  dfRef1 <- data.frame(ensGene = names(genomeInfo[[name1]]$ensGene[olRef1[,1]]), mcols(datRef1.gr[olRef1[,2]]) )
  
  datRef2.gr <- GRanges(seqnames = Rle(dat2$chr[dat2$genome == "ref"]),
                        ranges = IRanges(start = dat2$start[dat2$genome == "ref"], end = dat2$end[dat2$genome == "ref"]))
  mcols(datRef2.gr) <- dat2[dat2$genome == "ref", 4:ncol(dat2)]
  seqlevels(datRef2.gr) <- seqlevels(genomeInfo[[name1]]$ensGene)
  seqinfo(datRef2.gr) <- seqinfo(genomeInfo[[name1]]$ensGene)
  
  olRef2 <- as.matrix(findOverlaps(query = genomeInfo[[name1]]$ensGene, subject = datRef2.gr))
  dfRef2 <- data.frame(ensGene = names(genomeInfo[[name1]]$ensGene[olRef2[,1]]), mcols(datRef2.gr[olRef2[,2]]) )
  
  
  mRef1 <- merge(x = dfRef1, y = orthoDf, by.x = 1, by.y = (1:4)[colnames(orthoDf) == name1] )
  
  mRef2 <- merge(x = dfRef2, y = orthoDf, by.x = 1, by.y = (1:4)[colnames(orthoDf) == name1] )
  
  for(r in repGroups){
    repChoice = r
    mRef1.2 <- mRef1[ mRef1$repGroup == repChoice,]
    
    mRef2.2 <- mRef2[ mRef2$repGroup == repChoice,]
    
    
    layout(matrix(1:4,nrow = 4))
    par(mar = c(2,2,2,2), oma = c(5,5,5,5))
    
    #### other species 1
    if(length(genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef1.2[mRef1.2$repGroup == repChoice & mRef1.2$conState == "con", "ensGene"])]) > 2){
      plot(density(log10(width( genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef1.2[mRef1.2$repGroup == repChoice & mRef1.2$conState == "con", "ensGene"])])), 
      ), xlim = c(1,6), ylim= c(0,.8),lty = 1, lwd = 2, col = scols[name1],main = "")
      lines(density(log10(width( genomeInfo[[otherNames[1]]]$introns[mcols(genomeInfo[[otherNames[1]]]$introns)$name %in% unique(mRef1.2[mRef1.2$repGroup == repChoice & mRef1.2$conState == "con", otherNames[1] ] ) ])), 
      ), col = scols[otherNames[1]], lty = 2, lwd = 2)
      legend("topright", legend = c(name1, otherNames[1]), col = scols[c(name1, otherNames[1])], lwd = 2, bty = "n", lty = c(1,2))
      legend("topleft", bty = "n", legend = paste("n = ",  length(genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef1.2[mRef1.2$repGroup == repChoice & mRef1.2$conState == "con", "ensGene"]) ]) )  )
      title(main = "high in query")
    }else{
      plot.new()
    }
    
    if(length(genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef1.2[mRef1.2$repGroup == repChoice & mRef1.2$conState == "dif", "ensGene"])]) > 2){
      
      plot(density(log10(width( genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef1.2[mRef1.2$repGroup == repChoice & mRef1.2$conState == "dif", "ensGene"]) ])), 
      ),  xlim = c(1,6), ylim= c(0,.8),lty = 1, lwd = 2, col = scols[name1],main = "")
      lines(density(log10(width( genomeInfo[[otherNames[1]]]$introns[mcols(genomeInfo[[otherNames[1]]]$introns)$name %in% unique(mRef1.2[mRef1.2$repGroup == repChoice & mRef1.2$conState == "dif", otherNames[1] ])])), 
      ), col = scols[otherNames[1]], lty = 2, lwd = 2)
      legend("topright", legend = c(name1, otherNames[1]), col = scols[c(name1, otherNames[1])], lwd = 2, bty = "n", lty = c(1,2))
      legend("topleft", bty = "n", legend = paste("n = ",  length(genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef1.2[mRef1.2$repGroup == repChoice & mRef1.2$conState == "dif", "ensGene"]) ]) )  )
      title(main = "low in query")
    }else{
      plot.new()
    }
    
    
    #### other species 2
    
    if(length(genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef2.2[mRef2.2$repGroup == repChoice & mRef2.2$conState == "con", "ensGene"])]) > 2){
      
      plot(density(log10(width( genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef2.2[mRef2.2$repGroup == repChoice & mRef2.2$conState == "con", "ensGene"])])), 
      ), xlim = c(1,6), ylim= c(0,.8),lty = 1, lwd = 2, col = scols[name1],main = "")
      lines(density(log10(width( genomeInfo[[otherNames[2]]]$introns[mcols(genomeInfo[[otherNames[2]]]$introns)$name %in% unique(mRef2.2[mRef2.2$repGroup == repChoice & mRef2.2$conState == "con", otherNames[2] ] ) ])), 
      ), col = scols[otherNames[2]], lty = 2, lwd = 2)
      legend("topright", legend = c(name1, otherNames[2]), col = scols[c(name1, otherNames[2])], lwd = 2, bty = "n", lty = c(1,2))
      legend("topleft", bty = "n", legend = paste("n = ",  length(genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef2.2[mRef2.2$repGroup == repChoice & mRef2.2$conState == "con", "ensGene"]) ]) )  )
      title(main = "high in query")
    }else{plot.new()}
    
    if(length(genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef2.2[mRef2.2$repGroup == repChoice & mRef2.2$conState == "dif", "ensGene"])]) > 2){
      
      plot(density(log10(width( genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef2.2[mRef2.2$repGroup == repChoice & mRef2.2$conState == "dif", "ensGene"]) ])), 
      ),  xlim = c(1,6), ylim= c(0,.8),lty = 1, lwd = 2, col = scols[name1],main = "")
      lines(density(log10(width( genomeInfo[[otherNames[2]]]$introns[mcols(genomeInfo[[otherNames[2]]]$introns)$name %in% unique(mRef2.2[mRef2.2$repGroup == repChoice & mRef2.2$conState == "dif", otherNames[2] ])])), 
      ), col = scols[otherNames[2]], lty = 2, lwd = 2)
      legend("topright", legend = c(name1, otherNames[2]), col = scols[c(name1, otherNames[2])], lwd = 2, bty = "n", lty = c(1,2))
      legend("topleft", bty = "n", legend = paste("n = ",  length(genomeInfo[[name1]]$introns[mcols(genomeInfo[[name1]]$introns)$name %in% unique(mRef2.2[mRef2.2$repGroup == repChoice & mRef2.2$conState == "dif", "ensGene"]) ]) )  )
      title(main = "low in query")
    }else{plot.new()}
    
    title(main = paste(s,repChoice), outer = T,cex.main =2)
    mtext(text = "intron size (log10 bp)",side = 1,outer = TRUE, line = 1.25, cex = 1.25)
    mtext(text = "kernal density",side = 2,outer = TRUE, line = 1.25, cex = 1.25)
    
  }
}

dev.off()





# so we can get each conserved transcript
# average intro size of transcripts in each region








