

# we could select a few parameters and just run this as an Rscript across our species alignment

rm(list = ls())
setwd("Desktop/RTN_domains/")
refGenome = "hg19"
queGenome = "mm10"

datAll = NULL
for(i in c(10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000, 300000)){
  dat <- read.table(paste("data/chainAlignments/chainOut/",refGenome,"/",queGenome,"/merged.",
                          refGenome,".",queGenome,".",as.integer(i),".txt", sep = ""),
                    col.names = c("refChr",
                                  "refLen",
                                  "refStart",
                                  "refEnd",
                                  "refStrand",
                                  "refGap",
                                  "queChr",
                                  "queLen",
                                  "queStart",
                                  "quesEnd",
                                  "queStrand",
                                  "queGap",
                                  "chainID")
  )
  dat$group = i
  datAll <- rbind(datAll, dat)
}

datAll$group <- as.factor(datAll$group)


# maybe look at which minimum gap length gives us the most apropriate number of alignemnts

# the most number of bases above a certain length 
minRatio = .85
maxRatio = 1.15
blockLengths <- c(500000, 1e6,2*10^6,3e6,5e6)

pdf(file = paste("plots/alignmentChainStats/", refGenome, ".", queGenome,"_minRatio",minRatio,"_maxRatio",maxRatio ,".pdf", sep = ""), height = 8,width=5)

layout(mat= matrix(c(1,2,3,4), nrow = 4))
par(mar = c(0,0,0,0), oma = c(6,10,8,6))
for(i in 1:length(blockLengths)){
  datRefine <- datAll[datAll$refEnd - datAll$refStart > blockLengths[i],]
  datRefine <- datRefine[(datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) < maxRatio & 
                           (datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) > minRatio,]
  
  aggBP <- aggregate(as.numeric(datRefine$refEnd - datRefine$refStart), list(as.numeric(as.character(datRefine$group))), FUN = sum)
  
  if(i == 1){
    plot((aggBP[,1]),aggBP[,2], type = "b", xaxt = "n", ylim = c(0,max(aggBP[,2],na.rm = TRUE)), las = 2)
    mtext("aligned bp",side = 2,outer = TRUE,line = 4.5,at = .875)
    title(main = paste(refGenome,"chained alignment to",queGenome, "\nbetween",minRatio,"and",maxRatio,"ref length"),outer = TRUE, cex.main = 2)
  } else {
    lines((aggBP[,1]),aggBP[,2], type = "b", xaxt = "n", col = i)
  }
  
}
legend("bottomright",legend = paste(">" ,blockLengths/1e6, "mb"), col = 1:length(blockLengths), lwd = 1, title = "alignment blocks", bty = "n")


for(i in 1:length(blockLengths)){
  datRefine <- datAll[datAll$refEnd - datAll$refStart > blockLengths[i],]
  datRefine <- datRefine[(datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) < maxRatio & 
                           (datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) > minRatio,]
  
  aggGap <- aggregate(datRefine$refGap, list(as.numeric(as.character(datRefine$group))), FUN = sum)
  
  if(i == 1){
    plot((aggGap[,1]),aggGap[,2], type = "b", xaxt = "n", ylim = c(0,max(aggGap[,2])), las = 2)
    mtext("gap bp",side = 2,outer = TRUE,line = 4.5,at = .625)
    
  } else {
    lines((aggGap[,1]),aggGap[,2], type = "b", col = i )
  }
}

for(i in 1:length(blockLengths)){
  datRefine <- datAll[datAll$refEnd - datAll$refStart > blockLengths[i],]
  datRefine <- datRefine[(datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) < maxRatio & 
                           (datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) > minRatio,]
  
  aggBP <- aggregate(as.numeric(datRefine$refEnd - datRefine$refStart), list(as.numeric(as.character(datRefine$group))), FUN = sum)
  aggGap <- aggregate(datRefine$refGap, list(as.numeric(as.character(datRefine$group))), FUN = sum)
  
  if(i == 1){
    plot((aggGap[,1]),aggGap[,2]/aggBP[,2], type = "b", xaxt = "n", ylim = c(0,1), las = 2)
    mtext("gap rate",side = 2,outer = TRUE,line = 4.5,at = .375)
    
  } else {
    lines((aggGap[,1]),aggGap[,2]/aggBP[,2], type = "b", col = i )
  }
}
  
for(i in 1:length(blockLengths)){
  datRefine <- datAll[datAll$refEnd - datAll$refStart > blockLengths[i],]
  datRefine <- datRefine[(datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) < maxRatio & 
                           (datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) > minRatio,]
  
  aggNo <- aggregate(as.numeric(datRefine$refEnd - datRefine$refStart), list(as.numeric(as.character(datRefine$group))), FUN = length) 
  if(i == 1){
    plot((aggNo[,1]),aggNo[,2], type = "b", ylim = c(0,max(aggNo[,2])), las = 1)
    mtext("aligned blocks",side = 2,outer = TRUE,line = 4.5,at = .125)
    mtext("minimum gap size (bp)", side = 1, line = 3)
    
  } else {
    lines((aggNo[,1]),aggNo[,2], type = "b", col = i )
  }
}

dev.off()

# we could also look at the gap rate distribution
# only pick stuff with an appropriate gap rate



# what about the inheritance structure?
# when we look at differnt levels of alignments. 

