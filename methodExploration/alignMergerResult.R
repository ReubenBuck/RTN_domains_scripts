

# we could select a few parameters and just run this as an Rscript across our species alignment


setwd("Desktop/RTN_domains/")

datAll = NULL
for(i in c(10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000, 300000)){
  dat <- read.table(paste("data/chainAlignments/test/chainOut/mergedChain",as.integer(i),".txt", sep = ""),
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

layout(mat= matrix(c(1,2,3,4), nrow = 4))
par(mar = c(0,0,0,0), oma = c(6,10,8,6))
for(i in c(1000000, 2*10^6, 3*10^6, 4*10^6, 5*10^6)){
  datRefine <- datAll[datAll$refEnd - datAll$refStart > i,]
  datRefine <- datRefine[(datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) < 1.1 & (datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) > .9,]
  
  aggBP <- aggregate(datRefine$refEnd - datRefine$refStart, list(as.numeric(as.character(datRefine$group))), FUN = sum)
  
  if(i == 1000000){
    plot((aggBP[,1]),aggBP[,2], type = "b", xaxt = "n", ylim = c(0,11e8), las = 2)
    mtext("aligned bp",side = 2,outer = TRUE,line = 4.5,at = .875)
    title(main = "hg19 chained alignment to mm9",outer = TRUE, cex.main = 2)
  } else {
    lines((aggBP[,1]),aggBP[,2], type = "b", xaxt = "n", col = i / 1000000)
  }
  
}
legend("bottomright",legend = paste(">" ,1:5, "mb"), col = 1:5, lwd = 1, title = "alignment blocks", bty = "n")


for(i in c(1000000, 2*10^6, 3*10^6, 4*10^6, 5*10^6)){
  datRefine <- datAll[datAll$refEnd - datAll$refStart > i,]
  datRefine <- datRefine[(datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) < 1.1 & (datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) > .9,]
  
  aggGap <- aggregate(datRefine$refGap, list(as.numeric(as.character(datRefine$group))), FUN = sum)
  
  if(i == 1000000){
    plot((aggGap[,1]),aggGap[,2], type = "b", xaxt = "n", ylim = c(0,7e8), las = 2)
    mtext("gap bp",side = 2,outer = TRUE,line = 4.5,at = .625)
    
  } else {
    lines((aggGap[,1]),aggGap[,2], type = "b", col = i / 1000000)
  }
}

for(i in c(1000000, 2*10^6, 3*10^6, 4*10^6, 5*10^6)){
  datRefine <- datAll[datAll$refEnd - datAll$refStart > i,]
  datRefine <- datRefine[(datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) < 1.1 & (datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) > .9,]
  
  aggBP <- aggregate(datRefine$refEnd - datRefine$refStart, list(as.numeric(as.character(datRefine$group))), FUN = sum)
  aggGap <- aggregate(datRefine$refGap, list(as.numeric(as.character(datRefine$group))), FUN = sum)
  
  if(i == 1000000){
    plot((aggGap[,1]),aggGap[,2]/aggBP[,2], type = "b", xaxt = "n", ylim = c(.3,.7), las = 2)
    mtext("gap rate",side = 2,outer = TRUE,line = 4.5,at = .375)
    
  } else {
    lines((aggGap[,1]),aggGap[,2]/aggBP[,2], type = "b", col = i / 1000000)
  }
}
  
for(i in c(1000000, 2*10^6, 3*10^6, 4*10^6, 5*10^6)){
  datRefine <- datAll[datAll$refEnd - datAll$refStart > i,]
  datRefine <- datRefine[(datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) < 1.1 & (datRefine$refEnd-datRefine$refStart)/(datRefine$quesEnd - datRefine$queStart) > .9,]
  
  aggNo <- aggregate(datRefine$refEnd - datRefine$refStart, list(as.numeric(as.character(datRefine$group))), FUN = length) 
  if(i == 1000000){
    plot((aggNo[,1]),aggNo[,2], type = "b", ylim = c(0,300), las = 1)
    mtext("aligned blocks",side = 2,outer = TRUE,line = 4.5,at = .125)
    mtext("minimum gap size (bp)", side = 1, line = 3)
    
  } else {
    lines((aggNo[,1]),aggNo[,2], type = "b", col = i / 1000000)
  }
}



# we could also look at the gap rate distribution
# only pick stuff with an appropriate gap rate

