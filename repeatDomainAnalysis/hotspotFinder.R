### A script to extract bins as hotspots
### we could give the same ids to adjacent bins that make up the same regions
### just take the min coord and max coord for stuff with the same id

## needs to be converted to executable type script



### create a naming convention. 

setwd(dir = "~/Desktop/RTN_domains/")
rm(list = ls())

library(spdep)
library(GenomicRanges)


sizes <- 50e3

rGroups <- c("ancient", "new_SINE", "new_L1", "old_L1")

genome = "mm9"

Gstar = TRUE

sizesGstat <- NULL
noX = TRUE

load(file = paste("~/Desktop/RTN_domains/R_objects/rmskMapTables/binSizes/", genome, "/repData_", genome, "_",as.integer(sizes),"_size.RData", sep = ""))




if(noX){
  chrX <- (1:nrow(repDataList$bin))[ repDataList$bin$chr == "chrX"]
  repDataList$bin <- repDataList$bin[-chrX,]
  repDataList$repSummary <- repDataList$repSummary[-chrX,]
  repDataList$neighborMat <- repDataList$neighborMat[[1]][-chrX, -chrX]
}

mat2 <- repDataList$neighborMat
diag(mat2) <- 0
lw <- mat2listw(mat2)
if(Gstar)
  lwG <- nb2listw(include.self(lw$neighbours)) 

binGstat <- repDataList$bin
for(j in 1:length(rGroups)){
  if(Gstar)
    Gstat <- localG(x = repDataList$repSummary[,rGroups[j]]/repDataList$bin$Known * 1e6, listw = lwG)
  else
    Gstat <- localG(x = repDataList$repSummary[,rGroups[j]]/repDataList$bin$Known * 1e6, listw = lw)
  binGstat[,rGroups[j]] <- Gstat
}



#names(sizesGstat) <- paste("binSize", as.integer(sizes), sep = "_")
domainList = NULL
pAdjMethod = "fdr"
alpha <- .01


statG <- binGstat
domainAll <- NULL
for(i in 1:4){
  padj <- p.adjust(2*pnorm(-abs(statG[,rGroups[i]])),method = pAdjMethod)
  #plot((padj < .05)[1:500], type = "l")
  # so how do we get regions
  statG2 <- cbind(statG[,1:5], padj)
  section <- statG2[statG2$padj < alpha,]
  section$adjacent = c(1,(section$start[2:nrow(section)] - section$start[1:(nrow(section)-1)])/sizes - 1)
  
  #regions = data.frame(chr = section$chr[section$adjacent != 0], start = section$start[section$adjacent != 0], 
  #end = section$end[c(section$adjacent[2:nrow(section)] != 0, TRUE)])
  
  #domainList <- c(domainList, list(regions))
  a = 0
  for(j in 1:nrow(section)){
    if(section$adjacent[j] != 0){
      a = a + 1
    }
    section$group[j] <- a
  }
  domains <- data.frame(chr = section$chr, 
                        start = section$start, 
                        end = section$end, 
                        domainID = paste(rGroups[i], 1:nrow(section), section$group, sep="_")
  )
  domainAll <- rbind(domainAll,domains)
}


write.table(domainAll, paste("data/repeatHotspot/", genome, "Hotspots.bed"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)



#names(domainList) = rGroups



chrChoice = "chr1"
xlim = c(1e6, 35e6)

xVals <- repDataList$bin$start[repDataList$bin$chr == chrChoice] + (sizes/2)

layout(matrix(c(1,2,3,4,5,6,7,8), nrow = 8),heights = rep(c(3,1), 4))
par( oma = c(5,1,5,1))

for(i in 1:4){
  repVals <- repDataList$repSummary[repDataList$bin$chr == chrChoice,rGroups[i]]/max(repDataList$repSummary[,rGroups[i]])
  colChoice = c("darkblue", "aquamarine3", "purple", "orangered")[i]
  
  par(mar = c(0,5,.5,5))
  plot(1,xlim = xlim, ylim = c(0,1),type = "n", axes = FALSE, xlab = "", ylab = "")
  grid()
  
  lines(xVals,repVals, lwd =2, col =  colChoice)
  
  par(mar = c(.5,5,0,5))
  plot(1,xlim = xlim, ylim = c(0,1),type = "n", axes = FALSE, xlab = "", ylab = "")
  rect(xleft = domainList[[rGroups[i]]]$start[domainList[[rGroups[i]]]$chr == chrChoice], 
       xright = domainList[[rGroups[i]]]$end[domainList[[rGroups[i]]]$chr== chrChoice], 
       ytop = 1, ybottom = 0, col = colChoice, density = 20)
}
axis(side = 1, at = seq(xlim[1], xlim[2],length.out = 20))

title(main = sizes,outer = TRUE)


### maybe it is time to pull contact data


save(domainList,file = paste("R_objects/temp/",genome,"_domainList.RData", sep = ""))


nrow(domainList$ancient[domainList$ancient$chr == "chr1",])

