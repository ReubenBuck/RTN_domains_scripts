### A script to extract bins as hotspots
### we could give the same ids to adjacent bins that make up the same regions
### just take the min coord and max coord for stuff with the same id

## needs to be converted to executable type script


#turn this ito executable


### create a naming convention. 

# inputs
# 




setwd(dir = "~/Desktop/RTN_domains/")
rm(list = ls())

library(spdep)
library(GenomicRanges)
library(Matrix)


sizes <- 50e3

rGroups <- c("ancient", "new_SINE", "new_L1", "old_L1")

genome = "mm10"

Gstar = TRUE

sizesGstat <- NULL
noX = FALSE

load(file = paste("~/Desktop/RTN_domains/R_objects/rmskMapTables/binSizes/", genome, "/repData_", genome, "_",as.integer(sizes),"_size.RData", sep = ""))




if(noX){
  chrX <- (1:nrow(repDataList$bin))[ repDataList$bin$chr == "chrX"]
  repDataList$bin <- repDataList$bin[-chrX,]
  repDataList$repSummary <- repDataList$repSummary[-chrX,]
  repDataList$neighborMat <- repDataList$neighborMat[[1]][-chrX, -chrX]
} else {
  repDataList$neighborMat <- repDataList$neighborMat[[1]]
}

mat2 <- repDataList$neighborMat
Matrix::diag(mat2) <- 0
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


write.table(domainAll, paste("data/repeatHotspot/", genome, "Hotspots.bed", sep = ""),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)



save(domainList,file = paste("R_objects/domainList/",genome,"_domainList.RData", sep = ""))

