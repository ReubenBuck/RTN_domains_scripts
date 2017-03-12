rm(list = ls())

library(spdep)
library(GenomicRanges)


sizes <- c(20000, 50000, 100000,250000,500000,750000,1000000,1500000,2000000)

rGroups <- c("ancient", "new_SINE", "new_L1", "old_L1")


genome = "mm10"

mI <- NULL
noX = FALSE
for(i in 1:length(sizes)){
  
  load(file = paste("~/Desktop/RTN_domains/R_objects/rmskMapTables/binSizes/", genome, "/repData_", genome, "_",as.integer(sizes[i]),"_size",".RData", sep = ""))
  
  if(noX){
    chrX <- (1:nrow(repDataList$bin))[ repDataList$bin$chr == "chrX"]
    repDataList$bin <- repDataList$bin[-chrX,]
    repDataList$repSummary <- repDataList$repSummary[-chrX,]
    repDataList$neighborMat <- repDataList$neighborMat[[1]][-chrX, -chrX]
  }else{
    repDataList$neighborMat <- repDataList$neighborMat[[1]]
  }
  
  mat2 <- repDataList$neighborMat
  diag(mat2) <- 0
  lw <- mat2listw(mat2)
 
  mTmc <- NULL
  for(j in 1:length(rGroups)){
    mTmc <- c(mTmc, moran.mc(x = repDataList$repSummary[,rGroups[j]]/repDataList$bin$Known * 1e6,listw = lw,nsim = 1000)$statistic)
  }
  
  mI <- cbind(mI, mTmc)
  
  rm(repDataList)
}
rownames(mI) <- rGroups
colnames(mI) <- paste("binSize", as.integer(sizes), sep = "_")

pdf(file = paste("~/Desktop/RTN_domains/plots/binSizeExperiment/",genome,"/spatialClustering_",genome,"_NoX_",noX,".pdf", sep = ""), width = 6, height = 5)
layout(1)
par(mar = c(5,5,5,5),oma = c(0,0,0,0))
plot(sizes/1000,mI["ancient",], ylim = c(0,1), type = "b", 
     col = "darkblue", lwd = 3,
     main = paste("Spatial clustering of retrotransposons in", genome),
     xlab = "bin size (kb)", ylab = "Moran's I")
lines(sizes/1000,mI["new_SINE",], col = "aquamarine3", lwd = 3, type = "b")
lines(sizes/1000,mI["new_L1",], col = "purple", lwd = 3, type = "b")
lines(sizes/1000,mI["old_L1",], col = "orangered", lwd = 3, type = "b")
legend("bottomright", legend = rGroups, fill = c("darkblue", "aquamarine3", "purple", "orangered"), bty = "n") 
dev.off()

