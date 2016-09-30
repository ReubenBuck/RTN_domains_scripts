

# now we can begin to analyse the structure of repeats 

# simply read data in pull out hotspots and calculate clustering scores 

# get moran's I
rm(list = ls())

library(spdep)
library(GenomicRanges)


sizes <- c(50000, 100000,250000,500000,750000,1000000,1500000,2000000)

rGroups <- c("ancient", "new_SINE", "new_L1", "old_L1")

Gstar = FALSE

sizesGstat <- NULL
mI <- NULL
noX = TRUE
for(i in 1:length(sizes)){
  
  load(file = paste("Desktop/RTN_domain/R_objects/repMapData/repData_hg19_",as.integer(sizes[i]),".RData", sep = ""))
  
  if(noX){
    chrX <- (1:nrow(repDataList$bin))[ repDataList$bin$chr == "chrX"]
    repDataList$bin <- repDataList$bin[-chrX,]
    repDataList$repSummary <- repDataList$repSummary[-chrX,]
    repDataList$neighborMat <- repDataList$neighborMat[-chrX, -chrX]
  }
  
  mat2 <- repDataList$neighborMat
  diag(mat2) <- 0
  lw <- mat2listw(mat2)
  if(Gstar)
    lwG <- nb2listw(include.self(lw$neighbours)) 
  
  mTmc <- NULL
  binGstat <- repDataList$bin
  for(j in 1:length(rGroups)){
    mTmc <- c(mTmc, moran.mc(x = repDataList$repSummary[,rGroups[j]]/repDataList$bin$Known * 1e6,listw = lw,nsim = 100)$statistic)
    if(Gstar)
      Gstat <- localG(x = repDataList$repSummary[,rGroups[j]]/repDataList$bin$Known * 1e6, listw = lwG)
    else
      Gstat <- localG(x = repDataList$repSummary[,rGroups[j]]/repDataList$bin$Known * 1e6, listw = lw)
    binGstat[,rGroups[j]] <- Gstat
  }
  
  mI <- cbind(mI, mTmc)
  sizesGstat <- c(sizesGstat, list(binGstat))
  
  rm(repDataList)
}
rownames(mI) <- rGroups
colnames(mI) <- paste("binSize", as.integer(sizes), sep = "_")
names(sizesGstat) <- paste("binSize", as.integer(sizes), sep = "_")


# 
# layout(matrix(c(1,2), nrow = 2))
# chrSelect <- "chr20"
# a = 1
# b = 5
# c = 8
# plot(sizesGstat[[a]]$start[sizesGstat[[a]]$chr == chrSelect] + (sizes[a]/2), 
#      sizesGstat[[a]]$ancient[sizesGstat[[a]]$chr == chrSelect], 
#      type = "l", xlim = c(10000000, 60000000))
# lines(sizesGstat[[b]]$start[sizesGstat[[b]]$chr == chrSelect] + (sizes[b]/2), sizesGstat[[b]]$ancient[sizesGstat[[b]]$chr == chrSelect], col = 2)
# lines(sizesGstat[[c]]$start[sizesGstat[[c]]$chr == chrSelect] + (sizes[c]/2), sizesGstat[[c]]$ancient[sizesGstat[[c]]$chr == chrSelect], col = 3)
# 
# abline(h = 3)
# 
# plot(sizesGstat[[a]]$start[sizesGstat[[a]]$chr == chrSelect] + (sizes[a]/2), 
#      sizesGstat[[a]]$new_SINE[sizesGstat[[a]]$chr == chrSelect], 
#      type = "l", xlim = c(10000000, 60000000))
# lines(sizesGstat[[b]]$start[sizesGstat[[b]]$chr == chrSelect] + (sizes[b]/2), sizesGstat[[b]]$new_SINE[sizesGstat[[b]]$chr == chrSelect], col = 2)
# lines(sizesGstat[[c]]$start[sizesGstat[[c]]$chr == chrSelect] + (sizes[c]/2), sizesGstat[[c]]$new_SINE[sizesGstat[[c]]$chr == chrSelect], col = 3)
# 
# abline(h = 3)

  
chrChoice = "chr1"
pAdjMethod <- "fdr"
alpha <- .01
xlim = c(1, max(sizesGstat$binSize_50000$end[sizesGstat$binSize_50000$chr == chrChoice]))

fileName <- paste(c("Gi", "GStarI")[Gstar + 1],"_",pAdjMethod,"_" ,alpha,"_", chrChoice,":",as.integer(xlim[1]),"-",as.integer(xlim[2]),".pdf", sep = "")

pdf(paste("Desktop/RTN_domain/plots/binSizeExperiment/",fileName, sep = ""), width = 10,height = 5)
layout(matrix(c(1,2,3,4), nrow = 4))
par(mar = c(.5,5,.5,5), oma = c(5,1,5,1))

for(i in 1:4){
  repChoice = rGroups[i]
  colChoice = c("darkblue", "darkgreen", "purple", "red")[i]
    
  plot(1,xlim = xlim, ylim = c(0,1),type = "n", axes = FALSE, xlab = "", ylab = "")
  axis(side = 2,at = (1:length(sizes))/length(sizes) - (.5/length(sizes)), labels = sizes, las = 2)
  grid()
  
  for(j in 1:length(sizes)){
    
    sizeChoice = sizes[j]
    
    statG <- sizesGstat[[paste("binSize", as.integer(sizeChoice), sep ="_")]]
    padj <- p.adjust(2*pnorm(-abs(statG[,repChoice])),method = pAdjMethod)
    
    if(length(statG$start[padj < alpha & statG[,repChoice] > 0 & statG$chr == chrChoice]) > 0)
      rect(xleft = statG$start[padj < alpha & statG[,repChoice] > 0 & statG$chr == chrChoice], 
         xright = statG$start[padj < alpha & statG[,repChoice] > 0 & statG$chr == chrChoice] + sizeChoice, 
         ytop = (1/length(sizes)) * j, ybottom = ((1/length(sizes)) * j) - (1/length(sizes)), 
         density = -1, border = NA, col = colChoice)
  }
}
axis(seq(xlim[1], xlim[2],length.out = 20))
mtext(text = chrChoice,side = 1,outer = TRUE, line = 2.5)
title(main = paste(c("Gi,", "G*i,")[Gstar + 1], pAdjMethod, "adjusted P-value <", alpha), outer = TRUE)

dev.off()

# so it would good to make some comparisons 

# so fdr and bonferoni

# G and G*

pdf(file = "Desktop/RTN_domain/plots/binSizeExperiment/spatialClustering.pdf", width = 6, height = 5)
layout(1)
par(mar = c(5,5,5,5),oma = c(0,0,0,0))
plot(sizes/1000,mI["ancient",], ylim = c(0,1), type = "b", 
     col = "darkblue", lwd = 3,
     main = "Spatial clustering of retrotransposons",
     xlab = "bin size (kb)", ylab = "Moran's I")
lines(sizes/1000,mI["new_SINE",], col = "darkgreen", lwd = 3, type = "b")
lines(sizes/1000,mI["new_L1",], col = "purple", lwd = 3, type = "b")
lines(sizes/1000,mI["old_L1",], col = "red", lwd = 3, type = "b")
legend("bottomright", legend = rGroups, fill = c("darkblue", "darkgreen", "purple", "red"), bty = "n") 
dev.off()


