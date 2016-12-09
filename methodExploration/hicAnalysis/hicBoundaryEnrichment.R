### hic processing

rm(list = ls())
setwd("~/Desktop/RTN_domains/")

load(file = "R_objects/temp/rheMac3_domainList.RData")
load(file = "R_objects/rmskMapTables/binSizes/rheMac3/repData_rheMac3_50000_size.RData")

hicSum <- read.table(file = "data/hi-cSummary/GSE65126_HiC_macaque_liver_merged_50000.txt", header = TRUE)



hic <- hicSum


# I want to pull a square matrix out and look at the pattern


# if I could give a bunch of positions and identify the relevent contacts

repChoice = "new_SINE"


repDat <- data.frame(chr = as.character(gsub(pattern = "chr", replacement = "",x = repDataList$bin$chr)), 
                     start = repDataList$bin$start, 
                     rep =repDataList$repSummary[[repChoice]])


domains <- domainList[[repChoice]]

hist((domains$end - domains$start + 1), breaks= 100)

domains <- domains[(domains$end - domains$start + 1) >= 400000,]

positionsBoundary = data.frame(chr = gsub(pattern = "chr", 
                                          replacement = "", 
                                          x = c(as.character(domains$chr))) ,
                               position = c(domains$start - 1))
domainLength <- ((domains$end - domains$start + 1)/50000)
domainLength <- c(domainLength)

positionsRand = sample(x = 1:length(unique(paste(hic$chrom1, hic$start1))), size = nrow(positionsBoundary))
positionsRand = data.frame(chr = hic$chrom1[ positionsRand ], position = hic$start1[ positionsRand ])

# pull out regions where both contacts are completely within our range 

colourRamp <- colorRampPalette(c("white", "red"))
cols <- colourRamp(40)

for(pType in 1:2){
  
  positions = get(c("positionsBoundary","positionsRand")[pType])
  
  summedRegion <- matrix(0,nrow = 40, ncol = 40)
  summedRep <- NULL
  for(p in 1:nrow(positions)){
    pos = positions[p,]
    region = hic[hic$chrom1 == as.character(pos$chr), ]
    region = region[region$start1 >= pos$position - 1e6 & region$end1 <= pos$position + 1e6,]
    region = region[region$start2 >= pos$position - 1e6 & region$end2 <= pos$position + 1e6,]
    
    
    repRegion = repDat[repDat$chr == as.character(pos$chr),]
    repRegion = repRegion[repRegion$start >= pos$position - 1e6 & repRegion$start <= pos$position + 1e6,]
    repRegion$start = repRegion$start - (pos$position - 1e6)
    repRegion$start = (repRegion$start - 1)/50000 + 1
    
    region$start1 = region$start1 - (pos$position - 1e6)
    region$start2 = region$start2 - (pos$position - 1e6)
    region$end1 = region$end1 - (pos$position - 1e6)
    region$end2 = region$end2 - (pos$position - 1e6)
    
    regionIDs <- data.frame(ID1 = region$start1/50000 + 1,
                            ID2 = region$start2/50000 + 1,
                            value = region$observed_count/region$expected_count
    )
    
    regionMat <- matrix(0, nrow = 40, ncol =  40)
    for(i in 1:nrow(regionIDs)){
      regionMat[regionIDs$ID1[i], regionIDs$ID2[i]] <- regionIDs$value[i]
    }
    #image(sqrt(regionMat), col = cols)
#     if(p>(length(positions)/2) + 1){
#       regionMat[1:40,1:40] <- regionMat[40:1,40:1]
#     }
#      if(pType == 1 & p < (.5*length(domainLength)) + 1){
#       layout(c(1,2,3), heights = c(3,1,1))
#       par(mar = c(0,0,5,0), oma = c(5,5,5,5))
#       image((regionMat), col = cols, main = p,axes = FALSE)
#       par(mar=c(0,0,0,0))
#       repMat <- matrix( c(rep(0,20),rep(1,domainLength[p]),rep(0,20))[1:40],ncol = 1)
#       image(log10(repMat), col = cols, axes = FALSE)
#       plot(x = repRegion$start,y = repRegion$rep/max(repDat$rep), type = "l", ylim = c(0,1))
#       grid()
#     }
    summedRegion <- summedRegion + (regionMat)
    repPreSum <- rep(0,40)
    repPreSum[repRegion$start] <- repRegion$rep
    summedRep <- rbind(summedRep, repPreSum)
  }
  assign(x = c("boundaryRegion", "randomRegion")[pType],value = summedRegion)
  assign(x = c("boundaryRep", "randomRep")[pType],value = summedRep)
}


# the idea is that if our regions aren't actually hotspots, then of course we wont be able to see the association

colourRamp <- colorRampPalette(c("white", "red"))
cols <- colourRamp(40)


layout(matrix(c(1,2,3,4), nrow = 4), heights = c(3,1,3,1))
par(mar = c(0,0,0,0))
image((boundaryRegion), col = cols, axes = FALSE)
abline(v = .5)
plot(colMeans(boundaryRep), type = "l", ylim = c(0,80))
abline(v = 20.5)
lines(colMeans(boundaryRep)+(2*apply(X = boundaryRep,2,sd)), col = 2)
lines(colMeans(boundaryRep)-(2*apply(X = boundaryRep,2,sd)), col = 2)
image((randomRegion), col = cols, axes = FALSE)
abline(v = .5)
plot(colMeans(randomRep), type = "l", ylim = c(0,80))
abline(v = 20.5)
lines(colMeans(randomRep)+(2*apply(X = randomRep,2,sd)), col = 2)
lines(colMeans(randomRep)-(2*apply(X = randomRep,2,sd)), col = 2)



colourRamp <- colorRampPalette(c("blue","white", "red"))
cols <- colourRamp(40)
OtakeE <- (boundaryRegion)-(randomRegion)
#OtakeE <-OtakeE - 1
zMax <- max(abs(OtakeE))
image(OtakeE, col = cols, zlim = c(-zMax, zMax))
abline(h = .5)
abline(v = .5)
