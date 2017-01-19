
library(devtools)
setwd("~/Desktop/RTN_domains/")

rm(list = ls())

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")


g.X = "canFam3"

g.Y = "mm9"

genomes <- c(g.X,g.Y)

load(paste("R_objects/comparativeHiC/",g.X,"_",g.Y,".RData",sep = ""))

load(paste("R_objects/temp/",g.X,"_domainList.RData", sep = "")); assign(paste("g.X.domainList", sep = "_"), domainList)
load(paste("R_objects/temp/",g.Y,"_domainList.RData", sep = "")); assign(paste("g.Y.domainList", sep = "_"), domainList)
domainCols <- c("darkblue", "aquamarine3", "purple", "red")


squareSize = 50



# now we have the repeat information loaded in 
# this way we can look at comparative regions and their repeat enrichment
# might also be worth adding genes that aren't pairwise orhtos. 



head(combinedComplete)






# so it can be done, I've just pulled out the wrong triangle




# what kind of permutation cna we do to work out significant interaction differences



chr.X = "chr3"
chrChoiceX <- combinedComplete[combinedComplete$chr1.x == chr.X,]
chrChoiceX <- chrChoiceX[order(chrChoiceX$start1.x),]

chrCols <- as.factor(as.character(chrChoiceX$chr2.y))

layout(1);par(mar=c(0,0,0,0), oma = c(5,5,5,5))
plot(chrChoiceX$start1.x, chrChoiceX$start1.y, 
     col = rainbow(length(levels(chrCols)))[chrCols])
legend("topright", legend = levels(chrCols), fill = rainbow(length(levels(chrCols))), bty = "n")

abline(v=75e6)


xlimsX <- c(71.5e6,75e6)


genoMatrixX <- read.delim(paste("data/trialHicAnalysis/",g.X,"trial/",chr.X,"raw20kb.txt", sep = ""), 
                         header = TRUE, sep = "\t",row.names = 1)
genoMatrixX <- as.matrix(genoMatrixX[,2:ncol(genoMatrixX)])





xlimsY <- chrChoiceX[chrChoiceX$start1.x > xlimsX[1] & chrChoiceX$end2.x < xlimsX[2],]
xlimsY <- xlimsY[xlimsY$chr1.y == xlimsY$chr1.y[1],]
xlimsY <- range(c(xlimsY$start1.y, xlimsY$start2.y, xlimsY$end1.y, xlimsY$end2.y))
xlimsY[1] <- round((xlimsY[1]-3000000)/20e3) * 20e3
xlimsY[2] <- round((xlimsY[2]+3000000)/20e3) * 20e3

chrChoiceY <- chrChoiceX[chrChoiceX$start1.x > xlimsX[1]& chrChoiceX$end2.x < xlimsX[2],]
chrChoiceY <- chrChoiceY[chrChoiceY$chr1.y == chrChoiceY$chr1.y[1],]
if(chrChoiceY$start1.y[1] > chrChoiceY$start1.y[nrow(chrChoiceY)]){
  xlimsY <- rev(xlimsY)
}

xlimsX <- range(c(chrChoiceY$start1.x, chrChoiceY$start2.x, chrChoiceY$end1.x, chrChoiceY$end2.x))
xlimsX[1] <- round((xlimsX[1]-3000000)/20e3) * 20e3
xlimsX[2] <- round((xlimsX[2]+3000000)/20e3) * 20e3




chr.Y <- unique(chrChoiceY$chr1.y)
genoMatrixY <- read.delim(paste("data/trialHicAnalysis/",g.Y,"trial/",chr.Y,"raw20kb.txt", sep = ""), 
                          header = TRUE, sep = "\t",row.names = 1)
genoMatrixY <- as.matrix(genoMatrixY[,2:ncol(genoMatrixY)])








###### Fig Plot


ylims <- c(0, 0.01)
layout(c(1,2,3,4,5,6), height = c(2,.5,3,2,.5,3),width=c(2,2))
par(mar=c(0,0,0,0), oma=c(5,5,5,5))



plot(chrChoiceY$start1.x, y = (chrChoiceY$interaction.x/sum(chrChoiceY$interaction.x)), 
     xlim = xlimsX, type = "n", xaxt = "n", ylim = ylims, xaxs="i")
legend("topleft", legend = g.X,bty = "n", cex = 1.5)

for( i in 1:nrow(chrChoiceY)){
  polygon(x = c(chrChoiceY$start1.x[i], 
                (chrChoiceY$start1.x + ((chrChoiceY$end2.x - chrChoiceY$start1.x)/2))[i] ,
                chrChoiceY$end2.x[i]), 
          y = (c(0, 
                 (chrChoiceY$interaction.x[i]/sum(chrChoiceX$interaction.x, na.rm = T)),
                 0 )))
}

cols <- topo.colors(nrow(chrChoiceY))[1:nrow(chrChoiceY)] 
points((chrChoiceY$start1.x + ((chrChoiceY$end2.x - chrChoiceY$start1.x)/2)), 
       y = (chrChoiceY$interaction.x/sum(chrChoiceX$interaction.x)), 
       pch = 16, col = cols )

rect(xleft = chrChoiceY$start1.x, xright = chrChoiceY$end1.x, 
     ybottom = 0, ytop = .01* ylims[2], density = -1, col = 1)
rect(xleft = chrChoiceY$start2.x, xright = chrChoiceY$end2.x, 
     ybottom = 0, ytop = .01* ylims[2], density = -1, col = 1)

mat <- log10(genoMatrixX[ (xlimsX[1]/20e3):(xlimsX[2]/20e3),(xlimsX[1]/20e3):(xlimsX[2]/20e3) ])

par(new=TRUE)
plot(c(rep(NA,squareSize),insulationScore(mat,squareSize = squareSize),rep(NA,squareSize)), type = "l", col = 2, lwd = 2 , xaxs = "i", axes=F)


### plot RTN
plot( chrChoiceY$start1.x + ((chrChoiceY$start2.x - chrChoiceY$start1.x)/2), 
     y = log10( (chrChoiceY$interaction.x/sum(chrChoiceX$interaction.x)) / (chrChoiceY$interaction.y/sum(chrChoiceX$interaction.x)) ), 
     xlim = xlimsX, ylim = c(0,4),type = "n", xaxt = "n", xaxs="i")

for(i in 1:length(g.X.domainList)){
  domain <- g.X.domainList[[i]]
  domain <- domain[domain$chr == as.character(chr.X),]
  if(nrow(domain) > 0){
    rect(xleft = domain$start, ybottom = 4-i, ytop = 5-i,xright = domain$end, col = domainCols[i], density = -1)
  }
}





par(mar=c(7,0,0,0))
whiteRed <- colorRampPalette(c("white","red"))
image(triangulate(mat)[,(nrow(mat)/2):1], col = whiteRed(40), axes = F, ylim = c(.5,1), xlab = chr.X)

axis(side = 1,at = seq(0,1,length.out = 5),labels = seq(from = xlimsX[1], to = xlimsX[2], length.out = 5))



par(mar=c(0,0,0,0))


plot(chrChoiceY$start1.y, y = (chrChoiceY$interaction.y/sum(chrChoiceX$interaction.y)), 
     xlim = xlimsY, type = "n", xaxt = "n", ylim = ylims, xaxs="i")
legend("topleft", legend = g.Y,bty = "n", cex = 1.5)

if(xlimsY[1] < xlimsY[2]){
  for( i in 1:nrow(chrChoiceY)){
    polygon(x = c(chrChoiceY$start1.y[i], 
                  (chrChoiceY$start1.y + ((chrChoiceY$end2.y - chrChoiceY$start1.y)/2))[i] ,
                  chrChoiceY$end2.y[i]), 
            y = (c(0, 
                   (chrChoiceY$interaction.y[i]/sum(chrChoiceX$interaction.y, na.rm = T)),
                   0 )))
  }
  cols <- topo.colors(nrow(chrChoiceY))[1:nrow(chrChoiceY)] 
  points((chrChoiceY$start1.y + ((chrChoiceY$end2.y - chrChoiceY$start1.y)/2)), 
         y = (chrChoiceY$interaction.y/sum(chrChoiceX$interaction.y)), 
         pch = 16, col = cols )
}else{
  for( i in 1:nrow(chrChoiceY)){
    polygon(x = c(chrChoiceY$end1.y[i], 
                  (chrChoiceY$start2.y + ((chrChoiceY$end1.y - chrChoiceY$start2.y)/2))[i] ,
                  chrChoiceY$start2.y[i]), 
            y = (c(0, 
                   (chrChoiceY$interaction.y[i]/sum(chrChoiceX$interaction.y, na.rm = T)),
                   0 )))
  }
  cols <- topo.colors(nrow(chrChoiceY))[1:nrow(chrChoiceY)] 
  points((chrChoiceY$start2.y + ((chrChoiceY$end1.y - chrChoiceY$start2.y)/2)), 
         y = (chrChoiceY$interaction.y/sum(chrChoiceX$interaction.y)), 
         pch = 16, col = cols )
}




rect(xleft = chrChoiceY$start1.y, xright = chrChoiceY$end1.y, 
     ybottom = 0, ytop = .01* ylims[2], density = -1, col = 1)
rect(xleft = chrChoiceY$start2.y, xright = chrChoiceY$end2.y, 
     ybottom = 0, ytop = .01* ylims[2], density = -1, col = 1)

mat <- log10(genoMatrixY[ (xlimsY[1]/20e3):(xlimsY[2]/20e3),(xlimsY[1]/20e3):(xlimsY[2]/20e3) ])

par(new=TRUE)
plot(c(rep(NA,squareSize),insulationScore(mat, squareSize = squareSize),rep(NA,squareSize)), type = "l", col = 2, lwd = 2 , xaxs = "i")



# RTN plot 2
plot( chrChoiceY$start1.y + ((chrChoiceY$start2.y - chrChoiceY$start1.y)/2), 
      y = log10( (chrChoiceY$interaction.y/sum(chrChoiceX$interaction.y)) / (chrChoiceY$interaction.y/sum(chrChoiceX$interaction.y)) ), 
      xlim = xlimsY, ylim = c(0,4),type = "n", xaxt = "n", xaxs="i")

for(i in 1:length(g.Y.domainList)){
  domain <- g.Y.domainList[[i]]
  domain <- domain[domain$chr == as.character(chr.Y),]
  if(nrow(domain) > 0){
    rect(xleft = domain$start, ybottom = 4-i, ytop = 5-i,xright = domain$end, col = domainCols[i], density = -1)
  }
}




par(mar=c(7,0,0,0))

whiteRed <- colorRampPalette(c("white","red"))
image(triangulate(mat)[,(nrow(mat)/2):1], col = whiteRed(40), axes = F, ylim = c(.5,1), xlab = chr.Y)

axis(side = 1,at = seq(0,1,length.out = 5),labels = seq(from = xlimsY[1], to = xlimsY[2], length.out = 5))




# maybe we can get the score for the whole matrix and then see where we cna plot it
# change it so we cna see the mean


# we might be seeing changes in domain structure between old repeat hotspots and new repeat hotspots.
# It seems to be an over change in the line 
# the difficutly is now how to actually capture this
# find a few spots where we can see if this observation is persistant.
