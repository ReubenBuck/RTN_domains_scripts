rm(list = ls())


setwd("~/Desktop/RTN_domains/")

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")


g.X = "canFam3"

g.Y = "mm9"

genomes <- c(g.X,g.Y)

load(paste("R_objects/comparativeHiC/",g.X,"_",g.Y,".RData",sep = ""))


head(combinedComplete)



chr.X = "chr3"
chrChoiceX <- combinedComplete[combinedComplete$chr1.x == chr.X,]
chrChoiceX <- chrChoiceX[order(chrChoiceX$start1.x),]

chrCols <- as.factor(as.character(chrChoiceX$chr2.y))

layout(1);par(mar=c(0,0,0,0), oma = c(5,5,5,5))
plot(chrChoiceX$start1.x, chrChoiceX$start1.y, 
     col = rainbow(length(levels(chrCols)))[chrCols])
legend("topright", legend = levels(chrCols), fill = rainbow(length(levels(chrCols))), bty = "n")


head(chrChoiceX)

refDist <- chrChoiceX$start2.x - chrChoiceX$start1.x

queDist <- chrChoiceX$start2.y - chrChoiceX$start1.y



plot(chrChoiceX$start1.x,refDist)
lines(chrChoiceX$start1.x , abs(queDist))


# gene pairs with differences in space. 

load(paste("R_objects/temp/",g.X,"_domainList.RData", sep = "")); assign(paste("g.X.domainList", sep = "_"), domainList)
load(paste("R_objects/temp/",g.Y,"_domainList.RData", sep = "")); assign(paste("g.Y.domainList", sep = "_"), domainList)
domainCols <- c("darkblue", "aquamarine3", "purple", "red")


plot(chrChoiceX$start1.x,refDist - abs(queDist))
abline(h=0)

for(i in 1:length(g.X.domainList)){
  domain <- g.X.domainList[[i]]
  domain <- domain[domain$chr == as.character(chr.X),]
  if(nrow(domain) > 0){
    rect(xleft = domain$start, ybottom = -1500000, ytop = -1400000,xright = domain$end, col = domainCols[i], density = -1)
  }
}

# see if the difference in distance change correlates with changes 









