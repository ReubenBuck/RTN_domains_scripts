## what if we get the raw matrix for each chrmosome and look at the overlapping 50kb bins

## we can read in our data and go chromosome by chromosome with the hi-c matricies 


## maybe if we grab a few of the rows


rm(list = ls())

setwd("~/Desktop/RTN_domains/")

genomes <- c("canFam3", "rheMac8")
for(g in genomes){
  genoGenePairs <- read.table(paste("RTN_domains_scripts/orthoAnalysis/synPair/",g,"synPair.out", sep = ""), header = F, 
                              col.names = c("orthoPair",
                                            "chr1", "start1", "end1", "strand1", "orthoID1", "geneID1",
                                            "chr2", "start2", "end2", "strand2", "orthoID2", "geneID2"))
  
  genoChr <- list.files(paste("data/trialHicAnalysis/",g, "trial/",sep = ""))
  genoChr <- genoChr[grep("chr", genoChr)]
  genoChr <- gsub(pattern = "raw20kb.txt",replacement = "",x = genoChr)
  genoChr <- genoChr[genoChr != "chrM"] 
  genoChr <- genoChr[genoChr != "chrY"]
  
  genoComplete <- NULL
  for(chr in 1:length(genoChr)){
    
    genoMatrix <- read.delim(paste("data/trialHicAnalysis/",g,"trial/",genoChr[chr],"raw20kb.txt", sep = ""), 
                             header = TRUE, sep = "\t",row.names = 1)
    genoMatrix <- as.matrix(genoMatrix[,2:ncol(genoMatrix)])
    
    genoGenePairsChr <- genoGenePairs[genoGenePairs$chr1 == genoChr[chr],]
    coord <- cbind(round(genoGenePairsChr[,c(3,4,9,10)]/20000), genoGenePairsChr[,c(5,11)])
    
    genoGenePairsChr$interaction <- NA
    for(i in 1:nrow(coord)){
      interaction <- NA
    #  interaction <- genoMatrix[
     #   c(coord[i,1], coord[i,2])[(coord[i,5] == "-") + 1],
      #  c(coord[i,3], coord[i,4])[(coord[i,6] == "-") + 1]
       # ]
      
      interactionR <- coord[i,1]:coord[i,2]
      interactionC <- coord[i,3]:coord[i,4]
      if(any(interactionR %in% interactionC)){
        interaction = NA
      } else if ( any(c(interactionC, interactionR) == 0 )){
        interaction = NA
      }else {
        interaction <- sum(genoMatrix[interactionR,interactionC])
      }
      
      genoGenePairsChr$interaction[i] <- interaction
      
    }
    
    genoComplete <- rbind(genoComplete, genoGenePairsChr)
    print(genoChr[chr])
  }
  
  assign(x = paste(g,"Complete",sep = ""),value = genoComplete)
  
}
canFam3Complete$interaction = canFam3Complete$interaction / 
  ((round((canFam3Complete$end1 - canFam3Complete$start1)/20e3) + 1) * 
      (round((canFam3Complete$end2 - canFam3Complete$start2)/20e3) + 1))

rheMac8Complete$interaction = rheMac8Complete$interaction / 
  ((round((rheMac8Complete$end1 - rheMac8Complete$start1)/20e3) + 1) * 
     (round((rheMac8Complete$end2 - rheMac8Complete$start2)/20e3) + 1))

combined <- merge(x = canFam3Complete, y = rheMac8Complete, by.x = 1, by.y = 1)

combinedComplete <- combined[complete.cases(combined),]


plot(combinedComplete$interaction.y[order(combinedComplete$interaction.y)]/sum(combinedComplete$interaction.y, na.rm = T), type = "l")
lines(combinedComplete$interaction.x[order(combinedComplete$interaction.y)]/sum(combinedComplete$interaction.x, na.rm = T), col = 2, cex = .3)


chr1x <- combined[combined$chr1.x == "chr1",]
chr1x <- chr1x[order(chr1x$start1.x),]

plot(chr1x$start2.x, log10(chr1x$interaction.x), type = "l")
lines(chr1x$start2.x, log10(chr1x$interaction.y), col = 2)


chr1y <- combinedComplete[combinedComplete$chr1.y == "chr1",]
chr1y <- chr1y[order(chr1y$start1.y),]

chrCol <- as.integer(as.factor(chr1y$chr1.x))
chrCol <- topo.colors(39)[chrCol]

xlims <- c(30e6,50e6)
layout(c(1,2))
plot(chr1y$start2.y, (chr1y$interaction.y)/sum(chr1y$interaction.y,na.rm = TRUE), type = "p", xlim= xlims, pch =16 , cex = .3, ylim = c(0,0.01))
lines(chr1y$start2.y, (chr1y$interaction.x)/sum(chr1y$interaction.x,na.rm = TRUE), pch =16 , cex = .3, type = "h", col = chrCol)

plot(chr1y$start2.y, 
     log10(((chr1y$interaction.y)/sum(chr1y$interaction.y,na.rm = TRUE))/((chr1y$interaction.x)/sum(chr1y$interaction.x,na.rm = TRUE))),
     xlim = xlims, type = "h", col = chrCol)
abline(h= 1,lty = 2)
abline(h= -1, lty = 2)
#plot(chr1y$start2.y, 
 #    ((chr1y$interaction.y)/sum(chr1y$interaction.y,na.rm = TRUE))-((chr1y$interaction.x)/sum(chr1y$interaction.x,na.rm = TRUE)),
  #   xlim = xlims)

#image(log10(genoMatrix)[(xlims[1]/20e3):(xlims[2]/20e3),(xlims[1]/20e3):(xlims[2]/20e3)])


chr = "chr1"
chrChoiceX <- combinedComplete[combinedComplete$chr1.x == chr,]
chrChoiceX <- chrChoiceX[order(chrChoiceX$start1.x),]

g = "canFam3"
genoMatrix <- read.delim(paste("data/trialHicAnalysis/",g,"trial/",chr,"raw20kb.txt", sep = ""), 
                         header = TRUE, sep = "\t",row.names = 1)
genoMatrix <- as.matrix(genoMatrix[,2:ncol(genoMatrix)])


chrCol <- as.integer(as.factor(chrChoiceX$chr1.y))
chrCol <- topo.colors(39)[chrCol]

xlims <- c(10e6,105e6)
ylims <- c(0, 5)
layout(c(1,2), height = c(.5,2),width=c(2,2))
par(mar=c(0,0,0,0), oma=c(5,5,5,5))
plot(chrChoiceX$start1.x, y = (chrChoiceX$interaction.x/sum(chrChoiceX$interaction.x)), 
     xlim = xlims, type = "n", xaxt = "n", ylim = ylims)

for( i in 1:nrow(chrChoiceX)){
polygon(x = c(chrChoiceX$start1.x[i], 
              (chrChoiceX$start1.x + ((chrChoiceX$end2.x - chrChoiceX$start1.x)/2))[i] ,
              chrChoiceX$end2.x[i]), 
        y = (c(0, 
               -log10(chrChoiceX$interaction.x[i]/sum(chrChoiceX$interaction.x, na.rm = T)),
               0 )))
}
points((chrChoiceX$start1.x + ((chrChoiceX$end2.x - chrChoiceX$start1.x)/2)), 
       y = -log10(chrChoiceX$interaction.y/sum(chrChoiceX$interaction.y)), pch = 16, col = chrCol)
rect(xleft = chrChoiceX$start1.x, xright = chrChoiceX$end1.x, ybottom = 0, ytop = .05* ylims[2], density = -1, col = 1)
rect(xleft = chrChoiceX$start2.x, xright = chrChoiceX$end2.x, ybottom = 0, ytop = .05* ylims[2], density = -1, col = 1)

mat <- log10(genoMatrix[ (xlims[1]/20e3):(xlims[2]/20e3),(xlims[1]/20e3):(xlims[2]/20e3) ])
image(mat[upper.tri(mat)])


top <- (chrChoiceX$start1.x - xlims[1])/(xlims[2] - xlims[1])
bottom <- (chrChoiceX$end1.x - xlims[1])/(xlims[2] - xlims[1])
left <- (chrChoiceX$start2.x - xlims[1])/(xlims[2] - xlims[1])
right <- (chrChoiceX$end2.x - xlims[1])/(xlims[2] - xlims[1])

rect(xleft = left,xright = right, ytop = top, ybottom = bottom)
# interaction frequecy seems to be correlated with gene distance.



smoothScatter(log10(combinedComplete$end2.x - combinedComplete$start1.x), log10(combinedComplete$interaction.x))


layout(1)
smoothScatter(log10(combinedComplete$interaction.x/sum(combinedComplete$interaction.x, na.rm = T)), 
              log10(combinedComplete$interaction.y/sum(combinedComplete$interaction.y, na.rm = T)),nrpoints = Inf)
plot(log10(combinedComplete$interaction.x/sum(combinedComplete$interaction.x, na.rm = T)), 
     log10(combinedComplete$interaction.y/sum(combinedComplete$interaction.y, na.rm = T)), 
     pch = 16, cex = .1)

cor((combinedComplete$interaction.x),(combinedComplete$interaction.y))

plot((interactionPairs), ylim = c(0,100))



#image(log10(genoMatrix)[1:100,1:100])


hist(log10(combined$end1.x - combined$start1.x))



# stuff that had overlaps has been removed. 

# How do we transform the matrix


mat2 <- matrix(rnorm(100),nrow = 10)
mat2[lower.tri(mat2)] <- mat2[upper.tri(mat2)]








