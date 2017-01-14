## what if we get the raw matrix for each chrmosome and look at the overlapping 50kb bins

## we can read in our data and go chromosome by chromosome with the hi-c matricies 


## maybe if we grab a few of the rows

library(devtools)

rm(list = ls())


devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")


setwd("~/Desktop/RTN_domains/")

genomes <- c("canFam3", "mm9")
for(g in genomes){
  genoGenePairs <- read.table(paste("data/orthoAnalysis/canFam.vs.mm9.synPair/",g,".out", sep = ""), header = F, 
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

mm9Complete$interaction = mm9Complete$interaction / 
  ((round((mm9Complete$end1 - mm9Complete$start1)/20e3) + 1) * 
     (round((mm9Complete$end2 - mm9Complete$start2)/20e3) + 1))

combined <- merge(x = canFam3Complete, y = mm9Complete, by.x = 1, by.y = 1)

combinedComplete <- combined[complete.cases(combined),]

save(combinedComplete, file = paste("R_objects/comparativeHiC/", genomes[1],"_", genomes[2] ,".RData", sep = ""))
# big thing now is to do the compairson.
# it might also be worth smoothing over the contact map
# take the average contact from the surrounding area
# This way the contanct measurments would be less sensetive to regional issues.






