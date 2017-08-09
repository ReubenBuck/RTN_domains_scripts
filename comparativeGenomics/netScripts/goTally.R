


rm(list = ls())

options(stringsAsFactors = FALSE)
library(GO.db)



# terms unique to that level



getAllBPChildren <- function(goids) { 
  ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE)) 
  ans <- ans[!is.na(ans)]
}
getAllBPOffspring <- function(goids) { 
  ans <- unique(unlist(mget(goids, GOBPOFFSPRING), use.names=FALSE)) 
  ans <- ans[!is.na(ans)] 
} 
getAllBPAncestors <- function(goids) { 
  ans <- unique(unlist(mget(goids, GOBPANCESTOR), use.names=FALSE)) 
  ans <- ans[!is.na(ans)] 
} 
getAllBPOffspringList <- function(goids) { 
  ans <- mget(goids, GOBPOFFSPRING)
  ans <- ans[!is.na(ans)] 
} 

goterms = unlist(Term(GOTERM))


level1_BP_terms <- getAllBPChildren("GO:0008150") # 23 terms 
offsprings <- getAllBPOffspring(level1_BP_terms)


# use uniquely high terms at the third level
levelChoice <- level1_BP_terms
offsprings <- getAllBPOffspring(levelChoice)
termRameain <- levelChoice[!(levelChoice %in% offsprings)]
t <- goterms[termRameain]



offS <- getAllBPOffspringList(termRameain)

geneBasedGOTerms <- NULL
geno <- c("hg19", "mm10")
for(g in geno){
  
  if(g == "hg19"){
    fileChoices <- c("refIns", "refDel", "queIns", "queDel")
  }else(
    fileChoices <- c("queIns", "queDel", "refIns", "refDel")
  )
  
  
  dfResult <- data.frame(row.names = termRameain)
  for(f in fileChoices){
    
    
    geneNames <- read.table(paste("~/Desktop/RTN_domains/data/comparativeGenomics/hotspotsGenes/",g,"_",f,"_hotspotGenes.txt", sep = ""), 
                            header = TRUE, stringsAsFactors = FALSE)
    
    load(paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/goLists/",paste(g,f,sep = ""),"GoTermLists.RData", sep = ""))
    classic = TRUE
    
    
    cFisher$classicFisher <- gsub("< ","",cFisher$classicFisher)
    sigPC <- cFisher[p.adjust(cFisher$classicFisher,method = "fdr") < 0.05,"GO.ID"]
    fdrVals <- p.adjust(cFisher$classicFisher,method = "fdr")
    names(fdrVals) <- cFisher$GO.ID
    
    for(i in termRameain){
      dfResult[i,"term"] <- (t[i])
      dfResult[i,f] <- sum(offS[[i]] %in% sigPC)  / length(sigPC) / (sum( offS[[i]] %in% cFisher$GO.ID )/1000)
      print(length(offS[[i]]))
      
      if(any(offS[[i]] %in% sigPC)){
        
        childTerms <- offS[[i]][offS[[i]] %in% sigPC]
        
        if(any(geneNames$GO %in% childTerms)){
          childGeneName <- unique(geneNames[geneNames$GO %in% childTerms,])
          childInfoDF <- data.frame(genomicBackground = g,
                                    hotspotType = f,
                                    parentName = i,
                                    parentTerm = t[i],
                                    childName = childGeneName$GO,
                                    childTerm = goterms[childGeneName$GO],
                                    childFDR = fdrVals[childGeneName$GO],
                                    SYMBOL = childGeneName$SYMBOL,
                                    ENTREZID = childGeneName$ENTREZID)
          geneBasedGOTerms <- rbind(geneBasedGOTerms, childInfoDF)
        }
      }
    }
  }
  rownames(dfResult) <- dfResult$term
  dfResult <- dfResult[,fileChoices]
  dfResult[is.na(dfResult)] <- 0
  
  
  # make a stacked bar plot
 # pdf(paste("~/Desktop/goSum",g,".pdf", sep = ""), onefile = TRUE)
  par(mar= c(1,4,1,1), oma = c(20,3,4,2))
  layout(1:4)

  for(i in 1:4){
    barplot(dfResult[,i], beside = FALSE, col =  RColorBrewer::brewer.pal(4, "Dark2")[i], las = 2,
            ylim = c(0,.8), width = .9,space = .105, xaxs = "i")
    if(i == 1){title(main = g, line = 0, outer = T, cex.main = 1.5)}
    grid()
    legend("topright", c("hg19 gain", "hg19 loss", "mm10 gain", "mm10 loss")[i], bty = "n", cex = 1.5)
  }
  axis(side=1, at = .5:(length(t)-.5), labels  = t, las =2, tick = FALSE)
  mtext(text = "Significant child terms (proportion per 1000 annotations)", side = 2, outer = TRUE, cex = .7)

  mtext(text = "Parent terms", side = 1, outer = FALSE, cex = 1, line = 17)

 # dev.off()
  
  
  
  # the genes that contribute to each high level annotation
  
  
}
# how enriched terms are vs how many enriched terms there are.
# this is why classic works best for this method.

write.csv(geneBasedGOTerms,file = "Desktop/geneBasedGoTetmAnnotation.csv",quote = TRUE)


# report the relevant gene names for each of our GO terms
# simply go through the list of significant terms and pull them out



