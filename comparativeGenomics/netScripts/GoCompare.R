### now to do PCA biplot

rm(list = ls())




panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt, cex = sqrt(r*10))
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.1, txt2, cex = 1.2)
}


for(specRef in c("hg19", "mm10")){
  
  load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/goLists/",specRef,"refInsGoTermLists.RData", sep = ""))
  rownames(cFisher) <- cFisher$GO.ID
  refIns <- cFisher
  
  load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/goLists/",specRef,"refDelGoTermLists.RData", sep = ""))
  rownames(cFisher) <- cFisher$GO.ID
  cFisher <- cFisher[rownames(refIns),]
  refDel <- cFisher
  
  
  load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/goLists/",specRef,"queInsGoTermLists.RData", sep = ""))
  rownames(cFisher) <- cFisher$GO.ID
  cFisher <- cFisher[rownames(refIns),]
  queIns <- cFisher
  
  load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/goLists/",specRef,"queDelGoTermLists.RData", sep = ""))
  rownames(cFisher) <- cFisher$GO.ID
  cFisher <- cFisher[rownames(refIns),]
  queDel <- cFisher
  
  
  
  # what happens when we look at pvalues?
  colChoice <- "classicFisher"
  mat <- as.matrix(data.frame(refIns = as.numeric(gsub("< ","",refIns[,colChoice])), 
                              refDel =as.numeric(gsub("< ","",refDel[,colChoice])),
                              queIns =as.numeric(gsub("< ","",queIns[,colChoice])), 
                              queDel =as.numeric(gsub("< ","",queDel[,colChoice]))
  ))
  
  
  mat <- -log10(mat)
  rownames(mat) = rownames(refIns)
  
  mat <- mat[rownames(mat) != ("GO:0008150"),]
  
  
  if(specRef == "hg19"){
    colnames(mat) <- c("hg19 gain", "hg19 loss", "mm10 gain", "mm10 loss")
    hg19mat <- mat
    hg19stats = refIns
  }
  if(specRef == "mm10"){
    colnames(mat) <- c("mm10 gain", "mm10 loss", "hg19 gain", "hg19 loss")
    mm10mat <- mat
    mm10stats = refIns
  }
  
  
}

mm10mat <- mm10mat[,colnames(hg19mat)]


# make the plots pretty

gSize <- hg19stats$Annotated * 4/(max(hg19stats$Annotated))

pdf(file = "~/Desktop/RTN_domains/RTN_domain_plots/netGainLoss/goTables/hg19goComp.pdf")
pairs(hg19mat, upper.panel = panel.cor, cex = gSize, main = "hg19")
dev.off()

gSize <- mm10stats$Annotated * 4/(max(mm10stats$Annotated))

pdf(file = "~/Desktop/RTN_domains/RTN_domain_plots/netGainLoss/goTables/mm10goComp.pdf")
pairs(mm10mat, upper.panel = panel.cor, cex = gSize, main = "mm10")
dev.off()



# we can show both species genomes



