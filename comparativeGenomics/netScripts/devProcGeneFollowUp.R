

library(GenomicRanges)

rm(list = ls())

options(stringsAsFactors = FALSE)

specRef = "hg19"
specQue = "mm10"

load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/shiftData/",specRef,".expand.breaks.RData", sep = ""))
hg19refShift <- newSynthRefShift
load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/shiftData/",specQue,".expand.breaks.RData", sep = ""))
mm10refShift <- newSynthRefShift

load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/", specRef,"synthBinNorm.RData", sep = ""))
hg19SynthBinNorm.gr <- synthBinNorm.gr
load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/", specQue,"synthBinNorm.RData", sep = ""))
mm10SynthBinNorm.gr <- synthBinNorm.gr

load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/hotspots/", specRef, "repNoRep.RData", sep = ""))
hg19SigRange <- NULL
for(i in 1:4){
  if(names(repSigRanges)[i] == "refDel"){
    sigRange <- repSigRanges[[i]]
  }else{
    sigRange <- intersect(repSigRanges[[i]], noRepSigRanges[[i]])
  }
  hg19SigRange <- c(hg19SigRange, list(sigRange))
}
names(hg19SigRange) <- names(repSigRanges)

load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/hotspots/", specQue, "repNoRep.RData", sep = ""))
mm10SigRange <- NULL
for(i in 1:4){
  if(names(repSigRanges)[i] == "queDel"){
    sigRange <- repSigRanges[[i]]
  }else{
    sigRange <- intersect(repSigRanges[[i]], noRepSigRanges[[i]])
  }
  mm10SigRange <- c(mm10SigRange, list(sigRange))
}
names(mm10SigRange) <- names(repSigRanges)

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")



# import expression data and only select gene wich occur once in human
exprDiv <- read.csv("~/Desktop/expression_divergence.csv")
exprDiv <- exprDiv[exprDiv$Homology.Type == "ortholog_one2one",]
exprDiv <- exprDiv[!(exprDiv$Human.Entrez.ID %in% exprDiv$Human.Entrez.ID[duplicated(exprDiv$Human.Entrez.ID)]),]
exprDiv <- exprDiv[complete.cases(exprDiv),]
# read in human gene data for hg19
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
hg19gene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
mm10gene <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

exprDiv <- exprDiv[exprDiv$Human.Entrez.ID %in% names(hg19gene),]
exprDiv <- exprDiv[exprDiv$Mouse.Entrez.ID %in% names(mm10gene),]

hg19gene <- hg19gene[names(hg19gene) %in% exprDiv$Human.Entrez.ID]
hg19gene <- hg19gene[as.character(exprDiv$Human.Entrez.ID)]
mcols(hg19gene) <- exprDiv


hg19Stretch <- genoExpandStretch(hg19gene, synthGenome = hg19refShift, 
                expandedSeqlengths = seqlengths(hg19SynthBinNorm.gr))
strand(hg19Stretch) <- strand(hg19gene)

hg19Pro <- promoters(hg19Stretch)
outRange <- boxplot.stats(hg19Pro$Differential.connectivity)$stats[c(1,5)]
hg19Pro$divergent <- hg19Pro$Differential.connectivity > outRange[2] | hg19Pro$Differential.connectivity < outRange[1]
hg19Pro$bottomCCG <- rank(hg19Pro$Commonly.co.expressed.genes.ONE.to.ONE.homologs) < .1 * length(hg19Pro)
hg19Pro$topCCG <- rank(hg19Pro$Commonly.co.expressed.genes.ONE.to.ONE.homologs) > .9 * length(hg19Pro)


mm10gene <- mm10gene[names(mm10gene) %in% exprDiv$Mouse.Entrez.ID]
mm10gene <- mm10gene[as.character(exprDiv$Mouse.Entrez.ID)]
mcols(mm10gene) <- exprDiv

mm10Stretch <- genoExpandStretch(mm10gene, synthGenome = mm10refShift,
                               expandedSeqlengths = seqlengths(mm10SynthBinNorm.gr))
strand(mm10Stretch) <- strand(mm10gene)

mm10Pro <- promoters(mm10Stretch)
outRange <- boxplot.stats(mm10Pro$Differential.connectivity)$stats[c(1,5)]
mm10Pro$divergent <- mm10Pro$Differential.connectivity > outRange[2] | mm10Pro$Differential.connectivity < outRange[1]
mm10Pro$bottomCCG <- rank(mm10Pro$Commonly.co.expressed.genes.ONE.to.ONE.homologs) < .1 * length(mm10Pro)
mm10Pro$topCCG <- rank(mm10Pro$Commonly.co.expressed.genes.ONE.to.ONE.homologs) > .9 * length(mm10Pro)



goTermList <- read.csv("~/Documents/dna_turnover/workStationDesktop/RTN_domains/RTN_domain_plots/comparativeGenomics/geneBasedGoTetmAnnotation.csv", row.names = 1)
devProcGOID <- goTermList$parentName[goTermList$parentTerm == "developmental process"][1]


mcols(hg19Pro)$devTermLoss <- mcols(hg19Pro)$Human.Entrez.ID %in% goTermList$ENTREZID[goTermList$genomicBackground == "hg19" & 
                                                                          goTermList$hotspotType == "queDel" & 
                                                                          goTermList$parentName == devProcGOID]


mcols(mm10Pro)$devTermLoss <- exprDiv$Mouse.Entrez.ID %in% goTermList$ENTREZID[goTermList$genomicBackground == "mm10" & 
                                                                   goTermList$hotspotType == "refDel" & 
                                                                   goTermList$parentName == devProcGOID]


length(mm10Pro)
length(hg19Pro)


hg19orthoNo <- hg19PercentOL <- hg19MatPv <- hg19MatOr <- matrix(NA, nrow = 5, ncol = 3,
                  dimnames = list(c("refIns", "refDel", "queIns", "queDel", "devTermLoss"),
                                  c("divergent", "bottomCCG", "topCCG")))

for(i in c("refIns", "refDel", "queIns", "queDel", "devTermLoss")){
  for(j in c("divergent", "bottomCCG", "topCCG")){
    if(i == "devTermLoss"){
      ol <- hg19Pro$devTermLoss
    }else{
      ol <- overlapsAny(hg19Pro, hg19SigRange[[i]])
    }
    tab <- table(data.frame(ol, mcols(hg19Pro)[[j]]))
    tab <- tab[c("TRUE", "FALSE"), c("TRUE", "FALSE")]

    fTest <- fisher.test(tab)
    hg19MatOr[i,j] <- fTest$estimate
    hg19MatPv[i,j] <- fTest$p.value
    hg19PercentOL[i,j] <- tab[1,1]/sum(ol)
    hg19orthoNo[i,j] <- sum(ol)
  }
}





mm10PercentOL <- mm10orthoNo <- mm10MatPv <- mm10MatOr <- matrix(NA, nrow = 5, ncol = 3,
                                 dimnames = list(c("refIns", "refDel", "queIns", "queDel", "devTermLoss"),
                                                 c("divergent", "bottomCCG", "topCCG")))
for(i in c("refIns", "refDel", "queIns", "queDel", "devTermLoss")){
  for(j in c("divergent", "bottomCCG", "topCCG")){
    if(i == "devTermLoss"){
      ol <- mm10Pro$devTermLoss
    }else{
      ol <- overlapsAny(mm10Pro, mm10SigRange[[i]])
    }
    tab <- table(data.frame(ol, mcols(mm10Pro)[[j]]))
    tab <- tab[c("TRUE", "FALSE"), c("TRUE", "FALSE")]
    
    fTest <- fisher.test(tab)
    mm10MatOr[i,j] <- fTest$estimate
    mm10MatPv[i,j] <- fTest$p.value
    mm10PercentOL[i,j] <- tab[1,1]/sum(ol)
    mm10orthoNo[i,j] <- sum(ol)
  }
}

genoRate <- colSums(as.data.frame(mcols(mm10Pro))[,c("divergent", "bottomCCG", "topCCG")])
#colSums(as.data.frame(mcols(hg19Pro))[,c("divergent", "bottomCCG", "topCCG")])



# so lets print this off in tables 
# it appears we are seeing enrichemnt for genes with low numbers of commonly co-exoressed homologs between human and mouse


# next we classify genes as belonging to GRBs

# if genes that belong to certain are 

GRB <- read.table("~/Documents/dna_turnover/workStationDesktop/RTN_domains/RTN_domain_plots/comparativeGenomics/GRBs/41467_2017_524_MOESM2_ESM.txt", skip = 1)
GRB.gr <- GRanges(seqnames = Rle(GRB$V1), 
                  ranges = IRanges(start = GRB$V2, end = GRB$V3))

hg19gene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
hg19Pro <- promoters(hg19gene)
olGBR <- overlapsAny(hg19Pro, GRB.gr)

hg19ProStretch <-genoExpandStretch(hg19Pro, synthGenome = hg19refShift,
                                   expandedSeqlengths = seqlengths(hg19SynthBinNorm.gr))
hg19Pro$devTermLoss <- names(hg19Pro) %in% goTermList$ENTREZID[goTermList$genomicBackground == "hg19" & 
                                                                                 goTermList$hotspotType == "queDel" & 
                                                                                 goTermList$parentName == devProcGOID]


hg19MatOr <- data.frame(hg19MatOr)
hg19MatPv <- data.frame(hg19MatPv)
hg19orthoNo <- data.frame(hg19orthoNo)

hg19MatOr[,"GBR"] <- NA
hg19MatPv[,"GBR"] <- NA
hg19orthoNo[,"nonOrtho"] <- NA
for(i in c("refIns", "refDel", "queIns", "queDel", "devTermLoss")){
  if(i == "devTermLoss"){
    olRegion <- hg19Pro$devTermLoss
  }else{
    olRegion <- overlapsAny(hg19ProStretch, hg19SigRange[[i]])
  }
  
  tab <- table(data.frame(olGBR, olRegion))
  tab <- tab[c("TRUE", "FALSE"), c("TRUE", "FALSE")]
  
  fTest <- fisher.test(tab)
  
  hg19MatOr[i,"GBR"] <- fTest$estimate
  hg19MatPv[i, "GBR"] <- fTest$p.value
  hg19orthoNo[i,"nonOrtho"] <- sum(olRegion)
}



mm10gene <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
mm10Pro <- promoters(mm10gene)

mm10ProStretch <-genoExpandStretch(mm10Pro, synthGenome = mm10refShift,
                                   expandedSeqlengths = seqlengths(mm10SynthBinNorm.gr))
mm10Pro$devTermLoss <- names(mm10Pro) %in% goTermList$ENTREZID[goTermList$genomicBackground == "mm10" & 
                                                                 goTermList$hotspotType == "refDel" & 
                                                                 goTermList$parentName == devProcGOID]




mm10orthoNo <- data.frame(mm10orthoNo)
mm10orthoNo[,"nonOrtho"] <- NA
for(i in c("refIns", "refDel", "queIns", "queDel", "devTermLoss")){
  if(i == "devTermLoss"){
    olRegion <- mm10Pro$devTermLoss
  }else{
    olRegion <- overlapsAny(mm10ProStretch, mm10SigRange[[i]])
  }
  mm10orthoNo[i,"nonOrtho"] <- sum(olRegion)
}

# there is no significant overlap between genes in GBRs and genes in sigRanges

# interestingly development assocaited terms from increased loss regions are highly enriched in GBRs


# In reality this should be all the tables I need,
# we shoudl report the number used for each catacory




round(hg19PercentOL * 100, 1)

round(hg19MatOr, 2)

round(hg19MatPv, 3)

round(matrix(p.adjust(as.matrix(hg19MatPv), "fdr"), nrow = 5, dimnames = dimnames(hg19MatPv)),3)





round(mm10PercentOL * 100, 1)

round(mm10MatOr, 2)

round(mm10MatPv, 3)

round(matrix(p.adjust(as.matrix(mm10MatPv), "fdr"), nrow = 5, dimnames = dimnames(mm10MatPv)),3)




hg19orthoNo

mm10orthoNo

sum(olGBR)

genoRate

length(hg19gene)

length(mm10gene)

nrow(exprDiv)



