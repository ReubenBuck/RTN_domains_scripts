
rm(list = ls())
library(GenomicRanges)
options(stringsAsFactors = FALSE)

genomes <- c(ref = "mm10", que = "hg19")

if(genomes["ref"] == "hg19"){
  specClade <- c(ref = "primate", que = "rodent")
}else{
  specClade <- c(ref = "rodent", que = "primate")
}

inGroupSpecList <-list(primate = c("micMur", "tarSyr", "papHam", "panTro"),rodent = c("ochPri", "dipOrd", "rn"))

inDelGenomeList <- NULL

for(g in 1:2){
  
  filePath <- paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/data/comparativeGenomics/netAlignment/inGroupSpecies/processed/",genomes[3 -(g%%2 + 1)],"/",specClade[g], sep = "")
  newInsGenomes <- list.files(filePath, full.names = TRUE, recursive = TRUE)
  newInsGenomesNames <- list.files(filePath)
  
  filePath <- paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/data/comparativeGenomics/netAlignment/inGroupSpecies/processed/",genomes[3 -((g + 1)%%2 + 1)],"/",specClade[g], sep = "")
  newDelGenomes <- list.files(filePath, full.names = TRUE, recursive = TRUE)
  newDelGenomesNames <- list.files(filePath)
  
  # what order are we reading in the species?
  
  
  insGenomeList <- NULL
  combinedRanges <- GRanges()
  for(i in 1:length(inGroupSpecList[[specClade[g]]])){
    readPath <- newInsGenomes[grep(rev(inGroupSpecList[[specClade[g]]])[i], newInsGenomesNames)]
    inGroupSpec <- read.table(file = readPath, header = FALSE,
                              col.names = c("seqnames", "start", "end", 
                                            "queChr", "queStart", "queEnd", 
                                            "queStrand", "chainID"))
    inGroup.gr <- GRanges(inGroupSpec[,1:3])
    insGenomeList <- c(insGenomeList, list(inGroup.gr))
    if(i > 1){
      insGenomeList[[i]] <- setdiff(insGenomeList[[i]], combinedRanges)
    }
    combinedRanges <- reduce(c(combinedRanges, inGroup.gr))
  }
  names(insGenomeList) <- rev(inGroupSpecList[[specClade[g]]])
  
  
  delGenomeList <- NULL
  combinedRanges <- GRanges()
  for(i in 1:length(inGroupSpecList[[specClade[g]]])){
    readPath <- newDelGenomes[grep(inGroupSpecList[[specClade[g]]][i], newDelGenomesNames)]
    inGroupSpec <- read.table(file = readPath, header = FALSE,
                              col.names = c("seqnames", "start", "end", 
                                            "queChr", "queStart", "queEnd", 
                                            "queStrand", "chainID"))
    inGroup.gr <- GRanges(inGroupSpec[,1:3])
    delGenomeList <- c(delGenomeList, list(inGroup.gr))
    if(i > 1){
      delGenomeList[[i]] <- setdiff(delGenomeList[[i]], combinedRanges)
    }
    combinedRanges <- reduce(c(combinedRanges, inGroup.gr))
  }
  names(delGenomeList) <- inGroupSpecList[[specClade[g]]]
  
  
  inDelGenomeList <- c(inDelGenomeList, list(list(insGenomeList = insGenomeList, delGenomeList = delGenomeList)))
  
}

names(inDelGenomeList) <- c("ref", "que")


# we might as well readin the other stages
# this way we can get a complete picture of the timing

# the reuslt is a list 
# human 
#   human with insertion in human
#   mouse with deletion from human 
# mouse
#   mouse with insertion in mouse
#   human with deletion from mouse




filePath <- paste("Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/",
                  genomes["ref"],".stretch.rep.RData", sep = "")
load(filePath)


# ol <- findOverlaps(mcols(stretchedRef.gr)[["refRanges"]], minoverlap = 1)
# ol <- ol[!(isSelfHit(ol))]
# ol <- ol[(stretchedRef.gr[queryHits(ol)]$type == "refIns" | stretchedRef.gr[queryHits(ol)]$type == "queDel") ]
# ol <- ol[(stretchedRef.gr[subjectHits(ol)]$type == "refIns" | stretchedRef.gr[subjectHits(ol)]$type == "queDel") ]
# 
# 
# end(mcols(refS.gr)[[dataSelection$range[j]]][queryHits(ol)]) = start(mcols(refS.gr)[[dataSelection$range[j]]][queryHits(ol)]) + width(refS.gr)[queryHits(ol)]
# ol <- ol[!isRedundantHit(ol)]
# mcols(refS.gr)[[dataSelection$range[j]]][subjectHits(ol)] <- shift(mcols(refS.gr)[[dataSelection$range[j]]][subjectHits(ol)], shift = width(mcols(refS.gr)[[dataSelection$range[j]]][queryHits(ol)]))
# 
# 
# ol <- findOverlaps(mcols(stretchedRef.gr)[["queRanges"]], minoverlap = 1)
# ol <- ol[!(isSelfHit(ol))]
# ol <- ol[(stretchedRef.gr[queryHits(ol)]$type == "queIns" | stretchedRef.gr[queryHits(ol)]$type == "refDel") ]
# ol <- ol[(stretchedRef.gr[subjectHits(ol)]$type == "queIns" | stretchedRef.gr[subjectHits(ol)]$type == "refDel") ]




types = c("Ins", "Del")


dataSelection <- data.frame(type = c("refIns", "refDel", "queIns", "queDel"),
                            range = c("refRanges", "queRanges", "queRanges", "refRanges"),
                            lineage = c(rep(specClade["ref"], 2), rep(specClade["que"], 2)),
                            genome = c("ref", "ref","que", "que"),
                            listChoice = rep(c("insGenomeList", "delGenomeList"),2)
                            )


newInGroups <- GRanges()
for(j in 1:4){
  refS.gr <- stretchedRef.gr[stretchedRef.gr$type == dataSelection$type[j]]
  
  ol <- findOverlaps(mcols(refS.gr)[[dataSelection$range[j]]])
  ol <- ol[!(isSelfHit(ol))]
  end(mcols(refS.gr)[[dataSelection$range[j]]][queryHits(ol)]) = start(mcols(refS.gr)[[dataSelection$range[j]]][queryHits(ol)]) + width(refS.gr)[queryHits(ol)]
  ol <- ol[!isRedundantHit(ol)]
  mcols(refS.gr)[[dataSelection$range[j]]][subjectHits(ol)] <- shift(mcols(refS.gr)[[dataSelection$range[j]]][subjectHits(ol)], shift = width(mcols(refS.gr)[[dataSelection$range[j]]][queryHits(ol)]))
  
  
  for(i in names(inDelGenomeList[[dataSelection$genome[j]]][[dataSelection$listChoice[j]]])){
    
    ref.gr <- mcols(refS.gr)[[dataSelection$range[j]]]
    indel.gr <- inDelGenomeList[[dataSelection$genome[j]]][[dataSelection$listChoice[j]]][[i]]
    
    ol <- findOverlaps(ref.gr, indel.gr)
    
    pInt <- pintersect( ref.gr[queryHits(ol)], indel.gr[subjectHits(ol)]) 
    
    relStart <- start(pInt) - start(ref.gr[queryHits(ol)])
    if(dataSelection$range[j] == "queRanges" ){
      minusStrand <- refS.gr$sData[queryHits(ol)] == "-"
      relStart[minusStrand] <- (end(ref.gr[queryHits(ol)]) - end(pInt))[minusStrand]

    }
    
    newRange <- GRanges(seqnames = seqnames(refS.gr[queryHits(ol)]),
                        ranges = IRanges(start = start(refS.gr[queryHits(ol)]) + relStart,
                                         width = width(pInt)),
                        lineage = rep(dataSelection$lineage[j], length(pInt)),
                        species = rep(i, length(pInt)),
                        type = rep(dataSelection$type[j], length(pInt))
    )
    newInGroups <- c(newInGroups, newRange)
  }
  newRange <- setdiff(refS.gr, newInGroups)
  mcols(newRange)[["lineage"]] <- dataSelection$lineage[j]
  mcols(newRange)[["species"]] <- genomes[dataSelection$genome[j]]
  mcols(newRange)[["type"]] <- dataSelection$type[j]
  newInGroups <- c(newInGroups, newRange)
  
}


# found the bug, need to go through and edit the stretched ranges
newInGroups <- sort(sortSeqlevels(newInGroups))

seqinfo(newInGroups) <- seqinfo(stretchedRef.gr)

savePath <- paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/ingroupSpecies/",
                  genomes["ref"],".ingroup.RData", sep = "")

save(newInGroups, file = savePath)



