
specRef = "hg19"
queRef = "mm9"

load(paste("Desktop/RTN_domains/R_objects/mappedGaps/",specRef,".",specQue,".netData.RData",sep = ""))



resGenomeStat <- matrix(data = NA, ncol = 9, nrow = 2)
rownames(resGenomeStat) <- genomes
colnames(resGenomeStat) <- c("genome_size", "fills","ancestral" ,
                             "fills_percent_ancestral", "all_gaps", "placed_gaps", 
                             "placed_gaps_percentage","gain", "loss")
resGenomeStat <- as.data.frame(resGenomeStat)

# genome size
## mm10
resGenomeStat$genome_size[1] <- round( (sum(as.numeric(seqlengths(refSeqGaps.gr))) - sum(width(refSeqGaps.gr)))/1e6 )
## hg19
resGenomeStat$genome_size[2] <-round((sum(as.numeric(seqlengths(queSeqGaps.gr))) - sum(width(queSeqGaps.gr)))/1e6)

# all non gaps
## mm10
resGenomeStat$fills[1] <- round(sum(width(refFill.gr))/1e6)
## hg19
resGenomeStat$fills[2] <- round(sum(width(queFill.gr))/1e6)

#ancestral
resGenomeStat$ancestral[1] <- round(sum(width(refAncDna.gr))/1e6)
resGenomeStat$ancestral[2] <- round(sum(width(queAncDna.gr))/1e6)


# ancestral %
# mm10
resGenomeStat$fills_percent_ancestral[1] <- round(sum(width(GenomicRanges::intersect(refFill.gr, refAncDna.gr))) / sum(width(refFill.gr)) * 100, digits = 1)
# hg19
resGenomeStat$fills_percent_ancestral[2] <- round(sum(width(GenomicRanges::intersect(queFill.gr, queAncDna.gr))) / sum(width(queFill.gr)) * 100, digits = 1)


# all gaps
#mm10
resGenomeStat$all_gaps[1] <- round(sum(width(refFillGaps.gr))/1e6)
# hg19
resGenomeStat$all_gaps[2] <- round(sum(width(queFillGaps.gr))/1e6)


# placeable gaps 
#mm10
resGenomeStat$placed_gaps[1] <- round(sum(width(refGap.gr))/1e6)
#hg19
resGenomeStat$placed_gaps[2] <- round(sum(width(queGap.gr))/1e6)



# placeable gaps percentage
#mm10
resGenomeStat$placed_gaps_percentage[1] <- round(sum(width(refGap.gr))/sum(width(refFillGaps.gr))*100,digits = 1)
#hg19
resGenomeStat$placed_gaps_percentage[2] <- round(sum(width(queGap.gr))/sum(width(queFillGaps.gr))*100, digits = 1)


# gain
## mm10
resGenomeStat$gain[1] <- round(sum(width(refGapNonAnc.gr))/1e6)
# hg19
resGenomeStat$gain[2] <- round(sum(width(queGapNonAnc.gr))/1e6)

# loss
# mm10 
resGenomeStat$loss[1] <- round(sum(width(queGapAnc.gr))/1e6)
# hg19
resGenomeStat$loss[2] <- round(sum(width(refGapAnc.gr))/1e6)


