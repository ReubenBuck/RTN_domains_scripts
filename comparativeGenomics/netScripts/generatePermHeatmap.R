
rm(list = ls())


load("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/hotspots/hg19.permZ.RData")
hg19Z <- zAll

load("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/hotspots/mm10.permZkeepX.RData")
mm10Z <- zAll
load("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/hotspots/mm10.permZnoX.RData")
mm10ZnoX <- zAll

# replace recombination data in mouse because of missing X chromosome
mm10Z[,c("recHotSpot", "recombRate")] <- mm10ZnoX[,c("recHotSpot", "recombRate")]

rm(mm10ZnoX)

rownames(mm10Z) <- c("mm10.mm10.gain", "mm10.mm10.loss", "mm10.hg19.gain", "mm10.hg19.loss")
rownames(hg19Z) <- c("hg19.hg19.gain", "hg19.hg19.loss", "hg19.mm10.gain", "hg19.mm10.loss")


# generate col order based on clustering
hg19Z <- hg19Z[,!(colnames(hg19Z) == "dnaseActivity" | colnames(hg19Z) == "L1Motif" | colnames(hg19Z) == "prdm9Motif" | colnames(hg19Z) == "ctcfMotif")]
mm10Z <- mm10Z[,!(colnames(mm10Z) == "dnaseActivity" | colnames(mm10Z) == "L1Motif" | colnames(mm10Z) == "prdm9Motif" | colnames(mm10Z) == "ctcfMotif")]

# correct name and grouping for cols read in
colName <- c("DNase1 HS peaks", "LADs","Exons", "Introns", 
             "Ancient elements", "Lineage-specifc L1s", "Lineage-specifc SINEs",
             "Ancestral L1s","Recombination hotspots", "CpG islands",
             "GC content","Recombination rate")
colourNum <- c(1,2,1,1,3,4,4,3,5,1,1,5)

# arrange mm10 columns
mm10Z <- mm10Z[c("mm10.hg19.gain", "mm10.hg19.loss","mm10.mm10.gain", "mm10.mm10.loss"),]



zAll <- rbind(hg19Z,mm10Z)


# get genome distribtuion of feature indicators
# paste them together and get the corealtion
dfMcolAll <- NULL
for(specRef in c("hg19", "mm10")){
load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/",specRef,".synthBin.genome.variable.RData", sep = ""))

df<- data.frame(synthBin.gr)
colChoice <- c("dnasePeaks","ladCov","exon", "intron","ctcfMotif", 
               "L1Motif", "prdm9Motif", "ancient", "new_L1", "new_SINE", 
               "old_L1", "recHotSpot","cpgCov")

remainRefBases <- width(synthBin.gr)[1] - (df$queIns + df$refDel + df$seqGap)
df[remainRefBases > 0,colChoice] <- (df[remainRefBases > 0,colChoice]/remainRefBases[remainRefBases > 0]) * width(synthBin.gr)[1]
remainingBases <- (df$end - df$start + 1) - (df$missingGap + df$seqGap)
df <- df[remainingBases > 150001,]
df <- df[complete.cases(df),]
synthBinNormVar.gr <- GRanges(df)
seqinfo(synthBinNormVar.gr) <- seqinfo(synthBin.gr)
chrAll <- seqlevels(synthBinNormVar.gr)
chrAll <- chrAll[-grep("_", chrAll)]
chrAll <- chrAll[!(chrAll == "chrM" | chrAll == "chrY")]
synthBinNormVar.gr <- synthBinNormVar.gr[seqnames(synthBinNormVar.gr) %in% chrAll]
dfMcol <- data.frame(mcols(synthBinNormVar.gr))
dfMcol <- dfMcol[,c(colChoice,"gcContent", "dnaseActivity", "recombRate")]
dfMcolAll <- rbind(dfMcolAll, dfMcol)
}

dfMcolAll <- dfMcolAll[,!(colnames(dfMcolAll) == "dnaseActivity" | colnames(dfMcolAll) == "L1Motif" | colnames(dfMcolAll) == "prdm9Motif" | colnames(dfMcolAll) == "ctcfMotif")]
colnames(dfMcolAll) <- colName

# clustering pattern of genoic distribution of features
hClust <- hclust(dist(cor(dfMcolAll)))


# sort features according to their genomic locations along with their names
mm10Z <- mm10Z[,hClust$order]
hg19Z <- hg19Z[,hClust$order]
colName <- colName[hClust$order]
colourNum <- colourNum[hClust$order]


# need 5 colours for feature associations
GenomeCols <- RColorBrewer::brewer.pal(n = 5, name = "Set1")


# set the colour max
zAll <- 15

# plotting results

pdf(file = "~/Documents/dna_turnover/workStationDesktop/RTN_domains/RTN_domain_plots/netGainLoss/hotspotPlot/heatmapEasy.pdf")

#pdf(file = "~/Desktop/heatmap.pdf")

# set up plot space
layout(matrix(1:10, nrow = 2, byrow = TRUE), width = c(8,6,8,8,5), heights = c(5,1))
pal <- colorRampPalette(c("blue", "white", "red"))
par(mar = c(1,1,5,1), oma = c(0,1,0,0))


# plot the dendrogram
plot(as.dendrogram(hClust), horiz = T,yaxt = "n",leaflab = "none", xaxs = "i",ylim = c(1,ncol(hg19Z)))

# plot hg19 matrix and row names
plot.new()
title("Genomic features")

# set hg19 max values
hg19Z[hg19Z < 3 & hg19Z > -3] <- NA

hg19Z[hg19Z > 15 ] <- 15
hg19Z[hg19Z < -15 ] <- -15

# set mm10 max values
mm10Z[mm10Z < 3 & mm10Z > -3] <- NA

mm10Z[mm10Z > 15 ] <- 15
mm10Z[mm10Z < -15 ] <- -15

image(matrix(1, nrow = nrow(hg19Z), ncol = ncol(hg19Z)), col = "grey90", xaxt = "n", yaxt = "n", main = "DNA gain")
image(rbind(hg19Z[1,],mm10Z[1,],hg19Z[3,],mm10Z[3,]), 
      col = pal(40), zlim = c(-max(abs(zAll)), max(abs(zAll))), add = TRUE,
      xaxt = "n", yaxt = "n", main = "hg19")
axis(side = 1, at = seq(0,1,length.out = 4), labels = c("hg19", "mm10", "hg19", "mm10"), las = 2)
axis(side = 1, at = seq(0,1,length.out = 4), labels = c("hg19", "hg19","mm10" ,"mm10"), 
     las = 2, line = 3, tick = FALSE)
axis(side = 2, at = seq(0,1,length.out = ncol(hg19Z)), labels = NA, las = 1)
mtext(side = 2, at = seq(0,1,length.out = ncol(hg19Z)), text = colName, las = 1, col = GenomeCols[colourNum], cex = .7, line = 1)



# plot the mm10 matrix




image(matrix(1, nrow = nrow(hg19Z), ncol = ncol(hg19Z)), col = "grey90", xaxt = "n", yaxt = "n", main = "DNA loss")
image(rbind(hg19Z[2,],mm10Z[2,],hg19Z[4,],mm10Z[4,]),
      col = pal(40), zlim = c(-max(abs(zAll)), max(abs(zAll))), add = TRUE,
      xaxt = "n", yaxt = "n", main = "mm10")
axis(side = 1, at = seq(0,1,length.out = 4), labels = c("hg19", "mm10", "hg19", "mm10"), las = 2)
axis(side = 1, at = seq(0,1,length.out = 4), labels = c("hg19", "hg19","mm10" ,"mm10"), 
     las = 2, line = 3, tick = FALSE)


# plot colour legend
plot.new()
axis(side = 1, at = .5, labels = "Background", las = 1, line= 0.5, tick = FALSE)
axis(side = 1, at = .5, labels = "Genome", las = 1, line = 3.5, tick = FALSE)


par(new = TRUE, mar = c(30,3,5,3))

mat <- matrix(seq(-max(abs(zAll)), max(abs(zAll)), length.out = 19), nrow = 1)
image(mat, col = "grey90", xaxt = "n", yaxt = "n", main = "Z score")
mat[mat < 3 & mat > -3] <- NA
image(mat, col = pal(40), xaxt = "n", yaxt = "n", main = "Z score", add = TRUE)

# colour legend axis
# identify location of NA elements
upperCut <- min(which(as.integer(seq(-max(abs(zAll)), max(abs(zAll)), length.out = 1000)) == 3))/1000
axis(side = 2, at = c(0, 1-upperCut,upperCut, 1), labels = c(" < -15","-3","3"," > 15"), las = 2)

# plot feature legend
par(mar = c(0,1,0,0))
plot.new()
legend("topleft", title = "Feature indicators",
       legend = c("Gene-rich/Active", "Gene-poor/Silenced", "Ancestral DNA marker", "Source of DNA gain", "Genome instability"), 
       fill = GenomeCols , bty = "n")

plot.new()
plot.new()

dev.off()






