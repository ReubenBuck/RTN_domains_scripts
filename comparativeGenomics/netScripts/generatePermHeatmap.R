
rm(list = ls())


load("~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/hg19.permZ.RData")
hg19Z <- zAll

load("~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/mm10.permZkeepX.RData")
mm10Z <- zAll
load("~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/mm10.permZnoX.RData")
mm10ZnoX <- zAll

# replace recombination data in mouse because of missing X chromosome
mm10Z[,c("recHotSpot", "recombRate")] <- mm10ZnoX[,c("recHotSpot", "recombRate")]

rm(mm10ZnoX)

rownames(mm10Z) <- c("mm10.mm10.gain", "mm10.mm10.loss", "mm10.hg19.gain", "mm10.hg19.loss")
rownames(hg19Z) <- c("hg19.hg19.gain", "hg19.hg19.loss", "hg19.mm10.gain", "hg19.mm10.loss")


# generate col order based on clustering
hg19Z <- hg19Z[,!(colnames(hg19Z) == "dnaseActivity" | colnames(hg19Z) == "L1Motif" | colnames(hg19Z) == "prdm9Motif" | colnames(hg19Z) == "ctcfMotif")]
mm10Z <- mm10Z[,!(colnames(mm10Z) == "dnaseActivity" | colnames(mm10Z) == "L1Motif" | colnames(mm10Z) == "prdm9Motif" | colnames(mm10Z) == "ctcfMotif")]
colName <- c("DNase1 HS peaks", "LADs","Exons", "Introns", 
             "Ancient elements", "Lineage-specifc L1s", "Lineage-specifc SINEs",
             "Ancestral L1s","Recombination hotspots", "CpG islands",
             "GC content","Recombination rate")

zAll <- rbind(hg19Z,mm10Z)

# sorting
hClust <- hclust(dist(t(zAll)))
mm10Z <- mm10Z[,hClust$order]
hg19Z <- hg19Z[,hClust$order]

colName <- colName[hClust$order]


# get proper col names






# order mouse similar to human
mm10Z <- mm10Z[c("mm10.hg19.gain", "mm10.hg19.loss","mm10.mm10.gain", "mm10.mm10.loss"),]

pdf(file = "~/Desktop/RTN_domains/RTN_domain_plots/netGainLoss/hotspotPlot/heatmap.pdf")
layout(matrix(1:4, nrow = 1), width = c(8,8,8,5))
pal <- colorRampPalette(c("blue", "white", "red"))

par(mar = c(10,1,5,1))


#par(mar = c(10,10,5,1))
image(hg19Z, col = "white", zlim = c(-max(abs(zAll)), max(abs(zAll))),axes = FALSE)
#axis(side = 2, at = seq(0,1,length.out = ncol(hg19Z)), labels = colName, las = 2, line = -11)

hg19Z[hg19Z < 2 & hg19Z > -2] <- NA
image(hg19Z, col = pal(40), zlim = c(-max(abs(zAll)), max(abs(zAll))),xaxt = "n", yaxt = "n", main = "hg19")
axis(side = 1, at = seq(0,1,length.out = 4), labels = c("hg19 gain", "hg19 loss", "mm10 gain", "mm10 loss"), las = 2)
axis(side = 2, at = seq(0,1,length.out = ncol(hg19Z)), labels = colName, las = 2)



#par(mar = c(10,1,5,10))
mm10Z[mm10Z < 2 & mm10Z > -2] <- NA
image(mm10Z, col = pal(40), zlim = c(-max(abs(zAll)), max(abs(zAll))),xaxt = "n", yaxt = "n", main = "mm10")
axis(side = 1, at = seq(0,1,length.out = 4), labels = c("hg19 gain", "hg19 loss", "mm10 gain", "mm10 loss"), las = 2)

par(mar = c(35,3,5,3))

mat <- matrix(seq(-max(abs(zAll)), max(abs(zAll)), length.out = 19), nrow = 1)
mat[mat < 2 & mat > -2] <- NA
image(mat, col = pal(40), xaxt = "n", yaxt = "n", main = "Z score")

upperCut <- min(which(as.integer(seq(-max(abs(zAll)), max(abs(zAll)), length.out = 1000)) == 2))/1000

axis(side = 2, at = c(0, 1-upperCut,upperCut, 1), labels = c(-13,-2,2,13), las = 2)

dev.off()









