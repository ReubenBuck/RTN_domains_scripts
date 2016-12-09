setwd("~/Desktop/RTN_domains/")
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(expm)

rawMatrix <- read.delim("data/trialHicAnalysis/rawMatrix/rawMatrix.txt", header = TRUE, sep = "\t",row.names = 1)
rawMatrix <- as.matrix(rawMatrix[,2:ncol(rawMatrix)])

layout(c(1,2), heights = c(2,1))
par(mar = c(0,5,5,5))
image(log10(rawMatrix), xaxt = "n")
par(mar = c(5,5,0,5))
barplot(colSums(rawMatrix, na.rm = TRUE), space = 0)



OurTheme = theme(axis.line=element_blank(),
                 axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 legend.position="none",
                 panel.background=element_blank(),
                 panel.border=element_blank(),
                 panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 plot.background=element_blank())


data = reshape2::melt(rawMatrix)
g1 <- ggplot(data = data, aes(x = Var1, y = Var2,z = log10(value)))
g1 <- g1 + geom_raster(aes(fill = log10(value))) + guides(fill=FALSE)
g1 <- g1 + OurTheme + scale_fill_distiller(palette = "Spectral")

g2 <- ggplot(data = data, aes(Var1))
g2 <- g2 + geom_bar(aes(weight = value)) 
g2 <- g2 + OurTheme

grid.arrange(g1,g2,heights = c(3,1))


plot(apply(X = (rawMatrix),MARGIN = 2,FUN = sum, na.rm = TRUE))
filteredMatrix <- rawMatrix
filteredMatrix[colSums(rawMatrix) < 50000,] <- NA
filteredMatrix[,colSums(rawMatrix) < 50000] <- NA
layout(1)
par(mar = c(5,5,5,5))
image(log10(filteredMatrix))


testMat <- matrix(1:13,nrow=10,ncol = 10)
balancedMatrix <- balance(testMat, job = "B")

image(testMat)
image(scale(testMat,scale = TRUE,center = FALSE))

colSums(testMat)
colSums(scale(testMat))

library(HiTC)
library(Matrix)



GR <- GRanges(seqnames= Rle(rep("chr1",10)),
              ranges = IRanges(start = seq(1,by = 10,length.out = 10),width = rep(10,10)),
              names = 1:10)
htMatrix <- Matrix(data = as.matrix(dist(1:10)), dimnames = list(1:10,1:10))

HTob <- HTCexp(intdata = htMatrix,ygi = GR, xgi = GR)
isSymmetric(HTob)

normMat <- normICE(HTob)


plot(normMat)


##### lets do it to our matrix

GR <- GRanges(seqnames = Rle(gsub("\\..+", "",colnames(rawMatrix)))[1:23],
              ranges = IRanges(start = as.integer(gsub(".+\\.", "",colnames(rawMatrix)))[1:23]
                               , width = 10e6),
              names = colnames(rawMatrix)[1:23]
              )
GR <- GR[seqnames(GR) == "chr1"]


htMatrix <- Matrix(data = filteredMatrix[1:length(GR),1:length(GR)], 
                   dimnames = list(colnames(rawMatrix)[1:length(GR)], 
                                   colnames(rawMatrix)[1:length(GR)]
                                   )
                   )

HTob <- HTCexp(intdata = htMatrix,ygi = GR, xgi = GR)
isSymmetric(HTob)
normMat <- normICE(HTob)
layout(c(1,2))
par(mar=c(5,5,5,5))
image(as.matrix(log10(intdata(normMat))))
image(as.matrix(log10(htMatrix)))


hist(colSums(intdata(normMat)))
hist((colSums((htMatrix))))
