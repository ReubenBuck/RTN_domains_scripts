

rm(list = ls())


specRef <- "mm10"
specQue <- "hg19"

load(paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/",specRef,".stretch.RData", sep = ""))
stretchedRefAnc.gr <- stretchedRef.gr

load(paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/",specRef,".stretch.rep.RData", sep = ""))
stretchedRefRep.gr <- stretchedRef.gr

rm(stretchedRef.gr)

gapType = c("queIns", "refDel")
gapNames = c(paste(specQue,"gain"), paste(specRef, 'loss'))
names(gapNames) <- gapType

for(g in gapType){
  w <- width(stretchedRefAnc.gr[stretchedRefAnc.gr$type == g])
  assign(paste(g, "Anc", sep = ""), w)
  
  w <- width(stretchedRefRep.gr[stretchedRefRep.gr$type == g])
  assign(paste(g, "Rep", sep = ""), w)
  
  print(g)
}



for(g in gapType){
  pdf(file = paste("~/Desktop/RTN_domains/RTN_domain_plots/netGainLoss/gapAssigningDist/",
                   gsub(" ","",gapNames[g]), ".pdf", sep = "")
      )
  hist(log10(get(paste(g,"Rep",sep = ""))),
       breaks = seq(0,5.5,.05), 
       col = scales::alpha("black", .3),
       border = scales::alpha("black", .3),
       xlim = c(0,4), xaxt = "n",
       xlab = "Length (bp)",
       main = gapNames[g])
  hist(log10(get(paste(g,"Anc",sep = ""))), 
       breaks = seq(0,5.5,.05), add = TRUE, 
       col = scales::alpha("red", .3),
       border = scales::alpha("red", .3),
       xlim = c(0,4))
  
  aty <- axTicks(1)
  labels <- sapply(aty,function(i)
    as.expression(bquote(10^ .(i)))
  )
  axis(side = 1, at = aty,
       labels = labels)
  if(specRef == "hg19" & g == "refDel"){
  legend("topright", legend = c("Recent transposon","Ancestral element"), 
         fill = c(scales::alpha("black", .3), scales::alpha("red", .3)), bty = "n")
  }
  dev.off()
}




bitsAnc <- stretchedRefAnc.gr[stretchedRefAnc.gr$type == "queIns" | stretchedRefAnc.gr$type == "refDel"]
bitsRep <- stretchedRefRep.gr[stretchedRefRep.gr$type == "queIns" | stretchedRefRep.gr$type == "refDel"]

# the gaps stay 


ll <- list(a = expression(x^2 - 2*x + 1), b = as.name("Jim"),
           c = as.expression(exp(1)), d = call("sin", pi))
sapply(ll, typeof)


hist(log10(table(mcols(bitsAnc)$elementID )))

hist(log10(table(mcols(bitsRep)$elementID )), add = TRUE, col = 2, density = 0)







load("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/hg19.synthBin.rep.RData")
hg19synthBinRep.gr <- synthBin.gr

load("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/mm10.synthBin.rep.RData")
mm10synthBinRep.gr <- synthBin.gr

load("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/hg19synthBinNorm.RData")
hg19synthBinAnc.gr <- synthBinNorm.gr

load("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/mm10synthBinNorm.RData")
mm10synthBinAnc.gr <- synthBinNorm.gr


hg19synthBinRep.gr <- hg19synthBinRep.gr[as.integer(names(hg19synthBinAnc.gr))]
hg19dfRep <- data.frame(mcols(hg19synthBinRep.gr))[,c("refIns","refDel","queIns","queDel")] 
hg19dfRep <- hg19dfRep * 200000 / (200000 - hg19synthBinRep.gr$missingGap - hg19synthBinRep.gr$seqGap)

hg19dfAnc <- data.frame(mcols(hg19synthBinAnc.gr))[,c("refIns","refDel","queIns","queDel")] 


mm10synthBinRep.gr <- mm10synthBinRep.gr[as.integer(names(mm10synthBinAnc.gr))]
mm10dfRep <- data.frame(mcols(mm10synthBinRep.gr))[,c("refIns","refDel","queIns","queDel")] 
mm10dfRep <- mm10dfRep * 200000 / (200000 - mm10synthBinRep.gr$missingGap - mm10synthBinRep.gr$seqGap)

mm10dfAnc <- data.frame(mcols(mm10synthBinAnc.gr))[,c("refIns","refDel","queIns","queDel")] 




layout(matrix(1:4,nrow = 2))
hist(dfRep$refDel)
smoothScatter(dfRep$refDel, dfAnc$refDel)
plot.new()
hist(dfAnc$refDel)



# we can put all four on here
# maybe just do it for the one genome? and make 2 plots
# remeber to use the col brewer

library(MASS)




datHg19Del <- data.frame(rep = hg19dfRep$refDel, anc = hg19dfAnc$refDel)
datMm10Del <- data.frame(rep = mm10dfRep$refDel, anc = mm10dfAnc$refDel)
datHg19Ins <- data.frame(rep = hg19dfRep$refIns, anc = hg19dfAnc$refIns)
datMm10Ins <- data.frame(rep = mm10dfRep$refIns, anc = mm10dfAnc$refIns)

densHg19Del <- kde2d(datHg19Del$rep, datHg19Del$anc)
densMm10Del <- kde2d(datMm10Del$rep, datMm10Del$anc)
densHg19Ins <- kde2d(datHg19Ins$rep, datHg19Ins$anc)
densMm10Ins <- kde2d(datMm10Ins$rep, datMm10Ins$anc)

hg19DelRep <- hist(datHg19Del$rep, breaks = 100, plot = FALSE)
hg19InsRep <- hist(datHg19Ins$rep, breaks = 100, plot = FALSE)
Mm10DelRep <- hist(datMm10Del$rep, breaks = 100, plot = FALSE)
Mm10InsRep <- hist(datMm10Ins$rep, breaks = 100, plot = FALSE)

hg19DelAnc <- hist(datHg19Del$anc, breaks = 100, plot = FALSE)
hg19InsAnc <- hist(datHg19Ins$anc, breaks = 100, plot = FALSE)
Mm10DelAnc <- hist(datMm10Del$anc, breaks = 100, plot = FALSE)
Mm10InsAnc <- hist(datMm10Ins$anc, breaks = 100, plot = FALSE)



pal <- RColorBrewer::brewer.pal(4,"Dark2")

lims = c(0,90e3)


pdf(file = "~/Desktop/RTN_domains/RTN_domain_plots/netGainLoss/gapAssigningDist/jointDist.pdf")
layout(matrix(1:4,nrow = 2), heights = c(1,2), widths = c(2,1))

par(mar = c(0,5,5,0))
plot(hg19DelRep$mids, hg19DelRep$density, xlim = lims, axes = FALSE, col = 2, type = "n", ylab = "")
polygon(Mm10DelRep$mids, Mm10DelRep$density, col = scales::alpha(pal[4],.2), border = scales::alpha(pal[4],1))
polygon(hg19DelRep$mids, hg19DelRep$density, col = scales::alpha(pal[2],.2), border = scales::alpha(pal[2],1))
polygon(hg19InsRep$mids, hg19InsRep$density, col = scales::alpha(pal[1],.2), border = scales::alpha(pal[1],1))
polygon(Mm10InsRep$mids, Mm10InsRep$density, col = scales::alpha(pal[3],.2), border = scales::alpha(pal[3],1))

par(mar = c(5,5,0,0))
contour(densHg19Ins, add=FALSE, col = scales::alpha(pal[1],.5), 
        xlim = lims, ylim = lims,drawlabels = FALSE, lwd = 2,
        ylab = "Ancestral elements (bp)", xlab = "Recent transposons (bp)")
grid()
contour(densHg19Del, add = TRUE, col = scales::alpha(pal[2],.5),drawlabels = FALSE, lwd = 2)
contour(densMm10Ins, add = TRUE, col = scales::alpha(pal[3],.5),drawlabels = FALSE, lwd = 2)
contour(densMm10Del, add = TRUE, col = scales::alpha(pal[4],.5),drawlabels = FALSE, lwd = 2)

par(mar=c(1.5,1.5,2,2))
plot.new()
legend("bottomleft", legend = c("hg19 gain", "hg19 loss", "mm10 gain", "mm10 loss"),
      fill = c(scales::alpha(pal[1],.2),
               scales::alpha(pal[2],.2),
               scales::alpha(pal[3],.2),
               scales::alpha(pal[4],.2)),
      border = c(scales::alpha(pal[1],1),
                 scales::alpha(pal[2],1),
                 scales::alpha(pal[3],1),
                 scales::alpha(pal[4],1)),
      bty = "n")

par(mar = c(5,0,0,5))
plot(hg19DelAnc$density, hg19DelAnc$mids, type = "n", ylim = lims, axes = FALSE, col = 2, xlab = "")
polygon(Mm10DelAnc$density, Mm10DelAnc$mids, col = scales::alpha(pal[4],.2), border = scales::alpha(pal[4],1))
polygon(hg19InsAnc$density, hg19InsAnc$mids, col = scales::alpha(pal[1],.2), border = scales::alpha(pal[1],1))
polygon(c(0,hg19DelAnc$density), c(0,hg19DelAnc$mids), col = scales::alpha(pal[2],.2), border = scales::alpha(pal[2],1))
polygon(Mm10InsAnc$density, Mm10InsAnc$mids, col = scales::alpha(pal[3],.2), border = scales::alpha(pal[3],1))

dev.off()
# actually have a way of plotting it now to compare each speceis and the profiles





VennDiagramHg19 <- data.frame(repUniq = c(318, 1406, 1182, 1181), 
                          intersect = c(1098, 314, 2050, 1980), 
                          ancUniq = c(795, 1791, 687, 811))

VennDiagramMm10 <- data.frame(repUniq = c(380, 1460, 1271, 1110), 
                              intersect = c(1018, 350, 2057, 1956), 
                              ancUniq = c(723, 1756, 653, 858))


rownames(VennDiagramHg19) <- rownames(VennDiagramMm10) <- c("hg19 gain", "hg19 loss", "mm10 gain", "mm10 loss")

pdf(file = "Desktop/RTN_domains/RTN_domain_plots/netGainLoss/gapAssigningDist/hg19venn.pdf", width = 6.2)
layout(1)
par(mar = c(5,5,5,5))

plot.new()
points(c(.335, .665), c(.5,.5), cex = 38)
points(c(.335, .665), c(.5,.5), cex = 38, pch = 16, 
       col = c(scales::alpha(1, .1), scales::alpha(2, .1)))

for(i in 1:3){
text(x = rep(c(.18,.5,.82)[i],4),
     y = seq(.65,.35,length.out = 4),
     labels = VennDiagramHg19[,i]*200000/1e6, 
     cex = sqrt((VennDiagramHg19[,i]/max(as.matrix(VennDiagramHg19)))*5 ),
     col = pal)
}
text(c(.18, .82), c(.87,.87), c("Recent tansposon", "Ancestral element"))
title(main = "hg19")
dev.off()


pdf(file = "Desktop/RTN_domains/RTN_domain_plots/netGainLoss/gapAssigningDist/mm10venn.pdf", width = 6.2)
layout(1)
par(mar = c(5,5,5,5))

plot.new()
points(c(.335, .665), c(.5,.5), cex = 38)
points(c(.335, .665), c(.5,.5), cex = 38, pch = 16, 
       col = c(scales::alpha(1, .1), scales::alpha(2, .1)))

for(i in 1:3){
  text(x = rep(c(.18,.5,.82)[i],4),
       y = seq(.65,.35,length.out = 4),
       labels = VennDiagramMm10[,i]*200000/1e6, 
       cex = sqrt((VennDiagramMm10[,i]/max(as.matrix(VennDiagramMm10)))*5 ),
       col = pal)
}
text(c(.18, .82), c(.87,.87), c("Recent tansposon", "Ancestral element"))
title(main = "mm10")

dev.off()









