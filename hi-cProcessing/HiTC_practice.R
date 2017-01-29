require(HiTC)

require(HiCDataHumanIMR90)
data(Dixon2012_IMR90)
hox <- extractRegion(hic_imr90_40$chr6chr6, chr="chr6", from=30e6, to=58e6)
di<- directionalityIndex(hox)

layout(c(1,2))
mapC(hox)
barplot(di, col=ifelse(di>0,"darkred","darkgreen"))

# might be worth converting our results into Hitc pkg


