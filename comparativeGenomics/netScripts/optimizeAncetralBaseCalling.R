

hg19.fills <- read.table("Desktop/ancDistribution/hg19/hg19.ancFillsDistribution")
hg19.newRep <- read.table("Desktop/ancDistribution/hg19/hg19.ancNewRepDistribution")
sum(hg19.newRep$V2)

mm10.fills <- read.table("Desktop/ancDistribution/mm10/mm10.ancFillsDistribution")
mm10.newRep <- read.table("Desktop/ancDistribution/mm10/mm10.ancNewRepDistribution")


# select a threshold at which we reduce the impact of false positives and false negatives

plot( (cumsum(mm10.fills$V2)/sum(mm10.fills$V2) * 100), type = "l", 
     main = "mm10", xlab = "Coverage level", ylab = "Overlapping ancestral bases (cumulative %)")
grid()
lines( (cumsum(mm10.newRep$V2)/sum(mm10.newRep$V2)  * 100 ), col = 2)
legend("bottomright", col = c(1,2), legend = c("Fills", "Recent repeats"), lty = 1, bty = "n")

abline(v = 4)
abline(v = 5)

plot( (cumsum(hg19.fills$V2)/sum(hg19.fills$V2) * 100), type = "l", 
      main = "hg19", xlab = "Coverage level", ylab = "Overlapping ancestral bases (cumulative %)")
grid()
lines( (cumsum(hg19.newRep$V2)/sum(hg19.newRep$V2)  * 100 ), col = 2)
legend("bottomright", col = c(1,2), legend = c("Fills", "Recent repeats"), lty = 1, bty = "n")


# this is the trade off, so our best trade off is around 5

plot( (cumsum(hg19.newRep$V2)/sum(hg19.newRep$V2)) - (cumsum(hg19.fills$V2)/sum(hg19.fills$V2)) )

lines( (cumsum(mm10.newRep$V2)/sum(mm10.newRep$V2)) - (cumsum(mm10.fills$V2)/sum(mm10.fills$V2)) )

plot(mm10.fills$V2/sum(mm10.fills$V2), ylim = c(0,1))
lines(mm10.newRep$V2/sum(mm10.newRep$V2))

abline(v = 4)

abline(v = 6)

# all the falsely identified bases

layout(c(1,2))
plot(  (cumsum(hg19.fills$V2)) , cumsum(hg19.newRep$V2), type = "n" ,
       main = "hg19")
text((cumsum(hg19.fills$V2)) , cumsum(hg19.newRep$V2), hg19.fills$V1)

plot(  (cumsum(mm10.fills$V2)) , cumsum(mm10.newRep$V2), type = "n" ,
       main = "mm10")
text((cumsum(mm10.fills$V2)) , cumsum(mm10.newRep$V2), hg19.fills$V1, col =2)


plot(  (cumsum(hg19.fills$V2)/sum(hg19.fills$V2)) , cumsum(hg19.newRep$V2)/sum(hg19.newRep$V2) ,
       xlim = c(0,1), ylim = c(0,1), xaxs = "i", yaxs = "i", type = "b", pch = 16, cex = .5,
       ylab = "Ancestral bases overlaping recent repeats (cummulative proportion)", 
       xlab = "Ancestral bases overlaping net fills (cummulative proportion)",
       main = "Minimum ancestral coverage threshold")
lines((cumsum(mm10.fills$V2)/sum(mm10.fills$V2)) , 
      cumsum(mm10.newRep$V2)/sum(mm10.newRep$V2), type = "b", pch = 16, col = 2, cex = .5)
text((cumsum(hg19.fills$V2)/sum(hg19.fills$V2)) , 
     cumsum(hg19.newRep$V2)/sum(hg19.newRep$V2), hg19.fills$V1, pos = 1)
text((cumsum(mm10.fills$V2)/sum(mm10.fills$V2)) , 
     cumsum(mm10.newRep$V2)/sum(mm10.newRep$V2), mm10.fills$V1, col = 2, pos = 1)
abline(a=0,b=1, lty = 2)
grid()
legend("bottomright", c("hg19 coverage", "mm10 coverage"), col = c(1,2), lty = 1)

# the percentage of false hits vs the percentage of true hits


# type 1 error  
# type 2 error

pdf(file = "Desktop/RTN_domains/RTN_domain_plots/netGainLoss/gapAssigningDist/optimiseOutGroupCov.pdf")
y =  (cumsum(hg19.newRep$V2)/sum(hg19.newRep$V2))
x =  (cumsum(hg19.fills$V2)/sum(hg19.fills$V2))
plot( x, y ,type = "n",
      ylab  = "Type 1 errors (cumulative proportion)",
      xlab  = "Type 2 errors (proportional increase)",
      ylim = c(0,1),
      xlim = c(0,1))

text(x,y, cex = c(rep(.7,5), 1.3, rep(.7, 32 - 6))  )

y =  (cumsum(mm10.newRep$V2)/sum(mm10.newRep$V2))
x =  (cumsum(mm10.fills$V2)/sum(mm10.fills$V2))
text(x,y, cex = c(rep(.7,3), 1.3, rep(.7, 38 - 6)) , col = 2 )

legend("bottomright", legend = c("hg19","mm10"), fill = c(1,2), bty = "n")

dev.off()
# looking for the best trade of in type 2 errors 


# percentage of type 1 errors

# increase in type 2 error











plot(  (cumsum(mm10.fills$V2)) , cumsum(mm10.newRep$V2), type = "n" ,
       main = "mm10", xlim = c(0, 200e6))
text((cumsum(mm10.fills$V2)) , cumsum(mm10.newRep$V2), hg19.fills$V1)
