
# intron and intergenic is what is left over

# do we include gene number too 





df <- data.frame(synthBin.gr)
df$chrType = "autosome"
df[df$seqnames == "chrX","chrType"] <- "sexChr"

df[df$seqGap + df$missingGap > 20000,] <- NA
df <- df[complete.cases(df),]

df$gcContent <- df$gcContent * 100

modRefIns <- (lm(data = df, (refIns) ~ gcContent + dnasePeaks  + exon + chrType + distFromTelomere))
modRefDel <- (lm(data = df, (refDel) ~ gcContent + dnasePeaks  + exon  + chrType + distFromTelomere ))
modQueIns <- (lm(data = df, (queIns) ~ gcContent + dnasePeaks  + exon  + chrType + distFromTelomere ))
modQueDel <- (lm(data = df, (queDel) ~ gcContent + dnasePeaks  + exon  + chrType + distFromTelomere ))

summary(modRefIns)
summary(modRefDel)
summary(modQueIns)
summary(modQueDel)


modRefTurnover <- lm(data = df, (refIns + refDel) ~ gcContent + dnasePeaks + dnaseActivity + exon + intron + chrType + distFromTelomere)
modQueTurnover <- lm(data = df, (queIns + queDel) ~ gcContent + dnasePeaks + dnaseActivity + exon + intron + chrType+ distFromTelomere)
modRefNetGain <- lm(data = df, (refIns - refDel) ~ gcContent + dnasePeaks + dnaseActivity + exon + intron + chrType+ distFromTelomere)
modQueNetGain <- lm(data = df, (queIns - queDel) ~ gcContent + dnasePeaks + dnaseActivity + exon + intron + chrType+ distFromTelomere)

summary(modRefNetGain)
summary(modRefTurnover)
summary(modQueNetGain)
summary(modQueTurnover)

layout(matrix(1:4, nrow = 2))
plot(modQueDel)

h <- hist(residuals(modQueDel), breaks = 100)
plot(h$mids, sort(rnorm(h$mids)))



mod <- lm(data = df, queDel ~ gcContent)


s <- summary(mod)
layout(1)
plot(mod)

plot(mod$residuals)


# telomere distance and centromere distance are also important 


# relative contribution to variability explained

# distace to centromere and distance to telomere.

mod <- modRefNetGain
d <- drop1(mod)


for( i in 1:length(d$RSS)){
  print(d$`Sum of Sq`[i]/sum(anova(mod)[2][1:(length(coefficients(mod))-1),]))
}






modRefIns2 <- lm(data = df , (refIns) ~  gcContent + exon + chrType)


sum(anova(modRefIns)[2][1:4,]) - sum(anova(modRefIns2)[2][1:3,])

(summary(modRefIns)$r.squared - summary(modRefIns2)$r.squared)/summary(modRefIns)$r.squared


# not sure if im doing this right might play witg it at home

# DNase1 can explain 87% of the R squred value

# explained mostly by Dnase1 activity 


library(spdep)
library(igraph)



df <- data.frame(synthBin.gr)
df$chrType = "autosome"
df[df$seqnames == "chrX","chrType"] <- "sexChr"

df[df$seqGap + df$missingGap > 20000,] <- NA
df <- df[complete.cases(df),]

df$gcContent <- df$gcContent * 100

ol<-findOverlaps(synthBin.gr, maxgap = 30*width(synthBin.gr)[1])
ol <- ol[!(isRedundantHit(ol))]
# remove NA hits
ol <- ol[!is.na(df$refIns[queryHits(ol)])]
ol <- ol[!is.na(df$refIns[subjectHits(ol)])]
olMat <- data.frame(ol)
G <- graph.data.frame(d = olMat,directed=FALSE)
weight <- olMat$subjectHits - olMat$queryHits
weight <- -(weight - (max(weight) + 1)) / max(weight)
weight[isSelfHit(ol)] <- 0
E(G)$weight <- weight
A <- as_adjacency_matrix(G,type="both",names=FALSE,sparse=TRUE,edges = FALSE, attr = "weight")




dfChr <- df[df$seqnames == "chr2",]

A <- A[df$seqnames == "chr2"  ,df$seqnames == "chr2" ]
wMat <- mat2listw(A)


# con




spatMod <- errorsarlm(data = dfChr, refDel ~ gcContent + dnasePeaks + dnaseActivity + exon + intron, listw = wMat)

summary(spatMod)


qqnorm(residuals(spatMod))

acf(residuals(spatMod),lag.max = 40)

plot(fitted(spatMod),residuals(spatMod))


ESS <- sum((fitted(spatMod) - mean((dfChr$refDel)))^2)
RSS <- sum(residuals(spatMod)^2)

1 - (RSS / (RSS + ESS))

bptest.sarlm(spatMod)


spatModAuto <- spautolm(data = dfChr, queDel ~ gcContent + dnasePeaks + dnaseActivity + exon + intron, listw = wMat)

ESS <- sum((fitted(spatModAuto) - mean((dfChr$refDel)))^2)
RSS <- sum(residuals(spatModAuto)^2)

1 - (RSS / (RSS + ESS))

qqnorm(residuals(spatModAuto))
plot(residuals(spatModAuto))
bptest(spatModAuto)



# so now we have a modeling appraoch that might actually work


lm.morantest(spatMod, wMat)

moran.test(x = residuals(spatMod), wMat)
moran.test(x = residuals(modChr1), wMat)


# so now there is no autocorrelation in our residuals
# therefore we could succesfully remove the error

# next to get the R squered to determine our variance explained by the model.



modChr1 <- lm(data = dfChr, (queIns - queDel) ~ gcContent + dnasePeaks + dnaseActivity + exon + intron)

spCor <- sp.correlogram(neighbours = wMat$neighbours, var = dfChr$refIns, order = 10)

plot.spcor(spCor)


acf(residuals(modChr1),lag.max = 100)



h <- hist(scale(residuals(spatMod)), breaks = 1000)
plot( sort(rnorm(h$mids)), h$mids, xlim = c(-4,4), ylim = c(-4,4))


ESS <- sum((fitted(modChr1) - mean((df$queIns - df$queDel)))^2)
RSS <- sum((fitted(modChr1) - mean((df$queIns - df$queDel)))^2)

ts <- as.ts(modChr1)

acf(residuals(spatMod))
acf(residuals(modChr1))

# it is removing the auto correlation

# so we have a significant Morans I 
# and a significant autocorrelation
# when we apply the new spatial model we can factor that out.
lm.morantest(modQueIns, listw = wMat)

