
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
