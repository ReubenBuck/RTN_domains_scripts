A <- read.table("~/Desktop/RTN_domains/data/chainAlignments/test/longAlign.txt")

plot(A$V2, type = "l")


# essentially we can probably mask over some of these regions

# given taht most gaps are pretty small it's worth arguing that they can be part of the smae region



plot(density(log10(A$V2)), col =2)
lines(density(log10(A$V3)), col = 3)
lines(density(log10(A$V1)))

# we want to know how many intervals we need. 
# if we include every gap we have to create approximatly 100 regions, 
# set of suitable aligning regions is what we are aiming for here. 
# as we join these the distribution is likley to change 

# when both species are below the 25% level, we can merge the two regions and perform this iterativly. 
quantile(A$V1)[2]

l <- A[A$V2 < 13 & A$V3 < 13,]
# at the end we will end up with gaps that are big 


# chain 672911181 chr14 107349540 + 25157171 106332175 chr12 121257530 + 45737574 114669059 1
# 

# we just read our data in and find out our min gap size at level a
# begin merging alignmnets 




alignHead <- read.table(file = "~/Desktop/RTN_domains/data/chainAlignments/test/hg19.mm9.chainHeader.txt", nrows = 1000)


len = alignHead[,7] - alignHead[,6]

hist(len)
plot(ecdf(len))

plot(len, len/(2.8*10^9))


plot(sort(len/(2.8*10^9),decreasing = TRUE))


len2 <- sort(len, decreasing = TRUE)
summer <- NULL
for(i in 1:1000){
  summer <- c(summer, sum(as.numeric(len2[1:i])))
}

par(mar = c(7,7,7,7))
plot(log10(len2),(summer), xlim=c(log10(len2[1]), 3), type = "l",
     xlab = "minimum human chain length (log10 bp)", ylab = "total human chain length (bp)")

