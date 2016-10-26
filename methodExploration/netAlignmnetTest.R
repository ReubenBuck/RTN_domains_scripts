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

