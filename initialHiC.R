A <- data.frame(col1 = c("A","A","A","B", "B","B", "C", "C", "C"), col2 = c("A","B", "C","A","B", "C", "A","B", "C"), values = 1:9,stringsAsFactors = FALSE)

mat <- matrix(0, nrow =  3, ncol = 3,dimnames = list(c("A","B","C"), c("A", "B","C")))
mat[as.character(A$col1),as.character(A$col2)] = A$values 



# if we build the function this means that a new matrix we be built each time

matrixBuild <- function(df, matName){
  colnames(df) <- c("col1", "col2", "value")
  assign(matName,)
}


as.matrix.data.frame(A$values, rownames.force = A$col1)


mat1 <- matrix(NA)

as.table(A)
mat

hic <- read.table(file = "~/Downloads/GSE65126_HiC_dog_liver_merged_50000.txt", nrows = 200000, header = T)
hic <- hic[hic$chrom1 == 1,]
hic$OoverE <- hic$observed_count/hic$expected_count
hicIDs <- data.frame(ID1 = as.integer(as.factor(hic$start1[order(hic$start1)])),
                     ID2 = as.integer(as.factor(hic$start2[order(hic$start1)])),
                     value = hic$OoverE
                     )

hicMat <- matrix(0, nrow = max(c(hicIDs$ID1,hicIDs$ID2)), ncol =  max(c(hicIDs$ID1,hicIDs$ID2)))
for(i in 1:nrow(hicIDs)){
  hicMat[hicIDs$ID1[i], hicIDs$ID2[i]] <- hicIDs$value[i]
}

image(log10(hicMat[1000:1200,1000:1200]))

hicMat[is.na(hicMat)] <- 0
hicEig <- eigen(log10(hicMat+1))

plot(1:ncol(hicEig$vectors),smooth(as.numeric(hicEig$vectors[1,])), type = "l")


pca <- prcomp(log10(hicMat + 1))
plot(pca$x[,1])


layout(matrix(c(1,2), nrow = 2))
image(log10(hicMat[1000:1500,1000:1500]))
plot(pca$x[,2], type = "l")
par(new = TRUE)
plot(1:nrow(hicEig$vectors[1000:1500,]),(as.numeric(hicEig$vectors[1000:1500,2])), type = "l", col = 2)








