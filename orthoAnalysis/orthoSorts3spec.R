### Reciprical best clique

rm(list = ls())

setwd("Desktop/RTN_domains/")

options(stringsAsFactors = FALSE)

# potential user inputs
snames = c("hg19", "mm9", "canFam3")
ensGids = c("ENST","ENSMUST" ,"ENSCAFT")

names(snames) = ensGids
  
names(ensGids) = snames

colAsignment <- function(datObject, snames, ensGids){
  datCols <- snames[gsub(pattern = "[0-9].+",replacement = "",x = datObject[1,])]
  datObject <- rbind(datObject[grep(pattern = ensGids[datCols[1]], x = datObject[,1]),1:2], 
                     datObject[-grep(pattern = ensGids[datCols[1]], x = datObject[,1]),2:1])
  colnames(datObject) = datCols
  return(datObject)
}

# readin recip output for three seperate comparisons
dat1 <- read.table(file = "data/orthoAnalysis/orthoParisRecip/orthoMouseHumanDeci.txt", colClasses = "character")[,1:2]
dat1 <- colAsignment(datObject = as.matrix(dat1), snames = snames, ensGids = ensGids)

dat2 <- read.table(file = "data/orthoAnalysis/orthoParisRecip/orthoDogMouseDeci.txt", colClasses = "character")[,1:2]
dat2 <- colAsignment(datObject = as.matrix(dat2), snames = snames, ensGids = ensGids)

dat3 <- read.table(file = "data/orthoAnalysis/orthoParisRecip/orhtoDogHumanDeci.txt", colClasses = "character")[,1:2]
dat3 <- colAsignment(datObject = as.matrix(dat3), snames = snames, ensGids = ensGids)



# assign IDs to each col



m <- merge(dat1, dat2, 
      by.x = c(1:2,1:2)[c(colnames(dat1) == colnames(dat2)[1] , colnames(dat1) == colnames(dat2)[2]) ], 
      by.y = c(1:2,1:2)[c(colnames(dat2) == colnames(dat1)[1] , colnames(dat2) == colnames(dat1)[2]) ], all = FALSE)

m2 <- merge(m, dat3[,colnames(m[,2:3])], by.x = 2,by.y = 1, all= FALSE)

ortho = m2[m2[,3] == m2[,4],c(1,2,3)]
colnames(ortho)[3] = colnames(m)[3]


orthoDf <- data.frame(orthoID = paste("ortho", 1:nrow(ortho), sep = "_"), ortho)

sEnsGene <- list(NA,NA,NA)
names(sEnsGene) <- snames
for(s in snames){
  
  mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = s)
  ensGene <- dbGetQuery(mychannel, "SELECT * FROM ensGene;")
  
  
  orthoDf <- orthoDf[orthoDf[,s] %in% ensGene$name,]
  
  rownames(ensGene) <- ensGene$name
  
  ensGene <- data.frame(ensGene[orthoDf[,s],], orthoID = orthoDf$orthoID)
  
  rownames(ensGene) <- ensGene$orthoID
  
  
  sEnsGene[[s]] <- ensGene
}

# make sure we keep genes with the same exon count
for(s in snames){
  sEnsGene[[s]] <- sEnsGene[[s]][orthoDf$orthoID,]
}

orthoDf = orthoDf[sEnsGene[["hg19"]]$exonCount == sEnsGene[["mm9"]]$exonCount & 
      sEnsGene[["hg19"]]$exonCount == sEnsGene[["canFam3"]]$exonCount & 
      sEnsGene[["mm9"]]$exonCount == sEnsGene[["canFam3"]]$exonCount,]

write.table(x = orthoDf, file ="data/orthoAnalysis/orhtoList/ortholist.txt", quote = F,sep = "\t",row.names = F,col.names = T)


