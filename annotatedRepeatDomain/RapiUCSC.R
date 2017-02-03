
# bioconductor gene coordinates ar epreety bad 

# maybe using a UCSC api might be better 


library("RMySQL")


mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = "hg19")


query <- function(...) dbGetQuery(mychannel, ...)


A <- query("SELECT name, chrom, txStart, txEnd FROM mm9.knownGene;")



dbListTables(con = mychannel)

(dbObj = mychannel)
RMySQL::show(object = mychannel)


m <- dbListTables(mychannel)

a <- dbGetQuery(mychannel, "SELECT * FROM phastCons46wayPlacental LIMIT 1000;")


# so write a function that can extract introns from UCSC tables


# Rtrack layer is probably another effective way to analyze specifc sections of data.
s1name = "hg19"
s2name = "mm9"
s3name = "canFam3"

datas2 <- read.table(file = paste("data/repeatHotspot/", s1name, "/", s1name,"_",s2name,"_conDif.txt", sep = ""), header = T)
dat <- datas2[ datas2$genome == "ref" & datas2$conState == "dif" & datas2$repGroup == "new_L1", ]

library(rtracklayer)

seqInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo")

seqlevels(range) = seqInfo[,1]
seqlengths(range) = seqInfo[,2]


mySession <- browserSession()
genome(mySession) <- "hg19"
track.names <- trackNames(ucscTableQuery(mySession))
head(track.names)

track <- track.names[grep("Cons 46", names(track.names))]

table.names = tableNames(ucscTableQuery(mySession))

table.names[grep("46way",table.names )]

table = "phastCons46way"

range <- GRanges(seqnames = Rle(dat$chr), 
                 ranges = IRanges(start = dat$start, width = 1000))


cons <- ucscTableQuery(mySession, track=track,
 range=range[1:100], table=table)

a <- track(cons)

# seems to be a limit on bp number
# not sure how this works with a ranged list object but i can see it getting 
# quite memmory intensive quite quickly.
# maybe i could go home and look at removing regions

## Not run: 
session <- browserSession()
genome(session) <- "mm9"
trackNames(session) ## list the track names
## choose the Conservation track for a portion of mm9 chr1
query <- ucscTableQuery(session, "Conservation",
                        GRangesForUCSCGenome("mm9", "chr12",
                                             IRanges(57795963, 57815592)))
## list the table names
tableNames(query)
## get the phastCons30way track
tableName(query) <- "phastCons30way"
## retrieve the track data
track(query)  # a GRanges object
## get a data.frame summarizing the multiple alignment
tableName(query) <- "multiz30waySummary"
getTable(query)

s <- rtracklayer::score(track(query))













