rm(list=ls())

load("Desktop/RTN_domains/R_objects/mappedGaps/mm10.hg19.netData.RData")



extractClustersFromSelfHits <- function(hits)
{
  stopifnot(is(hits, "Hits"))
  N <- queryLength(hits)
  stopifnot(N == subjectLength(hits))
  h <- union(hits, t(hits))
  qh <- queryHits(h)
  sh <- subjectHits(h)
  cid <- cid0 <- seq_len(N)  # cluster ids
  while (TRUE) {
    cid2 <- pmin(cid, selectHits(h, "first"))
    if (identical(cid2, cid))
      break
    cid <- cid2
    h <- Hits(qh, cid[sh], N, N)
  }
  unname(splitAsList(cid0, cid))
}


library(GenomicRanges)

# an ajusted max gap, something that will change gieven the size of the elements

ol <- findOverlaps(queFill.gr$queRanges)
ol <- ol[queryHits(ol)!=subjectHits(ol)]
quePint <- pintersect(queFill.gr$queRanges[queryHits(ol)], queFill.gr$queRanges[subjectHits(ol)])
queSideDiff <- log10(width(queFill.gr$queRanges[queryHits(ol)])/width(quePint))
subSideDiff <- log10(width(queFill.gr$queRanges[subjectHits(ol)])/width(quePint))
ol <- ol[queSideDiff < .1 & subSideDiff < .1 & queSideDiff > -.1 & subSideDiff > -.1]

ol2 <- ol[!isRedundantHit(ol)]

# we can remoe redundant hits, this is good, this way we may be able to build pile
clusters <- extractClustersFromSelfHits(ol2)
clusterIDs <- rep(1:length(clusters), unlist(lapply(clusters,length)))
queSegDup <- queFill.gr[unlist(clusters)]
queSegDup$clusterID <- clusterIDs

sum(width(queSegDup))

# so there's ties in some of our connections
# exclusive clusters is what we are after

# edges that connect clusters



ol <- findOverlaps(refFill.gr$queRanges)
ol <- ol[queryHits(ol)!=subjectHits(ol)]
refPint <- pintersect(refFill.gr$queRanges[queryHits(ol)], refFill.gr$queRanges[subjectHits(ol)])
queSideDiff <- log10(width(refFill.gr$queRanges[queryHits(ol)])/width(refPint))
subSideDiff <- log10(width(refFill.gr$queRanges[subjectHits(ol)])/width(refPint))
ol <- ol[queSideDiff < .1 & subSideDiff < .1 & queSideDiff > -.1 & subSideDiff > -.1]
ol2 <- ol[!isRedundantHit(ol)]

# we can remoe redundant hits, this is good, this way we may be able to build pile
clusters <- extractClustersFromSelfHits(ol2)
clusterIDs <- rep(1:length(clusters), unlist(lapply(clusters,length)))
refSegDup <- refFill.gr[unlist(clusters)]
refSegDup$clusterID <- clusterIDs



round((sum(width(queSegDup)) - sum(width(refSegDup)))/1e6)
round((sum(width(queFill.gr)) - sum(width(refFill.gr)))/1e6)
# could rope these areas off?

# interestingly our differences are equall to the differences reported else where.
# This is consistant with repeated gap structure
# it would be cool if I could actually identify the same segmental duplications talked about in the liturature.


# alternativly we could fetch the chains of each one of these and look at the proportion of the chain that is involved in these events and remove it from the genome
# the problem is that we are counting events at specific locations more than once,
# especially when we attempt to map our deletions back

# it makes sense to remove offending chains

# there is two ways for this to happen, 

# 

ol <- findOverlaps(queSegDup$queRanges, refFill.gr)

df <- data.frame(queSegDup[queryHits(ol)],
                 refFill.gr[subjectHits(ol)])

olF <- poverlaps(ranges(queSegDup[queryHits(ol)]), ranges(refFill.gr$queRanges[subjectHits(ol)]))
olF[as.character(seqnames(queSegDup[queryHits(ol)])) != as.character(seqnames(refFill.gr$queRanges[subjectHits(ol)]))] <- FALSE
# we get our seed regions based on an overlap with the ref genome

# df[olF,][135:165,]
# segmental duplication seed regions


# not all of our original hits are here
seedRegionOL <- ol[olF]
nonSeedRegionOL <- ol[!olF]

queSegSeed<- (queSegDup[queryHits(seedRegionOL)])
queSegNonSeed<- (queSegDup[queryHits(nonSeedRegionOL)])

ol2 <- findOverlaps(queSegSeed, queSegNonSeed)

queSegSeed[c(51,52,53)]
queSegNonSeed[c(457,458,459,460,461,462,463)]


ol2[1:10]
df[olF,][133:142,]


dfChr2 <- df[olF,][df[olF,]$queRanges.seqnames =="chr2",]
dfChr2 <- dfChr2[order(dfChr2$start.1),]

a <- dfChr2$seqnames[1:(nrow(dfChr2)-1)] != dfChr2$seqnames[2:(nrow(dfChr2))]
(1:nrow(dfChr2))[a]
dfChr2[55:75,]

# we can use gaps and interect it with our gap data to get the total seq lenght of seg duplications.

# lets pull out the bad regions?





queSegChainSeed <- unique(queFill.gr[queFill.gr$chainID %in% unique(queSegDup[queryHits(seedRegionOL)]$chainID)])
queSegChainNonSeed <- unique(queFill.gr[queFill.gr$chainID %in% unique(queSegDup[queryHits(nonSeedOL)]$chainID)])

ol2 <- findOverlaps(unique(queSegDup[queryHits(seedRegionOL)]), unique(queSegDup[queryHits(nonSeedOL)]))
df2 <- data.frame(unique(queSegDup[queryHits(seedRegionOL)])[queryHits(ol2)],  unique(queSegDup[queryHits(nonSeedOL)])[subjectHits(ol2)])

aggSeed <- aggregate(width(queSegChainSeed), by = list(queSegChainSeed$chainID), FUN = sum)
aggNonSeed <- aggregate(width(queSegChainNonSeed), by = list(queSegChainNonSeed$chainID), FUN = sum)


# seed regions and non seed regions



hist(log10(aggNonSeed$x), breaks = 100)
hist(log10(aggSeed$x), add = TRUE, breaks = 100, col = 2, density = 0)


Int <- intersect(aggSeed$Group.1, aggNonSeed$Group.1)


smoothScatter(log10(aggSeg$x),log10(aggChain$x))


sum(aggChain$x)

dim(aggSeg)
dim(aggChain)



ol <- findOverlaps(refFill.gr$queRanges,queSegDup)

smoothScatter(log10(width(refFill.gr$queRanges[queryHits(ol)])), log10(width(queSegDup[subjectHits(ol)])))

pInt <- pintersect()


# give them Ids for each clique of duplicates

# one member from each segDup pile should be supproted by an overlap from another species. 


queSegChain2 <- queFill.gr[queFill.gr$chainID %in% unique(queSegDup$chainID[subjectHits(ol)])]

aggChain2 <- aggregate(width(queSegChain2), by = list(queSegChain2$chainID), FUN = sum)
aggSeg2 <- aggregate(width(queSegDup[subjectHits(ol)]), by = list(queSegDup$chainID[subjectHits(ol)]), FUN = sum)

smoothScatter(log10(aggSeg2$x),log10(aggChain2$x))


hist(log10(aggChain$x), breaks = 100, ylim = c(0,100))
hist(log10(aggChain2$x), add = TRUE, breaks = 100, col = 2, density = 0)







