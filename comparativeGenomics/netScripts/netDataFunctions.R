rmSeqGapsFromNetOutput <- function(netOutput, seqGaps){
  # removes seq gaps from a net output file by taking the setdiff
  sDiff <- GenomicRanges::setdiff(netOutput, seqGaps)
  if(!all(overlapsAny(sDiff,netOutput))){
    stop("some netOutput ranges have been completely removed, check if net output file is correct")
  }
  ol <- findOverlaps(sDiff, netOutput)
  # this intersection makes sure no regions have been unnecesarily joined
  sDiff <- pintersect(sDiff[queryHits(ol)], netOutput[subjectHits(ol)])
  ol <- findOverlaps(sDiff, netOutput)
  if(!all(sDiff == sDiff[queryHits(ol)])){
    stop("not a single layer net output and orger has not been maintained, cannot correctly assign mcols")
  }
  mcols(sDiff) <- mcols(netOutput[subjectHits(ol)])
  return(sDiff)
}

