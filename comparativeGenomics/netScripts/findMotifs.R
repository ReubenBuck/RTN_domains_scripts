
rm(list = ls())

specRef <- "hg19"


# load in sequence data
if(specRef == "hg19"){
  library(BSgenome.Hsapiens.UCSC.hg19); wholeGenoSeq <- BSgenome.Hsapiens.UCSC.hg19
  
}
if(specRef == "mm10"){
  library(BSgenome.Mmusculus.UCSC.mm10); wholeGenoSeq <- BSgenome.Mmusculus.UCSC.mm10
}

# conserved motifs
L1Motif <- DNAString("AATTTT")
ctcfMotif <- DNAString("CCGCGNGGNGGCAG")

# species specifc motifs
if(specRef == "hg19"){
  prdm9Motif <- DNAString("CCNCCNTNNCCNC")
}
if(specRef == "mm10"){
  prdm9Motif <- DNAString("GNTGCTNCT")
}



motifList <- c(L1Motif = L1Motif, ctcfMotif = ctcfMotif, prdm9Motif =  prdm9Motif)

dict0 <- DNAStringSet(motifList)





writeHits <- function(seqname, matches, strand, file="", append=FALSE)
{
  if (file.exists(file) && !append)
    warning("existing file ", file, " will be overwritten with 'append=FALSE'")
  if (!file.exists(file) && append)
    warning("new file ", file, " will have no header with 'append=TRUE'")
  hits <- data.frame(seqnames=rep.int(seqname, length(matches)),
                     start=start(matches),
                     end=end(matches),
                     strand=rep.int("*", length(matches)),
                     patternID=names(matches),
                     check.names=FALSE)
  # maybe I can convert to a GRanges object and remove regions I don't want, since it is matching Ns
  write.table(hits, file=file, append=append, quote=FALSE, sep="\t",
              row.names=FALSE, col.names=!append)
}

runAnalysis1 <- function(dict0, BSgenomeSeq ,outfile="")
{
  genome <- BSgenomeSeq
  seqnames <- seqnames(genome)
  seqnames_in1string <- paste(seqnames, collapse=", ")
  cat("Target:", providerVersion(genome),
      "chromosomes", seqnames_in1string, "\n")
  append <- FALSE
  for (seqname in seqnames) {
    subject <- genome[[seqname]]
    cat(">>> Finding all hits in chromosome", seqname, "...\n")
    for (i in seq_len(length(dict0))) {
      patternID <- names(dict0)[i]
      pattern <- dict0[[i]]
      plus_matches <- matchPattern(pattern, subject, fixed = "subject")
      names(plus_matches) <- rep.int(patternID, length(plus_matches))
      writeHits(seqname, plus_matches, "+", file=outfile, append=append)
      append <- TRUE
      rcpattern <- reverseComplement(pattern)
      minus_matches <- matchPattern(rcpattern, subject, fixed = "subject")
      names(minus_matches) <- rep.int(patternID, length(minus_matches))
      writeHits(seqname, minus_matches, "-", file=outfile, append=append)
    }
    cat(">>> DONE\n")
  }
}




runAnalysis1(dict0 = dict0,BSgenomeSeq = wholeGenoSeq, outfile = paste("~/Desktop/",specRef,".motifs", sep = ""))





