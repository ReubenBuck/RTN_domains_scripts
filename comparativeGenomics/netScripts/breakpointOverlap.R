

library(GenomicRanges)

rm(list = ls())

options(stringsAsFactors = FALSE)



library(TxDb.Hsapiens.UCSC.hg19.knownGene)
hg19gene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

bPoint <- read.table("~/Downloads/hglft_genome_7af_163030.bed", 
                     col.names = c("seqnames", "start", "end"))
bPoint.gr <- GRanges(bPoint)
seqinfo(bPoint.gr) <- seqinfo(hg19gene)
# overlap with gaps 
# genes
# look at bins that associate with these breakpoints and measure their content 


# we can also look at proximity to break points

specRef = "hg19"
specQue = "mm10"

devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")

load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/shiftData/",specRef,".expand.breaks.RData", sep = ""))
hg19refShift <- newSynthRefShift



pdf(file="Documents/RTN_domain_writing/manuscriptRound2/bioArxiv/sup/TexFigs/breakpointOL.pdf", height = 10, width = 10)
layout(1:2)

mains <- c("Ancestral element", "Recent transposon")
names(mains) <- c("", ".rep")

for( i in 1:2){
load(paste("~/Documents/dna_turnover/workStationDesktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/", specRef,".synthBin",names(mains)[i],".RData", sep = ""))
hg19SynthBinNorm.gr <- synthBin.gr
rmBases <- rowSums(data.frame(mcols(hg19SynthBinNorm.gr)[,c("missingGap", "seqGap")]))
mcols(hg19SynthBinNorm.gr) <- (data.frame(mcols(hg19SynthBinNorm.gr))[,c("refIns", "refDel", "queIns", "queDel")]/(width(hg19SynthBinNorm.gr) - rmBases)) * 2e5
hg19SynthBinNorm.gr <- hg19SynthBinNorm.gr[width(hg19SynthBinNorm.gr) - rmBases > 150000]

bStretch <- genoExpandStretch(sort(bPoint.gr), synthGenome = hg19refShift, 
                                expandedSeqlengths = seqlengths(hg19SynthBinNorm.gr))


bpBin <- hg19SynthBinNorm.gr[overlapsAny(hg19SynthBinNorm.gr, bStretch, maxgap = 2e5)]
genomeBin <-  hg19SynthBinNorm.gr[!overlapsAny(hg19SynthBinNorm.gr, bStretch, maxgap = 2e5)]

df <- c(data.frame(mcols(unique(bpBin)))[,c("refIns","refDel", "queIns","queDel")]/1000,
                 data.frame(mcols(unique(genomeBin)))[,c("refIns","refDel", "queIns","queDel")]/1000)
names(df) <- paste(names(df), c(rep(".bp", 4), rep(".genome", 4)), sep = "")

boxplot(df[order(names(df), decreasing = TRUE)], 
        col= c(scales::alpha(1, .2), scales::alpha(2, .2)), 
        outline = FALSE, notch = TRUE,
        las = 2, names = FALSE, xaxt = "n", xaxs = "i", xlim = c(.75,8.25), 
        ylab = "DNA turnover per bin (kb)", xlab = "Gap annotation",
        ylim = c(0,103), main = mains[i])
axis(side = 1, at = seq(1.5, by = 2, length.out = 4),
     labels = c("hg19 gain", "hg19 loss", "mm10 gain", "mm10 loss"))
abline(v = seq(2.5, by = 2, length.out = 3), lty = 2)
legend("topright", legend = c("Genome", "Break region"), 
       fill = c(scales::alpha(1, .2), scales::alpha(2, .2)), bty = "n")

}


dev.off()

# bins that past the threshold

# ovelapping with 200kb sourounding breakpoints as identified in Lemaitre et al 2008, BMC bioinformatics



