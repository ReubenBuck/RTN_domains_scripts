
hg19 <- read.table(file = "~/Desktop/hg19AncestralGenomeCov.txt")
hg19 <- hg19[2:nrow(hg19),]

mm10 <- read.table(file = "~/Desktop/mm10AncestralGenomeCov.txt")
mm10 <- mm10[2:nrow(mm10),]


### so this is whole genome, look at this interms of covered bases only
pdf(file = "~/Desktop/RTN_domains/RTN_domain_plots/comparativeGenomics/ancestralDNA/cumulative.pdf", onefile = TRUE)

plot((hg19$V2), (cumsum((hg19$V3))/sum(hg19$V3)) * 100, xlim = c(1,38), ylim = c(0,100),
     type = "l", cex = .5, lwd = 3,
     ylab = "Ancestral DNA (cumulative %)", xlab = "Overlapping outgroup genomes")
lines((mm10$V2), (cumsum((mm10$V3))/sum(mm10$V3))*100, col =2, type = "l", cex = .5, lwd = 3)
legend("bottomright", legend = c("hg19", "mm10"), lty = 1, col = c(1,2), lwd = 3)

plot((hg19$V2), (cumsum((hg19$V3))), xlim = c(1,38),
     type = "l", cex = .5, lwd = 3,
     ylab = "Ancestral DNA (cumulative sum)", xlab = "Overlapping outgroup genomes")
lines((mm10$V2), (cumsum((mm10$V3))), col =2, type = "l", cex = .5, lwd = 3)
legend("bottomright", legend = c("hg19", "mm10"), lty = 1, col = c(1,2), lwd = 3)

dev.off()
