#!/usr/bin/env Rscript

# lift over script

# supply chian path 

# supply bed file

rm(list = ls())

library(optparse)

option_list = list(
  make_option(c("-b", "--bedFile"), type="character", default=NA, 
              help="bedfile to go in", metavar="character"),
  make_option(c("-c", "--chain"), type="character", default=NA, 
              help="chain file", metavar="character"),
  make_option(c("-o", "--outBed"), type="character", default=NA, 
              help="output file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


library(rtracklayer)
library(dplyr)



bed.df <- read.table(opt$bedFile, header = FALSE, 
                      col.names = c("chr", "start", "end", "hotspotID"))

bed.gr <- GRanges(seqnames = Rle(bed.df$chr),
                  ranges = IRanges(start = bed.df$start, end = bed.df$end),
                  hotspotID = bed.df$hotspotID
)


ch = import.chain(opt$chain)


lift.gr <- liftOver(bed.gr, ch)


names(lift.gr) <- unique(as.data.frame(lift.gr)$hotspotID)

lift.gr <- reduce(lift.gr)


tDF <- as.data.frame(lift.gr)


aDF <- summarise(group_by(tDF, group_name, seqnames, group), width = sum(width), start = min(start), end = max(end)) %>% 
  arrange(group_name, desc(width)) 
aDF <- filter(group_by(aDF, group) ,width == max(width)) %>%
  mutate(queWidth = end - start + 1)


aDF <- aDF[aDF$width/aDF$queWidth > .1 & aDF$width/50000 > .1,]



aDF <- arrange(aDF, group)

outDF <- data.frame(chr = aDF$seqnames, start = aDF$start, end = aDF$end, hotspotID = aDF$group_name)

write.table(outDF, file = opt$outBed, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)




