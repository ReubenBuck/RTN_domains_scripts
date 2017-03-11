### so how is this going to look 

# read in repeats and look across different reolutions to look at bin size enrichment 

# if we had a script that made repeat tables into R objects

# one idea is to download the raw rmsk files and use the chaining
# gives us a method to count insertions 

#setwd("/Users/labadmin/Desktop/RTN_domains/")
devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/rtnDomainFunctions.R")

library(dplyr)
#library(magrittr)
library(optparse)

rm(list = ls())


# need to include an error about having / in path way names


option_list = list(
  make_option(c("-i", "--inPath"), type="character", default=getwd(), 
              help="path to genome, default = current dir", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NA, 
              help="genome name", metavar="character"),
  make_option(c("-o", "--outPathRobject"), type="character", default=getwd(), 
              help="output path for processed repeats, default = current dir", metavar="character"),
  make_option(c("-p", "--outPathRplots"), type="character", default=getwd(), 
              help="output path for repeats plots, default = current dir", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

genome <- opt$genome

print(paste("genome =", opt$genome))
print(paste("inPath =", opt$inPath))
print(paste("outPathRobject =", opt$outPathRobject))
print(paste("outPathRplots =", opt$outPathRplots))

if(is.na(genome)){stop("need to specify genome")}
if(rev(strsplit(opt$inPath,split = "")[[1]])[1] != "/"){
  opt$inPath = paste(opt$inPath,"/",sep = "")
}
if(rev(strsplit(opt$outPathRobject,split = "")[[1]])[1] != "/"){
  opt$outPathRobject = paste(opt$outPathRobject,"/",sep = "")
}
if(rev(strsplit(opt$outPathRplots,split = "")[[1]])[1] != "/"){
  opt$outPathRplots = paste(opt$outPathRplots,"/",sep = "")
}


# Downloading repeats

rep_name <- paste( opt$inPath, genome,".fa.out" , sep = "")
rep <- read.table(file = rep_name, header = FALSE, skip = 3, 
                  col.names = c("SWscore", "perDiv", "perIns", "perDel", "genoChr", "genoStart", "genoEnd", "genoLeft", 
                                "strand", "repName", "repClass", "repStart", "repEnd", "repLeft", "repID"),
                  colClasses = c("integer", "double", "double", "double", "character", "integer", "integer","character", 
                                 "character", "character", "character", "character", "integer", "character", "character"),
                  fill = TRUE
                  )
rep <- rep[rep$repID != "",]
rep[rep$strand == "C",c("repStart", "repLeft")] <- rep[rep$strand == "C",c("repLeft", "repStart")]
rep[rep$strand == "C", "strand"] = "-"
rep$genoLeft <- substring(text = rep$genoLeft, 2, nchar(rep$genoLeft) - 1)
rep$repLeft <- substring(text = rep$repLeft, 2, nchar(rep$repLeft) - 1)



# we should make width proportional to the total number of repeats in each family
rep$genoLeft = as.integer(rep$genoLeft)
rep$repStart = as.integer(rep$repStart)
rep$repLeft = as.integer(rep$repLeft)





rep$repFamily <- NA
rep$repGroup <- NA

rep$repFamily[rep$repClass =="SINE/MIR"] <- "MIR"
rep$repFamily[rep$repClass =="LINE/L2" ] <- "L2"
rep$repGroup[rep$repFamily == "L2" | rep$repFamily == "MIR"] <- "ancient"


rep$repFamily[grep("L1ME", rep$repName)] <- "L1ME"
rep$repFamily[grep("L1MD", rep$repName)] <- "L1MD"
rep$repFamily[grep("L1MC", rep$repName)] <- "L1MC"
rep$repFamily[grep("L1MB", rep$repName)] <- "L1MB"
rep$repGroup[rep$repFamily == "L1ME" | rep$repFamily == "L1MD" | rep$repFamily == "L1MC" | rep$repFamily == "L1MB"] <- "old_L1"



if(genome == "hg19"){
  rep$repFamily[grep("AluJ", rep$repName)] = "AluJ"
  rep$repFamily[grep("AluS", rep$repName)] = "AluS"
  rep$repFamily[grep("AluY", rep$repName)] = "AluY"
  rep$repGroup[rep$repFamily == "AluJ" | rep$repFamily == "AluS" | rep$repFamily == "AluY" ] <- "new_SINE"
  
  
  rep$repFamily[grep("L1PB", rep$repName)] = "L1PB"
  rep$repFamily[grep("L1PA", rep$repName)] = "L1PA"
  rep$repFamily[grep("L1HS", rep$repName)] = "L1HS"
  rep$repFamily[grep("L1MA", rep$repName)] <- "L1MA"
  rep$repGroup[rep$repFamily == "L1PB" | rep$repFamily == "L1PA" | rep$repFamily == "L1HS" | rep$repFamily == "L1MA"] <- "new_L1"
  
}else if(genome == "mm9" | genome == "mm10"){
  
  rep$repFamily[grep("PB", rep$repName)] = "PB"
  rep$repFamily[grep("B1_", rep$repName)] = "B1"
  rep$repFamily[grep("B1F", rep$repName)] = "B1"
  rep$repFamily[grep("B2", rep$repName)] = "B2"
  rep$repFamily[grep("B3", rep$repName)] = "B3"
  rep$repFamily[grep("B4", rep$repName)] = "B4"
  rep$repGroup[rep$repFamily == "PB" | rep$repFamily == "B1" | rep$repFamily == "B2" | rep$repFamily == "B3" | rep$repFamily == "B4"] <- "new_SINE"
  
  rep$repFamily[grep("Lx", rep$repName)] = "Lx"
  rep$repFamily[grep("L1Md", rep$repName)] = "L1Md"
  rep$repFamily[grep("L1_Mus", rep$repName)] = "L1_Mus"
  rep$repFamily[grep("L1_Mur", rep$repName)] <- "L1_Mur"
  rep$repFamily[grep("L1_Mm", rep$repName)] <- "L1_Mm"
  rep$repFamily[grep("L1MA", rep$repName)] <- "L1MA"
  rep$repGroup[rep$repFamily == "Lx" | rep$repFamily == "L1Md" | rep$repFamily == "L1_Mus" | rep$repFamily == "L1_Mur" | rep$repFamily == "L1_Mm" | rep$repFamily == "L1MA"] <- "new_L1"
  
} else if(genome == "oryCun2"){
  
  rep$repFamily[rep$repName == "CSINE1"] = "CSINE1"
  rep$repFamily[rep$repName == "CSINE2"] = "CSINE2"
  rep$repFamily[rep$repName == "CSINE2B"] = "CSINE2B"
  rep$repFamily[rep$repName == "CSINE3A"] = "CSINE3A"
  rep$repGroup[rep$repFamily == "CSINE1" | rep$repFamily == "CSINE2" | rep$repFamily == "CSINE2B" | rep$repFamily == "CSINE3A"] <- "new_SINE"
  
  rep$repFamily[grep("L1A_O", rep$repName)] = "L1A_OC"
  rep$repFamily[grep("L1A2_O", rep$repName)] = "L1A2_OC"
  rep$repFamily[grep("L1C_O", rep$repName)] = "L1C_OC"
  rep$repFamily[grep("L1MA", rep$repName)] <- "L1MA"
  rep$repGroup[rep$repFamily == "L1A_OC" | rep$repFamily == "L1A2_OC" | rep$repFamily == "L1C_OC" | rep$repFamily == "L1MA"] <- "new_L1"
  
  
}else if(genome == "canFam3"){
  
  rep$repFamily[grep("SINEC_c", rep$repName)] <- "SINEC_c"
  rep$repFamily[grep("SINEC_b", rep$repName)] <- "SINEC_b"
  rep$repFamily[grep("SINEC_a", rep$repName)] <- "SINEC_a"
  rep$repFamily[grep("SINEC_old", rep$repName)] <- "SINEC_old"
  rep$repFamily[grep("SINEC_Cf", rep$repName)] <- "SINEC_Cf"
  rep$repGroup[rep$repFamily == "SINEC_c" | rep$repFamily == "SINEC_b" | rep$repFamily == "SINEC_a" | rep$repFamily == "SINEC_old" | rep$repFamily == "SINEC_Cf"] <- "new_SINE"
  
  
  rep$repFamily[grep("L1_Carn", rep$repName)] <- "L1_Carn"
  rep$repFamily[grep("L1_Canid", rep$repName)] <- "L1_Canid"
  rep$repFamily[grep("L1_Canis", rep$repName)] <- "L1_Canis"
  rep$repFamily[grep("L1_Cf", rep$repName)] <- "L1_Cf"
  rep$repFamily[grep("L1MA", rep$repName)] <- "L1MA"
  rep$repGroup[rep$repFamily == "L1_Carn" | rep$repFamily == "L1_Canid" | rep$repFamily == "L1_Canis" | rep$repFamily == "L1_Cf" | rep$repFamily == "L1MA"] <- "new_L1"
  
  
}else if(genome == "rheMac3"){
  rep$repFamily[grep("AluJ", rep$repName)] = "AluJ"
  rep$repFamily[grep("AluS", rep$repName)] = "AluS"
  rep$repFamily[grep("AluY", rep$repName)] = "AluY"
  rep$repGroup[rep$repFamily == "AluJ" | rep$repFamily == "AluS" | rep$repFamily == "AluY" ] <- "new_SINE"
  
  
  rep$repFamily[grep("L1PB", rep$repName)] = "L1PB"
  rep$repFamily[grep("L1PA", rep$repName)] = "L1PA"
  rep$repFamily[grep("L1_RS", rep$repName)] = "L1_RS"
  rep$repFamily[grep("L1MA", rep$repName)] <- "L1MA"
  rep$repGroup[rep$repFamily == "L1PB" | rep$repFamily == "L1PA" | rep$repFamily == "L1_RS" | rep$repFamily == "L1MA"] <- "new_L1"
  
}else{
  print("genome unknown")
}


newDir <- paste(opt$outPathRplots,paste(genome,"Plots/",sep=""), sep = "")
dir.create(newDir)

# so what plots are we looking for?

L1s <- filter(rep, repClass == "LINE/L1")
L1sCov <- aggregate(L1s$genoEnd - L1s$genoStart, list(L1s$repName), sum)
L1statAll <- boxplot(L1s$perDiv ~ L1s$repName, plot = FALSE)
L1statAll <- data.frame(repName = L1statAll$names, t(L1statAll$stats))
L1statAll <- merge(L1statAll, L1sCov, by = 1)
L1statAll <- merge(L1statAll, distinct(L1s, repName, repFamily, repGroup)[,c("repName", "repFamily", "repGroup")], by = 1)
L1statAll <- L1statAll[order(L1statAll$X3),]
L1statAll$plotNames <- ""
L1statAll$plotNames[is.na(L1statAll$repGroup)] <- as.character(L1statAll$repName[is.na(L1statAll$repGroup)])

pdf(file = paste(newDir,"/L1classification.pdf", sep = ""))
par(mar = c(5,5,5,5))
xBarplot <- barplot(L1statAll$X3,width = L1statAll$x,las = 1, horiz = FALSE, ylim = c(0,35), 
                    names = "", las = 2, col = rainbow(length(unique(L1statAll$repFamily)))[as.factor(L1statAll$repFamily)],
                    density = c(-1,20)[as.integer(as.factor(L1statAll$repGroup))], xlab = "L1 names", 
                    ylab = "divergence from consensus (%)", main = paste(genome,"L1 classification"))
arrows(x0 = xBarplot, x1 = xBarplot,y0 = as.double(as.character(L1statAll$X2)),
       y1 = as.double(as.character(L1statAll$X4)), angle = 90, code = 3, length = .05)
legend("topleft", title = "L1 families", legend = levels(as.factor(L1statAll$repFamily)), fill = rainbow(length(unique(L1statAll$repFamily))), bty = "n")
legend("top", title = "L1 groups", legend = c("new L1", "old L1"), fill = c(1,1), density = c(-1,20), bty = "n")
dev.off()

pdf(file = paste(newDir,"/L1ungrouped.pdf", sep = ""))
par(mar = c(10,5,5,5))
ungrouped <- aggregate(L1statAll$x,list(L1statAll$plotNames), sum)
ungrouped$Group.1[ungrouped$Group.1 == ""] = "Grouped L1s"
barplot(ungrouped$x/sum(ungrouped$x), names = ungrouped$Group.1, las = 2,ylim = c(0,1), 
        ylab = "Total L1 coverage (%)", main = paste(genome,"ungrouped L1 elements"))
dev.off()
# so now that we have the plot 

chainRatio <- table(as.factor(rep$repClass[match(unique(rep$repID), rep$repID)]))/table(as.factor(rep$repClass))
repCov <- aggregate(x = rep$genoEnd - rep$genoStart, by = list(rep$repClass), FUN = sum)

pdf(file = paste(newDir,"/allRepChained.pdf", sep = ""))
par(mar = c(5,10,5,5))
barplot(chainRatio, horiz = TRUE, las = 1,width = repCov$x, xlab = "chained repeat instances / unchained repeat instances")
dev.off()

# repeat family coverage



# repeat group coverage

# do a stacked barplot where we put chained infront of unchained
# the chained and unchained stuff only makes sense with insertion numbers

grouped <- filter(rep, repGroup == "new_SINE" | repGroup == "new_L1" | repGroup == "old_L1" | repGroup == "ancient")
chainedGroup <- distinct(grouped,repID,.keep_all = TRUE)

groupedTab <- table(grouped$repGroup)[c("new_SINE", "new_L1", "old_L1", "ancient")]
chainedGroupTab <- table(chainedGroup$repGroup)[c("new_SINE", "new_L1", "old_L1", "ancient")]

pdf(file = paste(newDir,"/groupChained.pdf", sep = ""))
par(mar = c(10,5,5,5))
barplot(groupedTab, col = c("aquamarine3", "purple", "red", "darkblue"), density = 40, space = c(.2,.2,.2,.2), ylim = c(0,1.8*10^6), 
        ylab = "repeat masker hits", main = paste(genome, "retrotransposon groups"), las = 3)
barplot(chainedGroupTab,col = c("aquamarine3", "purple", "red", "darkblue"), width = .9,space = c(.28,.335,.335,.335),add = T, axes = FALSE, names = FALSE)
legend("topright", legend = c("unchained", "chained"), density = c(40,-1), bty = "n")
dev.off()

familied <- filter(rep, !is.na(repFamily))
chainedFamily <- distinct(familied,repID, .keep_all = TRUE)

familyNames <- distinct(familied, repFamily, repGroup, .keep_all = TRUE)[,c("repFamily", "repGroup")]
familyNames$repGroup <- factor(familyNames$repGroup, levels = c("new_SINE", "new_L1", "old_L1", "ancient"))

familiedTab <- table(familied$repFamily)[familyNames$repFamily[order(familyNames$repGroup)]]
chainedFamilyTab <- table(chainedFamily$repFamily)[familyNames$repFamily[order(familyNames$repGroup)]]

pdf(file = paste(newDir,"/familyChained.pdf", sep = ""))
par(mar = c(10,5,5,5))
barplot(familiedTab, density = 40,  ylim = c(0,.8*10^6), las = 3,
        ylab = "repeat masker hits", main = paste(genome, "retrotransposon families"), 
        col = c("aquamarine3", "purple", "red", "darkblue")[as.integer(sort(familyNames$repGroup))],
        width = 1)
barplot(chainedFamilyTab, add = T, axes = FALSE, names = FALSE,
        col = c("aquamarine3", "purple", "red", "darkblue")[as.integer(sort(familyNames$repGroup))],
        width = .9, space = c(.28,rep(.333, length(chainedFamilyTab) - 1)))
legend("topright", legend = c("unchained", "chained"), density = c(40,-1), bty = "n")
legend("top", legend = levels(familyNames$repGroup), fill  = c("aquamarine3", "purple", "red", "darkblue"), bty = "n")
dev.off()

## we could probably do our binning here too across different sizes

rep$repUID <- paste(rep$genoChr,":" ,rep$genoStart,"-" ,rep$genoEnd, sep = "")
save(rep, file=paste(opt$outPathRobject,genome,".RData", sep = ""))

