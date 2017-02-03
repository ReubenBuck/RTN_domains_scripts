### venn diagrams


rm(list=ls())

setwd(dir = "Desktop/RTN_domains/")


repGroups = c("ancient", "new_SINE", "new_L1", "old_L1")
repCols = c("darkblue", "aquamarine3", "purple", "red")
snames <- c(s1name = "hg19", s2name = "mm9", s3name = "canFam3")

genomeIntersect <- NULL

vennAll <- NULL
for(name1 in snames){
  otherNames <- snames[snames != name1]
  for(name2 in 1:length(otherNames)){
    dat <- read.table(file = paste("data/repeatHotspot/",name1,"/",name1,"_",otherNames[name2],"_conDif.txt", sep = ""), header = TRUE)
    assign(paste("data",name2, sep = "_"), dat)
  }
  
  for(conState in c("con", "dif")){
    for(rep in repGroups){
      x = unique(as.character(data_1[data_1$genome == "ref"  & data_1$conState == conState & data_1$repGroup == rep,
                                     "hotspotID"]))
      y = unique(as.character(data_2[data_2$genome == "ref" & data_2$conState == conState & data_2$repGroup == rep, 
                                     "hotspotID"]))
      xy.int <- intersect(x,y)
      x.len <- length(x) - length(xy.int)
      y.len <- length(y) - length(xy.int)
      # save this data to a list object so it cna be plotted later. 
      venn1 = data.frame(ref = name1,que1 = otherNames[1], que2 = otherNames[2], 
                conState = conState, repGroup = rep ,que1S = x.len, 
                que2S = y.len, que1_que2_int = length(xy.int))
      vennAll = rbind(vennAll, venn1)
      
      df = data.frame(genome = name1,conState = conState , repGroup = rep ,domains = as.character(xy.int))
      genomeIntersect <- rbind(genomeIntersect, df)
    }
  }
}

write.table(genomeIntersect, file = "data/repeatHotspot/intersect.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# we can just save the intersect IDs. 
# the hotspots in dog, 
# which ones are conserved between dog and human 
# which ones are conserved between dog and mouse

# ones that are only conserved in dog


# 

# we can probably put all four in 

pdf(file = paste("plots/hotspotOverlap/vennDiagram/",snames[1],"_",snames[2],"_",snames[3],".pdf", sep = ""), width = 10,height = 10)
layout(matrix(1:6, nrow = 3, byrow = TRUE))
par(mar = c(1,1,1,1),oma=c(1,5,9,5))
for(ref in snames){
  
  for(conState in c("con", "dif")){
    
    plot.new()
    points(c(.37, .63),c(.5, .5),pch = 1, cex = 45)
    
    venn1 <- vennAll[ vennAll$ref== ref & vennAll$conState == conState,]
    
    venn1[,c("que1S","que2S","que1_que2_int")] <- round(venn1[,c("que1S","que2S","que1_que2_int")]/20, digits = 1)
    
    textSize = venn1[,c("que1S","que2S","que1_que2_int")] /  rowSums(venn1[,c("que1S","que2S","que1_que2_int")])
    
    
    #leftText
    text(x = rep(.25, 4),
         y = seq(from = .725, by = -.15, length.out = 4), 
         labels = venn1$que1S, col = repCols, cex = sqrt(textSize$que1S)*3.5)
    
    #rightText
    text(x =  rep(.75, 4),
         y = seq(from = .725, by = -.15, length.out = 4), 
         labels = venn1$que2S, col = repCols, cex = sqrt(textSize$que2S)*3.5)
    
    #centerText
    text(x = rep(.5, 4),
         y = seq(from = .725, by = -.15, length.out = 4), 
         labels = venn1$que1_que2_int, col = repCols, cex = sqrt(textSize$que1_que2_int)*3.5)
    
    #labels
    text(c(.2, .8),c(.95, .95),labels = c(as.character(venn1$que1[1]), as.character(venn1$que2[1])) )
    
    title(main = venn1$ref[1])
    
    
  }
}

mtext("HH domains", side = 3, line = 2, outer = TRUE, at = .25, cex = 2)
mtext("HL domains", side = 3, line = 2, outer = TRUE, at = .75, cex = 2)

dev.off()

