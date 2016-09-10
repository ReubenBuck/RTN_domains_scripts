# so with this script will be used to bin genomes and create repeat data objects with binned information 

library(optparse)

#path = 
#genome = commandArgs(trailingOnly = TRUE)[1]
#print(genome)


#load(genome)



option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path to genome", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="out.txt", 
              help="genome name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


print(opt$genome)
print(opt$path)
print(opt$out)






