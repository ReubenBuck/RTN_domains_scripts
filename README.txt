#### Analysis scripts #####

-rtnDomainFunctions.R

various functions called by multiple other scripts

- dataGatering/rmskDownload.R

Input is a repeat masker table from http://repeatmasker.org/genomicDatasets/RMGenomicDatasets.html . 
Output is a series of plots describing various repeat statistics and a processed repeat masker dataset that includes repFamily and repGroup information. 
Encoded with in the script is the repeat grouping system for species analysed. 

- dataGathering/rmskGenomeSegmentation.R

Input is an R object repeat masker table from rmskDownload.R.
Output is an R object list containing a list of repeats and the genomic bins they are found in.

