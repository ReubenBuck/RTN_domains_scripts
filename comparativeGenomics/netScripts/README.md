## Detecting regional variation in DNA gain or loss

A collection of custom scripts designed to calculate regional levels of DNA gain and loss between two species.

##### getGapsNet
Takes UCSC net file as input and extracts gaps.
Returns reference and query start and end coordinates for each gap as output.
	
##### getFillsNet
Takes UCSC net file as input and extracts fills.
Returns reference and query start and end coordinates for each fill as output.

##### formatData.R


##### sortRepeatFamilies


##### syntheticGenomeBuilder.R
Annotates gaps between a reference and a query as gain or loss events.
Places all gaps within the reference species genomic background.

##### synthChrom.R
Captures chromosomal distribution of DNA gain and loss by binning at a window size of 200 kb

##### netDataFunctions.R
important functions are the genome stretching ones.