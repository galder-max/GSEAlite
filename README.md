# GSEAlite
Gene Set Enrichment Analysis in R

# Install

> R CMD build GSEAlite
> R CMD install GSEAlite_1.0.tar.gz


# Run 

## Load library

> library(GSEAlite)


## Single or multiple genesets analysis

Follow the example in GSEA.listofgenes

> ?GSEA.listofgenes

To read in GMT files from the MSigDB

> ?readgmt

## Single geneset analysis + plot

Follow the example in GSEAplot

> ?GSEAplot


## Run on expression, methylation or other types of data


You can use any scoring function to use within GSEA for the enrichment.

Two such functions are already coded within the R package environment: snr.SNR and snr.FC, for the t-test-like signal-to-noise ratio and the fold-change, repsectively.

To use a new scoring function, define a new function in the global environment with its name following the syntax:

> snr.NAMEOFSCORINGFUNCTION = function(m, o)

For example you could redefine a function that is a median fold-change:

> snr.MedianFC = function( m, o) apply(m,1,function(x) median(x[o==levels(o)[1]]) -median(x[o==levels(o)[2]]))

GSEA can then be run using snr.MedianFC as scoring function, by specifying as input parameters:

> metric = "MedianFC"



