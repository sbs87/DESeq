source("http://bioconductor.org/biocLite.R")
biocLite("DESeq")
biocLite("pasilla")
library("DESeq")
datafile =system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )
datafile
pasillaCountTable = read.table( datafile, header=TRUE, row.names=1 )


pasillaDesign = data.frame(
  row.names = colnames( pasillaCountTable ),
  condition = c( "untreated", "untreated", "untreated",
                 "untreated", "treated", "treated", "treated" ),
  libType = c( "single-end", "single-end", "paired-end",
               "paired-end", "single-end", "paired-end", "paired-end" ) )

pairedSamples = pasillaDesign$libType == "paired-end"
countTable = pasillaCountTable[ , pairedSamples ]
condition = pasillaDesign$condition[ pairedSamples ]
head(countTable)

cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )
head(counts(cds))
library("reshape")
library("ggplot2")
temp<-melt(counts(cds))
temp2<-melt(counts(cds, normalized=TRUE ))
head(temp)
temp[temp$X1=="FBgn0025692",]
ggplot(temp2)+ geom_bar(aes(x=value,fill=X2))+ facet_wrap(~ X2,ncol=1)


cds = estimateDispersions( cds )
str( fitInfo(cds) )
plotDispEsts( cds )
head(fitInfo(cds))

res = nbinomTest( cds, "untreated", "treated" )
head(res)
head(counts(cds, normalized=TRUE ))
mean(c(86.06763,62.8034302))
74.4355310/78.155755
plotMA(res)
