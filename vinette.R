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
nbinomTest
nbinomTestForMatrices
head(res)
head(counts(cds, normalized=TRUE ))
mean(c(86.06763,62.8034302))
74.4355310/78.155755
plotMA(res)

## Under the hood of nbinomTest and nbinomTestForMatrices...

condA="untreated"
condB="treated"
colA <- conditions(cds) == condA
colB <- conditions(cds) == condB
rawScvA <- fData(cds)[, paste("disp", dispTable(cds)[condA], 
                              sep = "_")]
rawScvB <- fData(cds)[, paste("disp", dispTable(cds)[condB], 
                              sep = "_")]

countsA=counts(cds)[, colA]
countsB=counts(cds)[, colB]
sizeFactorsA=sizeFactors(cds)[colA]
sizeFactorsB=sizeFactors(cds)[colB]
dispsA=rawScvA
dispsB=rawScvB

kAs <- rowSums(cbind(countsA))
kBs <- rowSums(cbind(countsB))
mus <- rowMeans(cbind(t(t(countsA)/sizeFactorsA), t(t(countsB)/sizeFactorsB)))
fullVarsA <- pmax(mus * sum(sizeFactorsA) + dispsA * mus^2 *sum(sizeFactorsA^2), mus * sum(sizeFactorsA) * (1 + 1e-08))
fullVarsB <- pmax(mus * sum(sizeFactorsB) + dispsB * mus^2 *sum(sizeFactorsB^2), mus * sum(sizeFactorsB) * (1 + 1e-08))
?pmax
sumDispsA <- (fullVarsA - mus * sum(sizeFactorsA))/(mus *sum(sizeFactorsA))^2
sumDispsB <- (fullVarsB - mus * sum(sizeFactorsB))/(mus *sum(sizeFactorsB))^2
sapply(seq(along = kAs), function(i) {
  if (kAs[i] == 0 & kBs[i] == 0) 
    return(NA)
  ks <- 0:(kAs[i] + kBs[i])
  ps <- dnbinom(ks, mu = mus[i] * sum(sizeFactorsA), size = 1/sumDispsA[i]) * 
    dnbinom(kAs[i] + kBs[i] - ks, mu = mus[i] * sum(sizeFactorsB), 
            size = 1/sumDispsB[i])
  pobs <- dnbinom(kAs[i], mu = mus[i] * sum(sizeFactorsA), 
                  size = 1/sumDispsA[i]) * dnbinom(kBs[i], mu = mus[i] * 
                                                     sum(sizeFactorsB), size = 1/sumDispsB[i])
  stopifnot(pobs == ps[kAs[i] + 1])
  if (kAs[i] * sum(sizeFactorsB) < kBs[i] * sum(sizeFactorsA)) 
    numer <- ps[1:(kAs[i] + 1)]
  else numer <- ps[(kAs[i] + 1):length(ps)]
  min(1, 2 * sum(numer)/sum(ps))
})


nbinomTestForMatrices(counts(cds)[, colA], counts(cds)[,colB], sizeFactors(cds)[colA], sizeFactors(cds)[colB],rawScvA, rawScvB)
