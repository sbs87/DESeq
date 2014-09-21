source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
#browseVignettes("GenomicRanges")
library(GenomicRanges)
library(GenomicAlignments)
#browseVignettes("DESeq2")
library("GenomicFeatures")

txdb<-loadDb("/Volumes/Solid_State_120GB/Alignment_Files_Temp/h_sapiens_cufflinks.sqlite")
#library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
exbygene <- exonsBy(txdb, "gene")

#library(pasillaBamSubset)
setwd("/Volumes/Solid_State_120GB/Alignment_Files_Temp")
un1 <- "MED_rep1_th2out_accepted_hits.bam.sort.bam" #untreated1_chr4() # pointer to file name
un2 <- "MED_rep2_th2out_accepted_hits.bam.sort.bam"
gal <- readGAlignments(un1)
ga2 <- readGAlignments(un2)
save("gal",file="MED_rep1_accepted_hits.R")
save("ga2",file="MED_rep2_accepted_hits.R")
se.1 <- summarizeOverlaps(exbygene, un1, mode="IntersectionNotEmpty")
se.2 <- summarizeOverlaps(exbygene, un2, mode="IntersectionNotEmpty")
#biocLite("GenomicFeatures")

## Load bam files locally + sql lite file, then play around here. Will need GB of space. 
d =  pipe( 'ssh stsmith@diag-gateway.igs.umaryland.edu 
                      "cat /diag/cloud/stsmith/HostReponse_LacticAcid/data/tophat_out/MED_CT_rep1_th2out/MED_CT_rep1_th2out_accepted_hits.bam.sort.bam"' ), 
                header =F))
file('ssh stsmith@diag-gateway.igs.umaryland.edu "cat sratch/test.txt"')
readGAlignments(pipe('ssh stsmith@diag-gateway.igs.umaryland.edu "samtools view /diag/cloud/stsmith/HostReponse_LacticAcid/data/tophat_out/MED_CT_rep1_th2out/MED_CT_rep1_th2out_accepted_hits.bam.sort.bam"'))
