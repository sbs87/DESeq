# This script takes alignment files (bam) and creates counts tables as experimentClass that can be subsequently used by DESeq
# Steven Smith
# Setember 21, 2014
# The summarizeOverlaps step takes up a tremendous amount of RAM and processing power. This needs to be run on the grid.

### Set up enviornment
source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
#biocLite("GenomicAlignments")
library(GenomicAlignments)
#biocLite("GenomicFeatures")
library(GenomicFeatures)

setwd("/Volumes/Solid_State_120GB/Alignment_Files_Temp")
#setwd("/Users/stsmith/Projects/Host_Response_Lactic_Acid/Analytical_Results/tophat_out")

### Transcript file/database
#### Create a sql lite file using GFF, then save it
#h_sapiens_cufflinks<-makeTranscriptDbFromGFF("genes_no_seleno.gtf",format="gtf")
#saveDb(h_sapiens_cufflinks, file="h_sapiens_cufflinks.sqlite")

#### Create genomic features using the transcript database. In this case, create exons based on gene
#txdb<-loadDb("/Volumes/Solid_State_120GB 3/Alignment_Files_Temp/h_sapiens_cufflinks.sqlite")
#exbygene <- exonsBy(txdb, "gene")
#save(exbygene,file="/Volumes/Solid_State_120GB 3/Alignment_Files_Temp/TxDb.Hsapiens.ensGene.exbygene")

#### Or skip all this and load the previously created transcripts by gene file
load("TxDb.Hsapiens.ensGene.exbygene")

### Now, read in alignment files and create a BamFileList based on these pointers
alignment_file_list<-as.character(read.table("HR_LA_accepted_hits_list.txt",sep="",header=F)$V1)

BFL<-list() # List to hold all accetped_hits.bam files
for(i in 1:length(alignment_file_list)){
  alignment_file<-alignment_file_list[i] 
  BFL[[i]]<-BamFile(file = alignment_file) ## Add each file to the BamFileList
}
reads=BamFileList(BFL) ## Officially create a BamFileList object

### And now, the heavy lifting (summarieOverlaps of 56 alignment files)
sumOl<-summarizeOverlaps(exbygene, reads, mode="Union") ## union is reccomended by DESeq
