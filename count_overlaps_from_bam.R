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
counts<-data.frame()
for(i in 1:length(alignment_file_list)){
  print(str(i))
  alignment_file<-alignment_file_list[i] 
  #BFL[[i]]<-BamFile(file = alignment_file) ## Add each file to the BamFileList
  #so<-summarizeOverlaps(exbygene, alignment_file, mode="Union")
  #save(so,file=paste(alignment_file,".overlap",sep=""))
  load(paste(alignment_file,".overlap",sep=""))
  if(i == 1){
    counts<-assays(so)[[1]]
  }else{
    counts<-cbind(counts,assays(so)[[1]])
  }
}
#reads=BamFileList(BFL) ## Officially create a BamFileList object
colnames(counts)<-gsub("_th2out_accepted_hits.bam.sort.bam","",colnames(counts))
### And now, the heavy lifting (summarieOverlaps of 56 alignment files)
#sumOl<-summarizeOverlaps(exbygene, reads, mode="Union") ## union is reccomended by DESeq
save(counts,file="counts.R")
library(DESeq2)
(colData<-data.frame(sample_name=colnames(counts),isoform=c(rep("D",times=6),rep("DL",times=4),rep("HCL",times=6),rep("L",times=6),rep("MED",times=4)),
      pH=c(rep(4,times=4),rep(7,times=2),rep(4,times=3),7,rep(4,times=4),7,7,4,4,4,4,7,7,rep("MED",times=4)),
      infection=c(rep("CT",times=2),rep("uninfected",times=4),rep("CT",times=2),rep("uninfected",times=2),rep("CT",times=2),rep("uninfected",times=4),
                  rep("CT",times=2),rep("uninfected",times=4),rep("CT",times=2),rep("uninfected",times=2)),
      replicate=1*grepl(pattern = "rep1",colnames(counts))+2*grepl(pattern = "rep2",colnames(counts)))
)

browseVignettes("DESeq2")
dds<-DESeqDataSetFromMatrix(countData=counts,colData=colData,design=~infection + isoform)
colData(dds)$infection<-factor(colData(dds)$infection,levels=c("uninfected","CT"))
colData(dds)$isoform<-factor(colData(dds)$isoform,levels=c("MED","HCL","DL","L","D"))
dds
dds<-DESeq(dds)
res<-results(dds)
(res<-res[order(res$padj),])
save(res,file="Infection_Isoform.DESeq2Results")
write.table(rownames(head(res,n=200)),file="Infection_Isoform.topgenes.DESeq2Results",sep="\t")
plotMA(dds,ylim=c(-2,2))

mcols(res,use.names=TRUE)
colData(dds)
resultsNames(dds)
results(dds,"isoform_vs_MED")

rld<-rlogTransformation(dds,blind=T)
print(plotPCA(rld,intgroup=c("replicate")))
plotDispEsts(dds)

library("RColorBrewer")
library("gplots")
select<-order(rowMeans(counts(dds,normalized=T)),decreasing = T)[1:30]
hmcol<-colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(assay(rld)[select,],col=hmcol,
          Rowv=F,Colv=F,scale="none",
          dendrogram="none",trace="none",margin=c(10,6))

distsRL<-dist(t(assay(rld)))
mat<-as.matrix(distsRL)
rownames(mat)<-colnames(mat)<-with(colData(dds),
                                   paste(pH,infection,isoform,sep=":"))
heatmap.2(mat,trace="none",col=rev(hmcol),margin=c(13,13))
