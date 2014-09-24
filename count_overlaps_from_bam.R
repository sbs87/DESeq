# This script takes alignment files (bam) and creates counts tables as experimentClass that can be subsequently used by DESeq
# Steven Smith
# Setember 21, 2014
# The summarizeOverlaps step takes up a tremendous amount of RAM and processing power. This needs to be run on the grid.

### Set up enviornment
rm(list=ls())
#source("http://bioconductor.org/biocLite.R")
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

#browseVignettes("DESeq2")

dds<-DESeqDataSetFromMatrix(countData=counts,colData=colData,design=~isoform +pH | infection)
colData(dds)$infection<-factor(colData(dds)$infection,levels=c("uninfected","CT"))
colData(dds)$isoform<-factor(colData(dds)$isoform,levels=c("MED","HCL","DL","L","D"))
dds
dds<-DESeq(dds)
res<-results(dds)
(res<-res[order(res$padj),])
save(res,file="Infection_Isoform.DESeq2Results")
load("Infection_Isoform.DESeq2Results")
(sig.genes<-res[(abs(res$log2FoldChange)>1.5 & res$padj<0.01 & !is.na(res$padj)),])
res[(abs(res$log2FoldChange)<1.5 & res$padj>0.01 & !is.na(res$padj)),]
write.table(row.names(sig.genes),quote=F,sep="\t",file="signifigant_genes.txt",row.names = F,col.names = F)
res.counts<-assays(dds)$counts
#(top.counts<-data.frame(res.counts[rownames(res.counts) %in% rownames(head(res)),]))
#names(top.counts)<-colData[,1]

library(ggplot2)
ENS_hugo_map<-read.table("ENS_hugo_map.txt",sep="\t",header=T)
for(i in 1:length(ENS_hugo_map$Gene)){
  Ens.id<-as.character(ENS_hugo_map$Ensembl_Gene[i])
  Ens.id<-"ENSG00000135414"
  Hugo.id<-as.character(ENS_hugo_map[ENS_hugo_map$Ensembl_Gene==Ens.id,]$Gene)
  current.title<-paste(Ens.id,"_",Hugo.id,".ps",sep="")
  pick_gene<-data.frame(counts=res.counts[rownames(res.counts)==Ens.id,],colData(dds))
  rownames(pick_gene)<-as.character(colData(dds)$sample_name)
  pick_gene$isoform<-factor(pick_gene$isoform,levels = c("HCL","DL","L","D","MED"),ordered=T)
  pick_gene$pH<-factor(pick_gene$pH,levels = c("4","7","MED"),ordered=T)
  postscript(file=current.title,horizontal = T)
  f<-ggplot(pick_gene)+geom_bar(aes(x=infection,y=counts,fill=isoform),stat="identity") +facet_wrap(~pH+isoform,ncol=4)+ggtitle(current.title)+theme_bw()
  plot(f)
  dev.off()

}
#top.counts
plot(res$log2FoldChange, -log(res$padj,base = 10),main="Volcano Plot of DESeq2 Results, isoform + infection",xlab="log2(FoldChange)",ylab="-log10(Adjusted p-value)")

abline(h = -log(0.01,base=10),col='red')
abline(v=-1.5,col='red')
abline(v=1.5,col='red')

plotMA(dds,ylim=c(-2,2))

mcols(res,use.names=TRUE)
colData(dds)
resultsNames(dds)


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

library(DESeq)
library(edgeR)
edgeRUsersGuide()

counts<-assays(dds)$counts
counts[counts=="ENSG00000006652"]

res
