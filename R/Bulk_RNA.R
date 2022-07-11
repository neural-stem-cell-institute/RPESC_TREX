library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)
library(foreach)
library(parallel)
library(doMC)
library(GenomicAlignments)
registerDoMC(cores=2)
###########

hg19<-makeTranscriptDbFromUCSC("hg19","knownGene")
hg19e<-exons(hg19)
hg19eg<-exonsBy(hg19, "gene")
load("lincRNAbase")
RNA<-c(hg19eg,lincRNA)

x<-dir()
x<-x[3:27]

test<-foreach(i=1:length(x)) %do% {
y<-readGAlignments(x[1],format="BAM")
exon.sample.ID<-summarizeOverlaps(hg19e,y,mode="IntersectionNotEmpty")
RNA.sample.ID<-summarizeOverlaps(RNA,y,mode="IntersectionNotEmpty")
counts.exon.sample.ID<-assay(exon.sample.ID)
counts.RNA.sample.ID<-assay(RNA.sample.ID)
list(counts.RNA.sample.ID,counts.exon.sample.ID)
}


RNA.cts<-foreach(i=1:length(test),.combine='cbind') %dopar% {
test[[i]][[1]]}
colnames(RNA.cts)<-unlist(strsplit(x,".bam"))

exon.cts<-foreach(i=1:length(test),.combine='cbind') %dopar% {
test[[i]][[2]]}

colnames(exon.cts)<-unlist(strsplit(x,".bam"))
rownames(exon.cts)<-hg19e$exon_id

