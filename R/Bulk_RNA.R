library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)
library(foreach)
library(parallel)
library(doMC)
library(GenomicAlignments)
library(hypeR)
library(goseq)
library(edgeR)
library(DESeq)
library(DESeq2)
library(ggplot2)
library(rgl)
library(sva)
library(factoextra)
library(Seurat)
library(GO.db)
library(UpSetR)
library(org.Hs.eg.db)
library(treemapify)
library(rrvgo)
library(LSAfun)
library(GOfuncR)

###functions

########### read bulk data in from Bam files

hg19<-makeTxDbFromUCSC("hg19","knownGene")
hg19e<-exons(hg19)
hg19eg<-exonsBy(hg19, "gene")
load("lincRNAbase")
RNA<-c(hg19eg,lincRNA)

x<-dir()

test<-foreach(i=1:length(x)) %do% {
y<-readGAlignments(x[1],format="BAM")
RNA.sample.ID<-summarizeOverlaps(RNA,y,mode="IntersectionNotEmpty")
counts.RNA.sample.ID<-assay(RNA.sample.ID)
list(counts.RNA.sample.ID,counts.exon.sample.ID)
}

RNA.cts<-foreach(i=1:length(test),.combine='cbind') %dopar% {
test[[i]][[1]]}
colnames(RNA.cts)<-unlist(strsplit(x,".bam"))

############annotation files
SYMBOLS<-unlist(as.list(org.Hs.egSYMBOL))

####address batch effects and generate normalized expression matrix
y<-RNA.cts[rowSums(RNA.cts)>10,]
apply(y,2,sum)/min(apply(y,2,sum))

# 228.2    229.2    230.2    233.2    228.3    229.3    230.3    233.3    228.4    229.4    230.4 
# 4.519919 3.895793 4.400003 3.757948 3.610122 3.204684 2.941372 2.070260 3.556956 3.573734 3.976852 
# 233.4    228.5    229.5    230.5    233.5  228.7.1    229.7    230.7    233.7    228.8    229.8 
# 1.000000 3.643296 2.920039 3.455075 2.584524 3.434788 2.222245 2.617154 3.737366 3.101630 4.465679 
# 230.8    233.8 
# 3.817528 3.317360

# sample 233.4 has much fewer reads than the other sample, q value=0.3, it is an outlier
#so will drop from analysis    

batch<-rep(1:4,6)
y<-y[,-12]
batch<-batch[-12]
dat<-ComBat_seq(y,batch=batch)

group<-c(rep(2,4),rep(3,4),rep(3,3),rep(5,4),rep(7,4),rep(8,4))
Sue<-newCountDataSet(dat,group)
x<-estimateSizeFactors(Sue)
norm.cts<-t(t(y)/sizeFactors(x))
colnames(norm.cts)<-c(colnames(norm.cts)[1:15],"228.7",colnames(norm.cts)[17:23])

#########get genes with highest variance
treat<-factor(group)

y<-DGEList(counts=dat,group=treat)
y<-calcNormFactors(y)

design<-model.matrix(~group)

y<-estimateGLMCommonDisp(y,design)
y<-estimateGLMTagwiseDisp(y,design,prior.df=45)
x<-rank(y$tagwise.dispersion)
x1<-rownames(y)[which(x>(max(x)-4000))]

######
x<-prcomp(dist(t(norm.cts[x1,])))
fviz_pca_var(x,habillage=batch,repel=T) + scale_color_manual(values=c("deeppink","deepskyblue2","gold1","green"))

######SVD of variant genes######

test<-norm.cts[intersect(rownames(dat),x1),]

x<-svd(test)
eigenassays <- x$u
eigenexpressions <- x$d
eigenfeatures <- t(x$v)
colnames(eigenfeatures) <- colnames(test)
rownames(eigenfeatures) <- c(1:min(nrow(test),ncol(test)))
fractions <- eigenexpressions^2/sum(eigenexpressions^2)
names(fractions) <- c(1:min(nrow(test),ncol(test)))
fractions
round(100*(-sum(fractions*log(fractions))/log(ncol(eigenassays))))/100

y<-eigenexpressions
y[1]<-0
y[15:23]<-0
test <- eigenassays %*% diag(y) %*% eigenfeatures
x<-svd(test)
eigenassays <- x$u
eigenexpressions <- x$d
eigenfeatures <- t(x$v)
colnames(eigenfeatures) <- colnames(test)
rownames(eigenfeatures) <- c(1:min(nrow(test),ncol(test)))
fractions <- eigenexpressions^2/sum(eigenexpressions^2)
names(fractions) <- c(1:min(nrow(test),ncol(test)))
fractions
round(100*(-sum(fractions*log(fractions))/log(ncol(eigenassays))))/100
barplot(fractions)
assaycorrelations<- diag(eigenexpressions) %*% eigenfeatures
featurecorrelations <- t(eigenassays %*% diag(eigenexpressions))
Time4000_svd<-x

#####plot svd in 3D (for figure 1C)

ac<-assaycorrelations
fc<-featurecorrelations
colnames(fc)<-x1
open3d()
z<-ac[1,]
x<-ac[4,]
y<-ac[3,]

color<-rep(c("deeppink","deepskyblue2","gold1","green"),6)
color<-color[-12]

plot3d(x,y,z,col=color,type='s')

rgl.postscript("SVD_4000_color_celllines.pdf",fmt="pdf")

color<-c(rep("red",4),rep("green3",4),rep("deepskyblue",3),rep("purple",4),
         rep("brown",4),rep("gold2",4))

plot3d(x,y,z,col=color,type='s')
identify3d(x,y,z,labels=colnames(norm.cts),adj=c(-0.65,0.5))
rgl.postscript("SVD_4000_color_time.pdf",fmt="pdf")

### No apparent correlation with time, so switch to a xplant based strategy. Due to avaiable xplant data going back to raw expression data 
### No xplant datas for line 233 so dropped all 233 data from analysis 

xplanted.cts<-cbind(RNA.cts[,1:3],RNA.cts[,20],RNA.cts[,24],RNA.cts[,6],RNA.cts[,9:11],RNA.cts[,13:14],RNA.cts[,17],
                    RNA.cts[,19],RNA.cts[,23])

colnames(xplanted.cts)<-c("228.2","229.2","230.2","230.7","230.8","229.3","228.4",
                          "229.4","230.4","228.5","229.5","228.7","229.7","229.8")

group<-c("W2","W2","W2","NE","NE","P90","P90","P90","P90","P90","P90","P90","P90","P90")

Sue<-newCountDataSet(xplanted.cts,group)
x<-estimateSizeFactors(Sue)
norm.xplant<-t(t(xplanted.cts)/sizeFactors(x))
biplot(princomp(dist(t(norm.xplant))))

ave.cts<-data.frame(apply(norm.xplant[,1:3],1,mean),apply(norm.xplant[,4:5],1,mean),apply(cbind(norm.xplant[,6:14]),1,mean))
colnames(ave.cts)<-c("W2","NonEff","Eff")

#########get genes with highest variance
treat<-factor(c("W2","W2","W2","NE","NE","P90","P90","P90","P90","P90","P90","P90","P90","P90"))

y<-DGEList(counts=xplanted.cts,group=treat)
y<-calcNormFactors(y)

group<-factor(c(1,1,1,2,2,3,3,3,3,3,3,3,3,3))
design<-model.matrix(~group)

y<-estimateGLMCommonDisp(y,design)
y<-estimateGLMTagwiseDisp(y,design,prior.df=45)
x<-rank(y$tagwise.dispersion)
x1<-rownames(y$counts)[which(x>(max(x)-4000))]

test<-norm.xplant[x1,]

######SVD of variant genes######


x<-svd(test)
eigenassays <- x$u
eigenexpressions <- x$d
eigenfeatures <- t(x$v)
colnames(eigenfeatures) <- colnames(test)
rownames(eigenfeatures) <- c(1:min(nrow(test),ncol(test)))
fractions <- eigenexpressions^2/sum(eigenexpressions^2)
names(fractions) <- c(1:min(nrow(test),ncol(test)))
fractions
round(100*(-sum(fractions*log(fractions))/log(ncol(eigenassays))))/100

y<-eigenexpressions
y[1]<-0
y[14]<-0
test <- eigenassays %*% diag(y) %*% eigenfeatures
x<-svd(test)
eigenassays <- x$u
eigenexpressions <- x$d
eigenfeatures <- t(x$v)
colnames(eigenfeatures) <- colnames(test)
rownames(eigenfeatures) <- c(1:min(nrow(test),ncol(test)))
fractions <- eigenexpressions^2/sum(eigenexpressions^2)
names(fractions) <- c(1:min(nrow(test),ncol(test)))
fractions
round(100*(-sum(fractions*log(fractions))/log(ncol(eigenassays))))/100
barplot(fractions)
assaycorrelations<- diag(eigenexpressions) %*% eigenfeatures
featurecorrelations <- t(eigenassays %*% diag(eigenexpressions))
xplant_svd<-x

#####plot svd in 3D (for figure 4D)

ac<-assaycorrelations
fc<-featurecorrelations
colnames(fc)<-x1
z<-ac[1,]
x<-ac[2,]
y<-ac[3,]

color<-c("deeppink","deeppink","deeppink","gold1","gold1","deepskyblue2","deepskyblue2",
         "deepskyblue2","deepskyblue2","deepskyblue2","deepskyblue2","deepskyblue2",
         "deepskyblue2","deepskyblue2")
open3d()
plot3d(x,y,z,col=color,type='s')
identify3d(x,y,z,labels=names(z),adj=c(-0.65,0.5))
rgl.postscript("xplant_4000.pdf",fmt="pdf")

### data does separate based on transplant status. Move forward with this reduced dataset
### identify differentially expressed genes using DESeq2 and EdgeR
### DESeq
group<-c("W2","W2","W2","NE","NE","P90","P90","P90","P90","P90","P90","P90","P90","P90")
group<-data.frame(group)
Sue<-DESeqDataSetFromMatrix(xplanted.cts,colData=group,design=~group)
colData(Sue)$group<-factor(colData(Sue)$group,levels=c("W2","NE","P90"))
colnames(Sue)<-colnames(xplanted.cts)
Sue<-estimateSizeFactors(Sue)
Sue<-estimateDispersions(Sue)
Sue<-nbinomLRT(Sue,reduced=~1)
x<-as.matrix(results(Sue))
results.table.DEseq<-x
x<-x[which(x[,6]<=0.1),]
DEseq.genes<-data.frame(xplanted.cts[rownames(x),],x)
colnames(DEseq.genes)<-c(colnames(xplanted.cts),colnames(x))

### edgeR
treat<-factor(c("W2","W2","W2","NE","NE","P90","P90","P90","P90","P90","P90","P90","P90","P90"))

y<-DGEList(counts=xplanted.cts,group=treat)
y<-calcNormFactors(y)

group<-factor(c(1,1,1,2,2,3,3,3,3,3,3,3,3,3))
design<-model.matrix(~group)

y<-estimateGLMCommonDisp(y,design)
y<-estimateGLMTagwiseDisp(y,design,prior.df=45)
fit<-glmFit(y,design)
lrt<-glmLRT(fit)
topTags(lrt)
summary(dt<-decideTestsDGE(lrt))
y<-as.matrix(lrt$table)
FDR<-p.adjust(y[,4],"fdr")
y<-cbind(y,FDR)
results.table.edgeR<-y
x<-which(FDR<=0.1)
z<-y[names(x),]
edgeR.genes<-cbind(norm.xplant[rownames(z),],z)

### identify genes significant in both approaches

x<-intersect(rownames(edgeR.genes),rownames(DEseq.genes))
sig.genes<-data.frame(norm.xplant[x,],edgeR.genes[x,18],DEseq.genes[x,19])
ave.cts<-data.frame(apply(norm.xplant[,1:3],1,mean),apply(norm.xplant[,4:5],1,mean),apply(cbind(norm.xplant[,6:14]),1,mean))
colnames(ave.cts)<-c("W2","NE","P90")
FC<-apply(ave.cts,1,max)/apply(ave.cts,1,min)
y<-FC[which(FC>=2)]
y<-intersect(names(y),x)
sig.genes<-data.frame(sig.genes[y,],FC[y])
colnames(sig.genes)<-c(colnames(norm.xplant),"edgeR FDR","DEseq FDR","Fold Change")

#### Time of Maximum (TOM) expression
x<-ave.cts[rownames(sig.genes),]
y<-apply(x,1,function(z) colnames(x)[which(z==max(z))])
table(y)
# Eff NonEff     W2 
# 232    291    942

### Lnc-RNA TOM
x<-x[grep("TCONS",rownames(x)),]
y<-apply(x,1,function(z) colnames(x)[which(z==max(z))])
table(y)
# Eff NonEff     W2 
# 62     73     70 

#######clean up####
rm(dt)
rm(design)
rm(FDR)
rm(fit)
rm(group)
rm(lrt)
rm(Sue)
rm(treat)
rm(x)
rm(y)
rm(z)
                 
### Enrichment analysis
### GO:BP
GObp<-getgo(names(hg19eg),"hg19","knownGene","GO:BP")
genlen<-genelengths[sort(names(hg19eg))]
x<-intersect(rownames(sig.genes),names(hg19eg))
z<-setdiff(names(hg19eg),x)
zs<-rep(0,length(z))
xs<-rep(1,length(x))
names(zs)<-z
names(xs)<-x
zt<-c(zs,xs)[sort(names(hg19eg))]

pwf<-nullp(zt,bias.data=genlen,plot.fit=FALSE)
blue<-goseq(pwf,gene2cat=GObp)
FDR<-p.adjust(blue[,2],"fdr")
RPE.go.bp<-data.frame(blue[which(FDR<=0.1),],FDR[which(FDR<=0.1)])
colnames(RPE.go.bp)<-c(colnames(blue),"FDR")

### GO:MF
GOmf<-getgo(names(hg19eg),"hg19","knownGene","GO:MF")
blue<-goseq(pwf,gene2cat=GOmf)
FDR<-p.adjust(blue[,2],"fdr")
RPE.go.mf<-data.frame(blue[which(FDR<=0.1),],FDR[which(FDR<=0.1)])
colnames(RPE.go.mf)<-c(colnames(blue),"FDR")
         
### REACTOME Pathways
REACTOME <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
REACTOME<-REACTOME$genesets
x<-rownames(pwf)
x<-SYMBOLS[x]
pwf<-pwf[!(is.na(x)),]
rownames(pwf)<-x[!(is.na(x))]
blue<-goseq(pwf,gene2cat=REACTOME)
FDR<-p.adjust(blue[,2],"fdr")
RPE_reactome<-data.frame(blue[which(FDR<=0.1),],FDR[which(FDR<=0.1)])
colnames(RPE_reactome)<-c(colnames(blue),"FDR")
                
### make enrichment figures (for Figure 1 H&I)
sizes<-(-(log(RPE.go.bp$FDR)))
names(sizes)<-RPE.go.bp$category
sm<-calculateSimMatrix(RPE.go.bp$category,orgdb="org.Hs.eg.db",ont="BP",method="Rel")
x<-intersect(y,rownames(sm))
sizes<-sizes[x]
o <- order(sizes, decreasing=T, na.last = FALSE)
sm <- sm[o, o]
cluster <- cutree(hclust(as.dist(1 - sm)), h = 0.9)
clusterRep <- tapply(rownames(sm), cluster, function(x) x[which.max(sizes[x])])
red<-data.frame(go = rownames(sm), 
                cluster = cluster, 
                parent = clusterRep[cluster],
                parentSimScore = unlist(Map(seq_len(nrow(sm)), 
                                            clusterRep[cluster], f = function(i, j) sm[i,j])), 
                size = sizes[match(rownames(sm), names(sizes))],
                term = rrvgo:::getGoTerm(rownames(sm)), 
                parentTerm = rrvgo:::getGoTerm(clusterRep[cluster]))

ncolors<-gg_color_hue(length(unique(red$parentTerm)))
mcolors<-red$parentTerm
for(i in 1:length(ncolors)){
  mcolors[which(mcolors==unique(red$parentTerm)[i])]<-ncolors[i]
}
red$mcolors<-mcolors

p<-ggplot(red, aes(area=size,subgroup=parentTerm,fill=mcolors))
p<-p + geom_treemap()
p<-p + geom_treemap_subgroup_border()
p<-p + geom_treemap_subgroup_text(place='center',grow=F,alpha=1,min.size=6,reflow=T)
p<-p + theme(legend.position = "none")
p

RPE_reactome$logFDR<-(-(log(RPE_reactome$FDR)))
RPE_reactome$percent<-RPE_reactome$numDEInCat/RPE_reactome$numInCat
RPE_reactome<-RPE_reactome[order(RPE_reactome$percent, decreasing=T),]
x<-RPE_reactome[RPE_reactome$numInCat>30,]
x<-x[order(x$logFDR,decreasing=T),]
y<-unlist(strsplit(x$category,"REACTOME_"))
y<-y[seq(2,length(y),by=2)]
x$pretty<-y
z<-x[c(4:5,7,10,20,21:23,27,42),]
p<-ggplot(z,aes(pretty,logFDR))
p<-p+geom_bar(stat="identity")
p<-p + theme_bw()
p
#######clean up and save
rm(test)
rm(sm)
rm(dt)
rm(design)
rm(FDR)
rm(fit)
rm(group)
rm(lrt)
rm(Sue)
rm(treat)
rm(x)
rm(y)
rm(z)
rm(p)
rm(ac)
rm(fc)
rm(cluster)
rm(clusterRep)
rm(color)
rm(mcolors)
rm(ncolors)
rm(o)
rm(i)
rm(sizes)
rm(red)
rm(x1)
rm(fractions)

         
         
         
