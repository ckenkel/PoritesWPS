#To download DESeq package (you can comment these lines out, they only need to be run once ever)
BiocManager::install("pathview")
# biocLite("DESeq2")

# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")

setwd("/Users/carlykenkel/Dropbox/AIMSpostdoc/PoritesPatch")

 #To upload packages - you need to do this every time you open R
library(DESeq)
library(gplots) # for venn diagram
library(ggplot2)
library(RColorBrewer)
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(plotrix)
library(reshape2)
library(factoextra)
library(KOGMWU)


counts=read.table("AllCountsHost.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

head(counts) 


length(counts[,1])  #25222 for gmapper on Plob+C15 combined

names(counts)


readsleft=c()
for (column in names(counts)) {
	val=sum(counts[,column])
	readsleft=append(readsleft,val)}

RLtable=data.frame(cbind(names(counts),readsleft))
RLtable$readsleft=as.numeric(as.character(RLtable$readsleft))
write.csv(RLtable,"readsleft.csv",quote=F) 

#######################Creating table of conditions for your experiment

patch=c(1:length(names(counts)))
patch[grep("B",names(counts))]="Patch"
patch[grep("N",names(counts))]="Normal"
geno=c("X1","X1","X2","X2","X3","X3","X4","X4","X5","X5","X6","X6","X7","X7","X8","X8")
pc1=c("-0.8745566","0.3957006","-1.1903046","2.8307584","-1.6972843","1.9503512","-1.9548308","1.7646005","-1.8041639","0.7985436","-0.9543884","3.5344855","-1.7413601","0.3277427","-1.8828452","0.4975513")


conditions=data.frame(patch,geno,pc1)
head(conditions)



############################################original DESeq methods

real=newCountDataSet(counts.nobad1,conditions.nobad1) 
real=estimateSizeFactors(real)

# ####all the data you ever wanted about quality control

cds=estimateDispersions(real,method="blind")
vsdBlind=varianceStabilizingTransformation(cds)

v="/Users/drcarl/Dropbox/AIMSpostdoc/PoritesPatch/AQM2"

arrayQualityMetrics(vsdBlind,outdir=v,intgroup=c("geno"),force=TRUE) #check .html output file in new folder

###############Remove outliers as detected; repeat arrayQualityMetrics above after regenerating newCountDataSet
###############to confirm that all outliers have been removed 



counts.nobad1=counts[,-c(7,8)] #According to heatmap AQMv1, Sample 4N is serious outlier; 4B questionable. Remove this geno rep and try re-running 
conditions.nobad1=conditions[-c(7,8),]
conditions.nobad1$geno<-factor(conditions.nobad1$geno)

#Remainder are fine

######################
### DESEQ HOST #######
######################

library(DESeq)

real=newCountDataSet(counts.nobad1,conditions.nobad1) 
real=estimateSizeFactors(real) #library size (~total mapped reads)

sizeFactors(real)

real=estimateDispersions(real,method="pooled-CR",sharingMode="gene-est-only",modelFormula=count~geno+patch)  
#can use gene-est-only, IF have >7 reps per treatment combo
# quartz()
plotDispEsts(real)
plotDispEsts(dds) #looks reasonable...but may not be same as regular deseq?

vsd=getVarianceStabilizedData(real)
# write.csv(vsd, file="VSD_allsymgenes_nobadsam_Jan2015.csv", quote=F)

# ######################Determining quality filtering cutoffs

fit0=fitNbinomGLMs(real, count ~ 1) # null model: expression does not depend on anything
fit1=fitNbinomGLMs(real, count ~ geno)
fit2=fitNbinomGLMs(real, count ~ patch)
fit3=fitNbinomGLMs(real, count ~ geno+patch)
fit4=fitNbinomGLMs(real, count ~ geno+pc1)

pvals.g<-nbinomGLMTest(fit1,fit0) #testing significance of geno
pvals.p<-nbinomGLMTest(fit2,fit0) #testing significance of patch
pvals.pg<-nbinomGLMTest(fit3,fit1) #testing significance of patch on top of geno
pvals.pc1<-nbinomGLMTest(fit4,fit1) #testing significance of patch on top of gen


pvalue=pvals.pc1 #change to each type of pval: i, t and o
theta=seq(from=0,to=0.8,by=0.02)

filterChoices=data.frame(`mean`=rowMeans(counts(real)),`median`=apply((counts(real)),1,median),`min`=rowMin(counts(real)),`max`=rowMax(counts(real)),`sd`=rowSds(counts(real)))
rejChoices=sapply(filterChoices,function(f) filtered_R(alpha=0.1,filter=f,test=pvalue,theta=theta,method="BH"))
library("RColorBrewer")
myColours=brewer.pal(ncol(filterChoices),"Set1")

# #quartz()
# #windows()
matplot(theta,rejChoices,type="l",lty=1,col=myColours,lwd=2,xlab=expression(theta),ylab="number of rejections")
legend("bottomleft",legend=colnames(filterChoices),fill=myColours)

# #look for peak in graph - corresponds to correct theta and best-fit line for which metric to use - pick best theta for all tests
# #p and g peak at 0.1, pc1 is much higher. lets go with it. Increased abundance of read data might be improving inferrence.

# #######################Quality Filtering Data based on theta - get rid of genes with low variance

# #FOR host
rs=rowSds(counts(real)) #using standard deviation as quality filtering metric based on analyses above
theta=0.1 
use=(rs>quantile(rs,probs=theta)) ###
table(use) 
# use
# FALSE  TRUE 
 # 2592 22630


realFilt=real[use,]
vsd=getVarianceStabilizedData(realFilt)

######################## Now for the real Model Testing

fit0=fitNbinomGLMs(realFilt, count ~ 1) # null model: expression does not depend on anything
fit1=fitNbinomGLMs(realFilt, count ~ geno)
fit2=fitNbinomGLMs(realFilt, count ~ patch)
fit3=fitNbinomGLMs(realFilt, count ~ geno+patch)
fit4=fitNbinomGLMs(realFilt, count ~ geno+pc1)

# testing section
pvals.p<-nbinomGLMTest(fit3,fit1) #testing significance of patch, after accounting for genotype
pvals.g<-nbinomGLMTest(fit1,fit0) #testing significance of genotype
pvals.pc1<-nbinomGLMTest(fit4,fit1) #testing significance of pc1

#making non-convergent model p-values NA's - only occur in pvals.p
pvals.p=data.frame(pvals.p)
rownames(fit3)->rownames(pvals.p)
badmods=subset(fit3,(!fit3[,11]))
for ( i in rownames(badmods)){pvals.p[i,1]<-NA}


summary(pvals.p)

pvals.g=data.frame(pvals.g)
rownames(fit1)->rownames(pvals.g)
badmods=subset(fit1,(!fit1[,10]))
for ( i in rownames(badmods)){pvals.g[i,1]<-NA}

summary(pvals.g)

pvals.pc1=data.frame(pvals.pc1)
rownames(fit4)->rownames(pvals.pc1)
badmods=subset(fit4,(!fit4[,24])) #18
for ( i in rownames(badmods)){pvals.pc1[i,1]<-NA}

summary(pvals.pc1)

#multiple test correction - adjust p-values using Benjamini-Hochburg

adjp.p<-p.adjust(pvals.p$pvals.p,method="BH")
adjp.g<-p.adjust(pvals.g$pvals.g,method="BH")
adjp.pc1<-p.adjust(pvals.pc1$pvals.pc1,method="BH")

do<-(cbind(vsd, "fitPatch" = fit3$patchPatch, "adjp.p" = adjp.p,"pval.p" = pvals.p$pvals.p,"adjp.g" = adjp.g,"pval.g" = pvals.g$pvals.g,"adjp.pc1" = adjp.pc1,"pval.pc1" = pvals.pc1$pvals.pc1)) #creating table of all multiple test corrected p-values with variance stabilized count data 


write.csv(do, file="hostVSDandPVALS_no_g4_deseq1_4jun_plusPC1.csv", quote=F) #writing an output file of vsd plus p-values

########################## counting, venn diagram:
p<-read.csv("hostVSDandPVALS_no_g4_deseq1_4jun_plusPC1.csv")

#p<-data.frame(do)
patch=row.names(p[p$adjp.p<=0.1 & !is.na(p$adjp.p),])
geno=row.names(p[p$adjp.g<=0.1 & !is.na(p$adjp.g),])
pc1=row.names(p[p$adjp.pc1<=0.1 & !is.na(p$adjp.pc1),])


candidates=list("Patch, P<=0.1"=patch,"Geno, P<=0.1"=geno,"PC1, P<=0.1"=pc1)
library(gplots)
quartz()
venn(candidates)

################################
###### now for other visualizations; sorting out of 'interesting genes'
#################################

d2<-read.csv("hostVSDandPVALS_no_g4_deseq1_4jun_plusPC1.csv")
rownames(d2)<-d2$X

head(d2)

summary(d2) #22630 genes host

d2_cc<-d2[!is.na(d2$adjp.p),] #22630 genes host



GEs<-d2 #rename the dataset

patch=head(GEs[order(GEs$adjp.p),],topnum)
geno=head(GEs[order(GEs$adjp.g),],topnum)
sig=data.frame(rbind(patch[!(patch$X %in% geno$X),],geno)) #remove redundant isogroups that are sig for both
length(sig[,1]) 

#Host #Top 100 (197 non-overlapping); 500 (957); 1000 (1889); 2000 (3643)


#Transpose expression values only SYM
tsig=t(sig[,2:13]) 
tsig[1:12,1:5]

#Transpose expression values only HOST
tsig=t(sig[,2:15]) 
tsig[1:14,1:5]


# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE.
#library(devtools)

 
pca <- prcomp(tsig,
      center = TRUE,
      scale. = TRUE) 
summary(pca) #First 2 PCs explain 61% of the variation
str(pca)

scores=pca$x

scores[,1:2]

#par(mfrow=c(1,4))

#FOR HOSTS
plot(scores[,1], scores[,2], xlab="PCA 1", ylab="PCA 2",type="n",main="Patch/Geno Top 500")
points(scores[c(1,3,5,7,9,11,13),1],scores[c(1,3,5,7,9,11,13),2],pch=1,col="lightgreen")
points(scores[c(2,4,6,8,10,12,14),1],scores[c(2,4,6,8,10,12,14),2],pch=19,col="darkgreen")

plot(scores[,1], scores[,2], xlab="PCA 1", ylab="PCA 2",type="n",main="Geno")
points(scores[c(1:2),1],scores[c(1:2),2],pch=1,col="pink")
points(scores[c(3:4),1],scores[c(3:4),2],pch=1,col="violet")
points(scores[c(5:6),1],scores[c(5:6),2],pch=1,col="green")
points(scores[c(7:8),1],scores[c(7:8),2],pch=1,col="blue")
points(scores[c(9:10),1],scores[c(9:10),2],pch=1,col="red")
points(scores[c(11:12),1],scores[c(11:12),2],pch=1,col="orange")
points(scores[c(13:14),1],scores[c(13:14),2],pch=1,col="black")

#Picking out top 15 genes in each of the patch and geno categories, shows that patch
#is stronger driver of expression in this instance


#########Generating directional GO output for DESeq results - signed log p-values

#patch
d2$direction=ifelse(d2$fitPatch>0,1,-1) #red=higher expression in patch
d2$pval<-(-log((d2$pval.p+0.0000000001),10))*(d2$direction)

patch<-cbind("gene"=rownames(d2),"pval"=d2$pval) 

head(patch)

patch=patch[complete.cases(patch[,2]),]

write.csv(patch,file="GOpatchHost.csv",quote=F,row.names=F)

##########################################Categorical GO files for gomwu scripts

#For all DEGs regardless of direction of expression
d2$binary=ifelse(d2$adjp.p<=0.1 & !is.na(d2$adjp.p),1,0)
d2$binary=as.factor(d2$binary)
summary(d2)
names(d2)
GObinary<-d2[,c(1,25)]
write.csv(GObinary,file="GObinaryHostPatch.csv",quote=F,row.names=F)

#And a vsd file for just this "interesting" gene subset
#for all DEGs
VSDsigOnly<-d2[d2$binary==1,]
nrow(VSDsigOnly)
head(VSDsigOnly)
VSDsigOnly$direction=as.factor(VSDsigOnly$direction)
summary(VSDsigOnly) #For genes sig by patch in host, 320 downregulated, 255 upregulated out of 575
#For genes sig in sym, 4097 downregulated, 133 upregulated...not all downregulated, but still most of DE could be due to bleaching?

write.csv(VSDsigOnly,file="VSDs_GObinaryHostPatch.csv",quote=F,row.names=F)

#For PC1
#For all DEGs regardless of direction of expression
d2$binary=ifelse(d2$adjp.pc1<=0.1 & !is.na(d2$adjp.pc1),1,0)
d2$binary=as.factor(d2$binary)
summary(d2)
names(d2)
GObinary<-d2[,c(1,25)]
write.csv(GObinary,file="GObinaryHostPC1.csv",quote=F,row.names=F)

#And a vsd file for just this "interesting" gene subset
#for all DEGs
VSDsigOnly<-d2[d2$binary==1,]
nrow(VSDsigOnly)
head(VSDsigOnly)
VSDsigOnly$direction=as.factor(VSDsigOnly$direction)
summary(VSDsigOnly) 
write.csv(VSDsigOnly,file="VSDs_GObinaryHostPC1.csv",quote=F,row.names=F)


#For upreg DEGs 
d2_cc$binary=ifelse(d2_cc$adjp.p<=0.1 & d2_cc$fitPatch>0,1,0)
d2_cc$binary=as.factor(d2_cc$binary)
summary(d2_cc)
names(d2_cc)
GObinary<-d2_cc[,c(1,25)]
write.csv(GObinary,file="GObinaryUpregHostPatch.csv",quote=F,row.names=F)

VSDsigOnly<-d2_cc[d2_cc$binary==1,]
nrow(VSDsigOnly)
head(VSDsigOnly)

write.csv(VSDsigOnly,file="VSDs_GObinaryUpregHostPatch.csv",quote=F,row.names=F)

#For downreg DEGs 
d2_cc$binary=ifelse(d2_cc$adjp.p<=0.1 & d2_cc$fitPatch<0,1,0)
d2_cc$binary=as.factor(d2_cc$binary)
summary(d2_cc)
names(d2_cc)
GObinary<-d2_cc[,c(1,25)]
write.csv(GObinary,file="GObinaryDownregHostPatch.csv",quote=F,row.names=F)



########################KOG enrichments.
#####There's no reason to plot hierarchically as in GO, as pathways don't have a nested structure

GEh<-read.csv("GOpatchHost.csv")

#gene2kegg<-read.table("C15_iso2kegg.tab",sep="\t")
gene2kogh<-read.table("plob_iso2kogClassNR.tab",sep="\t")

host_kog_mwu=kog.mwu(GEh,gene2kogh)
host_kog_mwu

#Sig terms for WPS
                                                            # term nseqs delta.rank         pval         padj
# 10                              Energy production and conversion   289      -1254 1.538142e-12 3.537726e-11
# 7                                                   Cytoskeleton   387        694 6.539071e-06 7.519931e-05
# 6                                RNA processing and modification   451        509 3.702805e-04 2.838817e-03
# 13    Cell cycle control, cell division, chromosome partitioning   181        700 1.686444e-03 9.697053e-03
# 9                Translation, ribosomal structure and biogenesis   370        443 4.845646e-03 2.228997e-02
# 19                              Chromatin structure and dynamics   121        698 1.021416e-02 3.915430e-02

#Try comparisson with short-term heat-stress induced bleaching
data(adults.3dHeat.logFoldChange)
data(larvae.longTerm)
gene2kog<-read.table("ComparativeGEstudies/amil_iso2kogClassNR.tab",sep="\t")
alfc.lth=kog.mwu(adults.3dHeat.logFoldChange,gene2kog)
alfc.lth
# aposymbiotic coral larvae response to 5-day heat stress:
l.lth=kog.mwu(larvae.longTerm,gene2kog)
l.lth

# Load in additional data 
load("ComparativeGEstudies/MetaAnalysisFiles.RData")
# adult response to 95-day heat stress:
sid.lth=kog.mwu(Ssid.tempStatus,Ssid_gene2kog)
sid.lth

# symbiotic state in Aiptasia
aip.stat=kog.mwu(Aiptasia.symStatus,Aip_gene2kog)
aip.stat

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Amil.adult.3d"=alfc.lth,"Amil.larvae.5d"=l.lth,"Ssid.adult.95d"=sid.lth,"Aip.sym.status"=aip.stat,"Plob.WPS"=host_kog_mwu))

col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.75)(100)

# Making a heatmap with hierarchical clustering trees:
pheatmap(as.matrix(ktable),color=col,clustering_distance_cols="correlation")
# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)
# plotting individual delta-rank correlations:
par(mfrow=c(3,1))
corrPlot(x="Amil.adult.3d",y="Amil.larvae.5d",ktable,ylim=c(-1000,500))
corrPlot(x="Amil.adult.3d",y="Plob.WPS",ktable,ylim=c(-1500,1000))
corrPlot(x="Plob.WPS",y="Amil.larvae.5d",ktable,ylim=c(-1000,500),xlim=c(-1500,1000))



####################Now for some heatmaps - lets look at what genes annotated with interesting terms are 
gg=read.table("plob_iso2gene.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)  #the iso2gene file for your species - plob_iso2gene.tab C15_iso2gene.tab


kog=read.table("plob_iso2kogClassNR.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE) #the iso2kog annotation file for your spp - plob_iso2kogClassNR.tab
rownames(kog)<-kog$V1

expr=read.csv("VSDs_GObinaryHostPatch.csv") #expression data for your subset #VSDs_GObinarySymPatch_NoOutlierGEO.csv
edata=c(2:15) #only columns with expression data - host 2:15 ; sym 2:17

dfull=read.csv("hostVSDandPVALS_no_g4_deseq1_4jun_plusPC1.csv") #expression data for full dataset (for sanity check)
# hostVSDandPVALS_no_g4_deseq1_4jun.csv symVSDandPVALS_deseq1_4jun.csv

term=subset(kog,V2=="Cell cycle control, cell division, chromosome partitioning") #write term of interest here from your sig list#
is=rownames(term)
#is<-pc1

#####################loop through genes matching term in vsd subset
sel=c();gnms=c()
for ( i in is){
	if (i %in% expr$X){
		sel=rbind(sel,expr[expr$X==i,])
		gnms=append(gnms,substr(paste(as.character(gg$V2[gg$V1==i]),collapse="."),1,75))
	}
}
row.names(sel)=paste(gnms,sel$X,sep=".")
nrow(sel)
rownames(sel)


exp=sel[,edata]
if (length(exp[,1])==1) { exp=rbind(exp,exp) }
nrow(exp)
#will return number of genes with this annotation in sig expr dataset, can compare to global list to recreate pval
rownames(exp)


#####################then through genes matching GO term in whole dataset
sel2=c();gnms2=c()
for ( i in is){
	if (i %in% dfull$X){
		sel2=rbind(sel2,dfull[dfull$X==i,])
		gnms2=append(gnms2,substr(paste(as.character(gg$V2[gg$V1==i]),collapse="."),1,50))
	}
}
row.names(sel2)=paste(gnms2,sel2$X,sep=".")

exp2=sel2[,edata]
if (length(exp2[,1])==1) { exp2=rbind(exp2,exp2) }

########return number of genes matching GO terms in each
a=nrow(exp) #genes with GO annotation in specified data subset
b=nrow(exp2) #genes with GO annotation in entire dataset
c=(nrow(expr)-a) #genes without GO annotation in specified data subset
f=(nrow(dfull)-b) #genes without GO annotation in entire dataset

###SANITY check!  do the fisher test just for fun
		sig=c(a,c) # number of significant genes belonging and not belonging to the tested GO category
		ns=c(b,f) # number of not-significant genes belonging and not belonging to the tested GO category
		mm=matrix(c(sig,ns),nrow=2,dimnames=list(ns=c("kog","notkog"),sig=c("kog","notkog")))
		ff=fisher.test(mm,alternative="greater")
ff #should be a significant value for categorical KOG runs...
b #number of terms in the whole dataset that are annotated with the KOG term (should match ktable above)
################################


means=apply(exp,1,mean) # means of rows
expc=exp-means #rescale expression data so it's up and down relative to mean

#reorder columns to match healthy/wps
expc<-expc[,c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)] #host
 
library(RColorBrewer)
library(pheatmap)

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=1.1)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

#quartz()
#pdf("HeatmapPatchUp.pdf",width=8,height=17)
pheatmap(expc,color=col,cluster_cols=F,clustering_distance_rows="correlation") #plot the heatmap
#dev.off()


