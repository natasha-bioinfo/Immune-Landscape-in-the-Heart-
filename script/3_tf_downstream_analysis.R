#basic setting
setwd("/Users/jiahuiji/Dropbox/project/t-cell/results/TF-analysis")

#load required packages
library(limma) #linear model fit 
library(VennDiagram) #venn plot
library(pheatmap) #heatmap
library(stringr) #deal with string
library(Seurat)

#read in scenic auc activity
type="BC"

dic="/Users/jiahuiji/Dropbox/project/t-cell/results/scenic_threecelltypes/"
file_name="_AUCMatrix.csv"
file_full=paste0(dic,type,file_name)

auc_mat=read.table(file=file_full,row.names="Cell",header=T,sep=",")





#-------step1: Distribution-------
dis=c()
for (i in 1:length(colnames(auc_mat)))
{
	dis=append(dis,mean(auc_mat[,i]))
}


dis_log=c()
for (i in 1:length(colnames(auc_mat)))
{
	dis_log=append(dis_log,mean(log(auc_mat)[,i]))
}

pdf(paste0(type,"_distribution.pdf"))
hist(dis, main="Histogram of mean AUC")
d=density(dis)
plot(d, main="Distribution of mean AUC")
dev.off()

pdf(paste0(type,"_log_distribution.pdf"))
hist(log(dis), main="Histogram of mean AUC")
d=density(log(dis))
plot(d, main="Distribution of mean AUC")
dev.off()





#-------step2: Differential expression analysis-------
#create continuous covariant
time_point=c()
for (i in 1:length(rownames(auc_mat)))
{
	time_point=append(time_point,strsplit(rownames(auc_mat)[i],"_")[[1]][1])
}

#design matrix for continuous covarient
#lmfit
design=model.matrix(~time_point)
mat=t(auc_mat)
fit=lmFit(mat,design)
fit2=eBayes(fit)
summary(decideTests(fit2))

for (i in 2:length(unique(group)$group))
{
    res=topTable(fit2,coef=i,n=Inf, adjust="fdr",sort.by="P")
    write.table(res,file=paste0(type,"_",unique(group)$group[i],"_factor_deg.csv"),sep="\t",quote=F)
}

sig_tf=c()
for (i in 2:length(unique(group)$group))
{
    res=topTable(fit2,coef=i,n=Inf, adjust="fdr",sort.by="P")
    sig_tf=append(sig_tf,rownames(res)[which(res$adj.P.Val<0.05)])
}
sig=unique(sig_tf)
save(sig,file=paste0(type,"_all_factor_deg.rds"))









