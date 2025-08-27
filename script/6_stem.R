setwd("/Users/jiahuiji/Dropbox/project/t-cell/tf_analysis")
type="BC"


#load required package
library(stringr)


#step 1 transfer AUC matrix for 
auc_matrix=read.csv(file=paste0(type,"_AUCMatrix.csv"),header=T,row.names="Cell")
mat=t(auc_matrix)

#create group for cells
group=c()
for (i in 1:length(colnames(mat)))
{
    group=append(group,strsplit(colnames(mat)[i],"_")[[1]][1])
}


#create matrix for stem analysis
stem=matrix(ncol=length(unique(group)),nrow=length(rownames(mat)))
colnames(stem)=unique(group)
rownames(stem)=rownames(mat)
for (i in 1:length(unique(group)))
{
	for (j in 1:length(rownames(mat)))
	{
		stem[j,i]=mean(mat[j,which(group==unique(group)[i])])
	}
}

write.table(stem,file=paste0(type,"_stem.txt"),sep="\t",quote=F)