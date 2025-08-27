#define path and load required packages
setwd("/Users/jiahuiji/Dropbox/project/t-cell/results/tf_analysis")

library(stringr)
library(pheatmap) #heatmap
library(Seurat) #seurat


#parameter
type="BC"
folder=paste0("./",type,"_against0_wil_deg/")




#-------Step 1: read in degs from all time points------- 
time=c("Day1","Day3","Day5","Day7","Day14","Day28")
deg_all=c()
for (i in 1:length(time))
{
	deg=read.table(file=paste0(folder,time[i],"vsDay0_wil_deg.csv"),sep="\t",header=T)
	deg_sig=rownames(deg)[which(deg$p_val_adj<0.05)]
	deg_all=append(deg_all,deg_sig)
}
select_deg=unique(deg_all)







#-------Step 2: read in AUC matrix-------
auc_matrix=read.csv(file=paste0(type,"_AUCMatrix.csv"),header=T,row.names="Cell")
mat=t(auc_matrix)

rownames(mat)=gsub("[...]","",colnames(auc_matrix))







#-------Step 3: calculate mean auc activity score of differential TF-------
mat_select=mat[select_deg,]

#create group for cells
group=c()
for (i in 1:length(colnames(mat_select)))
{
    group=append(group,strsplit(colnames(mat_select)[i],"_")[[1]][1])
}


#create matrix for stem analysis
stem=matrix(ncol=length(unique(group)),nrow=length(rownames(mat_select)))
colnames(stem)=unique(group)
rownames(stem)=rownames(mat_select)
for (i in 1:length(unique(group)))
{
	for (j in 1:length(rownames(mat_select)))
	{
		stem[j,i]=mean(mat_select[j,which(group==unique(group)[i])])
	}
}

mean_activity_mat=stem
write.table(mean_activity_mat,file=paste0(type,"_sig_tf_mean_auc.csv"),sep="\t",quote=F)









#-------Step 4: cluster and plot-------
wd="/Users/jiahuiji/Dropbox/project/t-cell/results/tf_analysis/heatmap_tf_v2"
setwd(wd)
mean_activity_mat=read.table(file="BC_sig_tf_mean_auc.csv",sep="\t",header=T)
type="BC"

library(pheatmap)
library(Seurat)

pdf(paste0(type,"_sig_tf_heatmap_v2.pdf"))
pheatmap(t(mean_activity_mat), 
         cluster_cols=T, 
         cluster_rows=F, 
         fontsize_row=10,
         fontsize_col=2,
         show_colnames=T,
         show_rownames=T,
         scale="column",
         cellheight=20,
         border_color=NA,
         col=PurpleAndYellow(100),
         cutree_cols=12,
         treeheight_row=0,
         treeheight_col=0)
dev.off()

pdf(paste0(type,"_sig_tf_heatmap_nolabel_v2.pdf"))
pheatmap(t(mean_activity_mat), 
         cluster_cols=T, 
         cluster_rows=F, 
         fontsize_row=10,
         fontsize_col=5,
         show_colnames=F,
         show_rownames=T,
         scale="column",
         cellheight=20,
         border_color=NA,
         col=PurpleAndYellow(100),
         cutree_cols=12,
         treeheight_row=0,
         treeheight_col=0)
dev.off()




#extract heatmap generated clusters
out=pheatmap(mean_activity_mat, 
         cluster_cols=F, 
         cluster_rows=TRUE, 
         fontsize_row=4,
         show_colnames=T,
         scale="row",
         cellheight=4,
         border_color=NA,
         col=PurpleAndYellow(100))


tf_cluster=sort(cutree(out$tree_row, k=12))
write.table(tf_cluster,file=paste0(type,"_sig_tf_clusters.csv"),sep="\t",quote=F)










