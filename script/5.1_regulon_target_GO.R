#define path and load required packages
setwd("/Users/jiahuiji/Dropbox/project/t-cell/tf_analysis")



#load packages
library(stringr)



#parameter
type="TC"



#-------Step 1: Load TF clusters-------
dic=paste0("./",type,"_tf_clusters/")
cluster_file=paste0(dic,type,"_sig_tf_clusters.csv")
tf_cluster=read.table(file=cluster_file,header=T,sep="\t")




#-------Step 2: load regulons-------
regulon_file=paste0(type,"_all_regulon.rds")
load(regulon_file)
regulon_name=names(regulons)
names(regulons)=gsub("[(+)]","",regulon_name)



#target regulon
group=unique(tf_cluster[,1])
for (i in 1:length(group))
{
	#extract TF clusters
	target=rownames(tf_cluster)[which(tf_cluster[,1]==group[i])]

	#extract regulon targets
	regulon_target=c()
	for (j in 1:length(target))
	{
		number=which(names(regulons)==target[j])
		regulon_target=append(regulon_target,regulons[number][[1]])
	}
	go_input=unique(regulon_target)

	write.table(go_input,file=paste0(type,"_cluster",group[i],"_stem_target.txt"),sep="\t",quote=F,row.names=F,col.names=F)
}




#-------Step 3: GO enrichment analysis-------
library(WebGestaltR)

#parameter
type="BC"
number_of_cluster=12

wd="/Users/jiahuiji/Dropbox/project/t-cell/tf_analysis/"
setwd(paste0(wd,type,"_tf_clusters"))


for (i in 1:number_of_cluster)
{
	gene_file=paste0("./genelist/",type,"_cluster",i,"_target.txt")
	backfile=paste0("./background/",type,"_background.txt")
	cluster=paste0(type,"_",i)

	gene=read.table(file=gene_file, header=F, sep="\t")
	gene_list=gene$V1

	print(gene_file)
	print(backfile)

	outputDirectory=getwd()
	WebGestaltR(enrichMethod="ORA",
            ORGANISM="mmusculus",
            enrichDatabase="geneontology_Biological_Process",
            interestGene=gene_list,
            interestGeneType="genesymbol",
            referenceGeneFile=backfile,
            referenceGeneType="genesymbol",
            minNum=1,
            maxNum=2000,
            sigMethod="fdr",
            isOutput=TRUE,
            outputDirectory=outputDirectory,
            projectName=paste0(cluster,'_bp'))


	WebGestaltR(enrichMethod="ORA",
	            ORGANISM="mmusculus",
	            enrichDatabase="disease_Disgenet",
	            interestGene=gene_list,
	            interestGeneType="genesymbol",
	            referenceGeneFile=backfile,
	            referenceGeneType="genesymbol",
	            minNum=1,
	            maxNum=2000,
	            sigMethod="fdr",
	            isOutput=TRUE,
	            outputDirectory=outputDirectory,
	            projectName=paste0(cluster,'_Disgenet'))


	WebGestaltR(enrichMethod="ORA",
	            ORGANISM="mmusculus",
	            enrichDatabase="disease_GLAD4U",
	            interestGene=gene_list,
	            interestGeneType="genesymbol",
	            referenceGeneFile=backfile,
	            referenceGeneType="genesymbol",
	            minNum=1,
	            maxNum=2000,
	            sigMethod="fdr",
	            isOutput=TRUE,
	            outputDirectory=outputDirectory,
	            projectName=paste0(cluster,'_GLAD4U'))


	WebGestaltR(enrichMethod="ORA",
	            ORGANISM="mmusculus",
	            enrichDatabase="disease_OMIM",
	            interestGene=gene_list,
	            interestGeneType="genesymbol",
	            referenceGeneFile=backfile,
	            referenceGeneType="genesymbol",
	            minNum=1,
	            maxNum=2000,
	            sigMethod="fdr",
	            isOutput=TRUE,
	            outputDirectory=outputDirectory,
	            projectName=paste0(cluster,'_OMIM'))


	WebGestaltR(enrichMethod="ORA",
	            ORGANISM="mmusculus",
	            enrichDatabase="pathway_Reactome",
	            interestGene=gene_list,
	            interestGeneType="genesymbol",
	            referenceGeneFile=backfile,
	            referenceGeneType="genesymbol",
	            minNum=1,
	            maxNum=2000,
	            sigMethod="fdr",
	            isOutput=TRUE,
	            outputDirectory=outputDirectory,
	            projectName=paste0(cluster,'_Reactome'))

}







#-------Step 4: Plot select regulon activity-------
#load mean auc for regulons
dic=paste0("./",type,"_tf_clusters/")
auc_file=paste0(dic,type,"_sig_tf_mean_auc.csv")
auc_score=read.table(file=auc_file,header=T,sep="\t")
mat=as.matrix(auc_score)


cluster=unique(tf_cluster[,1])

pdf(paste0(type,"_cluster_activity.pdf"))
par(mfrow=c(3,3)) 
for (i in 1:length(cluster))
{
	tf_name=rownames(tf_cluster)[which(tf_cluster[,1]==cluster[i])]
	
	range01=function(x){(x-min(x))/(max(x)-min(x))}
	select=t(apply(mat[tf_name,], 1, range01))

	plot(select[1,], type="o", pch=20, col="black", lty=1, ylab="Mean AUC", cex=0,
		xaxt = "n",
		xlab="Day",
		ylim=c(min(select),max(select)),
		main=paste0("cluster",cluster[i]))
	axis(1, at=1:7, labels=c(0,1,3,5,7,14,28))
	for (j in 1:length(tf_name))
	{
	lines(select[j,], type="o", pch=20, col = "black", cex=0)
	}
}
dev.off()















