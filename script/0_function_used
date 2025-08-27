#regulon_to_network
#regulon_to_network_json
#deg_module_frame_to_list
#overlap_regulon_module
#community_GO


####### Required pacakge
library(Seurat) #open seurat object
library(monocle3) #trajectory analysis
library(igraph) #network analysis 










#######
#######
####### Transfer regulon to network (function)
####### Return and save network from regulon (result)
####### regulons: regulon list file extracted from loom file
####### 
#######


regulon_to_network=function(regulons){
	#deal regulon list
	regulons.name=names(regulons)
	TF.names=c()
	for (i in 1:length(regulons.name))
	{
		TF.names.i=c()
		TF.names.i[1:length(regulons[[i]])]=c(regulons.name[i])
		TF.names=append(TF.names,TF.names.i)
	}

	#change (+) to ""
	q=strsplit(TF.names,split=" ")
	TF.names=c()
	for (i in 1:length(q))
	{
		TF.names=append(TF.names,q[[i]][1])
	}

	TF.target=c()
	for (i in 1:length(regulons.name))
	{
		TF.target.i=c()
		TF.target.i[1:length(regulons[[i]])]=c(regulons[[i]][1:length(regulons[[i]])])
		TF.target=append(TF.target,TF.target.i)
	}

	interaction=c()
	for (i in 1:length(TF.names))
	{
		interaction=append(interaction,TF.names[i])
		interaction=append(interaction,TF.target[i])
	}

	network=graph(interaction) #tf-gene network 

	#save file
	output=paste0(Tcell_day,"_scenic_network.gml")
	write_graph(network,output,format=c('gml'))
	return(network)

}

####### done!











#######
#######
####### Transfer regulon to network, json version (function)
####### Return and save network from regulon (result)
####### regulons: regulon list file extracted from loom file
####### 
#######


regulon_to_network_json=function(regulons){
	#deal regulon list
	regulons.name=names(regulons)
	TF.names=c()
	for (i in 1:length(regulons.name))
	{
		TF.names.i=c()
		TF.names.i[1:length(regulons[[i]])]=c(regulons.name[i])
		TF.names=append(TF.names,TF.names.i)
	}

	#change (+) to ""
	q=TF.names
	TF.names=c()
	for (i in 1:length(q))
	{
		num=length(strsplit(q[i],split="")[[1]])
		name=substr(q[i],1,num-3)
		TF.names=append(TF.names,name)
	}

	TF.target=c()
	for (i in 1:length(regulons.name))
	{
		TF.target.i=c()
		TF.target.i[1:length(regulons[[i]])]=c(regulons[[i]][1:length(regulons[[i]])])
		TF.target=append(TF.target,TF.target.i)
	}

	interaction=c()
	for (i in 1:length(TF.names))
	{
		interaction=append(interaction,TF.names[i])
		interaction=append(interaction,TF.target[i])
	}

	network=graph(interaction) #tf-gene network 

	#save file
	output=paste0(Tcell_day,"_scenic_network.gml")
	write_graph(network,output,format=c('gml'))
	return(network)

}

####### done!

















#######
#######
####### Transfer deg module data frame from monocle3 into list (function)
####### Return a list of module information (result)
####### gene_module_df: deg module data frame from the output of monocle3
####### 
#######


deg_module_frame_to_list=function(gene_module_df)
{
	module_max=max(as.numeric(levels(gene_module_df$module)))

	deg_module_list=c()
	for (i in 1:module_max)
	{
		module=names(gene_module_df$module)[which(gene_module_df$module==i)]
		deg_module_list[[i]]=module
	}
	names(deg_module_list)=levels(gene_module_df$module)
	return(deg_module_list)
}

####### done!

















#######
#######
####### Conduct fisher's exact test to overlap regulon and monocle DEG modules (function)
####### Return a matrix of adjusted p-values (result)
####### regulons: regulon list extracted from loom file
####### deg_module_list: deg module list from function deg_module_frame_to_list
####### tcell_background: background gene list
#######
#######


overlap_regulon_module=function(regulons,deg_module_list,tcell_background)
{
	overlap_mat=matrix(nrow=length(regulons),ncol=length(deg_module_list))
	rownames(overlap_mat)=names(regulons)
	colnames(overlap_mat)=names(deg_module_list)

	for (i in 1:length(regulons))
	{
		#fisher matrix
		mat=matrix(nrow=2,ncol=2)

		for (j in 1:length(deg_module_list))
		{
			mat[1,1]=length(intersect(regulons[[i]],deg_module_list[[j]]))
			mat[1,2]=length(setdiff(regulons[[i]],deg_module_list[[j]]))
	        mat[2,1]=length(setdiff(deg_module_list[[j]],regulons[[i]]))
	        mat[2,2]=length(tcell_background)-mat[1,1]-mat[1,2]-mat[2,1]

	        overlap_mat[i,j]=fisher.test(mat,alternative="greater")$p.value
		}
	}

	overlap_mat_adj=matrix(p.adjust(overlap_mat,method="BH"),nrow=length(regulons),ncol=length(deg_module_list))
	rownames(overlap_mat_adj)=names(regulons)
	colnames(overlap_mat_adj)=names(deg_module_list)
	return(overlap_mat_adj)
}

####### done!








#######
#######
####### Compare clustering, calculate Jaccard index (function)
####### Save jaccard index results (result)
####### clus_1, clus_2: matrix with row per gene, column per graph indicating cluster membership
####### name: file name, string              
#######
#######


compare_cluster=function(clus_1,clus_2,name)
{
	row_num=length(unique(clus_1[,1]))
	col_num=length(unique(clus_2[,1]))
	jar_mat=matrix(nrow=row_num,ncol=col_num)

	for (i in 1:row_num)
	{
		for (j in 1:col_num)
		{
			gene_1=rownames(clus_1)[which(clus_1[,1]==i)]
			gene_2=rownames(clus_2)[which(clus_2[,1]==j)]

			a=intersect(gene_1,gene_2)
	        b=union(gene_1,gene_2)

	        jar_mat[i,j] = length(a)/length(b)
		}
	}

	write.table(jar_mat,file=paste0(name,"_jaccard.csv"),sep="\t",quote=F)

}

####### done!
