#load packages
library(ggplot2)
library(ggpubr)
library(grid)
library(cowplot)



#parameter
cell="BC"
setwd(paste0("/Users/jiahuiji/Dropbox/project/t-cell/results/tf_analysis/",cell,"_tf_clusters"))






#grid.newpage()
frame=data.frame()

for (i in 1:7)
{
	#parameter
	cluster_num=i
	wd=paste0(cell,"_cluster",cluster_num)
	go_type="_bp"
	dic=paste0(getwd(),"/","Project","_",wd,go_type,"/")
	file_name=paste0("enrichment_results_",wd,go_type,".txt")




	#read in go table
	go_table=read.table(file=paste0(dic,file_name),sep="\t",header=T)
	#caculate -log10(FDR)
	go_table$bar=-log(go_table$FDR)
	table=go_table[order(go_table$bar),]
	#calculate number of bars over or smaller than 10
	len=length(rownames(go_table))
	define=10
	if (len > define)
		{num=define}
	else
		{num=len}




	#plot
	new=go_table[1:num,c("description","bar")]
	new$cluster=rep(wd,num)

	frame=rbind(frame,new)	
}



#modify Inf in the frame
frame$bar[18:24]=29.5
frame$bar[28:37]=19.5



reorder_within <- function(x, by, within, fun = mean, sep = " ", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}




pdf(paste0(cell,"_bar.pdf"))

My_Theme=theme(
axis.text.y = element_text(size=6)
)

p=ggplot(data=frame, 
aes(x=reorder_within(description,bar,cluster), y=bar)) +
geom_bar(stat="identity", width=0.5) +
labs(y="-log(FDR)", x="GO terms - biological process") +
facet_grid(cluster ~ ., scales="free_y", space="free_y", drop=TRUE) 

print(p+coord_flip()+My_Theme)

dev.off()



















for (i in 1:7)
{

	#parameter
	cluster_num=i
	wd=paste0("TC_cluster",cluster_num)
	go_type="_bp"
	dic=paste0(getwd(),"/","Project","_",wd,go_type,"/")
	file_name=paste0("enrichment_results_",wd,go_type,".txt")




	#read in go table
	go_table=read.table(file=paste0(dic,file_name),sep="\t",header=T)
	#caculate -log10(FDR)
	go_table$bar=-log(go_table$FDR)
	table=go_table[order(go_table$bar),]
	#calculate number of bars over or smaller than 10
	len=length(rownames(go_table))
	define=10
	if (len > define)
		{num=define}
	else
		{num=len}



	#plot
	p=ggplot(data=go_table[1:num,], 
	aes(x=reorder(description, bar) , y=bar)) +
	geom_bar(stat="identity", width=0.5) +
	theme_minimal() +
	labs(y="-log(FDR)", x=wd)




	
}









pdf("test2.pdf")
p+coord_flip()
dev.off()

define_region=function(row, col)
	{
  		viewport(layout.pos.row = row, layout.pos.col = col)
    } 


	print(p+coord_flip(), vp=define_region(row=cluster_num, col=1:2))






