#load Seurat
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('Seurat')
library(Seurat)
library(Seurat)
library(dplyr)
library(cowplot)
library("ggplot2")
library(tidygraph)
library(clustree)
##7DTP data 
# Create a Seurat object for all files
data10 <- Read10X(data.dir = "~/Desktop/7daytimepoint/10")
data10 <- Read10X(data.dir = "~/Desktop/7daytimepoint/10")
Day0 <- CreateSeuratObject(counts = data10, min.cells = 3, min.features  = 200, project = "D0", assay = "RNA")
data13 <- Read10X(data.dir = "~/Desktop/7daytimepoint/13")
Day1 <- CreateSeuratObject(counts = data13, min.cells = 3, min.features  = 200, project = "D1", assay = "RNA")
data14 <- Read10X(data.dir = "~/Desktop/7daytimepoint/14")
Day3 <- CreateSeuratObject(counts = data14, min.cells = 3, min.features  = 200, project = "D3", assay = "RNA")
data15 <- Read10X(data.dir = "~/Desktop/7daytimepoint/15")
Day5 <- CreateSeuratObject(counts = data15, min.cells = 3, min.features  = 200, project = "D5", assay = "RNA")
data16 <- Read10X(data.dir = "~/Desktop/7daytimepoint/16")
Day7 <- CreateSeuratObject(counts = data16, min.cells = 3, min.features  = 200, project = "D7", assay = "RNA")
data17 <- Read10X(data.dir = "~/Desktop/7daytimepoint/17")
Day14 <- CreateSeuratObject(counts = data17, min.cells = 3, min.features  = 200, project = "D14", assay = "RNA")
data18 <- Read10X(data.dir = "~/Desktop/7daytimepoint/18")
Day28 <- CreateSeuratObject(counts = data18, min.cells = 3, min.features  = 200, project = "D28", assay = "RNA")
all.data <- merge(Day0, y = c(Day1, Day3, Day5, Day7, Day14, Day28 ), add.cell.ids = c("Day0","Day1","Day3","Day5", "Day7", "Day14", "Day28"), project = "all.data")
# mitochondrial genes
mito.genes <- grep(pattern = "^mt-", x = rownames(all.data@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(all.data@assays[["RNA"]][mito.genes, ])/Matrix::colSums(all.data@assays[["RNA"]])
all.data <- AddMetaData(object = all.data, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = all.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
FeatureScatter(object = all.data, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = all.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
all.data <- subset(x = all.data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000 & percent.mito >  -Inf & percent.mito < 0.125)
#Normalizing
norm.data <- NormalizeData(object = all.data, normalization.method = "LogNormalize", scale.factor = 10000)
#Scaling
norm.data <- FindVariableFeatures(object = norm.data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
Scale.data <- ScaleData(object = norm.data, vars.to.regress = c("nCount_RNA", "percent.mito"))
#PCA
Scale.data <- RunPCA(object = Scale.data,  npcs = 20, verbose = FALSE)
Scale.data <- FindNeighbors(Scale.data, reduction = "pca", dims = 1:20)
Scale.data<- FindClusters(Scale.data, resolution = 0.5, algorithm = 1)
#UMAP
Scale.data <- RunUMAP(Scale.data, reduction = "pca", dims = 1:20)
DimPlot(Scale.data, reduction = "umap", label = TRUE)
DimPlot(Scale.data, reduction = "umap")
FeaturePlot(Scale.data, features = c("Ptprc", "Ms4a1", "Cd79a", "Cd3e","Cd3g", "Cd3d"))
cluster.averages <- AverageExpression(Scale.data, return.seurat=TRUE)
Bcell<- subset(x = Bcell, idents = c('6',"22"))
Tcell<- subset(x = Tcell, idents = c("13"))
#Bcellscollagen
collagen.marker <- c("Cd34", "Fbln1", "Fbln2", "Cd3e", "Cd3d", "Cd3g", "Cd79a", "Ms4a1", "Col3a1","Col1a2","Col1a1","Col4a1","Col5a1","Col5a2","Col8a1","Col6a3","Col15a1","Col14a1","Col5a3","Col4a2")
DoHeatmap(object = cluster.averages, features = collagen.marker, draw.lines=FALSE, size = 3.5, label = TRUE)
#Tcellsub
Tcellsub <- c("Cd8b1", "Cd8a","Dapl1", "Klrd1", "Il7r","Ly6c2","Klrg1","Cxcr3","Il2rb","Klf2","Jund", "Cd4","Ly6c1", "Tcf7","Actb","Ikzf2","Il2ra","Ctla4","Tnfrsf18","Tnfrsf4","Izumo1r","Nrp1","Itgae","Il10","Cxcr6", "Il17a","Il17f","Il1r1","Maf","Rorc","Csf2","S100a6","Klra8","Ncr1","Nkg7","Tyrobp","Klrk1","Klre1","Klrb1c","Ifng","Id2","Il4","Ifitm3","Ifitm2","Bcl6","Cxcr5","Id3","Ube2c","Ccnb2","Tk1","Cd79a","Ms4a1","Ly6d")
#Generate DEG for webgestalt 
##TCELL
##Day0 vs Day1 comparison 
TcellsD0D1<- subset(x = Tcell, idents = c("D0", "D1"))
DimPlot(TcellsD0D1, reduction = "umap", label = TRUE, pt.size = .1)
TcellsD0D1.markers <- FindAllMarkers(object = TcellsD0D1, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(TcellsD0D1.markers, sep = '\t', quote = FALSE, file = "TcellsD0D1.markers.txt")
##Day0 vs Day3 comparison 
TcellsD0D3<- subset(x = Tcell, idents = c("D0", "D3"))
DimPlot(TcellsD0D3, reduction = "umap", label = TRUE, pt.size = .1)
TcellsD0D3.markers <- FindAllMarkers(object = TcellsD0D3, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(TcellsD0D3.markers, sep = '\t', quote = FALSE, file = "TcellsD0D3.markers.txt")
##Day0 vs Day5 comparison 
TcellsD0D5<- subset(x = Tcell, idents = c("D0", "D5"))
DimPlot(TcellsD0D5, reduction = "umap", label = TRUE, pt.size = .1)
TcellsD0D5.markers <- FindAllMarkers(object = TcellsD0D5, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(TcellsD0D5.markers, sep = '\t', quote = FALSE, file = "TcellsD0D5.markers.txt")
##Day0 vs Day7 comparison 
TcellsD0D7<- subset(x = Tcell, idents = c("D0", "D7"))
DimPlot(TcellsD0D7, reduction = "umap", label = TRUE, pt.size = .1)
TcellsD0D7.markers <- FindAllMarkers(object = TcellsD0D7, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(TcellsD0D7.markers, sep = '\t', quote = FALSE, file = "TcellsD0D7.markers.txt")
##Day0 vs Day14 comparison 
TcellsD0D14<- subset(x = Tcell, idents = c("D0", "D14"))
DimPlot(TcellsD0D14, reduction = "umap", label = TRUE, pt.size = .1)
TcellsD0D14.markers <- FindAllMarkers(object = TcellsD0D14, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(TcellsD0D14.markers, sep = '\t', quote = FALSE, file = "TcellsD0D14.markers.txt")
##Day0 vs Day28 comparison 
TcellsD0D28<- subset(x = Tcell, idents = c("D0", "D28"))
DimPlot(TcellsD0D28, reduction = "umap", label = TRUE, pt.size = .1)
TcellsD0D28.markers <- FindAllMarkers(object = TcellsD0D28, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(TcellsD0D28.markers, sep = '\t', quote = FALSE, file = "TcellsD0D28.markers.txt")
##BCELL
##Day0 vs Day1 comparison 
BcellsD0D1<- subset(x = Bcell, idents = c("D0", "D1"))
DimPlot(BcellsD0D1, reduction = "umap", label = TRUE, pt.size = .1)
BcellsD0D1.markers <- FindAllMarkers(object = BcellsD0D1, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(BcellsD0D1.markers, sep = '\t', quote = FALSE, file = "BcellsD0D1.markers.txt")
##Day0 vs Day3 comparison 
BcellsD0D3<- subset(x = Bcell, idents = c("D0", "D3"))
DimPlot(BcellsD0D3, reduction = "umap", label = TRUE, pt.size = .1)
BcellsD0D3.markers <- FindAllMarkers(object = BcellsD0D3, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(BcellsD0D3.markers, sep = '\t', quote = FALSE, file = "BcellsD0D3.markers.txt")
##Day0 vs Day5 comparison 
BcellsD0D5<- subset(x = Bcell, idents = c("D0", "D5"))
DimPlot(BcellsD0D5, reduction = "umap", label = TRUE, pt.size = .1)
BcellsD0D5.markers <- FindAllMarkers(object = BcellsD0D5, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(BcellsD0D5.markers, sep = '\t', quote = FALSE, file = "BcellsD0D5.markers.txt")
##Day0 vs Day7 comparison 
BcellsD0D7<- subset(x = Bcell, idents = c("D0", "D7"))
DimPlot(BcellsD0D7, reduction = "umap", label = TRUE, pt.size = .1)
BcellsD0D7.markers <- FindAllMarkers(object = BcellsD0D7, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(BcellsD0D7.markers, sep = '\t', quote = FALSE, file = "BcellsD0D7.markers.txt")
##Day0 vs Day14 comparison 
BcellsD0D14<- subset(x = Bcell, idents = c("D0", "D14"))
DimPlot(BcellsD0D14, reduction = "umap", label = TRUE, pt.size = .1)
BcellsD0D14.markers <- FindAllMarkers(object = BcellsD0D14, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(BcellsD0D14.markers, sep = '\t', quote = FALSE, file = "BcellsD0D14.markers.txt")
##Day0 vs Day28 comparison 
BcellsD0D28<- subset(x = Bcell, idents = c("D0", "D28"))
DimPlot(BcellsD0D28, reduction = "umap", label = TRUE, pt.size = .1)
BcellsD0D28.markers <- FindAllMarkers(object = BcellsD0D28, only.pos = TRUE, min.pct = 0.20, thresh.use = 0.25)
write.table(BcellsD0D28.markers, sep = '\t', quote = FALSE, file = "BcellsD0D28.markers.txt")

##healthy human heart
load("~/Desktop/data/human_immune_b.rds")
all.data <- human_immune
norm.data <- NormalizeData(object = all.data, normalization.method = "LogNormalize", scale.factor = 10000)
#Scaling
norm.data <- FindVariableFeatures(object = norm.data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
Scale.data <- ScaleData(object = norm.data, vars.to.regress = c("nCount_RNA", "percent.mito"))
#PCA 
Scale.data <- RunPCA(object = Scale.data,  npcs = 20, verbose = FALSE)
Scale.data <- FindNeighbors(Scale.data, reduction = "pca", dims = 1:20)
Scale.data<- FindClusters(Scale.data, resolution = 0.5, algorithm = 1)
#UMAP
Scale.data <- RunUMAP(Scale.data, reduction = "pca", dims = 1:20)
DimPlot(Scale.data, reduction = "umap", label = TRUE)

DimPlot(human_immune_b, reduction = "umap", label = FALSE, pt.size = .1)
Tcellsub <- c("CD8A","CD8B","KLF2","JUND","IL32","GZMB","GZMA","GZMH","GZMK","IL2RG","CD69","CXCR4","CD4","IL7R","CCR7","TNFAIP3","FOXP3","IL2RA","CTLA4","TNFRSF18","TNFRSF4","CD52","SLC25A6","LGALS3","HMGB1","RORC","TYROBP","CD79A","MS4A1","SLC2A3","MT-CO2","MT-CO1")
bcellsub <- C("JUN","CD83","CD86","CD70","JUNB","GZMB","CCL4","CD3C","IL10","TGFB1","MS4A1","CD74","IGHD","CD79B","HVCN1","IGHG1","IGHG2","IGHG3","IGHG4","CD38","IGHM","XBP1","BACH2","PAX5","BCL2","CD79A")





