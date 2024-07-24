setwd("/home/pkf/gendata2/bioinformatica/rcarriero/corso")

#Upload libraries

library(Seurat)
library(dplyr)
library(clustree)
library(ggraph)

#Upload Seurat object

tcells<- readRDS('S:/Day3/Treg.rds')
ncol(tcells)
nrow(tcells) #1132 cells and 26k genes

#Find variable features
tcells <- FindVariableFeatures(tcells, selection.method = "vst", nfeatures = 2000) #it extract 2k HVG on normalized counts (so on the data layer)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(tcells), 10)

# plot variable features with and without labels
pdf("plot_variable_features.pdf", width = 15, height = 10)
plot1 <- VariableFeaturePlot(tcells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

#Scale data
all.genes <- rownames(tcells)
tcells <- ScaleData(tcells, features = all.genes) #you need to specify the features, because default is HVG

#Run PCA
tcells <- RunPCA(tcells, features = VariableFeatures(object = tcells))

# Examine and visualize PCA results a few different ways
print(tcells[["pca"]], dims = 1:5, nfeatures = 5)


pdf("vizdim_pca_plot.pdf")
VizDimLoadings(tcells, dims = 1:2, reduction = "pca")
DimPlot(tcells, reduction = "pca") 
DimHeatmap(tcells, dims = 1, cells = 500, balanced = TRUE) #subset of 500 cells, the expr of genes that correlate with PC1, how cells separate based on the expression of those genes. This can be useful for example if you have a KO, to check that that gene is among the TOP
DimHeatmap(tcells, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#at the moment the active identity is responder or not, therefore the plots are automatically groupped by that identity.
#but you can change it in this way, and the plot will be colored differently
Idents(tcells) <-tcells$orig.ident
DimPlot(tcells, reduction = "pca") #otherwise you can specify group_by. But if you want to avoid that, you can have it by default setting it as active identity
#don't forget to set it back: Idents(tcells) <-tcells$'tcells@active.ident'

tcells <- JackStraw(tcells, num.replicate = 100) #to compute the p-value for each dimension
tcells <- ScoreJackStraw(tcells, dims = 1:20) #score each straw

pdf("JackStrawPlot.pdf")
JackStrawPlot(tcells, dims = 1:20) #>> we can take up to 14 (because 15 is not significant), or all of them except the non significant one
ElbowPlot(tcells, ndims = 40) #>> after 9, the line reaches a plateu, so you could stop after 9
dev.off()

#Select significant components
tcells <- FindNeighbors(tcells, dims = 1:14) #to identify the neighbours. Here you specify the number of dimension you selected from before

#check which is the best resolution
cells_clust_tree <- FindClusters(tcells, resolution = seq(from=0, to=2, by=0.1), print.output = 0, save.SNN = T) #assign at each cell the cluster of belonging

pdf("CL7_tregs_Resolution_tree_up_to_2.0_pca_1_17.pdf",  width=17, height=17)
clustree(cells_clust_tree)
dev.off()

#choose one resolution
tcells <- FindClusters(tcells, resolution = 0.4)
tcells <- RunUMAP(tcells, dims = 1:14) #you need to indicate again the number of dimensions because it's a recalculation

tcells$condition <- tcells$`tcells@active.ident`#if you want to plot also the responders
pdf("umap.pdf")
DimPlot(tcells, reduction = "umap", label = T)|
  DimPlot(tcells, reduction = "umap", label = T, group.by = "condition")
dev.off()

#find differentially expressed features:
#FindAllMakers calculate all the markers for each cluster. Done taking in consideration one cluster per time and the genes of that particular cluster against all the other clusters as one big cluster
tcells.markers <- FindAllMarkers(tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#you can also calculate the one corrected based on condition
tcells.markers.cor <- FindAllMarkers(tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, latent.vars = "condition")

#look at top genes for each cluster
tcells.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
tcells.markers %>% group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

pdf("heatmap.pdf")
DoHeatmap(tcells, features = top10$gene) + NoLegend()
dev.off()

#export it
write.table(tcells.markers, "tcells.markers.txt", sep="\t")

#Plot of marker genes, to visualize the expression of a gene of interest:

pdf("VlnPlot_FOXP3.pdf")
VlnPlot(tcells, features = 'FOXP3') #useful to see in which clusters is expressed a certain gene of interest
dev.off()

pdf("FeaturePlot_CXCL13.pdf")
FeaturePlot(tcells, features = 'CXCL13', label = T) #to color the UMAP based on a certain gene. 
dev.off()


pdf("DotPlot_IL1R2.pdf")
DotPlot(tcells, features = 'IL1R2') 
dev.off()

#you can also calculate markers of a cluster comparing it to a certain cluster or to a few clusters of choice
tcell.makers2 <- FindMarkers(tcells, ident.1 = "0", ident.2 = c("1", "2"), only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)


#Responder vs Non responder --> all cells. If interested in the markers that distinguish responders from non responders
Idents(tcells)

levels(tcells@meta.data$`tcells@active.ident`)[1] <-"Non_responder"
levels(tcells@meta.data$`tcells@active.ident`)[2] <-"Responder"
table(tcells@meta.data$`tcells@active.ident`)

Idents(tcells)<- tcells@meta.data$`tcells@active.ident`
R_vs_NR<- FindMarkers(tcells, ident.1 = 'Responder', ident.2 = 'Non_responder')
head(R_vs_NR)
write.table(R_vs_NR, file='R_vs_NR.txt', sep="\t")

#Responder vs Non responder--> cluster 0
Idents(tcells)<- tcells$RNA_snn_res.0.4
tcells_c0<- subset(tcells, idents = '0') #select only cells from cluster 0
Idents(tcells_c0)<-tcells_c0@meta.data$`tcells@active.ident`
Idents(tcells_c0)
c0_R_vs_NR<- FindMarkers(tcells_c0, ident.1 = 'Responder', ident.2 = 'Non_responder')
head(c0_R_vs_NR)
write.table(c0_R_vs_NR, file='c0_R_vs_NR.txt', sep="\t")

#this can be useful to recluster a certain cluster to see subpopulations within it. So you can export one and reanalyse it.