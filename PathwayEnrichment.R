# Load necessary libraries. Uncomment and run these lines if the libraries are not already installed.
# if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat") # Install Seurat package if not already installed
# devtools::install_github('satijalab/seurat-data') # Install SeuratData package from GitHub if not already installed
# if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler") # Install clusterProfiler package if not already installed
# if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db") # Install org.Hs.eg.db package if not already installed
# if (!requireNamespace("DOSE", quietly = TRUE)) BiocManager::install("DOSE") # Install DOSE package if not already installed
# 
library(Seurat)          # Load Seurat library for single-cell RNA-seq analysis
library(SeuratData)      # Load SeuratData library for example datasets
library(clusterProfiler) # Load clusterProfiler for gene set enrichment analysis
library(org.Hs.eg.db)    # Load org.Hs.eg.db for human gene annotations
library(DOSE)            # Load DOSE for enrichment analysis visualization
library(dplyr)           # Load dplyr for data manipulation and transformation
# 
# # Load the PBMC 3k dataset
# pbmc <- readRDS("S:/Day4.rds")
# 
# # Preprocess and normalize the data
# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") # Calculate the percentage of mitochondrial gene expression. The main source of problem come from here. Mt genes are not all removed, only the one with high fraction, with a threshold. This will be enough, it will not give noise, while ribosomal genes always give noise, so you remove them all.
# ribosomal_genes <- grep("^RPS|^RPL", rownames(pbmc), value = TRUE)  # Identify ribosomal genes. ^ means start with.
# pbmc <- pbmc[!rownames(pbmc) %in% ribosomal_genes, ] # Remove ribosomal genes from the dataset. They are ALL filtered out.
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) # Filter cells based on gene count and mitochondrial percentage. 
# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) # Normalize data using log normalization
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) # Identify the top 2000 variable features
# all.genes <- rownames(pbmc) # Get all gene names
# pbmc <- ScaleData(pbmc, features = all.genes) # Scale the data for all features (genes)
# 
# # Perform PCA (Principal Component Analysis)
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) # Run PCA using the top variable features
# 
# # Find clusters of cells
# pbmc <- FindNeighbors(pbmc, graph.name ="test", dims = 1:10) # Build a nearest neighbor graph using the first 10 PCA dimensions
# pbmc <- FindClusters(pbmc, graph.name = "test", algorithm = 1, resolution = 0.5) # Cluster cells using Louvain algorithm with resolution of 0.5. If you wanted to use Leiden clustering, then you need to specify it.
# 
# # Run UMAP (Uniform Manifold Approximation and Projection) for visualization
# pbmc <- RunUMAP(pbmc, dims = 1:10) # Perform UMAP dimensionality reduction using the first 10 PCA dimensions
# DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 6) # Plot UMAP results with cluster labels
pbmc <- readRDS("S:/Day4_parsed.rds")
# # Identify marker genes for each cluster
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25) # Find marker genes for each cluster with specified thresholds
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) # For each cluster, select top 2 marker genes based on average log2 fold change: the top 2 most enriched genes for each cluster

##################################################################
# Subsetting top 100 markers with adjusted p-values lower than 0.05 #
##################################################################

# because if we do the pathway enrichment analysis on all genes it would be too heavy
top100 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) # Get top 100 markers per cluster based on average log2 fold change
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0) # Filter for markers with adjusted p-value < 0.05

#now we gt out of the Seraut object, by saving all info in a data.frame. Seraut and Scanpy don't perform gene enrichment analysis
# Prepare gene lists for enrichment analysis
df <- top100pval[,7:6] # Select columns containing gene names and cluster identifiers
dfsample <- split(df$gene, df$cluster) # Split gene list by cluster
length(dfsample) # Print the number of genes in each cluster

# Convert gene symbols to ENTREZ IDs for enrichment analysis
dfsample$`0` = bitr(dfsample$`0`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`1` = bitr(dfsample$`1`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`2` = bitr(dfsample$`2`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`3` = bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`4` = bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`5` = bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`6` = bitr(dfsample$`6`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`7` = bitr(dfsample$`7`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#A certain % of input genes map to fail. Sometimes a gene symbol is not associated to entrez ID, so you will lose this information. Generally they are pseudogenes or genes not well described.

# Create a list of gene sets for each cluster
genelist <- list("C0" = dfsample$`0`$ENTREZID, 
                 "C1" = dfsample$`1`$ENTREZID,
                 "C2" = dfsample$`2`$ENTREZID,
                 "C3" = dfsample$`3`$ENTREZID,
                 "C4" = dfsample$`4`$ENTREZID,
                 "C5" = dfsample$`5`$ENTREZID,
                 "C6" = dfsample$`6`$ENTREZID,
                 "C7" = dfsample$`7`$ENTREZID)

# until now we just prepared the list for the tool

## RUN ENRICHMENT ANALYSIS (NOT OVER REPRESENTATION ANALYSIS, there is no log2FC, just a list of genes) ()
# Perform Gene Ontology (GO) enrichment analysis
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.05) # Compare GO term enrichment across clusters
dotplot(GOclusterplot, font.size = 8) # Plot the GO enrichment results with adjusted font size

# Perform KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway enrichment analysis
KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", pvalueCutoff = 0.05) # Compare KEGG pathway enrichment across clusters
dotplot(KEGGclusterplot, font.size = 8) # Plot the KEGG enrichment results with adjusted font size

## RUN OVER REPRESENTATION ANALYSIS
# Perform differential expression analysis for cluster 0
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = FALSE) # Find differentially expressed genes in cluster 0
# or you could do FindAllMarkers and then take only the cluster of interest. In this case we also need the negative (downregulated) genes. So we will take all the genes here, while before we took only the positive.

# Prepare gene list for Gene Set Enrichment Analysis (GSEA)
rnk0 <- cluster0.markers[4] # Extract average log2 fold changes
gene_list <- bitr(row.names(rnk0), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Convert gene symbols to ENTREZ IDs
gene_list <- gene_list$ENTREZID # Extract ENTREZ IDs
original_gene_list <- rnk0$avg_log2FC # Extract log2 fold change values
names(original_gene_list) <- gene_list # Assign ENTREZ IDs as names to log2 fold change values
gene_list <- original_gene_list[!is.na(names(original_gene_list))] # Remove NA values
gene_list = sort(gene_list, decreasing = TRUE) # Sort genes by log2 fold change in decreasing order. If they are not sorted you will get an error

# Perform GSEA (Gene Set Enrichment Analysis)
gse <- gseGO(geneList = gene_list, 
             ont = "BP", # Use Biological Process ontology, this is the gene list
             keyType = "ENTREZID", 
             nPerm = 100, # Number of permutations for statistical testing. It would be good to hav 1k here. In this case there is 100 just to save time in the tutorial.
             minGSSize = 3, # Minimum size of gene sets
             maxGSSize = 800, # Maximum size of gene sets
             pvalueCutoff = 0.05, # P-value cutoff for significance
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none", # No multiple testing correction
             verbose = FALSE,# Do not print additional progress information
             seed = 42) #seed is important to reproduce the same results when you have permutations

# Convert GSEA results to readable format
gse_readble <- setReadable(gse, org.Hs.eg.db, 'ENTREZID') # Convert GSEA results for easy interpretation

# Plot GSEA results
require(DOSE) # Ensure DOSE library is loaded
dotplot(gse_readble, showCategory = 10, split = ".sign") + facet_grid(.~.sign) # Plot GSEA results with categories and split by sign


goplot(gse_readble) # Plot GSEA results as GO terms

# Plot GSEA results again (may be redundant)
#goplot(gse) # Plot GSEA results (may be redundant if gse_readble is already plotted)
