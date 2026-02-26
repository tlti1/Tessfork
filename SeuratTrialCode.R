# Seurat analysis to generate UMAP used for Shiny graphs as a template for when we have real analysis.
# Last edited: 26/02/26. Tom
# Seurat UMAP analysis code is taken from: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# Analysis is currently mostly identical to the tutorial code, using the sample data they provide. (Place in same directory as app.R)
# R Shiny code has now been moved out of this script and into two new scripts (ui.R and server.R, both in the ShinyApp directory)
# Run this code to generate the 'output.rds' file needed to be used by the shiny app. Once this has been generated, you don't need to run this code again.
# This code will eventually be replaced with the real project code that will also generate a .rds file, meaning all of the shiny code work just fine without changes.
# Future goals:
#   - Reduce Seurat analysis code to be only whats needed for Shiny (Done)
#   - Improve Shiny code to allow the user to search for individual cells
#   - Change layout of Shiny to be ready for displaying multiple graphs on each tab
#   - Consider splitting server and ui code into separate scripts (Done)
#   - See if Seurat code can be on separate script (Or find a way of avoiding script running multiple times) (Done)

library(rsconnect)
library(dplyr)
library(Seurat)
library(patchwork)

# This directory sometimes needs to be changed for no understood reason.
# Either "./ShinyApp/filtered_gene_bc_matrices/hg19"
# Or "./filtered_gene_bc_matrices/hg19"
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = 'pbmc3k', min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

SaveSeuratRds(pbmc, file="./ShinyApp/output.rds")